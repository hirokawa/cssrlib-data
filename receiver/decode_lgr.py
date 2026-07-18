#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Qascom LuGRE messages decoder

 [1] Quascom, LuGRE Receiver Interface Control Document Issue 2.0, 2025

@author Rui Hirokawa
"""

from skyfield.framelib import itrs
from skyfield.api import load
import argparse
from glob import glob
import multiprocessing as mp
import numpy as np
import os
from pathlib import Path
import struct as st
from crccheck.crc import Crc24LteA
from cssrlib.gnss import uGNSS, uTYP, prn2sat, Obs, rSigRnx, gpst2time, \
    uSIG, rCST, gtime_t, time2gpst
from cssrlib.rawnav import rcvDec, rcvOpt


class navRec():
    """ class for navigation record """

    def __init__(self):
        self.t = gtime_t()
        self.nsat = 0
        self.pos = np.zeros(3)
        self.vel = np.zeros(3)
        self.std = np.zeros(3)
        self.clkb = 0.0
        self.clkd = 0.0
        self.ggto = 0.0
        self.dops = np.zeros(5)


class lgr(rcvDec):
    """ class for LuGRE Binary Format decoder """
    tow = -1
    week = -1

    def __init__(self, opt=None, prefix='', gnss_t='GE'):
        super().__init__(opt, prefix, gnss_t)

        sig_tbl = {
            uGNSS.GPS: {0: uSIG.L1C, 1: uSIG.L5Q},
            uGNSS.GAL: {2: uSIG.L1C, 3: uSIG.L5Q, 4: uSIG.L7Q}
        }

        obs_tbl = {
            uGNSS.GPS: [uSIG.L1C, uSIG.L5Q],
            uGNSS.GAL: [uSIG.L1C, uSIG.L5Q, uSIG.L7Q]
        }

        self.sig_t = {}
        for sys in sig_tbl.keys():
            if sys not in self.sig_tab.keys():
                continue

            self.sig_t[sys] = {}
            for key in sig_tbl[sys].keys():
                sig = sig_tbl[sys][key]
                self.sig_t[sys][key] = {
                    uTYP.C: rSigRnx(sys, uTYP.C, sig),
                    uTYP.L: rSigRnx(sys, uTYP.L, sig),
                    uTYP.D: rSigRnx(sys, uTYP.D, sig),
                    uTYP.S: rSigRnx(sys, uTYP.S, sig),
                }

        for sys in sig_tbl.keys():
            if sys not in self.sig_tab.keys():
                continue
            self.sig_tab[sys][uTYP.C] = []
            self.sig_tab[sys][uTYP.L] = []
            self.sig_tab[sys][uTYP.D] = []
            self.sig_tab[sys][uTYP.S] = []
            for sig in obs_tbl[sys]:
                self.sig_tab[sys][uTYP.C].append(rSigRnx(sys, uTYP.C, sig))
                self.sig_tab[sys][uTYP.L].append(rSigRnx(sys, uTYP.L, sig))
                self.sig_tab[sys][uTYP.D].append(rSigRnx(sys, uTYP.D, sig))
                self.sig_tab[sys][uTYP.S].append(rSigRnx(sys, uTYP.S, sig))

        self.nav = navRec()

        self.fn = open('lgr-nav.txt', 'w')

        self.fn.write("# week, tow, nsat, x, y, z, vx, vy, vz, cb, cd, " +
                      "ggto, sigp, sigv, sigt, pdop, hdop, vdop\n")

    def sync(self, buff, k):
        return buff[k] == 0x71  # 'q'

    def msg_len(self, msg, k):
        self.len = st.unpack_from('<L', msg, k+6)[0]  # payload length
        self.dlen = self.len+13
        return self.dlen

    def checksum(self, msg, k, maxlen=0):
        """ check the checksum of the message """
        len_ = self.dlen
        if len_ < 23:
            return False
        if maxlen > 0 and k+len_ >= maxlen:
            return False
        cs = Crc24LteA.calc(msg[k:k+len_-3])
        cs_ = msg[len_-1] << 16 | msg[len_-2] << 8 | msg[len_-3]

        # if self.monlevel > 0 and cs != cs_:
        #    print(f"checksum error: len={len_}")
        self.len = len_
        self.dlen = len_
        # return cs == cs_
        return True

    def svid2prn(self, svid, sigid):
        """ convert from svid/sigid to sys/prn """
        sys = uGNSS.GPS if sigid <= 1 else uGNSS.GAL
        prn = svid
        return sys, prn

    def decode_iqs(self, buff, k=10):
        """ decode I/Q Sample message """
        trx = st.unpack_from('<d', buff, k)[0]
        k += 8
        ns, blk, stype, sinv, qb = st.unpack_from('<LLBHB', buff, k)
        k += 12
        sf, f0, if_, bw = st.unpack_from('<dddd', buff, k)
        k += 32

        batchID = blk & 0xf
        nblk = (blk >> 4) & 0x3fff
        bid = (blk >> 18) & 0x3fff

        # stype 0:1ch,real,1:1ch,complex,2:2ch,real,3:2ch,complex
        # ns*qb/4 is stype = complex, ceil(ns*qb/8) otherwise
        if stype == 1 or stype == 3:
            sz = ns*qb/4
        else:
            sz = int(np.ceil(ns*qb/8))

        iqSamples = buff[k:k+sz]

    def decode_nav(self, buff, k=10):
        """ decode NAV message """
        trx = st.unpack_from('<d', buff, k)[0]
        k += 8
        appname = buff[k:k+10]
        k += 10
        wn, tow, sec = st.unpack_from('<HLf', buff, k)
        k += 14
        tow = float(tow)+sec

        # wn = int(trx//rCST.WEEK_SEC)
        # tow = trx-wn*rCST.WEEK_SEC

        # obs.time = gpst2time(wn, tow)
        self.tow = tow
        self.week = wn

        nsat, x, y, z, vx, vy, vz = st.unpack_from('<Hdddddd', buff, k)
        k += 50
        stdp, stdv, stdt, clkb, clkd, ggto = st.unpack_from('<ffdddd', buff, k)
        k += 40
        gdop, pdop, hdop, vdop, tdop = st.unpack_from('<ddddd', buff, k)
        k += 40

        self.nav.t = gpst2time(wn, tow)
        self.nav.nsat = nsat
        self.nav.pos = np.array([x, y, z])
        self.nav.vel = np.array([vx, vy, vz])
        self.nav.clkb = clkb
        self.nav.clkd = clkd
        self.nav.ggto = ggto
        self.nav.std = np.array([stdp, stdv, stdt])
        self.nav.dops = np.array([gdop, pdop, hdop, vdop, tdop])

        return self.nav

    def output_nav(self, nav):
        f = self.fn
        week, tow = time2gpst(nav.t)
        f.write(f"{week:4d} {tow:6.1f} {nav.nsat:2d} ")
        f.write(f"{nav.pos[0]:11.1f} {nav.pos[1]:11.1f} {nav.pos[2]:11.1f} ")
        f.write(f"{nav.vel[0]:8.1f} {nav.vel[1]:8.1f} {nav.vel[2]:8.1f} ")
        f.write(f"{nav.clkb:10.1f} {nav.clkd:8.1f} {nav.ggto:4.1f} ")
        f.write(f"{nav.std[0]:8.1f} {nav.std[1]:8.1f} {nav.std[2]:8.1f} ")
        f.write(f"{nav.dops[1]:8.1f} {nav.dops[2]:8.1f} {nav.dops[3]:8.1f}\n")

    def decode_obs(self, buff, k=10):
        """ decode RAW message """
        obs = Obs()
        obs.sig = self.sig_tab

        nsig_max = 0
        for s in self.sig_tab:
            if len(self.sig_tab[s][uTYP.L]) > nsig_max:
                nsig_max = len(self.sig_tab[s][uTYP.L])

        self.nsig[uTYP.C] = nsig_max
        self.nsig[uTYP.L] = nsig_max
        self.nsig[uTYP.D] = nsig_max
        self.nsig[uTYP.S] = nsig_max

        obs.sat = np.empty(0, dtype=np.int32)

        trx, nm = st.unpack_from('<dH', buff, k)
        k += 10

        wn = int(trx//rCST.WEEK_SEC)
        tow = trx-wn*rCST.WEEK_SEC

        obs.time = gpst2time(wn, tow)
        self.tow = tow
        self.week = wn

        pr = {}
        cp = {}
        dp = {}
        lli = {}
        cn = {}

        for i in range(nm):
            sigid, svid, dop_, ddop_, cp_, pr_, cn_ = \
                st.unpack_from('<BHddddd', buff, k)
            k += 43

            if np.isnan(cn_):
                continue

            sys, prn = self.svid2prn(svid, sigid)

            if sys not in self.sig_tab:
                continue

            sat = prn2sat(sys, prn)
            sig = self.sig_t[sys][sigid][uTYP.L]

            if sig not in self.sig_tab[sys][sig.typ]:
                if self.monlevel > 1:
                    print("skip code={:}".format(sig))
                continue
            idx = self.sig_tab[sys][sig.typ].index(sig)
            wavelength = rCST.CLIGHT/sig.frequency()

            if sat not in pr.keys():
                pr[sat] = {}
                cp[sat] = {}
                dp[sat] = {}
                lli[sat] = {}
                cn[sat] = {}

            lli_ = 0

            pr[sat][idx] = pr_
            cp[sat][idx] = cp_/wavelength  # [m]->[cycle]
            dp[sat][idx] = dop_
            lli[sat][idx] = lli_
            cn[sat][idx] = cn_

            if sat not in obs.sat:
                obs.sat = np.append(obs.sat, sat)

            # print(f"OBS {gnss}:{prn:3d} {svid} {sigid:2d}")

        nsat = len(obs.sat)
        if nsat == 0:
            return None

        obs.sat.sort()
        obs.P = np.zeros((nsat, self.nsig[uTYP.C]), dtype=np.float64)
        obs.L = np.zeros((nsat, self.nsig[uTYP.L]), dtype=np.float64)
        obs.D = np.zeros((nsat, self.nsig[uTYP.D]), dtype=np.float64)
        obs.S = np.zeros((nsat, self.nsig[uTYP.S]), dtype=np.float64)
        obs.lli = np.zeros((nsat, self.nsig[uTYP.L]), dtype=np.int32)

        for k, sat in enumerate(obs.sat):
            for i in pr[sat].keys():
                obs.P[k][i] = pr[sat][i]
                obs.L[k][i] = cp[sat][i]
                obs.D[k][i] = dp[sat][i]
                obs.S[k][i] = cn[sat][i]
                obs.lli[k][i] = lli[sat][i]

        return obs

    def decode(self, buff, len_, sys=[], prn=[]):
        mt = buff[1:4]
        k = 4
        sender, dlen = st.unpack_from('<HL', buff, k)
        k += 6

        # print(f"cls={cls_:2x} id={id_:2x}")

        if mt == b'RAW':  # RAW
            obs = self.decode_obs(buff)
            if self.flg_rnxobs and obs is not None:
                self.re.rnx_obs_header(obs.time, self.fh_rnxobs)
                self.re.rnx_obs_body(obs, self.fh_rnxobs)
        elif mt == b'NAV':
            nav = self.decode_nav(buff)
            self.output_nav(nav)
        elif mt == b'IQS':
            self.decode_iqs(buff)
        return 0


def decode(f, opt, args):

    print("Decoding {}".format(f))

    bdir, fname = os.path.split(f)

    prefix = fname.removesuffix('.bin')[-4:]+'_'
    prefix = str(Path(bdir) / prefix) if bdir else prefix
    lgrdec = lgr(opt, prefix=prefix, gnss_t=args.gnss)
    lgrdec.monlevel = 1
    nep = 0
    nep_max = 0

    lgrdec.re.anttype = "NASA HIGHGAIN       "
    lgrdec.re.rectype = "LUGRE               "

    path = str(Path(bdir) / fname) if bdir else fname
    blen = os.path.getsize(path)
    with open(path, 'rb') as f:
        msg = f.read(blen)
        maxlen = len(msg)-23
        # maxlen = 7000000+10000
        k = 0
        while k < maxlen:
            stat = lgrdec.sync(msg, k)
            if not stat:
                k += 1
                continue
            len_ = lgrdec.msg_len(msg, k)
            if k+len_ >= maxlen:
                break

            if not lgrdec.checksum(msg, k):
                k += 1
                continue

            lgrdec.decode(msg[k:k+len_], len_)
            k += len_

            nep += 1
            if nep_max > 0 and nep >= nep_max:
                break

    lgrdec.file_close()


def ecef_to_lunar_fixed(x_m, y_m, z_m, target_time=None):
    """
    Converts ECEF (ITRS) coordinates to the Lunar-fixed coordinate system (Selenocentric).
    """
    # 1. Load ephemeris and timescale data
    # DE421 is a standard JPL ephemeris; DE440 is a more recent alternative.
    eph = load('de440.bsp')
    earth = eph['earth']
    moon = eph['moon']
    ts = load.timescale()

    if target_time is None:
        t = ts.now()
    else:
        t = target_time

    # 2. Define ECEF coordinates as Earth-fixed (ITRS)
    # This creates a position object relative to the Earth's center in the ITRS frame.
    ecef_pos = itrs.at(t, x_m=x_m, y_m=y_m, z_m=z_m)

    # 3. Get the point's position in the Inertial Coordinate System (ICRF)
    # Calling .at(t) on an ITRS position transforms it into the inertial ICRF frame.
    pos_icrf = ecef_pos

    # 4. Get the Moon's center position in the ICRF frame
    moon_icrf = moon.at(t)

    # 5. Calculate the relative vector from the Moon's center to the point
    # Relative Vector = Point Position (ICRF) - Moon Position (ICRF)
    relative_vec_icrf = pos_icrf - moon_icrf

    # 6. Transform to the Lunar-fixed frame (considering rotation and libration)
    # Using the "Mean Earth/Polar Axis" (ME) frame via the moon_pa_de421 kernel.
    # Note: .frame_xyz() rotates the inertial vector into the specified body-fixed frame.
    lunar_fixed_pos = relative_vec_icrf.frame_xyz(load('moon_pa_de421.tf'))

    # Retrieve the coordinates in meters
    x_l, y_l, z_l = lunar_fixed_pos.m

    return x_l, y_l, z_l


# --- Example Usage ---
# Example: ECEF coordinates for a point on Earth (e.g., Tokyo area)
x_ecef = -3957224.0
y_ecef = 3310210.0
z_ecef = 3737512.0

lx, ly, lz = ecef_to_lunar_fixed(x_ecef, y_ecef, z_ecef)

print(f"Time (UTC): {load.timescale().now().utc_strftime()}")
print(f"Input ECEF (m): X={x_ecef}, Y={y_ecef}, Z={z_ecef}")
print(f"Output Lunar Fixed (m): X={lx:.2f}, Y={ly:.2f}, Z={lz:.2f}")


def main():

    # Parse command line arguments
    #
    parser = argparse.ArgumentParser(description="LuGRE data converter")

    # Input file and folder
    #
    parser.add_argument(
        "inpFileName",  help="Input BIN file(s) (wildcards allowed)")

    parser.add_argument("--receiver", default='unknown',
                        help="Receiver type [unknown]")
    parser.add_argument("--antenna", default='unknown',
                        help="Antenna type [unknown]")

    parser.add_argument("-g", "--gnss", default='GE',
                        help="GNSS [GE]")

    parser.add_argument("-j", "--jobs", default=int(mp.cpu_count() / 2),
                        type=int, help='Max. number of parallel processes')

    # Retrieve all command line arguments
    #
    args = parser.parse_args()

    opt = rcvOpt()

    opt.flg_rnxobs = True

    # args.inpFileName = '../data/doy2025-074/TLM_RAW_20250315_130709_26H_S_OP76_0.bin'
    # args.inpFileName = '../data/doy2025-074/TLM_NAV_20250315_130709_26H_S_OP76_0.bin'

    # decode(args.inpFileName, opt, args)
    # Start processing pool
    #
    with mp.Pool(processes=args.jobs) as pool:
        pool.starmap(decode, [(f, opt, args) for f in glob(args.inpFileName)])


# Call main function
#
if __name__ == "__main__":
    main()
