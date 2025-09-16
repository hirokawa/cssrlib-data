#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Novatel Receiver messages decoder

 [1] NovAtel, OEM7 Commands and Logs Reference Manual v1A,
     November 2024

@author Rui Hirokawa
"""

from cssrlib.rawnav import rcvDec, rcvOpt
import argparse
from binascii import hexlify
import bitstruct.c as bs
from glob import glob
import multiprocessing as mp
import numpy as np
import os
from pathlib import Path
import struct as st

from cssrlib.gnss import uGNSS, uTYP, prn2sat, Obs, rSigRnx, gpst2time, uSIG, \
    timediff, gtime_t, copy_buff

CPSTD_VALID = 0.2           # stdev threshold of valid carrier-phase


class nov(rcvDec):
    """ class for Novatel Binary Format decoder """
    tow = -1
    week = -1

    def __init__(self, opt=None, prefix='', gnss_t='GECJ'):
        super().__init__(opt, prefix, gnss_t)

        self.sys_tbl = {
            0: uGNSS.GPS, 1: uGNSS.GLO, 2: uGNSS.SBS, 3: uGNSS.GAL,
            4: uGNSS.BDS, 5: uGNSS.QZS, 6: uGNSS.IRN
        }

        sig_tbl = {
            uGNSS.GPS: {0: uSIG.L1C, 5: uSIG.L2P, 9: uSIG.L2W, 14: uSIG.L5Q,
                        16: uSIG.L1L, 17: uSIG.L2S},
            uGNSS.SBS: {0: uSIG.L1C, 6: uSIG.L5I},
            uGNSS.GAL: {2: uSIG.L1C, 6: uSIG.L6B, 7: uSIG.L6C, 12: uSIG.L5Q,
                        17: uSIG.L7Q, 20: uSIG.L8Q},
            uGNSS.BDS: {0: uSIG.L2I, 1: uSIG.L7I, 2: uSIG.L6I,
                        4: uSIG.L2I, 5: uSIG.L7I, 6: uSIG.L6I, 7: uSIG.L1P,
                        9: uSIG.L5P, 11: uSIG.L7P},
            uGNSS.QZS: {0: uSIG.L1C, 14: uSIG.L5Q, 16: uSIG.L1L, 17: uSIG.L2S,
                        27: uSIG.L6L},
            uGNSS.GLO: {0: uSIG.L1C, 1: uSIG.L2C, 5: uSIG.L2P, 6: uSIG.L3Q},
            uGNSS.IRN: {0: uSIG.L5A}
        }

        obs_tbl = {
            uGNSS.GPS: [uSIG.L1C, uSIG.L2W, uSIG.L2S, uSIG.L5Q],
            uGNSS.SBS: [uSIG.L1C, uSIG.L5I],
            uGNSS.GAL: [uSIG.L1C, uSIG.L5Q, uSIG.L7Q, uSIG.L6B],
            uGNSS.BDS: [uSIG.L2I, uSIG.L1P, uSIG.L5P],
            uGNSS.QZS: [uSIG.L1C, uSIG.L2S, uSIG.L5Q],
            uGNSS.GLO: [uSIG.L1C, uSIG.L2C],
            uGNSS.IRN: [uSIG.L5A],
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

        self.bds_cnv1 = {}
        for k in range(uGNSS.BDSMAX):
            self.bds_cnv1[k] = bytearray(124)

        self.head_len = 28
        self.time_p = gtime_t()

    def sync(self, buff, k):
        return buff[k] == 0xAA and buff[k+1] == 0x44 and buff[k+2] == 0x12

    def msg_len(self, msg, k):
        return st.unpack_from('<H', msg, k+8)[0]+28

    def crc32(self, data, len_):
        poly_crc32 = 0xedb88320
        crc = 0

        for i in range(len_):
            crc ^= data[i]
            for _ in range(8):
                if crc & 1:
                    crc = (crc >> 1) ^ poly_crc32
                else:
                    crc >>= 1
        return crc

    def check_crc(self, msg, k):
        len_ = st.unpack_from('<H', msg, k+8)[0]+28
        crc_ = st.unpack_from('<I', msg, k+len_)[0]
        crc = self.crc32(msg[k:k+len_], len_)
        if crc_ != crc:
            if self.monlevel > 0:
                print("checksum error.")
            return False
        return True

    def decode_obs(self, buff):
        """ decode RANGEB """
        k = self.head_len

        obs = Obs()
        obs.time = self.time
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

        pr = {}
        cp = {}
        dp = {}
        lli = {}
        cn = {}

        nobs = st.unpack_from('<I', buff, k)[0]
        k += 4
        for _ in range(nobs):

            status = st.unpack_from('<I', buff, k+40)[0]

            sys_ = (status >> 16) & 0x7
            code = (status >> 21) & 0x1f
            ch = (status >> 5) & 0x1f
            clock = (status >> 12) & 0x1
            # parity = (status >> 11) & 0x1
            plock = (status >> 10) & 0x1
            halfc = 0x40 if (status >> 28) & 0x1 else 0

            sys = self.sys_tbl[sys_]
            prn = st.unpack_from('<H', buff, k)[0]

            if sys_ == 3:
                pass

            # if prn == 0:
            #    continue
            #    if prn == 0 or not parity:
            #        continue
            #    prn -= 37
            sat = prn2sat(sys, prn)

            gfrq, pr_, sig_pr, cp_, sig_cp, dop_, cn0, lockt = \
                st.unpack_from('<Hdfdffff', buff, k+2)
            k += 44

            if self.monlevel > 1:
                print(f"{sys_} {prn:3d} {code:2d} {ch:2d} {gfrq:2d}")

            if not clock:
                pr_ = 0.0
            if not plock:
                cp_ = dop_ = 0.0

            if sys not in self.sig_t.keys():
                continue

            sig = self.sig_t[sys][code][uTYP.L]

            if sig not in self.sig_tab[sys][sig.typ]:
                if self.monlevel > 1:
                    print("skip code={:}".format(sig))
                continue
            idx = self.sig_tab[sys][sig.typ].index(sig)

            if sat not in pr.keys():
                pr[sat] = {}
                cp[sat] = {}
                dp[sat] = {}
                lli[sat] = {}
                cn[sat] = {}

            slip = 0x01 if lockt == 0 else 0
            lli_ = slip + halfc*0

            pr[sat][idx] = pr_
            cp[sat][idx] = cp_
            dp[sat][idx] = dop_
            lli[sat][idx] = lli_
            cn[sat][idx] = cn0

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

    def decode_gps_lnav(self, buff, sys, prn, k):
        sat = prn2sat(sys, prn)
        msg = bytearray(40)
        i = k*8
        j = 2
        for _ in range(10):
            d = bs.unpack_from('u24', buff, i)[0]
            i += 24
            bs.pack_into('u24', msg, j, d)
            j += 32

        eph = self.rn.decode_gps_lnav(self.week, self.tow, sat, msg)
        if eph is not None:
            self.re.rnx_nav_body(eph, self.fh_rnxnav)

    def decode_gps_cnav(self, buff, sys, prn, src, k):
        sat = prn2sat(sys, prn)

        type_ = 0 if src in (69, 14787) else 1
        blen = 38

        fh_ = self.fh_qzscnav if sys == uGNSS.QZS else self.fh_gpscnav

        fh_.write("{:4d}\t{:6d}\t{:3d}\t{:1d}\t{:3d}\t".
                  format(self.week, int(self.tow), prn, type_, blen))
        msg = bytearray(40)
        for i in range(38):
            d = st.unpack_from('<B', buff, k)[0]
            fh_.write("{:02x}".format(d))
            msg[i] = d
            k += 1
        fh_.write("\n")

        eph = self.rn.decode_gps_cnav(
            self.week, self.tow, sat, msg)
        if eph is not None:
            self.re.rnx_nav_body(eph, self.fh_rnxnav)

    def decode(self, buff, len_, sys=[], prn=[]):

        len_h, id_ = st.unpack_from('<BH', buff, 3)
        sts, week, tow = st.unpack_from('<BHI', buff, 13)
        if sts == 20 or week == 0:
            return -1
        tow *= 1e-3
        self.head_len = len_h
        self.week = week
        self.tow = tow
        self.time = gpst2time(self.week, self.tow)
        if self.monlevel > 0:
            print(f"week={week} tow={tow:6.1f} id={id_:2d}")

        if id_ == 7:  # GPS L1 C/A ephemeris
            pass
        elif id_ == 25:  # RAWGPSSUBFRAME
            if self.flg_gpslnav:
                k = self.head_len
                _, prn, sid = st.unpack_from('<III', buff, k)
                k += 12
                self.decode_gps_lnav(buff, uGNSS.GPS, prn, k)

        elif id_ == 43:  # RANGEB
            obs = self.decode_obs(buff)
            if obs is not None:
                self.re.rnx_obs_header(obs.time, self.fh_rnxobs)
                self.re.rnx_obs_body(obs, self.fh_rnxobs)
        elif id_ == 93:  # RXSTATUS
            pass

        elif id_ == 722:  # GLORAWSTRING
            if self.flg_gloca:
                k = self.head_len
                slot, freq = st.unpack_from('<BB', buff, k)
                k += 2
                msg = bytearray(12)
                for i in range(11):
                    d = st.unpack_from('<B', buff, k)[0]
                    msg[i] = buff[k]
                    k += 1

                sat = prn2sat(uGNSS.GLO, slot)
                geph = self.rn.decode_glo_fdma(self.week, self.tow,
                                               sat, msg, freq)
                if geph is not None:
                    self.re.rnx_gnav_body(geph, self.fh_rnxnav)

        elif id_ == 723:  # GLONASS L1 C/A ephemeris
            pass
        elif id_ == 963:  # HWMONITOR
            pass
        elif id_ == 1122:  # GAL ephemeris
            pass
        elif id_ == 1330:  # QZSSRAWSUBFRAME
            if self.flg_qzslnav:
                k = self.head_len
                prn, sid = st.unpack_from('<II', buff, k)
                k += 8
                self.decode_gps_lnav(buff, uGNSS.QZS, prn, k)

        elif id_ == 1336:  # GLONASS Eph
            pass
        elif id_ == 1413:  # GALFNAVRAWPAGE
            k = self.head_len
            ch, prn = st.unpack_from('<II', buff, k)
            k += 8
            sat = prn2sat(uGNSS.GAL, prn)

            if self.flg_galfnav:
                type_ = 0
                self.fh_galfnav.write("{:4d}\t{:6d}\t{:3d}\t{:1d}\t{:3d}\t".
                                      format(self.week, int(self.tow), prn,
                                             type_, 27))

                msg = bytearray(31)
                for i in range(27):
                    d = st.unpack_from('<B', buff, k)[0]
                    self.fh_galfnav.write("{:02x}".format(d))
                    msg[i] = buff[k]
                    k += 1
                self.fh_galfnav.write("\n")

                eph = self.rn.decode_gal_fnav(self.week, self.tow, sat, 1, msg)
                if eph is not None:
                    self.re.rnx_nav_body(eph, self.fh_rnxnav)

        elif id_ == 1414:  # GALINAVRAWWORD
            k = self.head_len
            ch, prn, src = st.unpack_from('<III', buff, k)
            k += 12
            sat = prn2sat(uGNSS.GAL, prn)

            if self.flg_galinav:
                if src == 10433:  # E1B
                    type_ = 0
                elif src == 10499:  # E5b
                    type_ = 2
                else:
                    return -1

                self.fh_galinav.write("{:4d}\t{:6d}\t{:3d}\t{:1d}\t{:3d}\t".
                                      format(self.week, int(self.tow), prn,
                                             type_, 30))

                msg = bytearray(20)
                copy_buff(buff, msg, k*8, 2, 112)
                copy_buff(buff, msg, (k+14)*8, 122, 16)

                for i in range(16):
                    d = st.unpack_from('<B', buff, k)[0]
                    self.fh_galinav.write("{:02x}".format(d))
                    k += 1
                self.fh_galinav.write("\n")

                eph = self.rn.decode_gal_inav(self.week, self.tow,
                                              sat, 2, msg)
                if self.mode_galinav == 0 and eph is not None:
                    self.re.rnx_nav_body(eph, self.fh_rnxnav)
        elif id_ == 1695:  # BDSRAWNAVSUBFRAME
            if self.flg_bdsd12:
                k = self.head_len
                ch, prn, src, sid = st.unpack_from('<IIII', buff, k)
                k += 16
                sat = prn2sat(uGNSS.BDS, prn)
                msg = bytearray(40)
                ofst = k*8
                for i in range(10):
                    sz = 26 if i == 0 else 22
                    fmt = 'u'+str(sz)
                    d = bs.unpack_from(fmt, buff, ofst)[0]
                    ofst += sz
                    bs.pack_into(fmt, msg, 30*i, d)

                eph = None

                if prn > 5 and prn < 59:
                    eph = self.rn.decode_bds_d1(
                        self.week, self.tow, sat, msg)
                else:
                    eph = self.rn.decode_bds_d2(
                        self.week, self.tow, sat, msg)

                if eph is not None:
                    self.re.rnx_nav_body(eph, self.fh_rnxnav)

        elif id_ == 1696:  # BDS Eph
            pass
        elif id_ == 2105:  # NAVICRAWSUBFRAME
            if self.flg_irnnav:
                k = self.head_len
                ch, prn, sid = st.unpack_from('<III', buff, k)
                k += 12
                sat = prn2sat(uGNSS.IRN, prn)
                msg = bytearray(33)
                for i in range(33):
                    d = st.unpack_from('<B', buff, k)[0]
                    msg[i] = buff[k]
                    k += 1

                eph = self.rn.decode_irn_lnav(self.week, self.tow, sat, msg)
                if eph is not None:
                    self.re.rnx_nav_body(eph, self.fh_rnxnav)

        elif id_ == 2185:  # RAWSBASFRAME2
            if self.flg_sbas:
                if timediff(self.time, self.time_p) > 0:
                    self.sbas_frm = {}

                self.time_p = self.time
                k = self.head_len
                prn, ch, src, pre, _, frm = st.unpack_from('<IIBBHI', buff, k)
                k += 16

                itype = 0 if src == 1 else 1  # 0:L1, 1:L5

                if prn not in self.sbas_frm.keys():
                    self.sbas_frm[prn] = []
                if itype not in self.sbas_frm[prn]:
                    self.sbas_frm[prn].append(itype)

                    msg = bytearray(32)
                    for i in range(29):
                        d = st.unpack_from('<B', buff, k)[0]
                        msg[i] = buff[k]
                        k += 1

                    self.output_sbas(prn, msg, self.fh_sbas, itype)

                    sat = prn2sat(uGNSS.SBS, prn)
                    if self.flg_rnxnav:
                        seph = None
                        if itype == 0:
                            seph = self.rn.decode_sbs_l1(
                                self.week, self.tow, sat, msg)
                        if seph is not None:
                            self.re.rnx_snav_body(seph, self.fh_rnxnav)

        elif id_ == 2261:  # QZSSCNAVRAWMESSAGE
            if self.flg_qzscnav:
                k = self.head_len
                ch, prn, src, id_ = st.unpack_from('<IIII', buff, k)
                k += 16

                self.decode_gps_cnav(buff, uGNSS.QZS, prn, src, k)

        elif id_ == 2239:  # GALCNAVRAWPAGE
            k = self.head_len
            ch, prn, id_, page = st.unpack_from('<IIHH', buff, k)
            k += 12

            if self.flg_gale6:
                blen = 58
                type_ = 6
                self.fh_gale6.write("{:4d}\t{:6d}\t{:3d}\t{:1d}\t{:3d}\t".
                                    format(self.week, int(self.tow), prn,
                                           type_, blen))
                for i in range(58):
                    d = st.unpack_from('<B', buff, k)[0]
                    self.fh_gale6.write("{:02x}".format(d))
                    k += 1
                self.fh_gale6.write("\n")

        elif id_ == 2262:  # GPSCNAVRAWMESSAGE
            if self.flg_gpscnav:
                k = self.head_len
                ch, prn, src, id_ = st.unpack_from('<IIII', buff, k)
                k += 16

                self.decode_gps_cnav(buff, uGNSS.GPS, prn, src, k)

        elif id_ == 2373:  # BDSBCNAV1RAWMESSAGE
            if self.flg_bdsb1c:
                k = self.head_len
                ch, prn, page = st.unpack_from('<III', buff, k)
                k += 12
                blen = 110
                msg = bytearray(228)

                copy_buff(buff, msg, k*8, 608+264, 8)  # SF1
                k += 2
                copy_buff(buff, msg, k*8, 0, 576)  # SF1
                k += 72
                copy_buff(buff, msg, k*8, 608, 240)  # SF1
                k += 30

                sat = prn2sat(uGNSS.BDS, prn)
                eph = self.rn.decode_bds_b1c(self.week, self.tow, sat, msg)
                if eph is not None:
                    self.re.rnx_nav_body(eph, self.fh_rnxnav)

        elif id_ == 2374:  # BDSBCNAV2RAWMESSAGE
            if self.flg_bdsb2a:
                k = self.head_len
                ch, prn, page = st.unpack_from('<III', buff, k)
                k += 12
                blen = 36
                msg = bytearray(blen)

                pre = bs.unpack_from('u24', buff, k*8)[0]
                if pre != 0xe24de8:  # preamble
                    return -1
                k += 3
                copy_buff(buff, msg, k*8, 0, 264)

                sat = prn2sat(uGNSS.BDS, prn)
                eph = self.rn.decode_bds_b2a(self.week, self.tow, sat, msg)
                if eph is not None:
                    self.re.rnx_nav_body(eph, self.fh_rnxnav)

        elif id_ == 2411:  # BDSBCNAV3RAWMESSAGE
            if self.flg_bdsb2b:
                k = self.head_len
                ch, prn, page = st.unpack_from('<III', buff, k)
                k += 12
                blen = 62
                msg = bytearray(blen)

                pre = bs.unpack_from('u16', buff, k*8)[0]
                if pre != 0xeb90:  # preamble
                    return -1
                k += 2
                copy_buff(buff, msg, k*8, 0, 462)

                sat = prn2sat(uGNSS.BDS, prn)
                eph = self.rn.decode_bds_b2b(self.week, self.tow, sat, msg)
                if eph is not None:
                    self.re.rnx_nav_body(eph, self.fh_rnxnav)
        else:
            if self.monlevel > 0:
                print(f"missing {id_}")

        return 0


def decode(f, opt, args):

    print("Decoding {}".format(f))

    bdir, fname = os.path.split(f)

    prefix = fname[4:].removesuffix('.nvr')+'_'
    prefix = str(Path(bdir) / prefix) if bdir else prefix
    novdec = nov(opt, prefix=prefix, gnss_t=args.gnss)
    novdec.monlevel = 1
    nep = 0
    nep_max = 0

    if fname.startswith('nov'):
        novdec.re.anttype = "JAVRINGANT_DM   JVDM"
        novdec.re.rectype = "NOVATEL OEM 719     "
    else:
        novdec.re.anttype = args.antenna
        novdec.re.rectype = args.receiver

    path = str(Path(bdir) / fname) if bdir else fname
    blen = os.path.getsize(path)
    with open(path, 'rb') as f:
        msg = f.read(blen)
        maxlen = len(msg)-8
        # maxlen = 7000000+10000
        k = 0
        while k < maxlen:
            stat = novdec.sync(msg, k)
            if not stat:
                k += 1
                continue
            len_ = novdec.msg_len(msg, k)
            if k+len_ >= maxlen:
                break

            if not novdec.check_crc(msg, k):
                k += 1
                continue

            novdec.decode(msg[k:k+len_], len_)
            k += len_+4

            nep += 1
            if nep_max > 0 and nep >= nep_max:
                break

    novdec.file_close()


def main():

    # Parse command line arguments
    #
    parser = argparse.ArgumentParser(description="Novatel OEM converter")

    # Input file and folder
    #
    parser.add_argument(
        "inpFileName", help="Input OEM file(s) (wildcards allowed)")

    parser.add_argument("--receiver", default='unknown',
                        help="Receiver type [unknown]")
    parser.add_argument("--antenna", default='unknown',
                        help="Antenna type [unknown]")

    parser.add_argument("-g", "--gnss", default='GRECJ',
                        help="GNSS [GRECJ]")

    parser.add_argument("-j", "--jobs", default=int(mp.cpu_count() / 2),
                        type=int, help='Max. number of parallel processes')

    # Retrieve all command line arguments
    #
    args = parser.parse_args()

    opt = rcvOpt()

    opt.flg_rnxobs = True
    opt.flg_rnxnav = True

    opt.flg_gpslnav = True
    opt.flg_gpscnav = True

    opt.flg_qzsl6 = False
    opt.flg_qzslnav = True
    opt.flg_qzscnav = True
    opt.flg_qzsl1s = False

    opt.flg_gale6 = True
    opt.flg_galinav = True
    opt.flg_galfnav = True

    opt.flg_bdsb1c = True
    opt.flg_bdsb2a = True
    opt.flg_bdsb2b = True
    opt.flg_bdsd12 = True

    opt.flg_gloca = True

    opt.flg_irnnav = False
    opt.flg_sbas = True

    # Start processing pool
    #
    # with mp.Pool(processes=args.jobs) as pool:
    #    pool.starmap(decode, [(f, opt, args) for f in glob(args.inpFileName)])
    decode(args.inpFileName, opt, args)


# Call main function
#
if __name__ == "__main__":
    main()
