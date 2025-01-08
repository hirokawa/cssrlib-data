"""
u-blox Receiver UBX messages decoder

 [1] u-blox F9 HPG 1.51, u-blox F9 high precision GNSS receiver
     Interface description, 2024

@author Rui Hirokawa
"""

from glob import glob
import numpy as np
import os
import struct as st
import bitstruct.c as bs
from cssrlib.gnss import uGNSS, uTYP, prn2sat, Obs, rSigRnx, gpst2time, uSIG
from cssrlib.rawnav import rcvDec, rcvOpt
from binascii import hexlify

CPSTD_VALID = 0.2           # stdev threshold of valid carrier-phase


class ubx(rcvDec):
    """ class for u-blox Binary Format (UBX) decoder """
    tow = -1
    week = -1

    def __init__(self, opt=None, prefix='', gnss_t='GECJ'):
        super().__init__(opt, prefix, gnss_t)

        sig_tbl = {
            uGNSS.GPS: {0: uSIG.L1C, 3: uSIG.L2L, 4: uSIG.L2S, 6: uSIG.L5I,
                        7: uSIG.L5Q},
            uGNSS.SBS: {0: uSIG.L1C},
            uGNSS.GAL: {0: uSIG.L1C, 1: uSIG.L1B, 3: uSIG.L5I, 4: uSIG.L5Q,
                        5: uSIG.L7I, 6: uSIG.L7Q, 8: uSIG.L6B, 9: uSIG.L6C,
                        10: uSIG.L6A},
            uGNSS.BDS: {0: uSIG.L2I, 1: uSIG.L2I, 2: uSIG.L7I, 3: uSIG.L7I,
                        4: uSIG.L6I, 10: uSIG.L6I, 5: uSIG.L1P, 6: uSIG.L1D,
                        7: uSIG.L5P, 8: uSIG.L5D},
            uGNSS.QZS: {0: uSIG.L1C, 1: uSIG.L1Z, 4: uSIG.L2S, 5: uSIG.L2L,
                        8: uSIG.L5I, 9: uSIG.L5Q},
            uGNSS.GLO: {0: uSIG.L1C, 2: uSIG.L2C},
            uGNSS.IRN: {0: uSIG.L5A}
        }

        obs_tbl = {
            uGNSS.GPS: [uSIG.L1C, uSIG.L2L, uSIG.L5Q],
            uGNSS.SBS: [uSIG.L1C],
            uGNSS.GAL: [uSIG.L1C, uSIG.L5Q, uSIG.L7Q, uSIG.L6B],
            uGNSS.BDS: [uSIG.L1P, uSIG.L5P, uSIG.L6I],
            uGNSS.QZS: [uSIG.L1C, uSIG.L2L, uSIG.L5Q],
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

    def sync(self, buff, k):
        return buff[k] == 0xB5 and buff[k+1] == 0x62

    def msg_len(self, msg, k):
        return st.unpack_from('<H', msg, k+4)[0]+8

    def check_crc(self, msg, k):
        len_ = st.unpack_from('<H', msg, k+4)[0]+8
        i = k+len_-2
        cka_, ckb_ = msg[i:i+2]
        cka = ckb = 0
        for j in range(len_-4):
            cka += msg[k+j+2]
            cka &= 0xff
            ckb += cka
            ckb &= 0xff

        if cka_ != cka or ckb_ != ckb:
            if self.monlevel > 0:
                print("checksum error.")
            return False
        return True

    def svid2prn(self, gnss, svid):
        gnss_t = {0: uGNSS.GPS, 1: uGNSS.SBS, 2: uGNSS.GAL, 3: uGNSS.BDS,
                  5: uGNSS.QZS, 6: uGNSS.GLO, 7: uGNSS.IRN}
        sys = gnss_t[gnss]
        if sys == uGNSS.QZS:
            prn = svid + 192
        else:
            prn = svid
        return sys, prn

    def decode_obs(self, buff, k=6):
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

        tow, wn, leapS, nm, stat, ver = st.unpack_from('<dHbBBB', buff, k)
        k += 16

        obs.time = gpst2time(wn, tow)
        self.tow = tow
        self.week = wn

        pr = {}
        cp = {}
        dp = {}
        lli = {}
        cn = {}

        for i in range(nm):
            pr_, cp_, dop_, gnss, svid, sigid, freqid = st.unpack_from(
                '<ddfBBBB', buff, k)
            k += 24

            locktime, cn0, sig_pr, sig_cp, sig_dop, trk = st.unpack_from(
                '<HBBBBB', buff, k)
            k += 8

            sig_cp *= 0.004
            # tracking status:
            # b0: pr-valid
            # b1: cp-valid
            # b2: half-cycle valid
            # b3: half-cycle subtracted from phase
            if (trk & 1) == 0:
                pr_ = 0.0
            if (trk & 2) == 0 or cp_ == -0.5 or (sig_cp > CPSTD_VALID):
                cp_ = 0.0

            sys, prn = self.svid2prn(gnss, svid)

            if sys not in self.sig_tab:
                continue

            if sys == uGNSS.GLO and svid == 255:
                continue

            # ch = freqid-7 if sys == uGNSS.GLO else 0
            sat = prn2sat(sys, prn)
            sig = self.sig_t[sys][sigid][uTYP.L]

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

            slip = 0x01 if locktime == 0 else 0
            halfv = 0x02 if (trk & 4) == 0 else 0
            # halfc = 0x80 if (trk & 8) != 0 else 0

            lli_ = slip + halfv

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

    def decode_nav(self, buff, k=6):
        if self.week == -1:
            return None

        gnss, svid, sigid, freqid, len_, ch, ver = \
            st.unpack_from('<BBBBBBB', buff, k)
        k += 8

        sys, prn = self.svid2prn(gnss, svid)
        sat = prn2sat(sys, prn)

        msg = st.unpack_from('>'+len_*'L', buff, 14)
        b = bytes(np.array(msg, dtype='uint32'))
        fh_ = None
        eph = None
        geph = None
        seph = None
        type_ = sigid
        blen = len_*4
        if sys == uGNSS.GPS:
            if sigid == 0 and self.flg_gpslnav:  # L1C/A
                fh_ = self.fh_gpslnav
                type_ = 0
                eph = self.rn.decode_gps_lnav(self.week, self.tow, sat, b)
            elif sigid == 4 and self.flg_gpscnav:  # L2M
                fh_ = self.fh_gpscnav
                type_ = 1
                # eph = self.rn.decode_gps_cnav(self.week, self.tow, sat, b)
            elif sigid == 6 and self.flg_gpscnav:  # L5I
                fh_ = self.fh_gpscnav
                type_ = 2
                eph = self.rn.decode_gps_cnav(self.week, self.tow, sat, b)
        elif sys == uGNSS.GAL:
            if sigid == 1 or sigid == 3:  # I/NAV
                b_ = bytearray(b)
                b_[15:30] = b_[16:31]  # even 128bit, odd 128bit -> skip 8bit
                b_[30] = 0
                b = bytes(b_)

            if sigid == 1 and self.flg_galinav:  # E1B
                fh_ = self.fh_galinav
                type_ = 0
                eph = self.rn.decode_gal_inav(
                    self.week, self.tow, sat, type_, b)
            elif sigid == 3 and self.flg_galfnav:  # E5aI
                fh_ = self.fh_galfnav
                type_ = 1
                eph = self.rn.decode_gal_fnav(
                    self.week, self.tow, sat, type_, b)
            elif sigid == 5 and self.flg_galinav:  # E5bI
                fh_ = self.fh_galinav
                type_ = 2
                eph = self.rn.decode_gal_inav(
                    self.week, self.tow, sat, type_, b)
            elif sigid == 8 and self.flg_galcnav:  # E6B
                fh_ = self.fh_galcnav
                type_ = 6
        elif sys == uGNSS.BDS:
            if sigid == 0 and self.flg_bdsd12:  # B1I D1
                type_ = 0
                # eph = self.rn.decode_bds_d1(self.week, self.tow, sat, b)
            elif sigid == 1 and self.flg_bdsd12:  # B1I D2
                type_ = 0
            elif sigid == 2 and self.flg_bdsd12:  # B2I D1
                type_ = 1
            elif sigid == 3 and self.flg_bdsd12:  # B2I D2
                type_ = 1
            elif sigid == 4 and self.flg_bdsd12:  # B3I D1
                type_ = 2
                # eph = self.rn.decode_bds_d1(self.week, self.tow, sat, b, 2)
            elif sigid == 6 and self.flg_bdsb1c:  # B1Cd
                fh_ = self.fh_bdsb1c
                type_ = 3

                if blen == 36:
                    self.bds_cnv1[prn][76:109] = b[:33]
                    self.bds_cnv1[prn][123] = 1
                    blen = 0
                elif blen == 76:
                    self.bds_cnv1[prn][:76] = b[:76]
                    self.bds_cnv1[prn][123] |= 2
                    blen = 0
                elif blen == 12:
                    prn_, soh = bs.unpack_from('u6u8', b, 0)
                    self.bds_cnv1[prn][109] = soh
                    b = bytes(self.bds_cnv1[prn])
                    self.bds_cnv1[prn][123] |= 4
                    blen = 110
                    if self.bds_cnv1[prn][123] == 7:
                        eph = self.rn.decode_bds_b1c(
                            self.week, self.tow, sat, b)
            elif sigid == 8 and self.flg_bdsb2a:  # B2ad
                type_ = 5
            elif sigid == 10 and self.flg_bdsd12:  # B3I D2
                type_ = 2
        elif sys == uGNSS.GLO:
            freqid -= 7
            if sigid == 0 and self.flg_gloca:  # L1 OF
                type_ = 0
                geph = self.rn.decode_glo_fdma(self.week, self.tow,
                                               sat, b, freqid)
            elif sigid == 2 and self.flg_gloca:  # L2 OF
                type_ = 2
                geph = self.rn.decode_glo_fdma(self.week, self.tow,
                                               sat, b, freqid)
        elif sys == uGNSS.SBS:
            if sigid == 0 and self.flg_sbas:  # L1C/A
                fh_ = self.fh_sbas
                type_ = 0
                seph = self.rn.decode_sbs_l1(self.week, self.tow, sat, b)
        elif sys == uGNSS.QZS:
            if sigid == 0 and self.flg_qzslnav:  # L1C/A
                fh_ = self.fh_qzslnav
                type_ = 0
                eph = self.rn.decode_gps_lnav(self.week, self.tow, sat, b)
            elif sigid == 4 and self.flg_qzscnav:  # L2C
                fh_ = self.fh_qzscnav
                type_ = 1
                eph = self.rn.decode_gps_cnav(self.week, self.tow, sat, b)
            elif sigid == 8 and self.flg_qzscnav:  # L5I
                fh_ = self.fh_qzscnav
                type_ = 2
                eph = self.rn.decode_gps_cnav(self.week, self.tow, sat, b)

        if self.flg_rnxnav:
            if eph is not None:
                self.re.rnx_nav_body(eph, self.fh_rnxnav)
            if geph is not None:
                self.re.rnx_gnav_body(geph, self.fh_rnxnav)
            if seph is not None:
                self.re.rnx_snav_body(seph, self.fh_rnxnav)

        if fh_ is not None and blen > 0:
            fh_.write("{:4d}\t{:6d}\t{:3d}\t{:1d}\t{:3d}\t{:s}\n".
                      format(self.week, int(self.tow+0.01), prn, type_, blen,
                             hexlify(b[:blen]).decode()))
        if self.monlevel > 0:
            print(f"NAV gnss={gnss}:prn={svid:3d}({freqid:2d}):sig={sigid:2d}")

    def decode_l6msg(self, buff, k=6):
        ver, svid, cno, timetag, gd, nc, ch = st.unpack_from(
            'BBHLBBH', buff, k)
        k += 14
        b = buff[k:k+250]
        print(f"L6 {svid}:{ch}")
        if self.flg_qzsl6:
            self.fh_qzsl6.write("{:4d}\t{:6d}\t{:3d}\t{:1d}\t{:3d}\t{:s}\n".
                                format(self.week, self.tow, svid, ch, len_*4,
                                       hexlify(b).decode()))

    def decode_timegps(self, buff, k=6):
        """ decode NAV-TIMEGPS """
        itow, ftow, self.week, leaps, valid, tacc = \
            st.unpack_from('<LlhbBL', buff, k)
        k += 16
        self.tow = itow*1e-3+ftow*1e-9
        self.time = gpst2time(self.week, self.tow)
        self.leaps = leaps
        return 0

    def decode(self, buff, len_, sys=[], prn=[]):
        k = 2
        cls_, id_ = st.unpack_from('<BB', buff, k)
        k += 4

        # print(f"cls={cls_:2x} id={id_:2x}")

        if cls_ == 0x01:  # UBX-NAV
            if id_ == 0x20:  # UBX-NAV-TIMEGPS
                self.decode_timegps(buff)
        elif cls_ == 0x02:  # UBX-RXM
            if id_ == 0x15:  # UBX-RXM-RAWX
                obs = self.decode_obs(buff)
                self.obs = obs
                if self.flg_rnxobs and obs is not None:
                    self.re.rnx_obs_header(obs.time, self.fh_rnxobs)
                    self.re.rnx_obs_body(obs, self.fh_rnxobs)
            elif id_ == 0x13:  # UBX-RXM-SFRBX
                self.decode_nav(buff)
            elif id_ == 0x73:  # UBX-RXM-QZSSL6
                self.decode_l6msg(buff)
        return 0


if __name__ == "__main__":

    bdir = '../data/doy2025-004/'
    fnames = 'ubf9004c.ubx'

    gnss_t = 'GERCJ'

    opt = rcvOpt()
    opt.flg_qzsl6 = False
    opt.flg_qzslnav = True
    opt.flg_gpslnav = True
    opt.flg_qzscnav = True
    opt.flg_gpscnav = True
    opt.flg_qzsl1s = False
    opt.flg_gale6 = True
    opt.flg_galinav = True
    opt.flg_galfnav = True
    opt.flg_bdsb1c = True
    opt.flg_bdsb2a = False
    opt.flg_bdsb2b = False
    opt.flg_bdsd12 = True
    opt.flg_gloca = True
    opt.flg_irnnav = False
    opt.flg_sbas = True
    opt.flg_rnxnav = True
    opt.flg_rnxobs = True

    for f in glob(bdir+fnames):

        print("Decoding {}".format(f))

        bdir, fname = os.path.split(f)
        bdir += '/'

        prefix = bdir+fname[4:].removesuffix('.ubx')+'_'
        ubxdec = ubx(opt, prefix=prefix, gnss_t=gnss_t)
        ubxdec.monlevel = 1
        nep = 0
        nep_max = 0

        ubxdec.re.anttype = "JAVRINGANT_DM   JVDM"
        ubxdec.re.rectype = "UBLOX F9P           "

        blen = os.path.getsize(bdir+fname)
        with open(bdir+fname, 'rb') as f:
            msg = f.read(blen)
            maxlen = len(msg)-8
            # maxlen = 7000000+10000
            k = 0
            while k < maxlen:
                stat = ubxdec.sync(msg, k)
                if not stat:
                    k += 1
                    continue
                len_ = ubxdec.msg_len(msg, k)
                if k+len_ >= maxlen:
                    break

                if not ubxdec.check_crc(msg, k):
                    k += 1
                    continue

                ubxdec.decode(msg[k:k+len_], len_)
                k += len_

                nep += 1
                if nep_max > 0 and nep >= nep_max:
                    break

        ubxdec.file_close()
