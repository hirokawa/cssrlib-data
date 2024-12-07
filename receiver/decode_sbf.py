"""
Septentrio Receiver SBF messages decoder

 [1] mosaic-X5 Reference Guide, Applicable to version 4.14.0
     of the Firmware, 2023

 [2] PolaRX5 Reference Guide, Applicable to version 5.5.0
     of the Firmware, 2023

@author Rui Hirokawa
"""

from glob import glob
import numpy as np
import os
import struct as st

from cssrlib.gnss import uGNSS, uTYP, prn2sat, Eph, Obs, rSigRnx, gpst2time
from cssrlib.gnss import rCST, gst2time, pos2ecef, uSIG
from cssrlib.rawnav import RawNav, rcvDec, rcvOpt
from crccheck.crc import Crc16Xmodem
import bitstruct.c as bs


class sbf(rcvDec):
    """ class for Septentrio Binary Format (SBF) decoder """
    tow = -1
    week = -1

    def __init__(self, opt=None, prefix='', gnss_t='GECJ'):
        super().__init__(opt, prefix, gnss_t)

        # from 4.1.10 Signal Type [2]
        sig_tbl = {
            0: [uGNSS.GPS, uSIG.L1C],
            1: [uGNSS.GPS, uSIG.L1W],
            2: [uGNSS.GPS, uSIG.L2W],
            3: [uGNSS.GPS, uSIG.L2L],
            4: [uGNSS.GPS, uSIG.L5Q],
            5: [uGNSS.GPS, uSIG.L1L],
            6: [uGNSS.QZS, uSIG.L1C],
            7: [uGNSS.QZS, uSIG.L2L],
            8: [uGNSS.GLO, uSIG.L1C],
            9: [uGNSS.GLO, uSIG.L1P],
            10: [uGNSS.GLO, uSIG.L2P],
            11: [uGNSS.GLO, uSIG.L2C],
            12: [uGNSS.GLO, uSIG.L3Q],
            13: [uGNSS.BDS, uSIG.L1P],
            14: [uGNSS.BDS, uSIG.L5P],
            15: [uGNSS.IRN, uSIG.L5A],
            17: [uGNSS.GAL, uSIG.L1C],
            19: [uGNSS.GAL, uSIG.L6C],
            20: [uGNSS.GAL, uSIG.L5Q],
            21: [uGNSS.GAL, uSIG.L7Q],
            22: [uGNSS.GAL, uSIG.L8Q],
            24: [uGNSS.SBS, uSIG.L1C],
            25: [uGNSS.SBS, uSIG.L5I],
            26: [uGNSS.QZS, uSIG.L5Q],
            27: [uGNSS.QZS, uSIG.L6Z],
            28: [uGNSS.BDS, uSIG.L2I],
            29: [uGNSS.BDS, uSIG.L7I],
            30: [uGNSS.BDS, uSIG.L6I],
            32: [uGNSS.QZS, uSIG.L1L],
            33: [uGNSS.QZS, uSIG.L1Z],
            34: [uGNSS.BDS, uSIG.L7D],
            38: [uGNSS.QZS, uSIG.L1E],
            39: [uGNSS.QZS, uSIG.L5P],
        }

        self.sig_t = {}
        for key in sig_tbl.keys():
            s = sig_tbl[key]
            self.sig_t[key] = {
                uTYP.C: rSigRnx(s[0], uTYP.C, s[1]),
                uTYP.L: rSigRnx(s[0], uTYP.L, s[1]),
                uTYP.D: rSigRnx(s[0], uTYP.D, s[1]),
                uTYP.S: rSigRnx(s[0], uTYP.S, s[1]),
            }

        self.sig_t = {
            0: {uTYP.C: rSigRnx('GC1C'), uTYP.L: rSigRnx('GL1C'),
                uTYP.D: rSigRnx('GD1C'), uTYP.S: rSigRnx('GS1C')},
            1: {uTYP.C: rSigRnx('GC1W'), uTYP.L: rSigRnx('GL1W'),
                uTYP.D: rSigRnx('GD1W'), uTYP.S: rSigRnx('GS1W')},
            2: {uTYP.C: rSigRnx('GC2W'), uTYP.L: rSigRnx('GL2W'),
                uTYP.D: rSigRnx('GD2W'), uTYP.S: rSigRnx('GS2W')},
            3: {uTYP.C: rSigRnx('GC2L'), uTYP.L: rSigRnx('GL2L'),
                uTYP.D: rSigRnx('GD2L'), uTYP.S: rSigRnx('GS2L')},
            4: {uTYP.C: rSigRnx('GC5Q'), uTYP.L: rSigRnx('GL5Q'),
                uTYP.D: rSigRnx('GD5Q'), uTYP.S: rSigRnx('GS5Q')},
            5: {uTYP.C: rSigRnx('GC1L'), uTYP.L: rSigRnx('GL1L'),
                uTYP.D: rSigRnx('GD1L'), uTYP.S: rSigRnx('GS1L')},
            6: {uTYP.C: rSigRnx('JC1C'), uTYP.L: rSigRnx('JL1C'),
                uTYP.D: rSigRnx('JD1C'), uTYP.S: rSigRnx('JS1C')},
            7: {uTYP.C: rSigRnx('JC2L'), uTYP.L: rSigRnx('JL2L'),
                uTYP.D: rSigRnx('JD2L'), uTYP.S: rSigRnx('JS2L')},
            8: {uTYP.C: rSigRnx('RC1C'), uTYP.L: rSigRnx('RL1C'),
                uTYP.D: rSigRnx('RD1C'), uTYP.S: rSigRnx('RS1C')},
            9: {uTYP.C: rSigRnx('RC1P'), uTYP.L: rSigRnx('RL1P'),
                uTYP.D: rSigRnx('RD1P'), uTYP.S: rSigRnx('RS1P')},
            10: {uTYP.C: rSigRnx('RC2P'), uTYP.L: rSigRnx('RL2P'),
                 uTYP.D: rSigRnx('RD2P'), uTYP.S: rSigRnx('RS2P')},
            11: {uTYP.C: rSigRnx('RC2C'), uTYP.L: rSigRnx('RL2C'),
                 uTYP.D: rSigRnx('RD2C'), uTYP.S: rSigRnx('RS2C')},
            12: {uTYP.C: rSigRnx('RC3Q'), uTYP.L: rSigRnx('RL3Q'),
                 uTYP.D: rSigRnx('RD3Q'), uTYP.S: rSigRnx('RS3Q')},
            13: {uTYP.C: rSigRnx('CC1P'), uTYP.L: rSigRnx('CL1P'),
                 uTYP.D: rSigRnx('CD1P'), uTYP.S: rSigRnx('CS1P')},
            14: {uTYP.C: rSigRnx('CC5P'), uTYP.L: rSigRnx('CL5P'),
                 uTYP.D: rSigRnx('CD5P'), uTYP.S: rSigRnx('CS5P')},
            15: {uTYP.C: rSigRnx('IC5A'), uTYP.L: rSigRnx('IL5A'),
                 uTYP.D: rSigRnx('ID5A'), uTYP.S: rSigRnx('IS5A')},
            17: {uTYP.C: rSigRnx('EC1C'), uTYP.L: rSigRnx('EL1C'),
                 uTYP.D: rSigRnx('ED1C'), uTYP.S: rSigRnx('ES1C')},
            19: {uTYP.C: rSigRnx('EC6C'), uTYP.L: rSigRnx('EL6C'),
                 uTYP.D: rSigRnx('ED6C'), uTYP.S: rSigRnx('ES6C')},
            20: {uTYP.C: rSigRnx('EC5Q'), uTYP.L: rSigRnx('EL5Q'),
                 uTYP.D: rSigRnx('ED5Q'), uTYP.S: rSigRnx('ES5Q')},
            21: {uTYP.C: rSigRnx('EC7Q'), uTYP.L: rSigRnx('EL7Q'),
                 uTYP.D: rSigRnx('ED7Q'), uTYP.S: rSigRnx('ES7Q')},
            22: {uTYP.C: rSigRnx('EC8Q'), uTYP.L: rSigRnx('EL8Q'),
                 uTYP.D: rSigRnx('ED8Q'), uTYP.S: rSigRnx('ES8Q')},
            24: {uTYP.C: rSigRnx('SC1C'), uTYP.L: rSigRnx('SL1C'),
                 uTYP.D: rSigRnx('SD1C'), uTYP.S: rSigRnx('SS1C')},
            25: {uTYP.C: rSigRnx('SC5I'), uTYP.L: rSigRnx('SL5I'),
                 uTYP.D: rSigRnx('SD5I'), uTYP.S: rSigRnx('SS5I')},
            26: {uTYP.C: rSigRnx('JC5Q'), uTYP.L: rSigRnx('JL5Q'),
                 uTYP.D: rSigRnx('JD5Q'), uTYP.S: rSigRnx('JS5Q')},
            28: {uTYP.C: rSigRnx('CC2I'), uTYP.L: rSigRnx('CL2I'),
                 uTYP.D: rSigRnx('CD2I'), uTYP.S: rSigRnx('CS2I')},
            29: {uTYP.C: rSigRnx('CC7I'), uTYP.L: rSigRnx('CL7I'),
                 uTYP.D: rSigRnx('CD7I'), uTYP.S: rSigRnx('CS7I')},
            30: {uTYP.C: rSigRnx('CC6I'), uTYP.L: rSigRnx('CL6I'),
                 uTYP.D: rSigRnx('CD6I'), uTYP.S: rSigRnx('CS6I')},
            32: {uTYP.C: rSigRnx('JC1L'), uTYP.L: rSigRnx('JL1L'),
                 uTYP.D: rSigRnx('JD1L'), uTYP.S: rSigRnx('JS1L')},
            33: {uTYP.C: rSigRnx('JC1Z'), uTYP.L: rSigRnx('JL1Z'),
                 uTYP.D: rSigRnx('JD1Z'), uTYP.S: rSigRnx('JS1Z')},
            34: {uTYP.C: rSigRnx('CC7D'), uTYP.L: rSigRnx('CL7D'),
                 uTYP.D: rSigRnx('CD7D'), uTYP.S: rSigRnx('CS7D')},
            38: {uTYP.C: rSigRnx('JC1E'), uTYP.L: rSigRnx('JL1E'),
                 uTYP.D: rSigRnx('JD1E'), uTYP.S: rSigRnx('JS1E')},
            39: {uTYP.C: rSigRnx('JC5P'), uTYP.L: rSigRnx('JL5P'),
                 uTYP.D: rSigRnx('JD5P'), uTYP.S: rSigRnx('JS5P')},
        }

    def sync(self, buff, k):
        return buff[k] == 0x24 and buff[k+1] == 0x40

    def msg_len(self, msg, k):
        return st.unpack_from('<H', msg, k+6)[0]

    def check_crc(self, msg, k):
        crc_, id_, len_ = st.unpack_from('<HHH', msg, k+2)
        crc = Crc16Xmodem.calc(msg[k+4:k+len_])
        if crc_ != crc:
            if self.monlevel > 0:
                print("checksum error.")
            return False
        return True

    def svid2prn(self, svid):
        if svid == 0:
            sys = uGNSS.NONE
            return 0
        elif svid <= 37:
            sys = uGNSS.GPS
            prn = svid
        elif svid <= 61:
            sys = uGNSS.GLO
            prn = svid-37
        elif svid <= 62:
            sys = uGNSS.GLO
            prn = 0
        elif svid <= 68:
            sys = uGNSS.GLO
            prn = svid-38
        elif svid >= 71 and svid <= 106:
            sys = uGNSS.GAL
            prn = svid-70
        elif svid >= 120 and svid <= 140:
            sys = uGNSS.SBS
            prn = svid
        elif svid >= 141 and svid <= 180:
            sys = uGNSS.BDS
            prn = svid-140
        elif svid >= 181 and svid <= 190:
            sys = uGNSS.QZS
            prn = svid-180+192
        elif svid >= 191 and svid <= 197:
            sys = uGNSS.IRN
            prn = svid-190
        elif svid >= 198 and svid <= 215:
            sys = uGNSS.SBS
            prn = svid-57
        elif svid >= 216 and svid <= 222:
            sys = uGNSS.IRN
            prn = svid-208
        elif svid >= 223 and svid <= 245:
            sys = uGNSS.BDS
            prn = svid-182

        return sys, prn

    def decode_head(self, buff, k):
        tow, wn, svid = st.unpack_from('<LHB', buff, k)
        sys, prn = self.svid2prn(svid)
        self.tow = tow*1e-3
        self.week = wn
        return sys, prn

    def decode_obs(self, buff, k=8):
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

        obs.P = np.empty((0, self.nsig[uTYP.C]), dtype=np.float64)
        obs.L = np.empty((0, self.nsig[uTYP.L]), dtype=np.float64)
        obs.D = np.empty((0, self.nsig[uTYP.D]), dtype=np.float64)
        obs.S = np.empty((0, self.nsig[uTYP.S]), dtype=np.float64)
        obs.lli = np.empty((0, self.nsig[uTYP.L]), dtype=np.int32)
        obs.sat = np.empty(0, dtype=np.int32)

        tow, wn, nb1, sb1len, sb2len = st.unpack_from('<LHBBB', buff, k)
        obs.time = gpst2time(wn, tow*1e-3)

        k += 9
        self.tow = tow*1e-3
        self.week = wn
        flags, clkjump = st.unpack_from('<BB', buff, k)
        k += 3
        for i in range(nb1):
            ch, typ, svid, misc, P0, dop, cp0, cp1, cn0, lock, info, nb2 = \
                st.unpack_from('<BBBBLlHbBHBB', buff, k)
            k += sb1len

            sys, prn = self.svid2prn(svid)
            sat = prn2sat(sys, prn)

            if sys not in self.sig_tab:
                k += nb2*sb2len
                continue

            pr = np.zeros(nsig_max, dtype=np.float64)
            cp = np.zeros(nsig_max, dtype=np.float64)
            dp = np.zeros(nsig_max, dtype=np.float64)
            ll = np.zeros(nsig_max, dtype=np.int32)
            cn = np.zeros(nsig_max, dtype=np.float64)

            sig = typ & 0x1f
            # ant_id = (typ >> 5) & 0x7
            ch = 0
            if sig == 31:
                sig = (info >> 3)+32
            elif sig in (8, 9, 10, 11):  # GLONASS
                ch = (info >> 3) - 8

            code = self.sig_t[sig][uTYP.L]
            freq1 = code.frequency(ch)
            lam1 = rCST.CLIGHT/freq1

            P1, L1, D1, S1 = 0.0, 0.0, 0.0, 0.0
            lli = 0

            if misc & 0x1f != 0 or P0 != 0:
                P1 = (misc & 0xf)*4294967.296 + P0*1e-3
            if dop != -(1 << 31):
                D1 = dop*1e-4

            if P1 != 0.0 and freq1 > 0.0 and lock != 65535 and \
                    (cp1 != -128 or cp0 != 0):
                L1 = P1/lam1+(cp1*65.536+cp0*1e-3)
                if info & 4 != 0:
                    lli = 2

            if cn0 != 255:
                if sig in (1, 2):
                    S1 = cn0*0.25
                else:
                    S1 = cn0*0.25 + 10.0

            if code not in self.sig_tab[sys][code.typ]:
                if self.monlevel > 1:
                    print("skip code={:}".format(code))
                k += nb2*sb2len
                continue
            idx = self.sig_tab[sys][code.typ].index(code)

            pr[idx] = P1
            cp[idx] = L1
            dp[idx] = D1
            ll[idx] = lli
            cn[idx] = S1

            if self.monlevel >= 2:
                print("{:6d} {:3d} {:} {:14.3f} {:14.3f} {:10.3f} {:4.2f}".
                      format(int(self.tow), prn, code, P1, L1, D1, S1))

            for j in range(nb2):
                typ, ltime, cn0, ofst1, cp1, info, cofst0, cp0, dop0 = \
                    st.unpack_from('<BBBBbBHHH', buff, k)
                k += sb2len
                sig = typ & 0x1f

                dofst1, cofst1 = bs.unpack_from(
                    's5s3', ofst1.to_bytes(1, byteorder='big'), 0)

                if sig == 31:
                    sig = (info >> 3)+32

                if sig not in self.sig_t:
                    if self.monlevel > 1:
                        print("skip sig={:}".format(sig))
                    continue

                code = self.sig_t[sig][uTYP.L]
                freq2 = code.frequency(ch)
                lam2 = rCST.CLIGHT/freq2

                P2, L2, D2, S2 = 0.0, 0.0, 0.0, 0.0
                lli = 0

                if P1 != 0.0 and (cofst1 != -4 or cofst0 != 0):
                    P2 = P1 + cofst1*65.536+cofst0*1e-3

                if P2 != 0.0 and freq2 > 0.0 and (cp1 != -128 or cp0 != 0):
                    L2 = P2/lam2 + cp1*65.536+cp0*1e-3
                if lock != 255 and info & 4 != 0:
                    lli = 2

                if D1 != 0.0 and freq1 > 0.0 and freq2 > 0.0 and \
                        (dofst1 != -16 or dop0 != 0):
                    D2 = D1*(lam1/lam2)+dofst1*6.5536+dop0*1e-4

                if cn0 != 255:
                    if sig in (1, 2):
                        S2 = cn0*0.25
                    else:
                        S2 = cn0*0.25 + 10.0

                if code not in self.sig_tab[sys][code.typ]:
                    if self.monlevel > 1:
                        print("skip code={:}".format(code))
                    continue
                idx = self.sig_tab[sys][code.typ].index(code)

                pr[idx] = P2
                cp[idx] = L2
                dp[idx] = D2
                ll[idx] = lli
                cn[idx] = S2

                if self.monlevel >= 2:
                    print("{:6d} {:3d} {:} {:14.3f} {:14.3f} {:10.3f} {:4.2f}".
                          format(int(self.tow), prn, code, P2, L2, D2, S2))

            # Store prn and data
            #
            obs.P = np.append(obs.P, pr)
            obs.L = np.append(obs.L, cp)
            obs.D = np.append(obs.D, dp)
            obs.S = np.append(obs.S, cn)
            obs.lli = np.append(obs.lli, ll)
            obs.sat = np.append(obs.sat, sat)

        obs.P = obs.P.reshape(len(obs.sat), self.nsig[uTYP.C])
        obs.L = obs.L.reshape(len(obs.sat), self.nsig[uTYP.L])
        obs.D = obs.D.reshape(len(obs.sat), self.nsig[uTYP.D])
        obs.S = obs.S.reshape(len(obs.sat), self.nsig[uTYP.S])
        obs.lli = obs.lli.reshape(len(obs.sat), self.nsig[uTYP.L])

        return obs

    def decode_gpsnav(self, buff, k=8):
        sys, prn = self.decode_head(buff, k)
        sat = prn2sat(sys, prn)
        k += 8
        eph = Eph(sat)
        eph.week, eph.l2p, ura, eph.svh, l2data = st.unpack_from(
            '<HBBBB', buff, k)
        k += 6
        eph.iodc, eph.iode, iode3, eph.fit = st.unpack_from('<HBBB', buff, k)
        k += 6
        eph.tgd, toc, eph.af2, eph.af1, eph.af0 = st.unpack_from(
            '<fLfff', buff, k)
        k += 20
        eph.crs, deln, M0, eph.cuc, eph.e, eph.cus = st.unpack_from(
            '<ffdfdf', buff, k)
        k += 32
        eph.deln = deln*rCST.SC2RAD
        eph.M0 = M0*rCST.SC2RAD
        sqrta, eph.toes, eph.cic, OMG0, eph.cis, i0 = st.unpack_from(
            '<dLfdfd', buff, k)
        k += 36
        eph.OMG0 = OMG0*rCST.SC2RAD
        eph.i0 = i0*rCST.SC2RAD

        eph.crc, omg, OMGd, idot, wn_c, wn_e = st.unpack_from(
            '<fdffHH', buff, k)
        k += 28
        eph.omg = omg*rCST.SC2RAD
        eph.OMGd = OMGd*rCST.SC2RAD
        eph.idot = idot*rCST.SC2RAD

        eph.week += (self.week // 1024)*1024
        wn_c += (self.week // 1024)*1024
        wn_e += (self.week // 1024)*1024
        eph.A = sqrta*sqrta
        eph.toc = gpst2time(wn_c, toc)
        eph.toe = gpst2time(wn_e, eph.toes)
        eph.tot = gpst2time(self.week, self.tow)
        eph.mode = 0
        eph.sva = self.rn.urai2sva(ura)

        return eph

    def decode_galnav(self, buff, k=8):
        sys, prn = self.decode_head(buff, k)
        sat = prn2sat(sys, prn)
        k += 7
        eph = Eph(sat)
        src, sqrta, M0, e, i0, omg, OMG0, OMGd, idot, deln = st.unpack_from(
            '<Bddddddfff', buff, k)
        k += 61
        if src != 2:  # decode I/NAV only
            return None

        eph.M0 = M0*rCST.SC2RAD
        eph.i0 = i0**rCST.SC2RAD
        eph.omg = omg*rCST.SC2RAD
        eph.OMG0 = OMG0*rCST.SC2RAD
        eph.OMGd = OMGd*rCST.SC2RAD
        eph.idot = idot*rCST.SC2RAD
        eph.deln = deln*rCST.SC2RAD

        eph.cuc, eph.cus, eph.crc, eph.crs, eph.cic, eph.cis, eph.toes, toc = \
            st.unpack_from('<ffffffLL', buff, k)
        k += 32

        eph.af2, eph.af1, eph.af0, wn_e, wn_c, iodnav, svh = \
            st.unpack_from('<ffdHHHH', buff, k)
        k += 25
        sisa_f, sisa_i, _, eph.tgd, eph.tgd_b, _, cnavenc = \
            st.unpack_from('<BBBfffB', buff, k)
        k += 16

        eph.iode = iodnav
        eph.iodc = iodnav
        eph.week = wn_c
        eph.A = sqrta*sqrta
        eph.toc = gpst2time(wn_c, toc)
        eph.toe = gpst2time(wn_e, eph.toes)
        eph.tot = gpst2time(self.week, self.tow)
        eph.mode = 0
        eph.sva = self.rn.sisa2sva(sisa_i)

        eph.code = (1 << 9) | (1 << 2) | 1

        if svh & 1:
            e1b_dvs = (svh >> 1) & 0x1
            e1b_hs = (svh >> 2) & 0x3
        else:
            e1b_dvs, e1b_hs = 0, 0
        if (svh >> 4) & 1:
            e5b_dvs = (svh >> 5) & 0x1
            e5b_hs = (svh >> 6) & 0x3
        else:
            e5b_dvs, e5b_hs = 0, 0
        if (svh >> 8) & 1:
            e5a_dvs = (svh >> 9) & 0x1
            e5a_hs = (svh >> 10) & 0x3
        else:
            e5a_dvs, e5a_hs = 0, 0

        eph.svh = (e5b_hs << 7) | (e5b_dvs << 6) | (e1b_hs << 1) | (e1b_dvs)
        eph.svh |= (e5a_hs << 4) | (e5a_dvs << 3)

        return eph

    def decode(self, buff, len_, sys=[], prn=[]):
        k = 2
        _, id_, _ = st.unpack_from('<HHH', buff, k)
        k += 6
        blk_num = id_ & 0x1fff
        blk_rev = (id_ >> 13) & 0x7
        if self.monlevel > 1 and blk_num not in \
           (4002, 4004, 4006, 4007, 4017, 4018, 4019, 4020, 4021, 4022, 4023,
            4024, 4026, 4027, 4036, 4047, 4066, 4067, 4068, 4069, 4081, 4093,
                4095, 4218, 4219, 4228, 4242, 4246, 5891, 5894, 5896):
            print("block_num = {:d} rev={:d} len={:d}".format(
                blk_num, blk_rev, len_))

        if blk_num == 4002:  # GAL Decoded Message
            eph = self.decode_galnav(buff, k)
            if self.mode_galinav == 1 and eph is not None:
                self.re.rnx_nav_body(eph, self.fh_rnxnav)

        elif blk_num in (4006, 4007):  # PVT
            tow, wn, mode, err = st.unpack_from('<LHBB', buff, k)
            k += 8
            self.tow = tow*1e-3
            self.week = wn
            x, y, z, und, vx, vy, vz = st.unpack_from('<dddffff', buff, k)
            k += 40
            cog, cb, cd, tsys, datum, nrsv = st.unpack_from('<fdfBBB', buff, k)
            k += 19
            cinfo, refid, age, info, alert, nbase = st.unpack_from(
                '<BHHLBB', buff, k)
            k += 11
            ppp, latency, pacc, vacc, misc = st.unpack_from(
                '<HHHHB', buff, k)
            k += 9
            if blk_num == 4006:
                self.re.pos = np.array([x, y, z])
            else:
                self.re.pos = pos2ecef([x, y, z])

        elif blk_num in (4017, 4066):  # GPSRawCA, QZSRawCA
            sys, prn = self.decode_head(buff, k)
            sat = prn2sat(sys, prn)
            k += 7
            crcpass, _, src, _, ch = st.unpack_from('<BBBBB', buff, k)
            k += 5
            if (self.flg_gpslnav and sys == uGNSS.GPS) or \
                    (self.flg_qzslnav and sys == uGNSS.QZS):
                if crcpass != 1:
                    if self.monlevel > 0:
                        print("crc error in GPSRawCA/QZSRawCA " +
                              "{:6d}\t{:2d}\t{:1d}\t{:2d}".
                              format(int(self.tow), prn, crcpass, src))
                    return -1

                fh_ = self.fh_gpslnav if sys == uGNSS.GPS else self.fh_qzslnav

                blen = (300+7)//8
                fh_.write("{:4d}\t{:6d}\t{:3d}\t{:1d}\t{:3d}\t".
                          format(self.week, int(self.tow), prn, src, blen))

                msg = bytearray(40)
                for i in range(10):
                    d = st.unpack_from('<L', buff, k)[0]
                    fh_.write("{:08x}".format(d))
                    st.pack_into('>L', msg, i*4, d)
                    k += 4
                fh_.write("\n")
                msg = bytes(msg)
                eph = self.rn.decode_gps_lnav(self.week, self.tow, sat, msg)
                if eph is not None:
                    self.re.rnx_nav_body(eph, self.fh_rnxnav)

        elif blk_num in (4018, 4019, 4067, 4068):  # GPSRawL2C/L5, QZSRawL2C/L5
            sys, prn = self.decode_head(buff, k)
            k += 7
            crcpass, cnt, src, freq, ch = st.unpack_from('<BBBBB', buff, k)
            k += 5
            if (sys == uGNSS.GPS and self.flg_gpscnav) or \
               (sys == uGNSS.QZS and self.flg_qzscnav):
                if crcpass != 1:
                    if self.monlevel > 0:
                        print("crc error in GPSRawL2C/L5, QZSRawL2C/L5 " +
                              "{:6d}\t{:2d}\t{:1d}\t{:1d}\t{:2d}".
                              format(int(self.tow), prn, crcpass, cnt, src))
                    return -1

                fh_ = self.fh_gpscnav if sys == uGNSS.GPS else self.fh_qzscnav

                blen = (300+7)//8
                fh_.write("{:4d}\t{:6d}\t{:3d}\t{:1d}\t{:3d}\t".
                          format(self.week, int(self.tow), prn,
                                 src, blen))
                msg = bytearray(40)
                for i in range(10):
                    d = st.unpack_from('<L', buff, k)[0]
                    fh_.write("{:08x}".format(d))
                    st.pack_into('>L', msg, i*4, d)
                    k += 4
                msg = bytes(msg)
                fh_.write("\n")

                sat = prn2sat(sys, prn)
                eph = self.rn.decode_gps_cnav(self.week, self.tow, sat, msg)
                if eph is not None:
                    self.re.rnx_nav_body(eph, self.fh_rnxnav)

        elif blk_num in (4020, 4021):  # GEORawL1/GEORawL5
            sys, prn = self.decode_head(buff, k)
            k += 7
            crcpass, cnt, src, freq, ch = st.unpack_from('<BBBBB', buff, k)
            k += 5
            if self.flg_sbas:
                if crcpass != 1:
                    if self.monlevel > 0:
                        print("crc error in GEORawL1/5 " +
                              "{:6d}\t{:2d}\t{:1d}\t{:1d}\t{:2d}".
                              format(int(self.tow), prn, crcpass, cnt, src))
                    return -1
                self.fh_sbas.write(
                    "{:4d}\t{:6.1f}\t{:3d}\t{:2d}\t{:3d}\t".
                    format(self.week, self.tow, prn, src-24, 32))
                for i in range(8):
                    d = st.unpack_from('<L', buff, k)[0]
                    self.fh_sbas.write("{:08x}".format(d))
                    k += 4
                self.fh_sbas.write("\n")

        elif blk_num == 4022:  # GalRawFNAV
            sys, prn = self.decode_head(buff, k)
            sat = prn2sat(sys, prn)
            k += 7
            crcpass, cnt, src, freq, ch = st.unpack_from('<BBBBB', buff, k)
            k += 5
            if self.flg_galfnav:
                if src & 0x1f == 20:  # E5a
                    type_ = 1
                else:
                    return -1

                if crcpass != 1:
                    if self.monlevel > 0:
                        print("crc error in GALRawFNAV " +
                              "{:6d}\t{:2d}\t{:1d}\t{:1d}\t{:2d}".
                              format(int(self.tow), prn, crcpass,
                                     cnt, src & 0x1f))
                    return -1

                self.fh_galfnav.write("{:4d}\t{:6.1f}\t{:3d}\t{:1d}\t{:3d}\t".
                                      format(self.week, self.tow,
                                             prn, type_, 32))
                msg = bytearray(32)
                for i in range(8):
                    d = st.unpack_from('<L', buff, k)[0]
                    self.fh_galfnav.write("{:08x}".format(d))
                    st.pack_into('>L', msg, i*4, d)
                    k += 4
                msg = bytes(msg)

                self.fh_galfnav.write("\n")

                eph = self.rn.decode_gal_fnav(self.week, self.tow, sat, type_,
                                              msg)
                if eph is not None:
                    self.re.rnx_nav_body(eph, self.fh_rnxnav)

        elif blk_num == 4023:  # GALRawINAV
            sys, prn = self.decode_head(buff, k)
            sat = prn2sat(sys, prn)
            k += 7
            crcpass, cnt, src, freq, ch = st.unpack_from('<BBBBB', buff, k)
            k += 5
            if self.flg_galinav:
                if src & 0x1f == 17:  # E1B
                    type_ = 0
                elif src & 0x1f == 21:  # E5b
                    type_ = 2
                else:
                    return -1

                if crcpass != 1:
                    if self.monlevel > 0:
                        print("crc error in GALRawINAV " +
                              "{:6d}\t{:2d}\t{:1d}\t{:1d}\t{:2d}".
                              format(int(self.tow), prn, crcpass,
                                     cnt, src & 0x1f))
                    return -1

                self.fh_galinav.write("{:4d}\t{:6.1f}\t{:3d}\t{:1d}\t{:3d}\t".
                                      format(self.week, self.tow, prn,
                                             type_, 32))
                msg = bytearray(32)
                for i in range(8):
                    d = st.unpack_from('<L', buff, k)[0]
                    st.pack_into('>L', msg, i*4, d)
                    k += 4

                # GALRawINAV is missing tail bit (6) of even page
                # add 6 bits offset for odd page
                msg_ = bytearray(30)
                msg_[0:15] = msg[0:15]  # even page
                k = 114
                for i in range(15):
                    d = bs.unpack_from('u8', bytes(msg), k)[0]
                    bs.pack_into('u8', msg_, 120+i*8, d)
                    k += 8
                msg_ = bytes(msg_)

                for i in range(30):
                    self.fh_galinav.write("{:08x}".format(msg_[i]))
                self.fh_galinav.write("\n")

                eph = self.rn.decode_gal_inav(self.week, self.tow,
                                              sat, type_, msg_)
                if self.mode_galinav == 0 and eph is not None:
                    self.re.rnx_nav_body(eph, self.fh_rnxnav)

        elif blk_num == 4024:  # GALRawCNAV
            sys, prn = self.decode_head(buff, k)
            k += 7
            crcpass, cnt, src, freq, ch = st.unpack_from('<BBBBB', buff, k)
            k += 5
            if self.flg_gale6:
                if src & 0x1f == 19:
                    type_ = 6
                else:
                    return -1

                if crcpass != 1:
                    if self.monlevel > 0:
                        print("crc error in GALRawCNAV " +
                              "{:6d}\t{:2d}\t{:1d}\t{:1d}\t{:2d}".
                              format(int(self.tow), prn, crcpass, cnt, src))
                    return -1

                blen = (492+7)//8
                self.fh_gale6.write("{:4d}\t{:6d}\t{:3d}\t{:1d}\t{:3d}\t".
                                    format(self.week, int(self.tow), prn,
                                           src, blen))
                for i in range(16):
                    d = st.unpack_from('<L', buff, k)[0]
                    self.fh_gale6.write("{:08x}".format(d))
                    k += 4
                self.fh_gale6.write("\n")

        elif blk_num == 4026:  # GLORawCA
            sys, prn = self.decode_head(buff, k)
            k += 7
            crcpass, cnt, src, freq, ch = st.unpack_from('<BBBBB', buff, k)
            k += 5
            if self.flg_gloca:
                if crcpass != 1:
                    if self.monlevel > 0:
                        print("crc error in GLORawCA " +
                              "{:6d}\t{:2d}\t{:1d}\t{:1d}\t{:2d}".
                              format(int(self.tow), prn, crcpass, cnt, src))
                    return -1

                msg = bytearray(12)
                for i in range(3):
                    d = st.unpack_from('<L', buff, k)[0]
                    st.pack_into('>L', msg, i*4, d)
                    k += 4
                msg = bytes(msg)

                sat = prn2sat(sys, prn)
                geph = self.rn.decode_glo_fdma(self.week, self.tow,
                                               sat, msg, freq)
                if geph is not None:
                    self.re.rnx_gnav_body(geph, self.fh_rnxnav)

        elif blk_num == 4027:  # MeasEpoch
            obs = self.decode_obs(buff, k)
            if obs is not None:
                self.re.rnx_obs_header(obs.time, self.fh_rnxobs)
                self.re.rnx_obs_body(obs, self.fh_rnxobs)

        elif blk_num == 4047:  # BDSRaw
            sys, prn = self.decode_head(buff, k)
            k += 7
            # src 28: B1I (2I), 29: B2I (7I), 30: B3I (6I)
            crcpass, cnt, src, _, ch = st.unpack_from('<BBBBB', buff, k)
            k += 5
            if self.flg_bdsd12 and src == 28:  # only D1 is supported
                if crcpass != 1:
                    if self.monlevel > 0:
                        print("crc error in BDSRaw " +
                              "{:6d}\t{:2d}\t{:1d}\t{:1d}\t{:2d}".
                              format(int(self.tow), prn, crcpass, cnt, src))
                    return -1

                msg = bytearray(40)
                for i in range(10):
                    d = st.unpack_from('<L', buff, k)[0]
                    st.pack_into('>L', msg, i*4, d)
                    k += 4
                msg = bytes(msg)

                sat = prn2sat(sys, prn)
                eph = self.rn.decode_bds_d1(self.week, self.tow, sat, msg)
                if eph is not None:
                    self.re.rnx_nav_body(eph, self.fh_rnxnav)

        elif blk_num == 4069:  # QZSRawL6
            sys, prn = self.decode_head(buff, k)
            k += 7
            parity, rscnt, src, res, ch = st.unpack_from('<BBBBB', buff, k)
            k += 5
            if parity == 0:
                if self.monlevel > 0:
                    print("crc error in QZSRawL6 " +
                          "{:6d}\t{:2d}\t{:1d}\t{:1d}\t{:2d}".
                          format(int(self.tow), prn, parity, rscnt, src))
                return -1

            if self.flg_qzsl6:
                self.fh_qzsl6.write("{:4d}\t{:6.1f}\t{:3d}\t{:1d}\t{:3d}\t".
                                    format(self.week, self.tow, prn, src, 252))
                for i in range(63):
                    d = st.unpack_from('<L', buff, k)[0]
                    self.fh_qzsl6.write("{:08x}".format(d))
                    k += 4
                self.fh_qzsl6.write("\n")

        elif blk_num == 4093:  # NAVICRaw
            sys, prn = self.decode_head(buff, k)
            k += 7
            crcpass, cnt, src, _, ch = st.unpack_from('<BBBBB', buff, k)
            k += 5
            if self.flg_irnnav:
                if crcpass != 1:
                    if self.monlevel > 0:
                        print("crc error in NAVICRaw " +
                              "{:6d}\t{:2d}\t{:1d}\t{:1d}\t{:2d}".
                              format(int(self.tow), prn, crcpass, cnt, src))
                    return -1

                msg = bytearray(40)
                for i in range(10):
                    d = st.unpack_from('<L', buff, k)[0]
                    st.pack_into('>L', msg, i*4, d)
                    k += 4
                msg = bytes(msg)

                sat = prn2sat(sys, prn)
                eph = self.rn.decode_irn_lnav(self.week, self.tow, sat, msg)
                if eph is not None:
                    self.re.rnx_nav_body(eph, self.fh_rnxnav)

        elif blk_num == 4218:  # BDSRawB1C
            sys, prn = self.decode_head(buff, k)
            k += 7
            crcsf2, crcsf3, src, _, ch = st.unpack_from('<BBBBB', buff, k)
            k += 5
            if self.flg_bdsb1c and crcsf2 == 1 and crcsf3 == 1:
                # self.fh_bdsb1c.write("{:4d}\t{:6d}\t{:3d}\t{:1d}\t{:3d}\t".
                #                     format(self.week, int(self.tow), prn,
                #                            src, 225))
                # 1800 deinterleaved symbols of a BeiDou B1C
                # (B-CNAV1) navigation frame
                # 24 unused bits in NAVBits[56]
                # SF1 72sym, SF2 1200sym, SF3 528sym
                # BCH(21,6)+BCH(51,8)  64ary LDPC(200,100) 64ary LDPC(88,44)
                v = bytearray(228)
                for i in range(57):
                    d = st.unpack_from('<L', buff, k)[0]
                    k += 4
                    st.pack_into('>L', v, i*4, d)
                    # self.fh_bdsb1c.write("{:08x}".format(d))
                # self.fh_bdsb1c.write("\n")
                v = bytes(v)
                prn_ = bs.unpack_from('u6', v, 0)[0]
                if prn != prn_:
                    return

                # data2: 600b, errCorr2: 8b, data3: 264b, soh: 8b
                msg = bytearray(228)
                msg[0:75] = v[9:84]
                msg[75] = crcsf2 << 1 | crcsf3
                msg[76:109] = v[159:225]
                msg[109] = bs.unpack_from('u8', v, 21)[0]
                msg = bytes(msg)

                sat = prn2sat(sys, prn)
                eph = self.rn.decode_bds_b1c(self.week, self.tow, sat, msg)
                if eph is not None:
                    self.re.rnx_nav_body(eph, self.fh_rnxnav)

        elif blk_num == 4219:  # BDSRawB2a
            sys, prn = self.decode_head(buff, k)
            k += 7
            crcpass, cnt, src, _, ch = st.unpack_from('<BBBBB', buff, k)
            k += 5
            if self.flg_bdsb2a:
                if crcpass != 1:
                    if self.monlevel > 0:
                        print("crc error in BDSRawB2a " +
                              "{:6d}\t{:2d}\t{:1d}\t{:1d}\t{:2d}".
                              format(int(self.tow), prn, crcpass, cnt, src))
                    return -1

                msg = bytearray(40)
                for i in range(10):
                    d = st.unpack_from('<L', buff, k)[0]
                    st.pack_into('>L', msg, i*4, d)
                    k += 4
                msg = bytes(msg)

                sat = prn2sat(sys, prn)
                eph = self.rn.decode_bds_b2a(self.week, self.tow, sat, msg)
                if eph is not None:
                    self.re.rnx_nav_body(eph, self.fh_rnxnav)

        elif blk_num in (4221, 4227):  # GPSRawL1C, QZSRawL1C
            sys, prn = self.decode_head(buff, k)
            k += 7
            crcsf2, crcsf3, src, _, ch = st.unpack_from('<BBBBB', buff, k)
            k += 5
            if (sys == uGNSS.GPS and self.flg_gpscnav2) or \
               (sys == uGNSS.QZS and self.flg_qzscnav2):
                if crcsf2 != 1 or crcsf3 != 1:
                    if self.monlevel > 0:
                        print("crc error in GPSRawL1C, QZSRawL1C " +
                              "{:6d}\t{:2d}\t{:1d}\t{:1d}\t{:2d}".
                              format(int(self.tow), prn,
                                     crcsf2 << 1 | crcsf3, cnt, src))
                    return -1

                fh_ = self.fh_gpscnav2 if sys == uGNSS.GPS \
                    else self.fh_qzscnav2

                blen = (1800+7)//8
                fh_.write("{:4d}\t{:6d}\t{:3d}\t{:1d}\t{:3d}\t".
                          format(self.week, int(self.tow), prn, src, blen))
                msg = bytearray(228)
                for i in range(57):
                    d = st.unpack_from('<L', buff, k)[0]
                    fh_.write("{:08x}".format(d))
                    st.pack_into('>L', msg, i*4, d)
                    k += 4
                msg = bytes(msg)
                fh_.write("\n")

                sat = prn2sat(sys, prn)
                eph = self.rn.decode_gps_cnav2(self.week, self.tow, sat, msg)
                if eph is not None:
                    self.re.rnx_nav_body(eph, self.fh_rnxnav)

        elif blk_num in (4228, 4246):  # QZSRawL1S, QZSRawL5S
            src_t = {24: 0, 25: 1, 33: 2, 39: 3}  # L1C/A, L5, L1S, L5S
            sys, prn = self.decode_head(buff, k)
            k += 7
            crc, cnt, src, freq, ch = st.unpack_from('<BBBBB', buff, k)
            k += 5
            if sys == uGNSS.QZS and (self.flg_qzsl1s or self.flg_qzsl5s):
                if crc != 1:
                    if self.monlevel > 0:
                        print("crc error in QZSRawL1S, QZSRawL5S " +
                              "{:6d}\t{:2d}\t{:1d}\t{:1d}\t{:2d}".
                              format(int(self.tow), prn, crc, cnt, src))
                    return -1

                if src not in src_t.keys():
                    if self.monlevel > 0:
                        print("src not recgonised in QZSRawL1S/QZSRawL5S " +
                              "{:2d}".format(src))
                    return -1

                fh_ = self.fh_sbas

                blen = (250+7)//8
                fh_.write("{:4d}\t{:6.1f}\t{:3d}\t{:2d}\t{:3d}\t".
                          format(self.week, self.tow, prn, src_t[src], blen))

                msg = bytearray(32)
                for i in range(8):
                    d = st.unpack_from('<L', buff, k)[0]
                    fh_.write("{:08x}".format(d))
                    st.pack_into('>L', msg, i*4, d)
                    k += 4
                fh_.write("\n")

        elif blk_num == 4242:  # BDSRawB2b
            sys, prn = self.decode_head(buff, k)
            k += 7
            crcpass, _, src, _, ch = st.unpack_from('<BBBBB', buff, k)
            k += 5

            if self.flg_bdsb2b:
                flg_bdsppp = True if prn >= 59 else False

                if crcpass != 1:
                    if self.monlevel > 0:
                        print("crc error in BDSRawB2b " +
                              "{:6d}\t{:2d}\t{:1d}\t{:2d}".
                              format(int(self.tow), prn, crcpass, src))
                    return -1

                if flg_bdsppp:
                    self.fh_bdsb2b.write("{:4d}\t{:6d}\t{:3d}\t{:1d}\t{:3d}\t".
                                         format(self.week, int(self.tow), prn,
                                                src, 64))
                # 984 symbols of a BeiDou B2b navigation frame
                # 8 unused bits in NAVBits[30]
                msg = bytearray(64)
                for i in range(16):
                    d = st.unpack_from('<L', buff, k)[0]
                    st.pack_into('>L', msg, i*4, d)
                    k += 4

                    if flg_bdsppp:
                        if i == 0:
                            self.fh_bdsb2b.write("{:05x}".format(d & 0xfffff))
                        elif i == 15:
                            self.fh_bdsb2b.write(
                                "{:05x}{:06x}".format((d >> 12) & 0xffffc, 0))
                        else:
                            self.fh_bdsb2b.write("{:08x}".format(d))
                msg = bytes(msg)

                if flg_bdsppp:
                    self.fh_bdsb2b.write("\n")
                else:  # B2b-CNAV3
                    sat = prn2sat(sys, prn)
                    eph = self.rn.decode_bds_b2b(self.week, self.tow, sat, msg)
                    if eph is not None:
                        self.re.rnx_nav_body(eph, self.fh_rnxnav)

        elif blk_num in (4095, 5891):  # GPS/QZS Decoded Message
            eph = self.decode_gpsnav(buff, k)
            if eph is not None:
                self.re.rnx_nav_body(eph, self.fh_rnxnav)

        elif blk_num == 5894:  # GPS UTC Decoded Message
            None

        elif blk_num == 5896:  # SBAS L1 Decoded Message
            None

        return 0


if __name__ == "__main__":

    bdir = os.path.expanduser('~/Projects/CSSRlib/sbf/_sbf/')
    fnames = 'sep3238a.sbf'

    gnss_t = 'GERCJ'

    # bdir = '../data/doy244/'
    # fnames = 'sep3244*.sbf'

    # bdir = '../data/doy223/'
    # fnames = 'sept223v.sbf'

    # bdir = '../data/doy308/'
    # fnames = 'sept308b.sbf'

    opt = rcvOpt()
    opt.flg_qzsl6 = True
    opt.flg_qzslnav = True
    opt.flg_gpslnav = True
    opt.flg_qzscnav = True
    opt.flg_gpscnav = True
    opt.flg_qzscnav2 = True
    opt.flg_gpscnav2 = True
    opt.flg_qzsl1s = False
    opt.flg_qzsl5s = False
    opt.flg_gale6 = True
    opt.flg_galinav = True
    opt.flg_galfnav = True
    opt.flg_bdsb1c = True
    opt.flg_bdsb2a = True
    opt.flg_bdsb2b = True
    opt.flg_bdsd12 = True
    opt.flg_gloca = True
    opt.flg_irnnav = False
    opt.flg_sbas = False
    opt.flg_rnxnav = True
    opt.flg_rnxobs = True

    for f in glob(bdir+fnames):

        print("Decoding {}".format(f))

        bdir, fname = os.path.split(f)
        bdir += '/'

        prefix = bdir+fname[4:].removesuffix('.sbf')+'_'
        sbfdec = sbf(opt, prefix=prefix, gnss_t=gnss_t)
        sbfdec.monlevel = 1
        nep = 0
        nep_max = 0

        sbfdec.re.anttype = "JAVRINGANT_DM   JVDM"
        sbfdec.re.rectype = "SEPT POLARX5        "

        blen = os.path.getsize(bdir+fname)
        with open(bdir+fname, 'rb') as f:
            msg = f.read(blen)
            maxlen = len(msg)-5
            # maxlen = 400000
            k = 0
            while k < maxlen:
                stat = sbfdec.sync(msg, k)
                if not stat:
                    k += 1
                    continue
                if not sbfdec.check_crc(msg, k):
                    continue
                len_ = sbfdec.msg_len(msg, k)
                if k+len_ >= maxlen:
                    break

                sbfdec.decode(msg[k:k+len_], len_)
                k += len_

                nep += 1
                if nep_max > 0 and nep >= nep_max:
                    break

        sbfdec.file_close()
