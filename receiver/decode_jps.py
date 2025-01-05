"""
Javad Receiver GREIS messages decoder

[1] GREIS: GNSS Receiver External Interface Specification for 4.5.00,
     October, 2024

@author: Rui Hirokawa
"""

import os
import numpy as np
import struct as st
import bitstruct.c as bs
from cssrlib.gnss import epoch2time, time2gpst, prn2sat, uGNSS, uTYP, rSigRnx
from cssrlib.gnss import Obs, rCST, gpst2time, uSIG, copy_buff
from cssrlib.rawnav import rcvDec, rcvOpt
from glob import glob
from enum import IntEnum
from binascii import hexlify


def istxt(c):
    if '0' <= chr(c) and chr(c) <= '~':
        return True
    return False


def ishex(c):
    if '0' <= chr(c) and chr(c) <= '9':
        return True
    if 'A' <= chr(c) and chr(c) <= 'F':
        return True
    return False


class GNSS(IntEnum):
    GPS = 1
    GLO = 2
    SBS = 3
    GAL = 4
    QZS = 5
    BDS = 6
    IRN = 7
    LBD = 8
    GLO_C = 9


# from rtklib.h
class CODE(IntEnum):
    L1C = 1
    L1P = 2
    L1W = 3
    L1X = 12
    L1Z = 13
    L2C = 14
    L2X = 18
    L2P = 19
    L2W = 20
    L5X = 26
    L7I = 27
    L7X = 29
    L6X = 33
    L8X = 39
    L2I = 40
    L6I = 42
    L3X = 46
    L1I = 47
    L5A = 49
    L9A = 52


# from Table 3-8 in [1]
class SIG(IntEnum):
    QZS_L1C = 0
    QZS_L1Z = 1
    QZS_L6X = 2
    QZS_L2X = 3
    QZS_L5X = 4
    QZS_L1X = 5
    GPS_L1C = 0
    GPS_L1W = 1
    GPS_L2W = 2
    GPS_L2X = 3
    GPS_L5X = 4
    GPS_L1X = 5
    SBS_L1C = 0
    SBS_L5X = 4


class jps(rcvDec):
    tow = -1
    week = -1
    sat = []
    sys = []
    prn = []
    nsat = 0
    nmax = 96
    nsigmax = 7

    pr = []
    cp = []
    dp = []
    ch_t = {'X': 0, 'x': 0, '0': 0, 'C': 0, 'c': 0, '1': 1, '2': 2, '3': 3,
            '5': 4, 'l': 5}
    fn_t = [1.57542e9, 1.57542e9, 1.22760e9, 1.22760e9, 1.17645e9, 1.27875e9]
    data_L6D = bytearray(b'\x00'*(250))
    data_L6E = bytearray(b'\x00'*(250))

    types = [[uSIG.L1C, uSIG.L1W, uSIG.L2W, uSIG.L2X, uSIG.L5X, uSIG.L1X],
             [uSIG.L1C, uSIG.L1P, uSIG.L2P, uSIG.L2C, uSIG.L3X,        0],
             [uSIG.L1C,        0,        0,        0, uSIG.L5X,        0],
             [uSIG.L1X, uSIG.L8X, uSIG.L7X, uSIG.L6X, uSIG.L5X,        0],
             [uSIG.L1C, uSIG.L1Z, uSIG.L6X, uSIG.L2X, uSIG.L5X, uSIG.L1X],
             [uSIG.L2I, uSIG.L8X, uSIG.L7I, uSIG.L6I, uSIG.L5X, uSIG.L1X],
             [uSIG.L9A,        0,        0,        0, uSIG.L5A, uSIG.L1X],
             [0,        0,        0,        0,        0,        0],
             [uSIG.L4X,        0,        0, uSIG.L6X, uSIG.L3X,        0]]
    freqs = [[1, 1, 2, 2, 3, 1], [1, 1, 2, 2, 3, 0], [1, 0, 0, 0, 3, 0],
             [1, 6, 2, 4, 3, 0], [1, 1, 4, 2, 3, 1], [1, 6, 2, 3, 3, 1],
             [0, 0, 0, 0, 3, 1]]

    sys_t = {
        GNSS.GPS: uGNSS.GPS, GNSS.GLO: uGNSS.GLO, GNSS.GAL: uGNSS.GAL,
        GNSS.BDS: uGNSS.BDS, GNSS.QZS: uGNSS.QZS, GNSS.SBS: uGNSS.SBS,
        GNSS.IRN: uGNSS.IRN, GNSS.GLO_C: uGNSS.GLO,
    }

    rec = []
    mid_decoded = []

    def __init__(self, opt=None, prefix=''):
        super().__init__()

        self.pr_ref = np.zeros(self.nmax)
        self.PR_REF = np.zeros(self.nmax)
        self.pr = np.zeros((self.nmax, self.nsigmax))
        self.cp = np.zeros((self.nmax, self.nsigmax))
        self.dp = np.zeros((self.nmax, self.nsigmax))
        self.code = np.zeros((self.nmax, self.nsigmax), dtype=int)
        self.prc = np.zeros((self.nmax, self.nsigmax))
        self.cpc = np.zeros((self.nmax, self.nsigmax))
        self.CNO = np.zeros((self.nmax, self.nsigmax))
        self.CNOd = np.zeros((self.nmax, self.nsigmax))

        self.nsat = 0
        self.qzl6_time_p = -1
        self.tod = -1
        self.freqn = None

        self.sig_tab = {
            uGNSS.GPS: {
                uTYP.C: [rSigRnx('GC1C'), rSigRnx('GC1W'), rSigRnx('GC2W'),
                         rSigRnx('GC2X'), rSigRnx('GC5X')],
                uTYP.L: [rSigRnx('GL1C'), rSigRnx('GL1W'), rSigRnx('GL2W'),
                         rSigRnx('GL2X'), rSigRnx('GL5X')],
                uTYP.D: [rSigRnx('GD1C'), rSigRnx('GD1W'), rSigRnx('GD2W'),
                         rSigRnx('GD2X'), rSigRnx('GD5X')],
                uTYP.S: [rSigRnx('GS1C'), rSigRnx('GS1W'), rSigRnx('GS2W'),
                         rSigRnx('GS2X'), rSigRnx('GS5X')],
            },
            uGNSS.GLO: {
                uTYP.C: [rSigRnx('RC1C'), rSigRnx('RC1P'), rSigRnx('RC2C'),
                         rSigRnx('RC2P'), rSigRnx('RC3X')],
                uTYP.L: [rSigRnx('RL1C'), rSigRnx('RL1P'), rSigRnx('RL2C'),
                         rSigRnx('RL2P'), rSigRnx('RL3X')],
                uTYP.D: [rSigRnx('RD1C'), rSigRnx('RD1P'), rSigRnx('RD2C'),
                         rSigRnx('RD2P'), rSigRnx('RD3X')],
                uTYP.S: [rSigRnx('RS1C'), rSigRnx('RS1P'), rSigRnx('RS2C'),
                         rSigRnx('RS2P'), rSigRnx('RS3X')],
            },
            uGNSS.GAL: {
                uTYP.C: [rSigRnx('EC1X'), rSigRnx('EC5X'), rSigRnx('EC7X'),
                         rSigRnx('EC8X'), rSigRnx('EC6X')],
                uTYP.L: [rSigRnx('EL1X'), rSigRnx('EL5X'), rSigRnx('EL7X'),
                         rSigRnx('EL8X'), rSigRnx('EL6X')],
                uTYP.D: [rSigRnx('ED1X'), rSigRnx('ED5X'), rSigRnx('ED7X'),
                         rSigRnx('ED8X'), rSigRnx('ED6X')],
                uTYP.S: [rSigRnx('ES1X'), rSigRnx('ES5X'), rSigRnx('GS7X'),
                         rSigRnx('ES8X'), rSigRnx('ES6X')],
            },
            uGNSS.BDS: {
                uTYP.C: [rSigRnx('CC1X'), rSigRnx('CC2I'), rSigRnx('CC5X'),
                         rSigRnx('CC6I'), rSigRnx('CC8X'), rSigRnx('CC7I')],
                uTYP.L: [rSigRnx('CL1X'), rSigRnx('CL2I'), rSigRnx('CL5X'),
                         rSigRnx('CL6I'), rSigRnx('CL8X'), rSigRnx('CL7I')],
                uTYP.D: [rSigRnx('CD1X'), rSigRnx('CD2I'), rSigRnx('CD5X'),
                         rSigRnx('CD6I'), rSigRnx('CD8X'), rSigRnx('CD7I')],
                uTYP.S: [rSigRnx('CS1X'), rSigRnx('CS2I'), rSigRnx('CS5X'),
                         rSigRnx('CS6I'), rSigRnx('CS8X'), rSigRnx('CS7I')],
            },
            uGNSS.QZS: {
                uTYP.C: [rSigRnx('JC1C'), rSigRnx('JC1X'), rSigRnx('JC2X'),
                         rSigRnx('JC5X'), rSigRnx('JC6X')],
                uTYP.L: [rSigRnx('JL1C'), rSigRnx('JL1X'), rSigRnx('JL2X'),
                         rSigRnx('JL5X'), rSigRnx('JL6X')],
                uTYP.D: [rSigRnx('JD1C'), rSigRnx('JD1X'), rSigRnx('JD2X'),
                         rSigRnx('JD5X'), rSigRnx('JD6X')],
                uTYP.S: [rSigRnx('JS1C'), rSigRnx('JS1X'), rSigRnx('JS2X'),
                         rSigRnx('JS5X'), rSigRnx('JS6X')],
            },
            uGNSS.SBS: {
                uTYP.C: [rSigRnx('SC1C'), rSigRnx('SC5X')],
                uTYP.L: [rSigRnx('SL1C'), rSigRnx('SL5X')],
                uTYP.D: [rSigRnx('SD1C'), rSigRnx('SD5X')],
                uTYP.S: [rSigRnx('SS1C'), rSigRnx('SS5X')],
            },
            uGNSS.IRN: {
                uTYP.C: [rSigRnx('IC5A'), rSigRnx('IC1X')],
                uTYP.L: [rSigRnx('IL5A'), rSigRnx('IL1X')],
                uTYP.D: [rSigRnx('ID5A'), rSigRnx('ID1X')],
                uTYP.S: [rSigRnx('IS5A'), rSigRnx('IS1X')],
            },
        }

        if opt is not None:
            self.init_param(opt=opt, prefix=prefix)

    def crc8(self, src, cnt):
        res = 0
        for k in range(cnt):
            res = (((res << 2) | (res >> 6)) ^ src[k]) & 0xff
        return ((res << 2) | (res >> 6)) & 0xff

    def sync(self, buff, k):
        c = chr(buff[k])
        if (c == '\r' or c == '\n') and istxt(buff[k+1]) and istxt(buff[k+2]) \
           and ishex(buff[k+3]) and ishex(buff[k+4]) and ishex(buff[k+5]):
            return True
        return False

    def decode_eph(self, buff, sys, len_):
        j = 5
        tgd = [0, 0, 0]
        sv, tow, flags = st.unpack_from('<BLB', buff, j)
        j += 6
        iodc, toc, ura, health, wn, tgd[0], af2, af1, af0 = \
            st.unpack_from('<hlbBHffff', buff, j)
        j += 26
        toe, iode, rootA, ecc, m0, omega0, inc0, argPer = \
            st.unpack_from('<lhdddddd', buff, j)
        j += 54
        deln, omegaDot, idot, crc, crs, cuc, cus, cic, cis = \
            st.unpack_from('<fffffffff', buff, j)
        j += 36
        # print('sys={:d},eph prn={:3d} len={:d}'.format(sys,sv,len_))
        if sys == GNSS.GPS or sys == GNSS.QZS:
            if len_ <= 128:
                return
            navType, toe, toc, Adot, deln_dot, ura_oe, ura_oc, \
                ura_oc1, ura_oc2 = st.unpack_from('<Blldfbbbb', buff, j)
            j += 25
            if navType == 1 or navType == 2:  # CNAV
                isc_L1CA, isc_L2C, isc_L5I, isc_L5Q = \
                    st.unpack_from('<ffff', buff, j)
                j += 16
            elif navType == 3:  # CNAV2
                isc_L1CP, isc_L1CD = st.unpack_from('<ff', buff, j)
                j += 8
            daf0 = st.unpack_from('<f', buff, j)
            j += 4
        elif sys == GNSS.GAL:
            tgd[1], tgd[2], ai0, ai1, ai2, sfi, navType, daf0 = \
                st.unpack_from('<fffffBBf', buff, j)
            j += 26

    def tofreq(self, sig, sys):
        # GPS=1;GLO=2;SBS=3;GAL=4;QZS=5;BDS=6;IRN=7;LBD=8;GLO_C=9
        i = self.ch_t[sig]
        type_ = self.types[sys-1][i]
        freq = self.freqs[sys-1][i]-1
        return freq, type_

    def freq_sys(self, sys, freq, freqn):
        fn_t = [1.57542E9, 1.22760E9, 1.17645E9, 1.27875E9, 1.20714E9,
                1.191795E9, 2.492028E9]
        if sys == GNSS.GLO:
            if freq == 0:   # L1 FDMA
                fn = 1.60200E9+0.56250E6*freqn
            elif freq == 1:  # L2 FDMA
                fn = 1.24600E9+0.43750E6*freqn
            else:    # L3 CDMA
                fn = 1.202025E9
        elif sys == GNSS.BDS:
            if freq == 0:
                fn = 1.561098E9
            elif freq == 1:
                fn = 1.20714E9
            else:
                fn = 1.26852E9
        elif sys == GNSS.GAL and freq == 1:  # E5b
            fn = 1.20714E9
        else:
            fn = fn_t[freq]
        return fn

    def decode_obs(self):
        obs = Obs()
        obs.sig = self.sig_tab

        obs.time = gpst2time(self.week, self.tow)

        nsig_max = 0
        for s in self.sig_tab:
            if len(self.sig_tab[s][uTYP.L]) > nsig_max:
                nsig_max = len(self.sig_tab[s][uTYP.L])

        self.nsig[uTYP.C] = nsig_max
        self.nsig[uTYP.L] = nsig_max
        self.nsig[uTYP.D] = nsig_max
        self.nsig[uTYP.S] = nsig_max

        obs.sat = []
        nsat = len(self.prn)

        obs.P = np.zeros((nsat, self.nsig[uTYP.C]), dtype=np.float64)
        obs.L = np.zeros((nsat, self.nsig[uTYP.L]), dtype=np.float64)
        obs.D = np.zeros((nsat, self.nsig[uTYP.D]), dtype=np.float64)
        obs.S = np.zeros((nsat, self.nsig[uTYP.S]), dtype=np.float64)
        obs.lli = np.zeros((nsat, self.nsig[uTYP.L]), dtype=np.int32)
        j = 0
        kr = 0
        for k in range(nsat):
            if self.sys[k] == GNSS.GLO:
                prn = self.osn[kr]
                kr += 1
            else:
                prn = self.prn[k]

            sys = self.sys_t[self.sys[k]]
            if sys not in self.sig_tab.keys():
                continue
            if sys == uGNSS.SBS and self.prn[k] > 156:
                continue
            sat = prn2sat(sys, prn)
            if sat in obs.sat:
                jn = obs.sat.index(sat)
            else:
                jn = j

            for kk, sig_ in enumerate(self.sig_tab[sys][uTYP.L]):
                if sig_.sig in self.types[self.sys[k]-1]:
                    idx = self.types[self.sys[k]-1].index(sig_.sig)
                    obs.P[jn][kk] = self.pr[k][idx]*rCST.CLIGHT
                    obs.L[jn][kk] = self.cp[k][idx]
                    obs.D[jn][kk] = self.dp[k][idx]
                    obs.S[jn][kk] = self.CNO[k][idx]
                    if self.code[k][idx] & 0x0020:
                        obs.lli[jn][kk] += 1

            if sat not in obs.sat:
                obs.sat += [sat]
                j += 1

        nsat = len(obs.sat)
        obs.P = obs.P[:nsat]
        obs.L = obs.L[:nsat]
        obs.D = obs.D[:nsat]
        obs.S = obs.S[:nsat]
        obs.lli = obs.lli[:nsat]

        obs.P[np.isnan(obs.P)] = 0
        obs.L[np.isnan(obs.L)] = 0
        obs.D[np.isnan(obs.D)] = 0
        obs.S[np.isnan(obs.S)] = 0

        return obs

    def decode_nd(self, buff, sys=uGNSS.GPS):
        prn, time_, type_, len_ = st.unpack_from('<BLBB', buff, 5)
        sat = prn2sat(sys, prn)
        if self.monlevel >= 2:
            msg_ = 'gd' if sys == uGNSS.GPS else 'qd'
            print("[{:s}] prn={:2d} tow={:6d} type={:d} len={:d}".
                  format(msg_, prn, time_, type_, len_))
        msg = st.unpack_from('>'+len_*'L', buff, 12)
        b = bytes(np.array(msg, dtype='uint32'))
        # type_= 0: LNAV, 1:L2C CNAV, 2:L5 CNAV, 3: L1C CNAV2

        if self.flg_rnxnav:
            eph = None
            if type_ == 0:  # LNAV
                eph = self.rn.decode_gps_lnav(self.week, time_, sat, b)
            elif type_ == 1 or type_ == 2:  # CNAV
                eph = self.rn.decode_gps_cnav(self.week, time_, sat, b)
            elif type_ == 3:  # CNAV2
                msg = bytearray(228)  # recover original L1C msg structure
                # toi: 9b, data2: 600b, data3: 274b
                toi = bs.unpack_from('u9', b, 0)[0]
                bs.pack_into('u9', msg, 0, toi)
                for k in range(75):
                    b_ = bs.unpack_from('u8', b, k*8+9)[0]
                    bs.pack_into('u8', msg, k*8+52, b_)
                for k in range(35):
                    b_ = bs.unpack_from('u8', b, k*8+609)[0]
                    bs.pack_into('u8', msg, k*8+1252, b_)
                msg = bytes(msg)
                eph = self.rn.decode_gps_cnav2(self.week, time_, sat, msg)
            if eph is not None:
                self.re.rnx_nav_body(eph, self.fh_rnxnav)

        if sys == uGNSS.GPS and self.flg_gpslnav:
            self.fh_gpslnav.write(
                "{:4d}\t{:6d}\t{:3d}\t{:1d}\t{:3d}\t{:s}\n".
                format(self.week, time_, prn, type_, len_*4,
                       hexlify(b).decode()))

    def decode(self, buff, len_, sys=[], prn=[]):
        head = buff[0:2].decode()
        if head == 'RE':
            if self.monlevel > 1:
                print("[{:2s}] {:}".format(head, buff[5:5+len_].decode()))

        elif self.crc8(buff, len_-1) != buff[-1]:
            print("cs error")
            return -1

        if head == '~~':  # receiver time [RT] (epoch start)
            self.tod = st.unpack_from('<L', buff, 5)[0]*1e-3
            if self.monlevel > 1:
                print("[RT] tod={:.1f}".format(self.tod))
        elif head == 'RE':
            return
        elif head == '::':  # epoch time [ET] (epoch end)
            self.tod = st.unpack_from('<L', buff, 5)[0]*1e-3
            if self.monlevel > 1:
                print("[ET] tod={:.1f}".format(self.tod))

            obs = self.decode_obs()
            if self.flg_rnxobs and obs is not None:
                self.obs = obs
                self.re.rnx_obs_header(obs.time, self.fh_rnxobs)
                self.re.rnx_obs_body(obs, self.fh_rnxobs)

        elif head == 'GT':  # GPS time [GT]
            tow, wn, cycle = st.unpack_from('<LHB', buff, 5)  # ms
            self.week = wn+cycle*1024
            self.tow = tow*1e-3
            if self.monlevel > 1:
                print("[GT] tow={:.1f} week={:4d}".format(self.tow, self.week))
        elif head == 'RD':  # Receiver Date and Receiver Time
            # base 0:GPS,1:UTC USNO,2:GLO,3:UTC SU
            year, month, day, base = st.unpack_from('<HBBB', buff, 5)
            if self.tod >= 0:
                h = self.tod // 3600
                tmp = (self.tod % 3600)
                m = tmp//60
                s = tmp-m*60
                ep = [year, month, day, h, m, s]
                self.week, self.tow = time2gpst(epoch2time(ep))

            if self.monlevel > 1:
                print("[RD] {:d}/{:d}/{:d} {:d}".
                      format(year, month, day, base))
        elif head == 'SX':  # satellite index
            sys = []
            prn = []
            nsat = (len_-6)//2
            esi = st.unpack_from('<'+'B'*2*nsat, buff, 5)
            self.nsat = nsat
            self.freqn = np.zeros(nsat, dtype=int)
            for k in range(nsat):
                ssid = esi[k*2]
                svid = esi[k*2+1]
                if ssid == GNSS.GLO:
                    if svid > 127:
                        self.freqn[k] = svid-256
                    else:
                        self.freqn[k] = svid
                sys += [ssid]
                prn += [svid]
                #  print("{:d} {:3d}".format(sys[k],prn[k]))
            self.sys = sys
            self.prn = prn
        elif head == 'SI':  # satellite index (obsolete)
            nsat = (len_-6)//1
            # si = st.unpack_from('b'*nsat, buff, 5)
        elif head == 'EL':  # elevation
            nsat = (len_-6)//1
            # el = st.unpack_from('b'*nsat, buff, 5)
        elif head == 'AZ':  # azimuth
            nsat = (len_-6)//1
            # az = st.unpack_from('b'*nsat, buff, 5)
        elif head == 'DP':  # DOP
            nsat = (len_-6)//1
            hdop, vdop, tdop, soltype, edop = st.unpack_from('<fffBf', buff, 5)
        elif head == 'TO':  # Reference Time to Receiver Time Offset
            nsat = (len_-6)//1
            val, sval = st.unpack_from('dd', buff, 5)
        elif head == 'DO':  # Derivative of Receiver Time Offset
            nsat = (len_-6)//1
            val, sval = st.unpack_from('ff', buff, 5)
        elif head == 'PM':  # parameters
            None
            # print(buff[5:len-4])
        elif head == 'MF':  # messages format
            if self.monlevel >= 1:
                print(buff[5:len_-2])
        elif head == 'NN':  # GLONASS Satellite system number
            self.nsat_glo = (len_-6)//1
            self.osn = st.unpack_from('B'*self.nsat_glo, buff, 5)

        elif head == 'qd':  # QZSS Raw Navigation Data
            self.decode_nd(buff, sys=uGNSS.QZS)

        elif head == 'xd':  # QZSS L6 Message Data
            prn, time_, type_, len_ = st.unpack_from('<BLBB', buff, 5)
            if self.monlevel > 1:
                print("[xd] prn={:d} tow={:d} type={:d} len={:d}".
                      format(prn, time_, type_, len_))
            if self.week >= 0:
                if self.flg_qzsl6:
                    if prn_ref > 0 and prn != prn_ref:
                        return
                    msg_l6 = buff[12:12+len_]
                    self.fh_qzsl6.write(
                        "{:4d}\t{:6d}\t{:3d}\t{:1d}\t{:3d}\t{:s}\n".
                        format(self.week, time_, prn, type_, len_,
                               hexlify(msg_l6).decode()))

            # errCorr = st.unpack_from('<B', buff, 12+len_)
        elif head == 'cd':  # BeiDou Navigation data
            prn, time_, type_, len_ = st.unpack_from('<BLBB', buff, 5)
            time_ = (time_ + 14) % 604800  # BDST -> GPST
            ch = type_ & 0x3f
            B2bq = (type_ >> 6) & 1
            D2 = (type_ >> 7) & 1
            if self.monlevel >= 2:
                s = "{:2d} B2bq:{:1d} D2-GEO:{:1d}".format(ch, B2bq, D2)
                print("[cd] prn={:2d} tow={:6d} type={:s} len={:d}".
                      format(prn, time_, s, len_))

            # type_[0:5] 0:B1,1:B2,2:B3,3:B1C,5:B2a,6:B2b
            msg = st.unpack_from('>'+len_*'L', buff, 12)
            b = bytes(np.array(msg, dtype='uint32'))

            sat = prn2sat(uGNSS.BDS, prn)
            if self.week < 0:
                return

            if self.flg_rnxnav:
                eph = None
                if ch == 0:  # B1 (D1/D2)
                    if D2 == 0:
                        eph = self.rn.decode_bds_d1(self.week, time_, sat, b)
                    else:
                        eph = self.rn.decode_bds_d2(self.week, time_, sat, b)
                elif ch == 3:  # B1C
                    eph = self.rn.decode_bds_b1c(self.week, time_, sat, b)
                elif ch == 5:  # B2a
                    eph = self.rn.decode_bds_b2a(self.week, time_, sat, b)
                elif ch == 6 and prn < 59 and B2bq == 0:  # B2b
                    eph = self.rn.decode_bds_b2b(self.week, time_, sat, b, 0)

                if eph is not None:
                    self.re.rnx_nav_body(eph, self.fh_rnxnav)

            if ch == 6 and self.flg_bdsb2b and prn >= 59:  # B2b: BDS PPP
                self.fh_bdsb2b.write(
                    "{:4d}\t{:6d}\t{:3d}\t{:1d}\t{:3d}\t{:s}\n".
                    format(self.week, time_, prn, type_, len_*4,
                           hexlify(b).decode()))

        elif head == 'gd':  # GPS Navigation data
            self.decode_nd(buff, sys=uGNSS.GPS)

        elif head == 'id':  # NavIC Navigation data
            prn, time_, type_, len_ = st.unpack_from('<BLBB', buff, 5)
            sat = prn2sat(uGNSS.IRN, prn)
            # type 0 - L5, 1 - S, 2 - L1
            msg = st.unpack_from('>'+len_*'L', buff, 12)
            b = bytes(np.array(msg, dtype='uint32'))

            if self.flg_rnxnav:
                eph = None
                if type_ == 0:
                    eph = self.rn.decode_irn_lnav(self.week, time_, sat, b)
                elif type_ == 2:
                    # for L1
                    # data[0] – subframe 1 (toi)
                    # data[1…19] – subframe 2
                    # data[20…28] – subframe 3

                    msg = bytearray(228)  # recover original L1C msg structure
                    # toi: 9b, data2: 600b, data3: 274b
                    toi = bs.unpack_from('u32', b, 0)[0]
                    bs.pack_into('u9', msg, 0, toi)
                    copy_buff(b, msg, 32, 52, 600)
                    copy_buff(b, msg, 640, 1252, 274)
                    msg = bytes(msg)
                    eph = self.rn.decode_irn_l1nav(self.week, time_, sat, msg)

                if eph is not None:
                    self.re.rnx_nav_body(eph, self.fh_rnxnav)

            if self.monlevel >= 2:
                print(f"[id] time={time_:6d} prn={prn:2d} type={type_}")

        elif head == 'lD':  # Glonass Raw Navigation data
            svn, fcn, time_, type_, len_ = st.unpack_from('<BbLBB', buff, 5)
            # type 0 - L1, 2 - L2C, 3 - P1, 4 - P2
            msg = st.unpack_from('>'+len_*'L', buff, 13)
            b = bytes(np.array(msg, dtype='uint32'))
            sat = prn2sat(uGNSS.GLO, svn)

            # get 77 bit (25x3+2) in frame without hamming and time mark
            if type_ == 0:
                buff = bytearray(12)
                for k in range(4):
                    d = bs.unpack_from('u32', b, 32*k)[0]
                    if k < 3:
                        bs.pack_into('u25', buff, 25*k, d & 0x1ffffff)
                    else:
                        bs.pack_into('u2', buff, 25*k, (d >> 23) & 0x3)

            if self.flg_rnxnav and type_ == 0:
                geph = self.rn.decode_glo_fdma(
                    self.week, self.tow, sat, buff, fcn)

                if geph is not None:
                    self.re.rnx_gnav_body(geph, self.fh_rnxnav)

            if self.monlevel >= 2:
                print(f"[lD] time={time_:6d} svn={svn:2d} fcn{fcn:2d} " +
                      f"type={type_}")
        elif head == 'ud':  # Glonass CDMA Raw Navigation data
            prn, time_, type_, len_ = st.unpack_from('<BLBB', buff, 5)
            # type: 0 - L1, 1 - L2, 3 - L3
            sat = prn2sat(uGNSS.GLO, prn)
            msg = st.unpack_from('>'+len_*'L', buff, 12)
            b = bytes(np.array(msg, dtype='uint32'))

            if self.flg_rnxnav:
                geph = None
                if type_ == 0:  # L1OC
                    geph = self.rn.decode_glo_l1oc(self.week, self.tow, sat, b)
                elif type_ == 1:  # L2CSI
                    None
                elif type_ == 2:  # L3OC
                    geph = self.rn.decode_glo_l3oc(self.week, self.tow, sat, b)

                if geph is not None:
                    self.re.rnx_gnav_body(geph, self.fh_rnxnav)

            if self.monlevel >= 2:
                print(f"[ud] time={time_:6d} prn={prn:2d} type={type_}")

        elif head == 'ED':  # Galileo Raw Navigation data
            prn, time_, type_, len_ = st.unpack_from('<BLBB', buff, 5)
            sat = prn2sat(uGNSS.GAL, prn)
            if self.monlevel >= 2:
                print(f"[ED] time={time_:6d} prn={prn:2d} type={type_}")

            # [I/NAV]
            # even/odd,page-type,data(1/2),tail => 1,1,112,6
            # evan/odd,page-type,data(2/2),field1,crc,field2,tail
            # => 1,1,16,64,24,8,6
            # For E1B: field1(64) = OSNMA(40)+SAR(22)+spare(2), field2(8)=SSP
            # For E5B: field1(64) = resv, field2(8) = resv

            # [FNAV]
            # page-type(6), nav(208), crc(24), tail(6)

            # type_ = 0:E1B(INAV), 1:E5a(FNAV), 2:E5b(INAV), 6:E6(CNAV)
            if type_ == 0 or type_ == 2:  # INAV
                b = buff[12:]
                if self.flg_rnxnav:
                    eph = self.rn.decode_gal_inav(
                        self.week, time_, sat, type_, b)
                    if eph is not None:
                        self.re.rnx_nav_body(eph, self.fh_rnxnav)
                if self.flg_galinav and self.week >= 0:
                    self.fh_galinav.write(
                        "{:4d}\t{:6d}\t{:3d}\t{:1d}\t{:3d}\t{:s}\n"
                        .format(self.week, time_, prn, type_, len_,
                                hexlify(b).decode()))
            elif type_ == 1:  # FNAV
                b = buff[12:]
                if self.flg_rnxnav:
                    eph = self.rn.decode_gal_fnav(
                        self.week, time_, sat, type_, b)
                    if eph is not None:
                        self.re.rnx_nav_body(eph, self.fh_rnxnav)
                if self.flg_galfnav and self.week >= 0:
                    self.fh_galfnav.write(
                        "{:4d}\t{:6d}\t{:3d}\t{:1d}\t{:3d}\t{:s}\n"
                        .format(self.week, time_, prn, type_, len_,
                                hexlify(b).decode()))
            elif type_ == 6:  # CNAV
                if self.flg_gale6 and self.week >= 0:
                    self.fh_gale6.write(
                        "{:4d}\t{:6d}\t{:3d}\t{:1d}\t{:3d}\t{:s}\n"
                        .format(self.week, time_, prn, type_, len_,
                                hexlify(buff[12:]).decode()))

        elif head == 'WD':  # SBAS Navigation data
            prn, time_, type_, len_ = st.unpack_from('<BLBB', buff, 5)
            if self.monlevel >= 2:
                print("[WD] prn={:d} tow={:d} type={:d}".
                      format(prn, time_, type_))

            if self.flg_sbas and self.week >= 0:
                if sbs_ref > 0 and prn != sbs_ref:
                    return
                self.fh_sbas.write("{:4d}\t{:6d}\t{:3d}\t{:1d}\t{:3d}\t{:s}\n".
                                   format(self.week, time_, prn, type_, len_,
                                          hexlify(buff[12:]).decode()))

        elif head[0] == 'r' and head[1] in self.ch_t.keys():
            # Integer Pseudo-ranges
            ch = self.ch_t[head[1]]
            nsat = (len_-6)//4
            Ksys_t = {GNSS.GPS: 1e-11, GNSS.GLO: 1e-11, GNSS.GAL: 2e-11,
                      GNSS.SBS: 1e-11, GNSS.QZS: 2e-11, GNSS.BDS: 2e-11,
                      GNSS.IRN: 2e-11}
            Asys_t = {GNSS.GPS: 0.075, GNSS.GLO: 0.075, GNSS.GAL: 0.085,
                      GNSS.SBS: 0.125, GNSS.QZS: 0.125, GNSS.BDS: 0.105,
                      GNSS.IRN: 0.105}

            if self.sys == []:
                return
            spr = st.unpack_from('<'+'i'*nsat, buff, 5)
            for k in range(nsat):
                Ksys = Ksys_t[self.sys[k]]
                Asys = Asys_t[self.sys[k]]
                pr_ = (spr[k]*Ksys+Asys)*rCST.CLIGHT
                if head[1] == 'x' or head[1] == 'c':
                    self.pr_ref[k] = pr_
                else:
                    self.pr[k, ch] = pr_
                if self.monlevel >= 2:
                    print("{:2s} {:d} {:3d} {:14.3f}".
                          format(head, self.sys[k], self.prn[k], pr_))

        elif head[1] == 'p' and head[0] in self.ch_t.keys():
            # Integer Carrier-Phase
            ch = self.ch_t[head[0]]
            nsat = (len_-6)//4
            if self.sys == []:
                return
            rcp = st.unpack_from('<'+'i'*nsat, buff, 5)
            for k in range(nsat):
                if rcp[k] == 0x7fffffff or self.pr_ref[k] == 0:
                    continue
                ref = self.pr_ref[k]/rCST.CLIGHT
                freq, code = self.tofreq(head[0], self.sys[k])
                fn = self.freq_sys(self.sys[k], freq, self.freqn[k])
                # self.fn_t[ch]
                self.cp[k, ch] = (rcp[k]*rCST.P2_40+ref)*fn
                self.code[k, ch] = code
                if self.sys[k] == GNSS.GPS and self.prn[k] == 24:
                    sys = self.sys[k]
                # if self.sys[k] == GNSS.QZS:
                #    cp = self.cp[k, ch]
                if self.monlevel >= 2:
                    print("{:2s} {:d} {:3d} {:14.3f}".
                          format(head, self.sys[k], self.prn[k],
                                 self.cp[k, ch]))

        elif head[1] == 'r' and head[0] in self.ch_t.keys():
            # Integer Relative Pseudo-ranges
            ch = self.ch_t[head[0]]
            nsat = (len_-6)//2
            if self.sys == []:
                return
            srpr = st.unpack_from('<'+'h'*nsat, buff, 5)
            for k in range(nsat):
                if srpr[k] == 0x7fff or self.pr_ref[k] == 0:
                    continue
                ref = self.pr_ref[k]
                self.pr[k, ch] = (1e-11*srpr[k]+2e-7)*rCST.CLIGHT+ref
                # if self.sys[k] == GNSS.QZS:
                #    pr = self.pr[k, ch]
                if self.monlevel >= 2:
                    print("{:2s} {:d} {:3d} {:14.3f}".
                          format(head, self.sys[k], self.prn[k],
                                 self.pr[k, ch]))

        elif head[0] == 's' and head[1] in self.ch_t.keys():  # CNO x 256
            ch = self.ch_t[head[1]]
            nsat = (len_-6)//2
            self.CNO[0:nsat, ch] = st.unpack_from('<'+'h'*nsat, buff, 5)
            if self.monlevel >= 2:
                print("[{:2s}] nsat={:d}".format(head, nsat))

        elif head[0] == 'j' and head[1] in self.ch_t.keys():  # CNO x 256 data
            ch = self.ch_t[head[1]]
            nsat = (len_-6)//2
            self.CNO[0:nsat, ch] = st.unpack_from('<'+'h'*nsat, buff, 5)
            if self.monlevel >= 1:
                print("[{:2s}] nsat={:d}".format(head, nsat))

        elif head == 'ID':  # ionosphere delay
            nsat = (len_-6)//4
            # iono = st.unpack_from('<'+'f'*nsat, buff, 5)

        elif head[1] == 'm' and head[0] in self.ch_t.keys():
            # Pseudo-ranges correction
            ch = self.ch_t[head[0]]
            nsat = (len_-6)//2
            self.prc[0:nsat, ch] = st.unpack_from('<'+'h'*nsat, buff, 5)
            mode = st.unpack_from('b', buff, 5+nsat*2)[0]
            if self.monlevel >= 1:
                print("[{:2s}] nsat={:d} mode={:d}".format(head, nsat, mode))

        elif head[1] == 'f' and head[0] in self.ch_t.keys():  # CP correction
            ch = self.ch_t[head[0]]
            nsat = (len_-6)//2
            self.cpc[0:nsat, ch] = st.unpack_from('<'+'h'*nsat, buff, 5)
            mode = st.unpack_from('b', buff, 5+nsat*2)[0]
            if self.monlevel >= 1:
                print("[{:2s}] nsat={:d} mode={:d}".format(head, nsat, mode))

        elif head[0] == 'R' and head[1].lower() in self.ch_t.keys():  # PR
            ch = self.ch_t[head[1].lower()]
            nsat = (len_-6)//8
            self.pr[:nsat, ch] = st.unpack_from('d'*nsat, buff, 5)

        elif head[0] == 'P' and head[1].lower() in self.ch_t.keys():
            # Carrier-Phase
            ch = self.ch_t[head[1].lower()]
            nsat = (len_-6)//8
            self.cp[:nsat, ch] = st.unpack_from('d'*nsat, buff, 5)

        elif head[0] == 'c' and head[1] in self.ch_t.keys():
            # smoothing corrections
            ch = self.ch_t[head[1]]
            nsat = (len_-6)//2
            # smooth = st.unpack_from('h'*nsat, buff, 5)
        elif head[0] == 'D' and head[1].lower() in self.ch_t.keys():
            # doppler [Hz*1e-4]
            ch = self.ch_t[head[1]]
            nsat = (len_-6)//4
            dp = st.unpack_from('<'+'l'*nsat, buff, 5)
            for k in range(nsat):
                self.dp[k, ch] = dp[k]*1e-4
                if dp[k] == 2147483647:
                    self.dp[k, ch] = 0.0

        elif head[0] == 'E' and head[1].lower() in self.ch_t.keys():
            # C/N [dB-Hz]
            ch = self.ch_t[head[1]]
            nsat = (len_-6)//1
            cnr = st.unpack_from('b'*nsat, buff, 5)
            for k in range(nsat):
                if cnr[k] == -1:
                    continue
                if self.monlevel >= 2:
                    print("{:2s} {:d} {:3d} {:3d}".
                          format(head, self.sys[k], self.prn[k], cnr[k]))

        elif head[1] == 'E' and head[0].lower() in self.ch_t.keys():
            # C/Nx4 [dB-Hz]
            ch = self.ch_t[head[0]]
            nsat = (len_-6)//1
            cnr = st.unpack_from('B'*nsat, buff, 5)
            for k in range(nsat):
                if cnr[k] == 255:
                    continue
                self.CNO[k, ch] = cnr[k]*0.25
        # elif head == 'SE':  # security (skip)
        #    data = st.unpack_from('bbbbb', buff, 5)
        elif head == 'GE':  # GPS ephemeris
            self.decode_eph(buff, GNSS.GPS, len_)
        elif head == 'QE':  # QZS ephemeris
            self.decode_eph(buff, GNSS.QZS, len_)
        elif head == 'EN':  # Galileo ephemeris
            self.decode_eph(buff, GNSS.GAL, len_)
        elif head == 'CN':  # BDS ephemeris
            self.decode_eph(buff, GNSS.BDS, len_)
        elif head == 'NE':  # GLO ephemeris
            sv, frqnum, dne, tk, tb = st.unpack_from('<Bbhll', buff, 5)
        elif head == 'WE':  # SBAS ephemeris
            sbasprn, gpsprn, iod, acc, tod = st.unpack_from('<BBBBL', buff, 5)
        elif head[0] == 'F' and head[1].lower() in self.ch_t.keys():
            # signal lock loop flags
            ch = self.ch_t[head[1]]
            nsat = (len_-6)//2
            flags = st.unpack_from('<'+'H'*nsat, buff, 5)
            self.code[:nsat, ch] = flags

        elif head[1] == 'd' and head[0] in self.ch_t.keys():
            # relative doppler [Hz*1e-4]
            ch = self.ch_t[head[0]]
            nsat = (len_-6)//2
            # srdp = st.unpack_from('h'*nsat, buff, 5)
        elif head in ('ST', 'SP', 'PV', 'PG', 'IE', 'UO'):
            # [ST] Solution Time-Tag
            # [SP] Position Covariance Matrix
            # [PV] Cartesian Position and Velocity
            # [PG] Geodetic Position
            # [IE] IRNSS Ephemeris
            # [UO] GPS UTC Time Parameters
            None
        else:
            print("[{:s}] undef".format(head))
        return 0


if __name__ == "__main__":

    # locl mode
    # turn off tracking of all GPS SVs but SVs #8
    # set,/par/lock/sat/qzss,n
    # set,/par/lock/sat/qzss/193,y

    # turn off GALILEO E6 signal tracking for all SVs
    # set,/par/lock/sig/qzss/l6,n
    # set,/par/lock/sig/qzss/l6/193,y

    year = 2023

    bdir = '../data/doy223/'
    fnames = 'jav3223v.jps'

    opt = rcvOpt()
    opt.flg_qzsl6 = True
    opt.flg_gale6 = True
    opt.flg_galinav = True
    opt.flg_galfnav = True
    opt.flg_bdsb2b = True
    opt.flg_bdsb1c = True
    opt.flg_sbas = False
    opt.flg_rnxnav = True

    # prn_ref = 199
    prn_ref = -1
    sbs_ref = -1

    for f in glob(bdir+fnames):

        print("Decoding {}".format(f))
        bdir, fname = os.path.split(f)
        bdir += '/'

        prefix = bdir+fname[4:].removesuffix('.jps')+'_'
        jpsdec = jps(opt=opt, prefix=bdir+fname[4:].removesuffix('.jps')+'_')
        jpsdec.monlevel = 1

        jpsdec.re.anttype = "JAVRINGANT_DM   JVDM"
        jpsdec.re.rectype = "JAVAD DELTA-3"

        blen = os.path.getsize(bdir+fname)
        with open(bdir+fname, 'rb') as f:
            msg = f.read(blen)
            maxlen = len(msg)-5
            # maxlen = 400000
            for k in range(maxlen):
                stat = jpsdec.sync(msg, k)
                if not stat:
                    continue
                k += 1
                len_ = int(msg[k+2:k+5], 16)+5
                jpsdec.decode(msg[k:k+len_], len_)
                k += len_

        jpsdec.file_close()
