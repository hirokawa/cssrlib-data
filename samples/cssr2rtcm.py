"""
Compact SSR messages to RTCM 3 messages converter

[1] RTCM Standard 10403.4 Differential GNSS Services -
    Version 3 with Amendment 1, November, 2024

@author: Rui Hirokawa
"""

import numpy as np
from binascii import unhexlify

from cssrlib.rtcm import rtcme
from cssrlib.cssrlib import cssr, sCType
from cssrlib.gnss import load_config

config = load_config('config_ppprtk.yml')

# msg, msg_e = decode_msg(v, tow, prn_ref, l6_ch, prn_ref_ext, l6_ch_ext)


def decode_msg(v, tow, prn_ref, l6_ch=0, prn_ref_ext=0, l6_ch_ext=0):
    """ find valid correction message """

    msg, msg_e = None, None

    vi_ = v[v['tow'] == tow]
    vi = vi_[(vi_['type'] == l6_ch) & (vi_['prn'] == prn_ref)]
    if len(vi) > 0:
        msg = unhexlify(vi['nav'][0])

    # load regional STEC info (experimental)
    if prn_ref_ext > 0:
        vi = vi_[(vi_['type'] == l6_ch_ext) &
                 (vi_['prn'] == prn_ref_ext)]
        if len(vi) > 0:
            msg_e = unhexlify(vi['nav'][0])

    return msg, msg_e


def encode_msg(cs, re, msgtype=1057, maxlen=1024):

    k = 0
    buff = bytearray(maxlen)
    msg = bytearray(maxlen)

    re.msgtype = msgtype  # Orbit GPS
    re.udi = cs.udi
    re.datum = cs.datum
    re.nsat_n = cs.nsat_n
    re.sat_n = cs.sat_n
    re.lc = cs.lc
    re.grid = cs.grid
    i = re.encode(buff)

    len_ = (i+7)//8
    if k+len_+6 >= maxlen:
        return -1

    re.set_body(msg, buff, k, len_)
    re.set_sync(msg, k)
    re.set_len(msg, k, len_)
    re.set_checksum(msg, k)
    k += re.dlen

    return msg[:k]


if __name__ == "__main__":

    file_l6 = '../data/doy2025-233/233h_qzsl6.txt'
    file_rtcm = '../data/doy2025-233/233h_qzsl6.rtc'

    dtype = [('wn', 'int'), ('tow', 'int'), ('prn', 'int'),
             ('type', 'int'), ('len', 'int'), ('nav', 'S500')]
    v = np.genfromtxt(file_l6, dtype=dtype)

    griddef = config['griddef']

    nep = 10

    cs = cssr()
    cs.monlevel = 2
    cs.read_griddef(griddef)

    fc = open(file_rtcm, 'wb')
    if not fc:
        print("RTCM message file cannot open.")

    re = rtcme()
    re.gid = 1  # TBD
    re.gtype = 0
    re.ofst = 0
    re.nm = 0  # number of metadata model

    prn_ref = 199
    l6_ch = 0  # L6D
    tow = 370800-1

    maxlen = len(cs.buff)

    # re.encode(cs.buff)
    k = 0

    for ne in range(nep):
        tow += 1
        msg_, _ = decode_msg(v, tow, prn_ref, l6_ch)
        re.tow = tow
        if msg_ is not None:
            cs.decode_l6msg(msg_, 0)
            if cs.fcnt == 5:  # end of sub-frame
                cs.decode_cssr(bytes(cs.buff), 0)
                k = 0

        if ne == 0:
            msg = encode_msg(cs, re, 60)  # metadata
            fc.write(msg)
            msg = encode_msg(cs, re, 61)  # grid
            fc.write(msg)

        if cs.lc[0].cstat & (1 << sCType.MASK) == 0:
            continue

        if cs.lc[0].cstat & (1 << sCType.CLOCK):
            msg = encode_msg(cs, re, 1058)  # GPS
            fc.write(msg)
            msg = encode_msg(cs, re,   65)  # GAL
            fc.write(msg)
            msg = encode_msg(cs, re,   67)  # QZS
            fc.write(msg)
            cs.lc[0].cstat ^= (1 << sCType.CLOCK)

        if cs.lc[0].cstat & (1 << sCType.ORBIT):
            msg = encode_msg(cs, re, 1057)  # GPS
            fc.write(msg)
            msg = encode_msg(cs, re,   62)  # GAL
            fc.write(msg)
            msg = encode_msg(cs, re,   64)  # QZS
            fc.write(msg)
            cs.lc[0].cstat ^= (1 << sCType.ORBIT)

        if cs.lc[0].cstat & (1 << sCType.URA):
            msg = encode_msg(cs, re, 1061)  # GPS
            fc.write(msg)
            msg = encode_msg(cs, re,   74)  # GAL
            fc.write(msg)
            msg = encode_msg(cs, re,   76)  # QZS
            fc.write(msg)
            cs.lc[0].cstat ^= (1 << sCType.URA)

        if cs.lc[0].cstat & (1 << sCType.CBIAS):
            msg = encode_msg(cs, re, 1059)  # GPS
            fc.write(msg)
            msg = encode_msg(cs, re,   68)  # GAL
            fc.write(msg)
            msg = encode_msg(cs, re,   70)  # QZS
            fc.write(msg)
            cs.lc[0].cstat ^= (1 << sCType.CBIAS)

        if cs.lc[0].cstat & (1 << sCType.PBIAS):
            msg = encode_msg(cs, re,   85)  # GPS
            fc.write(msg)
            msg = encode_msg(cs, re,   87)  # GAL
            fc.write(msg)
            msg = encode_msg(cs, re,   89)  # QZS
            fc.write(msg)
            cs.lc[0].cstat ^= (1 << sCType.PBIAS)

    fc.close()
