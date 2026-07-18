"""
Compact SSR messages to RTCM 3 messages converter

[1] RTCM Standard 10403.4 Differential GNSS Services -
    Version 3 with Amendment 1, November, 2024

@author: Rui Hirokawa
"""

import argparse
from glob import glob
import multiprocessing as mp
import numpy as np
from binascii import unhexlify

from cssrlib.rtcm import rtcme
from cssrlib.cssrlib import cssr, sCType
from cssrlib.gnss import load_config, char2sys

config = load_config('config_ppprtk.yml')


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


def encode_msg(cs, re, msgtype, maxlen=1024):
    """ encode RTCM message """

    k = 0
    buff = bytearray(maxlen)
    msg = bytearray(maxlen)

    re.msgtype = msgtype
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


def out_ssr(cs, re, fc, sct, sys_t=None, inet=0):
    """ output SSR messages """

    if cs.lc[inet].cstat & (1 << sct) == 0:
        if sct not in [sCType.META, sCType.GRID]:
            return

    if sys_t is None:
        mt = re.sct2mt(sct)
        msg = encode_msg(cs, re, mt)
        fc.write(msg)
    else:
        for sys in sys_t:
            mt = re.sct2mt(sct, sys)
            if mt < 0:
                return
            msg = encode_msg(cs, re, mt)
            fc.write(msg)
    cs.lc[inet].cstat ^= (1 << sct)


def process(infile, outfile, args):
    dtype = [('wn', 'int'), ('tow', 'int'), ('prn', 'int'),
             ('type', 'int'), ('len', 'int'), ('nav', 'S500')]
    v = np.genfromtxt(args.inpFileName, dtype=dtype)

    griddef = config['griddef']

    cs = cssr()
    cs.monlevel = config['cs']['monlevel']
    cs.read_griddef(griddef)

    fc = open(outfile, 'wb')
    if not fc:
        print("RTCM message file cannot open.")

    re = rtcme()
    re.gtype = 1

    re.gid = args.gid
    prn_ref = args.prnref
    l6_ch = args.l6ch  # L6D
    tow = v[0]['tow']-1
    nep = 3600

    re.inet = re.gid
    # maxlen = len(cs.buff)

    sys_t = char2sys(args.gnss)

    for ne in range(nep):
        tow += 1
        msg_, _ = decode_msg(v, tow, prn_ref, l6_ch)
        re.tow = tow
        if msg_ is not None:
            cs.decode_l6msg(msg_, 0)
            if cs.fcnt == 5:  # end of sub-frame
                cs.decode_cssr(bytes(cs.buff), 0)

        if ne == 0:
            out_ssr(cs, re, fc, sCType.META)
            out_ssr(cs, re, fc, sCType.GRID)

        if cs.lc[0].cstat & (1 << sCType.MASK) == 0:
            continue

        out_ssr(cs, re, fc, sCType.CLOCK, sys_t)
        out_ssr(cs, re, fc, sCType.ORBIT, sys_t)
        out_ssr(cs, re, fc, sCType.URA, sys_t)
        out_ssr(cs, re, fc, sCType.CBIAS, sys_t)
        out_ssr(cs, re, fc, sCType.PBIAS, sys_t, inet=re.inet)
        out_ssr(cs, re, fc, sCType.TROP, inet=re.inet)
        out_ssr(cs, re, fc, sCType.STEC, sys_t, inet=re.inet)

    fc.close()


if __name__ == "__main__":

    # Parse command line arguments
    #
    parser = argparse.ArgumentParser(
        description="QZS L6 (CSSR) to RTCM SSR converter")

    parser.add_argument("inpFileName",
                        help="Input QZS L6 file(s) (wildcards allowed)")
    parser.add_argument("-g", "--gnss", default='GEJ',
                        help="GNSS [GEJ]")
    parser.add_argument("--prnref", type=int, default=199,
                        help="QZS satellite PRN [193-210] (default:199)")
    parser.add_argument("--l6ch", type=int, default=0,
                        help="QZS satellite L6 channel [0|1] (default:0)")
    parser.add_argument("--gid", type=int, default=7,
                        help="Network ID [1-12] (default:7)")
    parser.add_argument("-j", "--jobs", default=int(mp.cpu_count() / 2),
                        type=int, help='Max. number of parallel processes')

    args = parser.parse_args()

    # args.inpFileName = '../data/doy2025-233/233h_qzsl6.txt'
    foutname = args.inpFileName.removesuffix('.txt')+'.rtcm3'

    # process(args.inpFileName, foutname, args)
    # Start processing pool
    #
    with mp.Pool(processes=args.jobs) as pool:
        pool.starmap(process, [(f, foutname, args)
                     for f in glob(args.inpFileName)])
