# -*- coding: utf-8 -*-
"""
Galileo HAS decoding using the live data set

"""

import numpy as np
import bitstruct as bs
from binascii import unhexlify

from cssrlib.cssr_has import cssr_has

# receiver log
dataset = 0

if dataset == 0:

    file_has = '../data/gale6.txt'
    dtype = [('wn', 'int'), ('tow', 'int'), ('prn', 'int'),
             ('type', 'int'), ('len', 'int'), ('nav', 'S124')]
    v = np.genfromtxt(file_has, dtype=dtype)

else:

    file_has = '../data/SBF_GALRawCNAV.txt'
    dtype = [('tow', 'float64'), ('wn', 'int'),  ('prn', 'S3'), ('validity', 'str'),
             ('num1', 'int'), ('signal', 'str'), ('num2', 'int'), ('num3', 'int'),
             ('nav', 'S144')]
    v = np.genfromtxt(file_has, dtype=dtype, delimiter=',')
    v = v[v['validity'] == b'Passed']

    # Eliminate whitespace
    for i, nav in enumerate(v['nav']):
        v[i]['nav'] = (b''.join(nav.split()))

i = 0
tow = np.unique(v['tow'])


mid_ = -1
mid_decoded = []
rec = []
has_pages = np.zeros((255, 53), dtype=int)

# available from Galileo HAS ICD
file_gm = "Galileo-HAS-SIS-ICD_1.0_Annex_B_Reed_Solomon_Generator_Matrix.txt"
gMat = np.genfromtxt(file_gm, dtype="u1", delimiter=",")
dec = cssr_has()
dec.mon_level = 2

for i, t in enumerate(tow):
    vi = v[v['tow'] == t]
    for vn in vi:
        prn = (vn['prn'])
        buff = unhexlify(vn['nav'])

        i = 14
        if bs.unpack_from('u24', buff, i)[0] == 0xaf3bc3:
            continue
        hass, res, mt, mid, ms, pid = bs.unpack_from('u2u2u2u5u5u8', buff, i)
        ms += 1
        i += 24
        if hass >= 2:  # 0:test,1:operational,2:res,3:dnu
            continue
        print("tow={:6d} prn={} hass={:1d} mt={:1d} mid={:2d} ms={:2d} pid={:3d}"
              .format(int(t), prn, hass, mt, mid, ms, pid))

        if mid_ == -1 and mid not in mid_decoded:
            mid_ = mid
        if mid == mid_ and pid-1 not in rec:
            page = bs.unpack_from('u8'*53, buff, i)
            rec += [pid-1]
            has_pages[pid-1, :] = page

        if len(rec) >= ms:
            print("data collected mid={:2d} ms={:2d}".format(mid_, ms))
            HASmsg = dec.decode_has_page(rec, has_pages, gMat, ms)
            dec.decode_cssr(HASmsg)
            rec = []
            has_pages = np.zeros((255, 53), dtype=int)
            mid_decoded += [mid_]
            if len(mid_decoded) > 32:
                mid_decoded = mid_decoded[-32:]
            mid_ = -1
