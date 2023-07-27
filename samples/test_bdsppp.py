"""
BDS PPP decoding using the live data set

[1] BeiDou Navigation Satellite System Signal In Space
Interface Control Document Precise Point Positioning Service Signal PPP-B2b
(Version 1.0) , July 2020

"""
from binascii import unhexlify
import numpy as np

from cssrlib.cssr_bds import cssr_bds

# receiver log
file_has = '../data/bdsb2b.txt'

dtype = [('wn', 'int'), ('tow', 'int'), ('prn', 'int'),
         ('type', 'int'), ('len', 'int'), ('nav', 'S124')]

v = np.genfromtxt(file_has, dtype=dtype)

i = 0
tow = np.unique(v['tow'])

dec = cssr_bds()
dec.mon_level = 2

prn_a = 60

for i, t in enumerate(tow):
    vi = v[v['tow'] == t]
    for vn in vi:
        prn = int(vn['prn'])
        buff = unhexlify(vn['nav'])
        if prn == prn_a:
            dec.decode_cssr(buff, 0)
