"""
sample for QZNMA-CLAS

[1] Quasi-Zenith Satellite System Interface Specification
    Signal Authentication Service (IS-QZSS-SAS-002 Draft),
    October, 2026

Note:
    to use the package for QZNMA-CLAS, the user needs to
    install the public keys provided by QSS.

@author Rui Hirokawa

"""

import os
from sys import exit as sys_exit
from binascii import unhexlify
import numpy as np
from cssrlib.cssrlib import cssr
from cssrlib.gnss import epoch2time, time2gpst, time2doy
from cssrlib.qznma import qznma
import matplotlib.pyplot as plt

dtype = [('wn', 'int'), ('tow', 'int'), ('prn', 'int'),
         ('type', 'int'), ('len', 'int'), ('nav', 'S512')]

l6_mode = 1  # 0: from receiver log, 1: from archive (.l6d/.l6e)
ep = [2023, 2, 13, 9, 0, 0]
prn_ref = -1
nep = 300
flg_plot = True

time = epoch2time(ep)
week, tow0 = time2gpst(time)
year = ep[0]
doy = int(time2doy(time))
session = chr(ord('a') + ep[3])

bdir = f'../data/doy{year}-{doy:03d}/'
if not os.path.exists(bdir):
    bdir = f'C:/work/gps/doy{year}-{doy:03d}/'

pubk_bdir = '../data/pubkey/qznma'
if not os.path.exists(os.path.join(pubk_bdir, 'CMA_001.der')):
    if os.path.exists('C:/work/gps/clas/CMA_001.der'):
        pubk_bdir = 'C:/work/gps/clas'
    else:
        print('please install CLAS public key file (CMA_001.der) from QSS.')
        sys_exit(0)

cs = cssr()
cs.week = week
cs.tow0 = (int(tow0) // 3600) * 3600
cs.pubk_bdir = pubk_bdir
cs.monlevel = 1
cs.read_griddef('../data/clas_grid.def')

qz = qznma()
qz.pubk_bdir = pubk_bdir
qz.cs = cs
cs.nma = qz

if l6_mode == 1:
    file_l6d = bdir + f'{year}{doy:03d}{session.upper()}.l6d'
    file_l6e = bdir + f'{year}{doy:03d}{session.upper()}.l6e'
    fd = open(file_l6d, 'rb')
    fe = open(file_l6e, 'rb')
    tow_ = np.arange(int(tow0), int(tow0) + nep)
else:
    navfile_n = bdir + f'{doy:03d}{session}_qzsl6.txt'
    vn = np.genfromtxt(navfile_n, dtype=dtype)
    tow_ = np.arange(int(vn['tow'][0]), int(vn['tow'][-1]) + 1)
    nep = min(nep, len(tow_))

vstatus = np.zeros((nep, 7), dtype=int)

for k in range(nep):
    tow = int(tow_[k])

    if l6_mode == 1:
        msg_d = fd.read(250)
        msg_e = fe.read(250)
        if not msg_d or not msg_e:
            nep = k
            break

        cs.decode_l6msg(msg_d, 0)
        if cs.fcnt == 5:
            cs.decode_cssr(cs.buff, 0)

        qz.decode(tow, msg_n=bytearray(msg_e))
    else:
        vid_ = vn[(vn['tow'] == tow) & (vn['type'] == 0)]
        vie_ = vn[(vn['tow'] == tow) & (vn['type'] == 1)]
        if prn_ref > 0:
            vid_ = vid_[vid_['prn'] == prn_ref]
            vie_ = vie_[vie_['prn'] == prn_ref]

        for vid in vid_:
            msg_d = unhexlify(vid['nav'])
            cs.decode_l6msg(msg_d, 0)
            if cs.fcnt == 5:
                cs.decode_cssr(bytes(cs.buff), 0)

        for vie in vie_:
            msg_e = unhexlify(vie['nav'])
            qz.decode(tow, msg_n=bytearray(msg_e))

    mac_tow = tow - 10
    if mac_tow in cs.clas_vstatus:
        for ctr, ok in cs.clas_vstatus[mac_tow].items():
            vstatus[k, ctr - 1] = 1 if ok else -1

if l6_mode == 1:
    fd.close()
    fe.close()

if flg_plot:
    fsz = 11
    t = tow_[:nep] - tow_[0]

    fig, ax = plt.subplots()
    for ctr in range(1, 8):
        idx = np.where(vstatus[:nep, ctr - 1] != 0)[0]
        if len(idx) > 0:
            plt.plot(t[idx], vstatus[idx, ctr - 1] * ctr, '.',
                     label=f'CTR{ctr}')
    plt.grid()
    plt.legend(loc='lower right', ncol=4)
    plt.xlim([0, t[-1] if len(t) > 0 else nep])
    plt.ylim([0, 8])
    ax.set_yticks(range(1, 8))
    ax.set_yticklabels([f'CTR{c}' for c in range(1, 8)])
    plt.ylabel('authenticated CTR', fontsize=fsz)
    plt.xlabel('time [s]', fontsize=fsz)
    plt.savefig(f'qznma-clas-{doy:03d}-status.png')
    plt.show()

