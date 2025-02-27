"""
test script for Galileo OSNMA

[1] Galileo Open Service Navigation Message Authentication (OSNMA)
    Signal-in-Space Interface Control Document (SIS ICD), October, 2023.

[2] Galileo Open Service Navigation Message Authentication (OSNMA)
    Receiver Guidelines Issue 1.3, January, 2024.

Note:
    to use the package for OSNMA, the user needs to
    install the public keys provided by EUSPA.

@author Rui Hirokawa

"""
import os
import numpy as np
import cssrlib.osnma as om
from sys import exit as sys_exit
from binascii import unhexlify, hexlify
import matplotlib.pyplot as plt

tofst = -2  # time offset to synchronize tow
mt_file = 'OSNMA_MerkleTree_20240115100000_newPKID_1.xml'

if not os.path.exists('../data/pubkey/osnma/'+mt_file):
    print('please install OSNMA_MerkleTree*.xml from EUSPA.')
    sys_exit(0)

nma = om.osnma(mt_file)

nma.flg_slowmac = False

year = 2025
doy = 46
session = 'r'

file_galinav = f'../data/doy{year}-{doy:03d}/{doy:03d}{session}_galinav.txt'

dtype_ = [('tow', 'i8'), ('wn', 'i8'), ('prn', 'i8'),
          ('mt', 'i8'), ('k', 'i8'), ('nma', 'S10'),
          ('wt', 'i8'), ('nav', 'S32')]

dtype_ = [('wn', 'int'), ('tow', 'float'), ('prn', 'int'),
          ('type', 'int'), ('len', 'int'), ('nav', 'S512')]

v = np.genfromtxt(file_galinav, dtype=dtype_)

i = 0

v = v[v['type'] == 0]  # E1 only
tow = np.unique(v['tow'])
ntow = len(tow)
nsat = np.zeros((ntow, 3), dtype=int)
vstatus = np.zeros(ntow, dtype=int)

# nep = 90
# nep = 180
nep = 300
nep = 1799

for i, t in enumerate(tow[0:nep]):
    vi = v[v['tow'] == t]
    for vn in vi:
        tow_ = int(vn['tow'])+tofst
        prn = int(vn['prn'])
        nma.prn_a = prn
        msg = unhexlify(vn['nav'])  # I/NAV (120bit+120bit)
        nav, nma_b = nma.load_gal_inav(msg)
        nma.save_gal_inav(nav, prn, tow_)
        if nma_b[0] != 0:  # for connected satellite
            nma.decode(nma_b, int(vn['wn']), tow_, prn)
            nsat[i, 1] += 1

    nsat[i, 0] = len(vi)
    nsat[i, 2] = nma.nsat  # authenticated sat
    vstatus[i] = nma.status


if True:

    tmax = 300

    fig, ax = plt.subplots()
    plt.plot(tow-tow[0], nsat[:, 0], label='tracked')
    plt.plot(tow-tow[0], nsat[:, 1], label='connected')
    plt.plot(tow-tow[0], nsat[:, 2], label='authenticated')
    plt.grid()
    plt.legend()
    plt.xlim([0, tmax])
    ax.set_xticks(np.arange(0, 300, 30))
    plt.ylabel('number of satellites')
    plt.xlabel('time [s]')
    plt.savefig('osnma-{0:d}-nsat-{1:d}.png'.format(doy, tmax))
    plt.show()

    y = np.ones(ntow)
    lbl_t = ['rootkey-loaded', 'rootkey-verified', 'keychain-verified',
             'utc-verified', 'auth-position']
    fig, ax = plt.subplots()
    for k in range(5):
        idx = np.where(vstatus & (1 << k))
        plt.plot(tow[idx]-tow[0], y[idx]*(k+1), '.', label=lbl_t[k])
    plt.grid()
    ax.set_yticks(np.arange(0, 6))
    ax.set_xticks(np.arange(0, 300, 30))
    plt.legend()
    plt.ylim([0, 6])
    plt.xlim([0, tmax])
    plt.ylabel('status')
    plt.xlabel('time [s]')
    plt.savefig('osnma-{0:d}-status.png'.format(doy))
    plt.show()
