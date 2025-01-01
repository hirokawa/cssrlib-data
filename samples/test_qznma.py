"""
sample for QZNMA

@author Rui Hirokawa

"""

from binascii import unhexlify
import numpy as np
from cssrlib.qznma import qznma
import matplotlib.pyplot as plt

prn_ref = 199
navmode = 4  # 1:LNAV, 2:CNAV, 3:CNAV2, 4:L6
doy = 305

qz = qznma(prn_ref=prn_ref)

if navmode == 1:
    navfile = '../data/doy2024-305/305a_qzslnav.txt'
elif navmode == 2:
    navfile = '../data/doy2024-305/305a_qzscnav.txt'
elif navmode == 3:
    navfile = '../data/doy2024-305/305a_qzscnav2.txt'
elif navmode == 4:
    navfile = '../data/doy2024-305/305a_qzsl6.txt'
    navfile_gpslnav = '../data/doy2024-305/305a_gpslnav.txt'
    navfile_gpscnav = '../data/doy2024-305/305a_gpscnav.txt'
    # navfile_gpscnav2 = '../data/doy2024-305/305a_gpscnav2.txt'
    navfile_galinav = '../data/doy2024-305/305a_galinav.txt'
    navfile_galfnav = '../data/doy2024-305/305a_galfnav.txt'

    # load navigation message
    qz.load_navmsg_lnav(navfile_gpslnav)
    qz.load_navmsg_cnav(navfile_gpscnav)
    qz.load_navmsg_inav(navfile_galinav)
    qz.load_navmsg_fnav(navfile_galfnav)

dtype = [('wn', 'int'), ('tow', 'float'), ('prn', 'int'),
         ('type', 'int'), ('len', 'int'), ('nav', 'S512')]
v = np.genfromtxt(navfile, dtype=dtype)

tow_ = np.unique(v['tow'])
nep = len(tow_)
# nep = 1200

nsat = np.zeros((nep, 2), dtype=int)
vstatus = np.zeros(nep, dtype=int)

for k in range(nep):
    if navmode != 4:
        vi = v[(v['tow'] == tow_[k]) & (v['prn'] == prn_ref)]
    else:  # L6E
        vi = v[(v['tow'] == tow_[k]) &
               (v['prn'] == prn_ref) & (v['type'] == 1)]
    if len(vi) == 0:
        qz.flag_e = 0
        continue
    msg = unhexlify(vi['nav'][0])
    prn = vi['prn'][0]
    tow = vi['tow'][0]

    qz.decode(tow, msg, navmode)

    nsat[k, 0] = qz.count_tracked_sat(tow)
    nsat[k, 1] = qz.nsat

if True:

    tmax = 900

    fig, ax = plt.subplots()
    plt.plot(tow_-tow_[0], nsat[:, 0], label='tracked')
    plt.plot(tow_-tow_[0], nsat[:, 1], label='authenticated')
    plt.grid()
    plt.legend()
    plt.xlim([0, tmax])
    # ax.set_xticks(np.arange(0, 300, 30))
    plt.ylabel('number of satellites')
    plt.xlabel('time [s]')
    # plt.savefig('qznma-{0:d}-nsat-{1:d}.png'.format(doy, tmax))
    plt.show()
