"""
 static test for PPP (IGS)
"""
import matplotlib.pyplot as plt
import numpy as np
from os.path import expanduser
import sys

import cssrlib.gnss as gn
from cssrlib.cssrlib import cssr
from cssrlib.gnss import ecef2pos, Nav, time2gpst, time2doy, timediff, epoch2time
from cssrlib.peph import peph, readpcv, searchpcv, biasdec
from cssrlib.pppigs import rtkinit, pppigspos
from cssrlib.rinex import rnxdec

l6file = '../data/2021078M.l6'
griddef = '../data/clas_grid.def'
"""
navfile = '../data/SEPT078M.21P'
obsfile = '../data/SEPT078M.21O'

# based on GSI F5 solution
xyz_ref = [-3962108.673,   3381309.574,   3668678.638]
pos_ref = ecef2pos(xyz_ref)
"""

# Start epoch and number of epochs
#
ep = [2021, 3, 19, 12, 0, 0]
time = epoch2time(ep)
doy = int(time2doy(time))
nep = 720  # 600

# Files
#
atxfile = expanduser('~/GNSS_DAT/IGS/ANTEX/igs14.atx')

orbfile = expanduser(
    '~/GNSS_DAT/COD0IGSRAP/{:4d}/COD0IGSRAP_{:4d}{:03d}0000_01D_15M_ORB.SP3')\
    .format(ep[0], ep[0], doy)
clkfile = expanduser(
    '~/GNSS_DAT/COD0IGSRAP/{:4d}/COD0IGSRAP_{:4d}{:03d}0000_01D_30S_CLK.CLK')\
    .format(ep[0], ep[0], doy)

bsxfile = expanduser('~/GNSS_DAT/COD0IGSRAP/{:4d}/COD0IGSRAP_{:4d}{:03d}0000_01D_01D_OSB.BIA')\
    .format(ep[0], ep[0], doy)

navfile = expanduser('~/GNSS_NAV/IGS/{:4d}/BRDC00IGS_R_{:4d}{:03d}0000_01D_MN.rnx')\
    .format(ep[0], ep[0], doy)

obsfile = expanduser('~/GNSS_OBS/IGS/HIGHRATE/{:4d}/{:03d}/CHOF00JPN_S_{:4d}{:03d}{:02d}{:02d}_15M_01S_MO.rnx')\
    .format(ep[0], doy, ep[0], doy, ep[3], ep[4])


xyz_ref = [-3946217.2224, 3366689.3786, 3698971.7536]
pos_ref = ecef2pos(xyz_ref)

cs = cssr()
cs.monlevel = 1
cs.week = 2149
cs.read_griddef(griddef)

rnx = rnxdec()
nav = Nav()
orb = peph()

# Decode RINEX NAV data
#
nav = rnx.decode_nav(navfile, nav)

# Load precise orbits and clock offsets
#
nav = orb.parse_sp3(orbfile, nav)
nav = rnx.decode_clk(clkfile, nav)

# Load code and phase biases from Bias-SINEX
#
bsx = biasdec()
bsx.parse(bsxfile)

# Load ANTEX data for satellites and stations
#
nav.pcvs, pcvr = readpcv(atxfile)

# Intialize data structures for results
#
t = np.zeros(nep)
tc = np.zeros(nep)
enu = np.ones((nep, 3))*np.nan
sol = np.zeros((nep, 4))
dop = np.zeros((nep, 4))
ztd = np.zeros((nep, 1))
smode = np.zeros(nep, dtype=int)

# Load RINEX OBS file header
#
if rnx.decode_obsh(obsfile) >= 0:

    # Get equipment information
    #
    rcv = rnx.rcv
    ant = rnx.ant
    print("Receiver:", rcv)
    print("Antenna :", ant)

    if 'UNKNOWN' in ant or ant.strip() == "":
        print("ERROR: missing antenna type in RINEX OBS header!")
        sys.exit(-1)

    # Set PCO/PCV information
    #
    pcv = searchpcv(ant, time, pcvr)
    if pcv is None:
        print("ERROR: missing antenna type <{}> in ANTEX file!".format(ant))
        sys.exit(-1)

    # TODO: set the PCO,PCV for receiver antenna properly!!
    #nav.ant_pco = pcv.off

    # Position
    #
    rr = rnx.pos
    rtkinit(nav, rnx.pos)
    pos = ecef2pos(rr)

    inet = cs.find_grid_index(pos)

    fc = open(l6file, 'rb')
    if not fc:
        print("L6 messsage file cannot open.")
        sys.exit(-1)

    for ne in range(nep):

        # print(ne)

        obs = rnx.decode_obs()
        week, tow = time2gpst(obs.t)

        cs.decode_l6msg(fc.read(250), 0)
        if cs.fcnt == 5:  # end of sub-frame
            cs.week = week
            cs.decode_cssr(cs.buff, 0)

        if ne == 0:
            t0 = nav.t = obs.t
            t0.time = t0.time//30*30
            cs.time = obs.t
            nav.time_p = t0

        cstat = cs.chk_stat()
        if cstat:
            pppigspos(nav, obs, orb, bsx, cs)

        t[ne] = timediff(nav.t, t0)
        tc[ne] = timediff(cs.time, t0)

        sol = nav.xa[0:3] if nav.smode == 4 else nav.x[0:3]
        enu[ne, :] = gn.ecef2enu(pos_ref, sol-xyz_ref)
        idx = nav.idx_ztd
        ztd[ne] = nav.xa[idx] if nav.smode == 4 else nav.x[idx]
        smode[ne] = nav.smode

    fc.close()
    rnx.fobs.close()

fig_type = 1
ylim = 0.4

idx4 = np.where(smode == 4)[0]
idx5 = np.where(smode == 5)[0]
idx0 = np.where(smode == 0)[0]

fig = plt.figure(figsize=[7, 9])

if fig_type == 1:

    lbl_t = ['east[m]', 'north[m]', 'up[m]']
    for k in range(3):
        plt.subplot(4, 1, k+1)
        plt.plot(t[idx0], enu[idx0, k], 'r.')
        plt.plot(t[idx5], enu[idx5, k], 'y.')
        plt.plot(t[idx4], enu[idx4, k], 'g.')

        plt.xticks(np.arange(0, nep+1, step=30))
        plt.ylabel(lbl_t[k])
        plt.grid()
        #plt.axis([0, ne, -ylim, ylim])

    plt.subplot(4, 1, 4)
    plt.plot(t[idx0], ztd[idx0]*1e2, 'r.', markersize=8, label='none')
    plt.plot(t[idx5], ztd[idx5]*1e2, 'y.', markersize=8, label='float')
    plt.plot(t[idx4], ztd[idx4]*1e2, 'g.', markersize=8, label='fix')
    plt.xticks(np.arange(0, nep+1, step=30))
    plt.ylabel('ztd [cm]')
    plt.grid()
    plt.xlabel('time[s]')
    plt.legend()

elif fig_type == 2:

    ax = fig.add_subplot(111)

    #plt.plot(enu[idx0, 0], enu[idx0, 1], 'r.', label='stdpos')
    plt.plot(enu[idx5, 0], enu[idx5, 1], 'y.', label='float')
    plt.plot(enu[idx4, 0], enu[idx4, 1], 'g.', label='fix')

    plt.xlabel('easting [m]')
    plt.ylabel('northing [m]')
    plt.grid()
    plt.axis('equal')
    plt.legend()
    #ax.set(xlim=(-ylim, ylim), ylim=(-ylim, ylim))


plotFileFormat = 'eps'
plotFileName = '.'.join(('test_pppigs', plotFileFormat))

plt.savefig(plotFileName, format=plotFileFormat, bbox_inches='tight')

# plt.show()

if nav.fout is not None:
    nav.fout.close()
