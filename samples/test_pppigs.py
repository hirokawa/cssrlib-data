"""
 static test for PPP (IGS)
"""
import matplotlib.pyplot as plt
import numpy as np
from os.path import expanduser, exists
import cssrlib.gnss as gn
from cssrlib.gnss import ecef2pos, Nav
from cssrlib.gnss import time2gpst, time2doy, timediff, epoch2time
from cssrlib.gnss import rSigRnx
from cssrlib.gnss import sys2str
from cssrlib.peph import atxdec, searchpcv
from cssrlib.peph import peph, biasdec
from cssrlib.pppigs import rtkinit, pppigspos
from cssrlib.rinex import rnxdec
import sys

"""
navfile = '../data/SEPT078M.21P'
obsfile = '../data/SEPT078M.21O'

# based on GSI F5 solution
xyz_ref = [-3962108.673,   3381309.574,   3668678.638]
pos_ref = ecef2pos(xyz_ref)
"""

# Start epoch and number of epochs
#
#ep = [2021, 3, 19, 12, 0, 0]
ep = [2022, 4, 1, 12, 0, 0]
time = epoch2time(ep)
year = ep[0]
doy = int(time2doy(time))
nep = 3600

# Files
#
atxfile = expanduser('~/GNSS_DAT/IGS/ANTEX/igs14.atx')

orbfile = expanduser(
    '~/GNSS_DAT/COD0IGSRAP/{:4d}/COD0IGSRAP_{:4d}{:03d}0000_01D_05M_ORB.SP3')\
    .format(year, year, doy)

clkfile = expanduser(
    '~/GNSS_DAT/COD0IGSRAP/{:4d}/COD0IGSRAP_{:4d}{:03d}0000_01D_30S_CLK.CLK')\
    .format(year, year, doy)

bsxfile = expanduser(
    '~/GNSS_DAT/COD0IGSRAP/{:4d}/COD0IGSRAP_{:4d}{:03d}0000_01D_01D_OSB.BIA')\
    .format(year, year, doy)

navfile = expanduser(
    '~/GNSS_NAV/IGS/{:4d}/BRDC00IGS_R_{:4d}{:03d}0000_01D_MN.rnx')\
    .format(year, year, doy)

obsfile = expanduser(
    '~/GNSS_OBS/VGS/HOURLY/{:4d}/{:03d}/CHOF00JPN_S_{:4d}{:03d}{:02d}{:02d}_01H_01S_MO.rnx')\
    .format(year, doy, year, doy, ep[3], ep[4])

if not exists(orbfile):
    orbfile = orbfile.replace('05M_ORB', '15M_ORB')

xyz_ref = [-3946217.1932, 3366689.4557, 3698971.7703]
pos_ref = ecef2pos(xyz_ref)

# Define signals to be processed
#
sigs = [rSigRnx("GC1C"), rSigRnx("GC2W"),
        rSigRnx("GL1C"), rSigRnx("GL2W"),
        rSigRnx("GS1C"), rSigRnx("GS2W"),
        rSigRnx("EC1X"), rSigRnx("EC5X"),
        rSigRnx("EL1X"), rSigRnx("EL5X"),
        rSigRnx("ES1X"), rSigRnx("ES5X")]

rnx = rnxdec()
rnx.setSignals(sigs)

# Set rover antenna
#
#rnx.ant = "{:16s}{:4s}".format("JAVRINGANT_DM", "SCIS")

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
atx = atxdec()
atx.readpcv(atxfile)

# Set satelite antenna PCO/PCV data
#
nav.sat_ant = atx.pcvs

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
    print("Receiver:", rnx.rcv)
    print("Antenna :", rnx.ant)
    print()

    if 'UNKNOWN' in rnx.ant or rnx.ant.strip() == "":
        print("ERROR: missing antenna type in RINEX OBS header!")
        sys.exit(-1)

    # Set PCO/PCV information
    #
    nav.rcv_ant = searchpcv(atx.pcvr, rnx.ant,  rnx.ts)
    if nav.rcv_ant is None:
        print("ERROR: missing antenna type <{}> in ANTEX file!".format(rnx.ant))
        sys.exit(-1)

    # Print available signals
    #
    print("Available signals")
    for sys, sigs in rnx.sig_map.items():
        txt = "{:7s} {}".format(sys2str(sys),
                                ' '.join([sig.str() for sig in sigs.values()]))
        print(txt)
    print()

    print("Selected signals")
    for sys, tmp in rnx.sig_tab.items():
        txt = "{:7s} ".format(sys2str(sys))
        for _, sigs in tmp.items():
            txt += "{} ".format(' '.join([sig.str() for sig in sigs]))
        print(txt)
    print()

    # Position
    #
    rr = rnx.pos
    rtkinit(nav, rnx.pos)
    pos = ecef2pos(rr)

    # TODO: disabled for testing!
    #
    nav.tidecorr = False

    # Loop over number of epoch from file start
    #
    for ne in range(nep):

        obs = rnx.decode_obs()
        week, tow = time2gpst(obs.t)

        # Set intial epoch
        #
        if ne == 0:
            t0 = nav.t = obs.t
            t0.time = t0.time//30*30
            nav.time_p = t0

        # Call PPP module with IGS products
        #
        pppigspos(nav, obs, orb, bsx)

        # Save output
        #
        t[ne] = timediff(nav.t, t0)

        sol = nav.xa[0:3] if nav.smode == 4 else nav.x[0:3]
        enu[ne, :] = gn.ecef2enu(pos_ref, sol-xyz_ref)
        idx = nav.idx_ztd
        ztd[ne] = nav.xa[idx] if nav.smode == 4 else nav.x[idx]
        smode[ne] = nav.smode

    rnx.fobs.close()

fig_type = 1
ylim = 0.4

idx4 = np.where(smode == 4)[0]
idx5 = np.where(smode == 5)[0]
idx0 = np.where(smode == 0)[0]

fig = plt.figure(figsize=[7, 9])

if fig_type == 1:

    lbl_t = ['East [m]', 'North [m]', 'Up [m]']
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
    plt.ylabel('ZTD [cm]')
    plt.grid()
    plt.xlabel('Time [s]')
    plt.legend()

elif fig_type == 2:

    ax = fig.add_subplot(111)

    #plt.plot(enu[idx0, 0], enu[idx0, 1], 'r.', label='stdpos')
    plt.plot(enu[idx5, 0], enu[idx5, 1], 'y.', label='float')
    plt.plot(enu[idx4, 0], enu[idx4, 1], 'g.', label='fix')

    plt.xlabel('easting [m]')
    plt.ylabel('Northing [m]')
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
