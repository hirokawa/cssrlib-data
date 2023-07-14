"""
 static test for PPP (IGS)
"""
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np

import cssrlib.gnss as gn
from cssrlib.gnss import ecef2pos, Nav
from cssrlib.gnss import time2doy, time2str, timediff, epoch2time
from cssrlib.gnss import rSigRnx
from cssrlib.gnss import sys2str
from cssrlib.peph import atxdec, searchpcv
from cssrlib.peph import peph, biasdec
from cssrlib.pppigs import rtkinit, ppppos, IT
from cssrlib.rinex import rnxdec

# Start epoch and number of epochs
#
ep = [2021, 9, 22, 6, 30, 0]

time = epoch2time(ep)
year = ep[0]
doy = int(time2doy(time))
nep = 360

navfile = '../data/SEPT2650.21P'
obsfile = '../data/SEPT265G.21O'

orbfile = '../data/COD0IGSRAP_{:4d}{:03d}0000_01D_15M_ORB.SP3'\
    .format(year, doy)

clkfile = '../data/COD0IGSRAP_{:4d}{:03d}0000_01D_30S_CLK.CLK'\
    .format(year, doy)

bsxfile = '../data/COD0IGSRAP_{:4d}{:03d}0000_01D_01D_OSB.BIA'\
    .format(year, doy)

# based on pnt solution
xyz_ref = [-3961951.1326752,  3381198.11019757,  3668916.0417232]
pos_ref = ecef2pos(xyz_ref)

# Define signals to be processed
#
sigs = [rSigRnx("GC1C"), rSigRnx("GC2W"),
        rSigRnx("GL1C"), rSigRnx("GL2W"),
        rSigRnx("GS1C"), rSigRnx("GS2W"),
        rSigRnx("EC1C"), rSigRnx("EC5Q"),
        rSigRnx("EL1C"), rSigRnx("EL5Q"),
        rSigRnx("ES1C"), rSigRnx("ES5Q")]

atxfile = '../data/igs14.atx'

rnx = rnxdec()
rnx.setSignals(sigs)

nav = Nav()
orb = peph()

# Positioning mode
# 0:static, 1:kinematic
#
nav.pmode = 1

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

# Logging level
#
nav.monlevel = 2  # TODO: enabled for testing!

# Load RINEX OBS file header
#
if rnx.decode_obsh(obsfile) >= 0:

    # Position
    #
    rr = rnx.pos
    rtkinit(nav, rnx.pos, 'test_pppigs2.log')
    pos = ecef2pos(rr)

    if 'UNKNOWN' in rnx.ant or rnx.ant.strip() == '':
        rnx.ant = "{:16s}{:4s}".format("JAVRINGANT_DM", "SCIS")

    # Get equipment information
    #
    nav.fout.write("Receiver: {}\n".format(rnx.rcv))
    nav.fout.write("Antenna : {}\n".format(rnx.ant))
    nav.fout.write("\n")

    if 'UNKNOWN' in rnx.ant or rnx.ant.strip() == "":
        nav.fout.write("ERROR: missing antenna type in RINEX OBS header!\n")

    # Set PCO/PCV information
    #
    nav.rcv_ant = searchpcv(atx.pcvr, rnx.ant,  rnx.ts)
    if nav.rcv_ant is None:
        nav.fout.write("ERROR: missing antenna type <{}> in ANTEX file!\n"
                       .format(rnx.ant))

    # nav.fout.write available signals
    #
    nav.fout.write("Available signals\n")
    for sys, sigs in rnx.sig_map.items():
        txt = "{:7s} {}\n".format(sys2str(sys),
                                  ' '.join([sig.str() for sig in sigs.values()]))
        nav.fout.write(txt)
    nav.fout.write("\n")

    nav.fout.write("Selected signals\n")

    for sys, tmp in rnx.sig_tab.items():
        txt = "{:7s} ".format(sys2str(sys))
        for _, sigs in tmp.items():
            txt += "{} ".format(' '.join([sig.str() for sig in sigs]))
        nav.fout.write(txt+"\n")
    nav.fout.write("\n")

    # Loop over number of epoch from file start
    #
    for ne in range(nep):

        # Get new epoch, exit after last epoch
        #
        obs = rnx.decode_obs()
        if obs.t.time == 0:
            break

        # Set initial epoch
        #
        if ne == 0:
            nav.t = deepcopy(obs.t)
            t0 = deepcopy(obs.t)
            t0.time = t0.time//30*30
            nav.time_p = t0

        # Call PPP module with IGS products
        #
        ppppos(nav, obs, orb, bsx)

        # Save output
        #
        t[ne] = timediff(nav.t, t0)/60

        sol = nav.xa[0:3] if nav.smode == 4 else nav.x[0:3]
        enu[ne, :] = gn.ecef2enu(pos_ref, sol-xyz_ref)

        ztd[ne] = nav.xa[IT(nav.na)] if nav.smode == 4 else nav.x[IT(nav.na)]
        smode[ne] = nav.smode

        nav.fout.write("{} {:14.4f} {:14.4f} {:14.4f} {:14.4f} {:14.4f} {:14.4f} {:2d}\n"
                       .format(time2str(obs.t),
                               sol[0], sol[1], sol[2],
                               enu[ne, 0], enu[ne, 1], enu[ne, 2],
                               smode[ne]))

    rnx.fobs.close()

fig_type = 1
ylim = 0.4

idx4 = np.where(smode == 4)[0]
idx5 = np.where(smode == 5)[0]
idx0 = np.where(smode == 0)[0]

fig = plt.figure(figsize=[7, 9])

if fig_type == 1:

    lbl_t = ['East [m]', 'North [m]', 'Up [m]']
    x_ticks = np.arange(0, nep/60+1, step=1)

    for k in range(3):
        plt.subplot(4, 1, k+1)
        plt.plot(t[idx0], enu[idx0, k], 'r.')
        plt.plot(t[idx5], enu[idx5, k], 'y.')
        plt.plot(t[idx4], enu[idx4, k], 'g.')

        plt.xticks(x_ticks)
        plt.ylabel(lbl_t[k])
        plt.grid()
        #plt.axis([0, ne, -ylim, ylim])

    plt.subplot(4, 1, 4)
    plt.plot(t[idx0], ztd[idx0]*1e2, 'r.', markersize=8, label='none')
    plt.plot(t[idx5], ztd[idx5]*1e2, 'y.', markersize=8, label='float')
    plt.plot(t[idx4], ztd[idx4]*1e2, 'g.', markersize=8, label='fix')
    plt.xticks(x_ticks)
    plt.ylabel('ZTD [cm]')
    plt.grid()
    plt.xlabel('Time [min]')
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
plotFileName = '.'.join(('test_pppigs2', plotFileFormat))

plt.savefig(plotFileName, format=plotFileFormat, bbox_inches='tight')

# plt.show()
