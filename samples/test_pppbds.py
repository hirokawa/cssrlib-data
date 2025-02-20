"""
 static test for PPP (BeiDou PPP)
"""
from binascii import unhexlify
from copy import deepcopy
import matplotlib.pyplot as plt
import matplotlib.dates as md
import numpy as np
from sys import exit as sys_exit
from sys import stdout

import cssrlib.gnss as gn
from cssrlib.gnss import ecef2pos, Nav
from cssrlib.gnss import time2gpst, time2doy, time2str, timediff, epoch2time
from cssrlib.gnss import rSigRnx
from cssrlib.gnss import sys2str
from cssrlib.peph import atxdec, searchpcv
from cssrlib.cssr_bds import cssr_bds
from cssrlib.pppssr import pppos
from cssrlib.rinex import rnxdec


# Select test case
#
dataset = 2

# Start epoch and number of epochs
#
if dataset == 0:
    ep = [2023, 7, 8, 4, 0, 0]
    navfile = '../data/doy2023-189/BRD400DLR_S_20231890000_01D_MN.rnx'
    obsfile = '../data/doy2023-189/SEPT1890.23O'
    file_bds = '../data/doy2023-189/bdsb2b_189e.txt'
elif dataset == 1:
    ep = [2023, 8, 11, 21, 0, 0]
    navfile = '../data/doy2023-223/BRD400DLR_S_20232230000_01D_MN.rnx'
    # navfile = '../data/doy2023-223/NAV223.23p'
    # obsfile = '../data/doy2023-223/SEPT223Z.23O'  # MOSAIC-CLAS
    obsfile = '../data/doy2023-223/SEPT223Y.23O'  # PolaRX5
    file_bds = '../data/doy2023-223/223v_bdsb2b.txt'
elif dataset == 2:
    ep = [2025, 2, 15, 13, 0, 0]
    navfile = '../data/doy2025-046/046n_rnx.nav'
    obsfile = '../data/doy2025-046/046n_rnx.obs'  # PolaRX5
    file_bds = '../data/doy2025-046/046n_bdsb2b.txt'
    xyz_ref = [-3962108.6726, 3381309.4719, 3668678.6264]


time = epoch2time(ep)
year = ep[0]
doy = int(time2doy(time))

nep = 900*4

dtype = [('wn', 'int'), ('tow', 'int'), ('prn', 'int'),
         ('type', 'int'), ('len', 'int'), ('nav', 'S124')]
v = np.genfromtxt(file_bds, dtype=dtype)

prn_ref = 59  # satellite PRN to receive BDS PPP collection

xyz_ref = [-3962108.7007, 3381309.5532, 3668678.6648]
pos_ref = ecef2pos(xyz_ref)

# Define signals to be processed
#
gnss = "GC"
sigs = []
if 'G' in gnss:
    sigs.extend([rSigRnx("GC1C"), rSigRnx("GC2W"),
                 rSigRnx("GL1C"), rSigRnx("GL2W"),
                 rSigRnx("GS1C"), rSigRnx("GS2W")])
if 'C' in gnss:
    sigs.extend([rSigRnx("CC1P"), rSigRnx("CC5P"),
                 rSigRnx("CL1P"), rSigRnx("CL5P"),
                 rSigRnx("CS1P"), rSigRnx("CS5P")])

rnx = rnxdec()
rnx.setSignals(sigs)

nav = Nav()

# Positioning mode
# 0:static, 1:kinematic
#
nav.pmode = 0

# Decode RINEX NAV data
#
nav = rnx.decode_nav(navfile, nav)

cs = cssr_bds()
cs.monlevel = 0
"""
cs = cssr_bds('test_pppbds_ssr.log')
cs.monlevel = 2
"""

# Load ANTEX data for satellites and stations
#

if time > epoch2time([2022, 11, 27, 0, 0, 0]):
    atxfile = '../data/antex/igs20.atx'
else:
    atxfile = '../data/antex/igs14.atx'

atx = atxdec()
atx.readpcv(atxfile)

# Initialize data structures for results
#
t = np.zeros(nep)
enu = np.ones((nep, 3))*np.nan
sol = np.zeros((nep, 4))
ztd = np.zeros((nep, 1))
smode = np.zeros(nep, dtype=int)

# Logging level
#
nav.monlevel = 1  # TODO: enabled for testing!

# Load RINEX OBS file header
#
if rnx.decode_obsh(obsfile) >= 0:

    # Auto-substitute signals
    #
    rnx.autoSubstituteSignals()

    # Initialize position
    #
    ppp = pppos(nav, rnx.pos, 'test_pppbds.log')
    nav.elmin = np.deg2rad(5.0)

    # Get equipment information
    #
    nav.fout.write("FileName: {}\n".format(obsfile))
    nav.fout.write("Start   : {}\n".format(time2str(rnx.ts)))
    if rnx.te is not None:
        nav.fout.write("End     : {}\n".format(time2str(rnx.te)))
    nav.fout.write("Receiver: {}\n".format(rnx.rcv))
    nav.fout.write("Antenna : {}\n".format(rnx.ant))
    nav.fout.write("\n")

    # Set satellite PCO/PCV information
    #
    nav.sat_ant = atx.pcvs

    # Set receiver PCO/PCV information, check antenna name and exit if unknown
    #
    # NOTE: comment out the line with 'sys_exit(1)' to continue with zero
    #       receiver antenna corrections!
    #
    if 'UNKNOWN' in rnx.ant or rnx.ant.strip() == "":
        nav.fout.write("ERROR: missing antenna type in RINEX OBS header!\n")
        sys_exit(1)
    else:
        nav.rcv_ant = searchpcv(atx.pcvr, rnx.ant,  rnx.ts)
        if nav.rcv_ant is None:
            nav.fout.write("ERROR: missing antenna type <{}> in ANTEX file!\n"
                           .format(rnx.ant))
            sys_exit(1)

    if nav.rcv_ant is None:
        nav.fout.write("WARNING: no receiver antenna corrections applied!\n")
        nav.fout.write("\n")

    # Print available signals
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

    # Skip epochs until start time
    #
    obs = rnx.decode_obs()
    while time > obs.t and obs.t.time != 0:
        obs = rnx.decode_obs()

    # Loop over number of epoch from file start
    #
    for ne in range(nep):

        week, tow = time2gpst(obs.t)
        cs.week = week
        cs.tow0 = tow//86400*86400

        # Set initial epoch
        #
        if ne == 0:
            nav.t = deepcopy(obs.t)
            t0 = deepcopy(obs.t)
            t0.time = t0.time//30*30
            nav.time_p = t0

        vi = v[(v['tow'] == tow) & (v['prn'] == prn_ref)]
        if len(vi) > 0:
            buff = unhexlify(vi['nav'][0])
            # prn, rev = bs.unpack_from('u6u6', buff, 0)
            cs.decode_cssr(buff, 0)

        # Call PPP module with BDS-PPP corrections
        #
        if (cs.lc[0].cstat & 0xf) == 0xf:
            ppp.process(obs, cs=cs)

        # Save output
        #
        t[ne] = timediff(nav.t, t0)/86400.0

        sol = nav.xa[0:3] if nav.smode == 4 else nav.x[0:3]
        enu[ne, :] = gn.ecef2enu(pos_ref, sol-xyz_ref)

        ztd[ne] = nav.xa[ppp.IT(nav.na)] \
            if nav.smode == 4 else nav.x[ppp.IT(nav.na)]
        smode[ne] = nav.smode

        nav.fout.write("{} {:14.4f} {:14.4f} {:14.4f} "
                       "ENU {:7.3f} {:7.3f} {:7.3f}, 2D {:6.3f}, mode {:1d}\n"
                       .format(time2str(obs.t),
                               sol[0], sol[1], sol[2],
                               enu[ne, 0], enu[ne, 1], enu[ne, 2],
                               np.sqrt(enu[ne, 0]**2+enu[ne, 1]**2),
                               smode[ne]))

        # Log to standard output
        #
        stdout.write('\r {} ENU {:7.3f} {:7.3f} {:7.3f}, 2D {:6.3f}, mode {:1d}'
                     .format(time2str(obs.t),
                             enu[ne, 0], enu[ne, 1], enu[ne, 2],
                             np.sqrt(enu[ne, 0]**2+enu[ne, 1]**2),
                             smode[ne]))

        # Get new epoch, exit after last epoch
        #
        obs = rnx.decode_obs()
        if obs.t.time == 0:
            break

    # Send line-break to stdout
    #
    stdout.write('\n')

    # Close RINEX observation file
    #
    rnx.fobs.close()

    # Close output file
    #
    if nav.fout is not None:
        nav.fout.close()

fig_type = 1
ylim = 1.0

idx4 = np.where(smode == 4)[0]
idx5 = np.where(smode == 5)[0]
idx0 = np.where(smode == 0)[0]

fig = plt.figure(figsize=[7, 9])
fig.set_rasterized(True)

fmt = '%H:%M'

if fig_type == 1:

    lbl_t = ['East [m]', 'North [m]', 'Up [m]']

    for k in range(3):
        plt.subplot(4, 1, k+1)
        plt.plot(t[idx0], enu[idx0, k], 'r.')
        plt.plot(t[idx5], enu[idx5, k], 'y.')
        plt.plot(t[idx4], enu[idx4, k], 'g.')

        plt.ylabel(lbl_t[k])
        plt.grid()
        plt.ylim([-ylim, ylim])
        plt.gca().xaxis.set_major_formatter(md.DateFormatter(fmt))

    plt.subplot(4, 1, 4)
    plt.plot(t[idx0], ztd[idx0]*1e2, 'r.', markersize=8, label='none')
    plt.plot(t[idx5], ztd[idx5]*1e2, 'y.', markersize=8, label='float')
    plt.plot(t[idx4], ztd[idx4]*1e2, 'g.', markersize=8, label='fix')
    plt.ylabel('ZTD [cm]')
    plt.grid()
    plt.gca().xaxis.set_major_formatter(md.DateFormatter(fmt))

    plt.xlabel('Time [HH:MM]')
    plt.legend()

elif fig_type == 2:

    ax = fig.add_subplot(111)

    plt.plot(enu[idx0, 0], enu[idx0, 1], 'r.', label='none')
    plt.plot(enu[idx5, 0], enu[idx5, 1], 'y.', label='float')
    plt.plot(enu[idx4, 0], enu[idx4, 1], 'g.', label='fix')

    plt.xlabel('Easting [m]')
    plt.ylabel('Northing [m]')
    plt.grid()
    plt.axis('equal')
    plt.legend()
    # ax.set(xlim=(-ylim, ylim), ylim=(-ylim, ylim))

plotFileFormat = 'eps'
plotFileName = '.'.join(('test_pppbds', plotFileFormat))

plt.savefig(plotFileName, format=plotFileFormat, bbox_inches='tight', dpi=300)
# plt.show()
