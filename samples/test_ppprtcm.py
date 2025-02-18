"""
 static test for PPP (Galileo HAS IDD)
"""

import os
from copy import deepcopy
import matplotlib.pyplot as plt
import matplotlib.dates as md
import numpy as np
import sys
from sys import stdout

import cssrlib.gnss as gn
from cssrlib.gnss import ecef2pos, Nav
from cssrlib.gnss import time2gpst, time2doy, time2str, timediff, epoch2time
from cssrlib.gnss import rSigRnx, sys2str
from cssrlib.cssrlib import sCSSRTYPE
from cssrlib.peph import atxdec, searchpcv
from cssrlib.rtcm import rtcm
from cssrlib.pppssr import pppos
from cssrlib.rinex import rnxdec


# Select test case
#
icase = 2

# Start epoch and number of epochs
#
if icase == 1:  # Galileo HAS IDD
    ep = [2023, 8, 17, 2, 0, 0]
    navfile = '../data/doy229/OBE42023229c.nav'
    # navfile = '../data/doy229/BRD400DLR_S_20232290000_01D_MN.rnx'
    obsfile = '../data/doy229/OBE42023229c.obs'
    xyz_ref = [4186704.2262, 834903.7677, 4723664.9337]
    file_rtcm = '../data/doy229/idd2023229c.rtc'
    file_rtcm_log = '../data/doy229/idd2023229c.log'
elif icase == 2:  # JPL GDGPS  Mosaic-X5
    ep = [2024, 2, 12, 7, 0, 0]
    navfile = '../data/doy2024-043/043h_rnx.nav'
    # navfile = '../data/doy2024-043/BRD400DLR_S_20240430000_01D_MN.rnx'
    obsfile = '../data/doy2024-043/043h_rnx.obs'
    xyz_ref = [-3962108.7007, 3381309.5532, 3668678.6648]
    file_rtcm = '../data/doy2024-043/JPL32T2043h.rtcm3'
    file_rtcm_log = '../data/doy2024-043/JPL32T2043h.log'

time = epoch2time(ep)
year = ep[0]
doy = int(time2doy(time))

nep = 900*4


# Set user reference position
#
pos_ref = ecef2pos(xyz_ref)

# Define signals to be processed
#

if icase == 1:

    gnss = "GE"
    sigs = []
    if 'G' in gnss:
        sigs.extend([rSigRnx("GC1C"), rSigRnx("GC2W"),
                     rSigRnx("GL1C"), rSigRnx("GL2W"),
                     rSigRnx("GS1C"), rSigRnx("GS2W")])
    if 'E' in gnss:
        sigs.extend([rSigRnx("EC1C"), rSigRnx("EC7Q"),
                     rSigRnx("EL1C"), rSigRnx("EL7Q"),
                     rSigRnx("ES1C"), rSigRnx("ES7Q")])

elif icase == 2:

    gnss = "GE"
    sigs = []
    if 'G' in gnss:
        sigs.extend([rSigRnx("GC1C"), rSigRnx("GC2W"),
                     rSigRnx("GL1C"), rSigRnx("GL2W"),
                     rSigRnx("GS1C"), rSigRnx("GS2W")])
    if 'E' in gnss:
        sigs.extend([rSigRnx("EC1C"), rSigRnx("EC7Q"),
                     rSigRnx("EL1C"), rSigRnx("EL7Q"),
                     rSigRnx("ES1C"), rSigRnx("ES7Q")])


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

cs = rtcm(file_rtcm_log)
cs.monlevel = 1
cs.cssrmode = sCSSRTYPE.RTCM3_SSR
cs.inet = 0

if icase == 2:  # mask phase-bias for JPL GDGPS
    cs.mask_pbias = True

if True:
    fc = open(file_rtcm, 'rb')
    if not fc:
        print("RTCM message file cannot open.")

    blen = os.path.getsize(file_rtcm)
    msg = fc.read(blen)
    maxlen = len(msg)-5
    fc.close()

# Load ANTEX data for satellites and stations
#
atxfile = '../data/igs20.atx'
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
    ppp = pppos(nav, rnx.pos, 'test_ppprtcm.log')
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
    # NOTE: comment out the line with 'sys.exit(1)' to continue with zero
    #       receiver antenna corrections!
    #
    if 'UNKNOWN' in rnx.ant or rnx.ant.strip() == "":
        nav.fout.write("ERROR: missing antenna type in RINEX OBS header!\n")
        sys.exit(1)
    else:
        nav.rcv_ant = searchpcv(atx.pcvr, rnx.ant,  rnx.ts)
        if nav.rcv_ant is None:
            nav.fout.write("ERROR: missing antenna type <{}> in ANTEX file!\n"
                           .format(rnx.ant))
            sys.exit(1)

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

    k = 0
    # Loop over number of epoch from file start
    #
    for ne in range(nep):

        week, tow = time2gpst(obs.t)
        cs.week = week
        cs.tow0 = tow//3600*3600

        # Set initial epoch
        #
        if ne == 0:
            nav.t = deepcopy(obs.t)
            t0 = deepcopy(obs.t)
            t0.time = t0.time//30*30
            nav.time_p = t0

        while True:
            stat = cs.sync(msg, k)
            if stat is False:
                k += 1
                continue
            if not cs.checksum(msg, k, maxlen):
                k += 1
                continue

            tc = cs.decode_time(msg[k:k+cs.len+3])
            if (tc is not False) and timediff(tc, obs.t) > 0:
                break

            _, _, eph = cs.decode(msg[k:k+cs.len+3])
            k += cs.dlen

            if cs.msgtype in cs.eph_t.values():
                nav.eph.append(eph)

        # Call PPP module with HAS corrections
        #
        if (cs.lc[0].cstat & 0xe) == 0xe:
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
col_t = ['#d62728', '#1f77b4', '#2ca02c']  # tab:red, tab:blue, tab:green

if fig_type == 1:

    lbl_t = ['East [m]', 'North [m]', 'Up [m]']
    # nm = 4
    nm = 3

    for k in range(3):
        plt.subplot(nm, 1, k+1)
        plt.plot(t[idx0], enu[idx0, k], color=col_t[0], marker='.')
        plt.plot(t[idx5], enu[idx5, k], color=col_t[1], marker='.')
        plt.plot(t[idx4], enu[idx4, k], color=col_t[2], marker='.')

        plt.ylabel(lbl_t[k])
        plt.grid()
        plt.ylim([-ylim, ylim])
        plt.gca().xaxis.set_major_formatter(md.DateFormatter(fmt))

    if nm > 3:
        plt.subplot(nm, 1, 4)
        plt.plot(t[idx0], ztd[idx0]*1e2, color=col_t[0],
                 marker='.', markersize=8, label='none')
        plt.plot(t[idx5], ztd[idx5]*1e2, color=col_t[1],
                 marker='.',  markersize=8, label='float')
        plt.plot(t[idx4], ztd[idx4]*1e2, color=col_t[2],
                 marker='.', markersize=8, label='fix')
        plt.ylabel('ZTD [cm]')
        plt.grid()
        plt.gca().xaxis.set_major_formatter(md.DateFormatter(fmt))

    plt.xlabel('Time [HH:MM]')
    plt.legend()

elif fig_type == 2:

    ax = fig.add_subplot(111)

    plt.plot(enu[idx0, 0], enu[idx0, 1],
             color=col_t[0], marker='.', label='none')
    plt.plot(enu[idx5, 0], enu[idx5, 1],
             color=col_t[1], marker='.', label='float')
    plt.plot(enu[idx4, 0], enu[idx4, 1],
             color=col_t[2], marker='.', label='fix')

    plt.xlabel('Easting [m]')
    plt.ylabel('Northing [m]')
    plt.grid()
    plt.axis('equal')
    plt.legend()
    # ax.set(xlim=(-ylim, ylim), ylim=(-ylim, ylim))

plotFileFormat = 'png'
plotFileName = '.'.join(('test_ppprtcm', plotFileFormat))

plt.savefig(plotFileName, format=plotFileFormat, bbox_inches='tight', dpi=300)
# plt.show()
