"""
 static test for SBAS (L1 or DFMC)
"""
from binascii import unhexlify
from copy import deepcopy
import matplotlib.pyplot as plt
import matplotlib.dates as md
import numpy as np
from sys import stdout

from cssrlib.gnss import ecef2pos, Nav, ecef2enu
from cssrlib.gnss import time2gpst, time2doy, time2str, timediff, epoch2time
from cssrlib.gnss import rSigRnx, sys2str, uIonoModel
from cssrlib.peph import atxdec, searchpcv
from cssrlib.pntpos import stdpos
from cssrlib.sbas import sbasDec
from cssrlib.rinex import rnxdec
# from cssrlib.cssr_pvs import decode_sinca_line


# Select test case
#
dataset = 1

# Start epoch and number of epochs
#
if dataset == 0:  # MSAS, L1 SBAS
    ep = [2023, 8, 11, 21, 0, 0]
    navfile = '../data/brdc/BRD400DLR_S_20232230000_01D_MN.rnx'
    obsfile = '../data/2023-doy223/SEPT223Y.23O'  # PolaRX5
    file_sbas = '../data/2023-doy223/223v_sbas.txt'
    xyz_ref = [-3962108.6726, 3381309.4719, 3668678.6264]
    prn_ref = 137  # satellite PRN for SBAS correction
    sbas_type = 0  # L1: 0, L5: 1
    nf = 1

elif dataset == 1:  # MSAS, L1 SBAS
    ep = [2025, 2, 15, 17, 0, 0]
    navfile = '../data/doy2025-046/046r_rnx.nav'
    obsfile = '../data/doy2025-046/046r_rnx.obs'  # PolaRX5
    file_sbas = '../data/doy2025-046/046r_sbas.txt'
    xyz_ref = [-3962108.6726, 3381309.4719, 3668678.6264]
    prn_ref = 137  # satellite PRN for SBAS correction
    sbas_type = 0  # L1: 0, L5: 1
    nf = 1

elif dataset == 2:  # QZSS, L5 DFMC
    ep = [2023, 8, 11, 21, 0, 0]
    navfile = '../data/brdc/BRD400DLR_S_20232230000_01D_MN.rnx'
    obsfile = '../data/2023-doy223/SEPT223Y.23O'  # PolaRX5
    file_sbas = '../data/2023-doy223/223v_sbas.txt'
    xyz_ref = [-3962108.6726, 3381309.4719, 3668678.6264]
    prn_ref = 189  # satellite PRN for SBAS correction
    sbas_type = 1  # L1: 0, L5: 1
    nf = 2

time = epoch2time(ep)
year = ep[0]
doy = int(time2doy(time))

nep = 900*4
# nep = 360


pos_ref = ecef2pos(xyz_ref)

# Define signals to be processed
#
sigs = []

if sbas_type == 0:  # single frequency L1 SBAS
    nf = 1
    gnss = "G"
    if 'G' in gnss:
        sigs.extend([rSigRnx("GC1C"), rSigRnx("GL1C"), rSigRnx("GS1C")])

else:  # dual frequency
    nf = 2
    gnss = "GE"
    if 'G' in gnss:
        sigs.extend([rSigRnx("GC1C"), rSigRnx("GC5Q"),
                     rSigRnx("GL1C"), rSigRnx("GL5Q"),
                     rSigRnx("GS1C"), rSigRnx("GS5Q")])
    if 'E' in gnss:
        sigs.extend([rSigRnx("EC1C"), rSigRnx("EC5Q"),
                     rSigRnx("EL1C"), rSigRnx("EL5Q"),
                     rSigRnx("ES1C"), rSigRnx("ES5Q")])
    if 'J' in gnss:
        sigs.extend([rSigRnx("JC1C"), rSigRnx("JC5Q"),
                     rSigRnx("JL1C"), rSigRnx("JL5Q"),
                     rSigRnx("JS1C"), rSigRnx("JS5Q")])

rnx = rnxdec()
rnx.setSignals(sigs)

nav = Nav(nf=nf)

# Positioning mode
# 0:static, 1:kinematic
#
nav.pmode = 1

# Decode RINEX NAV data
#
nav = rnx.decode_nav(navfile, nav)

cs = sbasDec('test_sbas_cs.log')
cs.monlevel = 0

# Load ANTEX data for satellites and stations
#
atx = atxdec()
atx.readpcv('../data/antex/igs20.atx')

# Initialize data structures for results
#
t = np.zeros(nep)
enu = np.ones((nep, 3))*np.nan
sol = np.zeros((nep, 4))
ztd = np.zeros((nep, 1))
smode = np.zeros(nep, dtype=int)
nsat = np.zeros(nep, dtype=int)

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
    std = stdpos(nav, rnx.pos, 'test_sbas.log')
    std.monlevel = 1
    nav.elmin = np.deg2rad(5.0)

    std.ionoModel = uIonoModel.SBAS

    nav.nf = nf
    nav.rmode = 2 if nav.nf == 2 else 0  # L1/L5 iono-free combination
    nav.csmooth = True

    # Get equipment information
    #
    nav.fout.write("FileName: {}\n".format(obsfile))
    nav.fout.write("Start   : {}\n".format(time2str(rnx.ts)))
    if rnx.te is not None:
        nav.fout.write("End     : {}\n".format(time2str(rnx.te)))
    nav.fout.write("Receiver: {}\n".format(rnx.rcv))
    nav.fout.write("Antenna : {}\n".format(rnx.ant))
    nav.fout.write("\n")

    if 'UNKNOWN' in rnx.ant or rnx.ant.strip() == "":
        nav.fout.write("ERROR: missing antenna type in RINEX OBS header!\n")

    # Set PCO/PCV information
    #
    nav.sat_ant = atx.pcvs
    nav.rcv_ant = searchpcv(atx.pcvr, rnx.ant,  rnx.ts)
    if nav.rcv_ant is None:
        nav.fout.write("ERROR: missing antenna type <{}> in ANTEX file!\n"
                       .format(rnx.ant))

    # Print available signals
    #
    nav.fout.write("Available signals\n")
    for sys, sigs in rnx.sig_map.items():
        txt = "{:7s} {}\n".format(sys2str(sys),
                                  ' '.join([s.str() for s in sigs.values()]))
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

    dtype = [('wn', 'int'), ('tow', 'float'), ('prn', 'int'),
             ('type', 'int'), ('marker', 'S2'), ('nav', 'S124')]
    v = np.genfromtxt(file_sbas, dtype=dtype)

    # Loop over number of epoch from file start
    #
    for ne in range(nep):

        week, tow = time2gpst(obs.t)
        cs.week = week
        cs.tow0 = tow//86400*86400
        cs.time = obs.t

        # Set initial epoch
        #
        if ne == 0:
            nav.t = deepcopy(obs.t)
            t0 = deepcopy(obs.t)
            t0.time = t0.time//30*30
            nav.time_p = t0

        vi = v[(v['tow'] == tow) & (v['prn'] == prn_ref)]
        if sbas_type == 0:  # L1
            vi = vi[vi['type'] <= 30]
        else:  # DFMC L5
            vi = vi[vi['type'] > 30]

        #       & (v['type'] == sbas_type)]
        if len(vi) > 0:
            buff = unhexlify(vi['nav'][0])
            cs.decode_cssr(buff, 0, src=sbas_type, prn=prn_ref)

        # cs.check_validity(obs.t)

        # Call PPP module with PVS corrections
        #
        std.process(obs, cs=cs)
        # std.process(obs)

        # Save output
        #
        t[ne] = timediff(nav.t, t0)/86400.0

        sol = nav.xa[0:3] if nav.smode == 4 else nav.x[0:3]
        enu[ne, :] = ecef2enu(pos_ref, sol-xyz_ref)

        # ztd[ne] = nav.xa[std.IT(nav.na)] \
        #    if nav.smode == 4 else nav.x[std.IT(nav.na)]
        smode[ne] = nav.smode
        # nsat[ne] = std.nsat

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
ylim_h = 10.0
ylim_v = 12.0

idx2 = np.where(smode == 2)[0]
idx1 = np.where(smode == 1)[0]
idx0 = np.where(smode == 0)[0]

fig = plt.figure(figsize=[7, 9])
fig.set_rasterized(True)

fmt = '%H:%M'

if fig_type == 1:

    lbl_t = ['East [m]', 'North [m]', 'Up [m]']

    for k in range(3):
        ylim = ylim_h if k < 2 else ylim_v

        plt.subplot(3, 1, k+1)
        plt.plot(t[idx0], enu[idx0, k], 'r.', label='none')
        plt.plot(t[idx2], enu[idx2, k], 'y.', label='SBAS/DGPS')
        plt.plot(t[idx1], enu[idx1, k], 'g.', label='standalone')

        plt.ylabel(lbl_t[k])
        plt.grid()
        plt.ylim([-ylim, ylim])
        plt.gca().xaxis.set_major_formatter(md.DateFormatter(fmt))

    plt.xlabel('Time [HH:MM]')
    plt.legend()

elif fig_type == 2:

    ax = fig.add_subplot(111)

    plt.plot(enu[idx0, 0], enu[idx0, 1], 'r.', label='none')
    plt.plot(enu[idx2, 0], enu[idx2, 1], 'y.', label='SBAS/DGPS')
    plt.plot(enu[idx1, 0], enu[idx1, 1], 'g.', label='standalone')

    plt.xlabel('Easting [m]')
    plt.ylabel('Northing [m]')
    plt.grid()
    plt.axis('equal')
    plt.legend()
    # ax.set(xlim=(-ylim, ylim), ylim=(-ylim, ylim))

plotFileFormat = 'eps'
plotFileName = '.'.join(('test_sbas', plotFileFormat))

plt.savefig(plotFileName, format=plotFileFormat, bbox_inches='tight', dpi=300)
# plt.show()
