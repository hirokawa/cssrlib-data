"""
 static test for stand-alone positioning
"""
import matplotlib.pyplot as plt
import matplotlib.dates as md
import numpy as np
from sys import stdout

from cssrlib.rinex import rnxdec
from cssrlib.gnss import ecef2pos, timediff, ecef2enu
from cssrlib.gnss import rSigRnx, epoch2time, Nav, time2str
from cssrlib.pntpos import stdpos

if False:
    xyz_ref = [-3962108.673,   3381309.574,   3668678.638]
    navfile = '../data/doy2021-078/SEPT078M.21P'
    obsfile = '../data/doy2021-078/SEPT078M.21O'
    ep = [2021, 3, 19, 12, 0, 0]
else:
    xyz_ref = [-3962108.6726, 3381309.4719, 3668678.6264]
    ep = [2023, 8, 11, 21, 0, 0]
    navfile = '../data/doy2023-223/BRD400DLR_S_20232230000_01D_MN.rnx'
    # navfile = '../data/doy2023-223/BRDC00IGS_R_20232230000_01D_MN.rnx'
    # navfile = '../data/doy2023-223/NAV223.23p'
    # obsfile = '../data/doy2023-223/SEPT223Z.23O'  # MOSAIC-CLAS
    obsfile = '../data/doy2023-223/SEPT223Y.23O'  # PolaRX5

pos_ref = ecef2pos(xyz_ref)
nep = 360

time = epoch2time(ep)

# Define signals to be processed
#
sigs = [rSigRnx("GC1C"), rSigRnx("EC1C"), rSigRnx("JC1C"),
        rSigRnx("GL1C"), rSigRnx("EL1C"), rSigRnx("JL1C"),
        rSigRnx("GS1C"), rSigRnx("ES1C"), rSigRnx("JS1C")]

# sigs = [rSigRnx("GC1C"), rSigRnx("GL1C"), rSigRnx("GS1C")]

rnx = rnxdec()
rnx.setSignals(sigs)


nav = Nav(nf=1)
nav.pmode = 1  # 0: static, 1: kinematic
nav = rnx.decode_nav(navfile, nav)

# cs = sbasDec('test_stdpos.log')
# cs.monlevel = 0

t = np.zeros(nep)
enu = np.zeros((nep, 3))
dop = np.zeros((nep, 4))
smode = np.zeros(nep, dtype=int)

if rnx.decode_obsh(obsfile) >= 0:

    # Auto-substitute signals
    #
    rnx.autoSubstituteSignals()

    # Initialize position
    #
    std = stdpos(nav, rnx.pos, 'test_stdpos.log')
    nav.elmin = np.deg2rad(5.0)

    sol = np.zeros((nep, nav.nx))

    # Skip epochs until start time
    #
    obs = rnx.decode_obs()
    while time > obs.t and obs.t.time != 0:
        obs = rnx.decode_obs()

    for ne in range(nep):

        if ne == 0:
            t0 = nav.t = obs.t
        t[ne] = timediff(obs.t, t0)/86400.0

        std.process(obs, cs=None)

        sol[ne, :] = nav.x
        enu[ne, :] = ecef2enu(pos_ref, sol[ne, 0:3]-xyz_ref)
        dop[ne, :] = std.dop

        smode[ne] = nav.smode

        nav.fout.write("{:s} {:14.4f} {:14.4f} {:14.4f} "
                       "ENU {:7.3f} {:7.3f} {:7.3f}, 2D {:6.3f}, mode {:1d}\n"
                       .format(time2str(obs.t),
                               sol[ne, 0], sol[ne, 1], sol[ne, 2],
                               enu[ne, 0], enu[ne, 1], enu[ne, 2],
                               np.sqrt(enu[ne, 0]**2+enu[ne, 1]**2),
                               smode[ne]))

        # Log to standard output
        #
        stdout.write('\r {:s} ENU {:7.3f} {:7.3f} {:7.3f}, 2D {:6.3f}, '
                     'mode {:1d}'
                     .format(time2str(obs.t),
                             enu[ne, 0], enu[ne, 1], enu[ne, 2],
                             np.sqrt(enu[ne, 0]**2+enu[ne, 1]**2),
                             smode[ne]))

        obs = rnx.decode_obs()
        if obs.t.time == 0:
            break

    rnx.fobs.close()


fig_type = 1
ylim_h = 2.0
ylim_v = 6.0

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
plotFileName = '.'.join(('test_stdpos', plotFileFormat))

plt.savefig(plotFileName, format=plotFileFormat, bbox_inches='tight', dpi=300)
