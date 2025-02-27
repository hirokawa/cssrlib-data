"""
 static test for stand-alone positioning
"""
import matplotlib.pyplot as plt
import matplotlib.dates as md
import numpy as np
from sys import stdout

from cssrlib.rinex import rnxdec
from cssrlib.gnss import Nav, rSigRnx, ecef2pos, ecef2enu
from cssrlib.gnss import epoch2time, time2doy, time2str, timediff
from cssrlib.pntpos import stdpos

# Start epoch and number of epochs
#
dataset = 3

if dataset == 0:  # SETP078M.21O
    ep = [2021, 3, 19, 12, 0, 0]
    xyz_ref = [-3962108.6617, 3381309.5232, 3668678.6410]
elif dataset == 1:  # SETP1890.23O
    ep = [2023, 7, 8, 4, 0, 0]
    xyz_ref = [-3962108.7063, 3381309.5703, 3668678.6690]
elif dataset == 2:  # SETP223Z.23O
    ep = [2023, 8, 11, 21, 0, 0]
    xyz_ref = [-3962108.7063, 3381309.5703, 3668678.6690]
elif dataset == 3:  # 046r_rnx.obs
    ep = [2025, 2, 15, 17, 0, 0]
    xyz_ref = [-3962108.6819, 3381309.5707, 3668678.6750]
else:
    print("ERROR: no RINEX data set selected!")
    exit(1)

time = epoch2time(ep)
year = ep[0]
doy = int(time2doy(time))

nep = 900*4

bdir = '../data/doy{:04d}-{:03d}/'.format(year, doy)

if dataset == 0:
    let = chr(ord('A')+ep[3])
    navfile = bdir+'SEPT{:03d}{}.{:02d}P'.format(doy, let, year % 2000)
    obsfile = bdir+'SEPT{:03d}{}.{:02d}O'.format(doy, let, year % 2000)
if dataset == 1:
    let = '0'
    navfile = bdir+'SEPT{:03d}{}.{:02d}P'.format(doy, let, year % 2000)
    obsfile = bdir+'SEPT{:03d}{}.{:02d}O'.format(doy, let, year % 2000)
elif dataset == 2:
    # let = 'Z' # MOSAIC-CLAS
    let = 'Y'  # PolaRX5
    navfile = '../data/brdc/' + \
        'BRD400DLR_S_{:04d}{:03d}0000_01D_MN.rnx'.format(year, doy)
    obsfile = bdir+'SEPT{:03d}{}.{:02d}O'.format(doy, let, year % 2000)
else:
    let = chr(ord('a')+ep[3])
    navfile = bdir+'{:03d}{}_rnx.nav'.format(doy, let)
    obsfile = bdir+'{:03d}{}_rnx.obs'.format(doy, let)

pos_ref = ecef2pos(xyz_ref)

# Define signals to be processed
#
gnss = "GEJ"
sigs = []
if 'G' in gnss:
    sigs.extend([rSigRnx("GC1C"), rSigRnx("GL1C"), rSigRnx("GS1C")])
if 'E' in gnss:
    sigs.extend([rSigRnx("EC1C"), rSigRnx("EL1C"), rSigRnx("ES1C")])
if 'C' in gnss:
    sigs.extend([rSigRnx("CC2I"), rSigRnx("CL2I"), rSigRnx("CS2I")])
if 'C' in gnss:
    sigs.extend([rSigRnx("JC1C"), rSigRnx("JL1C"), rSigRnx("JS1C")])

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
