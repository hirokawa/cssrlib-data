"""
 static test for PPP-RTK (QZSS CLAS)
"""
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
from sys import exit as sys_exit
from sys import stdout

import cssrlib.gnss as gn
from cssrlib.cssrlib import cssr
from cssrlib.gnss import ecef2pos, Nav, time2gpst, timediff, time2str, time2doy
from cssrlib.gnss import rSigRnx, sys2str, epoch2time
from cssrlib.peph import atxdec, searchpcv
from cssrlib.ppprtk import ppprtkpos
from cssrlib.rinex import rnxdec
from binascii import unhexlify

l6_mode = 0  # 0: from receiver log, 1: from archive on QZSS
dataset = 2

if l6_mode == 1:

    ep = [2021, 3, 19, 12, 0, 0]
    xyz_ref = [-3962108.673, 3381309.574, 3668678.638]
    navfile = '../data/doy2021-078/SEPT078M.21P'
    obsfile = '../data/doy2021-078/SEPT078M.21O'
    l6file = '../data/doy2021-078/2021078M.l6'

else:

    if dataset == 0:

        ep = [2023, 8, 11, 21, 0, 0]
        xyz_ref = [-3962108.7007, 3381309.5532, 3668678.6648]
        navfile = '../data/doy2023-223/NAV223.23p'
        obsfile = '../data/doy2023-223/SEPT223Y.23O'  # PolaRX5
        file_l6 = '../data/doy2023-223/223v_qzsl6.txt'

    elif dataset == 2:

        ep = [2025, 2, 15, 14, 0, 0]
        xyz_ref = [-3962108.6726, 3381309.4719, 3668678.6264]

        time = epoch2time(ep)
        year = ep[0]
        doy = int(time2doy(time))
        let = chr(ord('a')+ep[3])

        bdir = '../data/doy{:04d}-{:03d}/'.format(year, doy)

        navfile = bdir+'{:03d}{}_rnx.nav'.format(doy, let)
        obsfile = bdir+'{:03d}{}_rnx.obs'.format(doy, let)  # SEPT MOSAIC-X5
        file_l6 = bdir+'{:03d}{}_qzsl6.txt'.format(doy, let)

    prn_ref = 199  # QZSS PRN
    l6_ch = 0  # 0:L6D, 1:L6E

time = epoch2time(ep)

atxfile = '../data/antex/igs14.atx'
griddef = '../data/clas_grid.def'

pos_ref = ecef2pos(xyz_ref)

cs = cssr()
cs.monlevel = 1
cs.week = time2gpst(time)[0]
cs.read_griddef(griddef)

# Define signals to be processed
#
gnss = "GEJ"  # "GEJ"
sigs = []
if 'G' in gnss:
    sigs.extend([rSigRnx("GC1C"), rSigRnx("GC2W"),
                 rSigRnx("GL1C"), rSigRnx("GL2W"),
                 rSigRnx("GS1C"), rSigRnx("GS2W")])
if 'E' in gnss:
    sigs.extend([rSigRnx("EC1C"), rSigRnx("EC5Q"),
                 rSigRnx("EL1C"), rSigRnx("EL5Q"),
                 rSigRnx("ES1C"), rSigRnx("ES5Q")])
if 'J' in gnss:
    sigs.extend([rSigRnx("JC1C"), rSigRnx("JC2L"),
                 rSigRnx("JL1C"), rSigRnx("JL2L"),
                 rSigRnx("JS1C"), rSigRnx("JS2L")])

rnx = rnxdec()
rnx.setSignals(sigs)

nav = Nav()
nav = rnx.decode_nav(navfile, nav)
nep = 900*4

# Load ANTEX data for satellites and stations
#
atx = atxdec()
atx.readpcv(atxfile)

t = np.zeros(nep)
enu = np.ones((nep, 3))*np.nan
sol = np.zeros((nep, 4))
smode = np.zeros(nep, dtype=int)

if rnx.decode_obsh(obsfile) >= 0:

    # Auto-substitute signals
    #
    rnx.autoSubstituteSignals()

    # Initialize position
    #
    ppprtk = ppprtkpos(nav, rnx.pos, 'test_ppprtk.log')

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
        txt = "{:7s} {}\n".format(sys2str(sys), ' '.
                                  join([sig.str() for sig in sigs.values()]))
        nav.fout.write(txt)
    nav.fout.write("\n")

    nav.fout.write("Selected signals\n")
    for sys, tmp in rnx.sig_tab.items():
        txt = "{:7s} ".format(sys2str(sys))
        for _, sigs in tmp.items():
            txt += "{} ".format(' '.join([sig.str() for sig in sigs]))
        nav.fout.write(txt+"\n")
    nav.fout.write("\n")

    # Get grid location
    #
    pos = ecef2pos(rnx.pos)
    inet = cs.find_grid_index(pos)

    if l6_mode == 1:
        fc = open(l6file, 'rb')
        if not fc:
            nav.fout.write("ERROR: cannot open L6 message file {}!"
                           .format(l6file))
            sys_exit(-1)
    else:
        dtype = [('wn', 'int'), ('tow', 'int'), ('prn', 'int'),
                 ('type', 'int'), ('len', 'int'), ('nav', 'S500')]
        v = np.genfromtxt(file_l6, dtype=dtype)

    # Skip epoch until start time
    #
    obs = rnx.decode_obs()
    while time > obs.t and obs.t.time != 0:
        obs = rnx.decode_obs()

    for ne in range(nep):

        week, tow = time2gpst(obs.t)

        if l6_mode == 1:
            cs.decode_l6msg(fc.read(250), 0)
            if cs.fcnt == 5:  # end of sub-frame
                cs.week = week
                cs.decode_cssr(cs.buff, 0)
        else:
            vi = v[(v['tow'] == tow) & (v['type'] == l6_ch)
                   & (v['prn'] == prn_ref)]
            if len(vi) > 0:
                cs.decode_l6msg(unhexlify(vi['nav'][0]), 0)
                if cs.fcnt == 5:  # end of sub-frame
                    cs.decode_cssr(bytes(cs.buff), 0)

        if ne == 0:
            nav.t = deepcopy(obs.t)
            t0 = deepcopy(obs.t)
            t0.time = t0.time//30*30
            cs.time = obs.t
            nav.time_p = t0

        cstat = cs.chk_stat()
        if cstat:
            ppprtk.process(obs, cs=cs)

        t[ne] = timediff(nav.t, t0)/60

        sol = nav.xa[0:3] if nav.smode == 4 else nav.x[0:3]
        enu[ne, :] = gn.ecef2enu(pos_ref, sol-xyz_ref)
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

    # Close RINEX observation and CLAS correction file
    #
    if l6_mode == 1:
        fc.close()
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

if fig_type == 1:

    lbl_t = ['East [m]', 'North [m]', 'Up [m]']
    for k in range(3):
        plt.subplot(3, 1, k+1)
        plt.plot(t[idx0], enu[idx0, k], 'r.', label='none')
        plt.plot(t[idx5], enu[idx5, k], 'y.', label='float')
        plt.plot(t[idx4], enu[idx4, k], 'g.', label='fix')

        if k == 2:
            plt.xlabel('Time [min]')
            plt.legend()
        plt.ylabel(lbl_t[k])
        plt.grid()
        plt.ylim([-ylim, ylim])

elif fig_type == 2:

    ax = fig.add_subplot(111)

    plt.plot(enu[idx0, 0], enu[idx0, 1], 'r.', label='none')
    plt.plot(enu[idx5, 0], enu[idx5, 1], 'y.', label='float')
    plt.plot(enu[idx4, 0], enu[idx4, 1], 'g.', label='fix')

    plt.xlabel('Easting [m]')
    plt.ylabel('Northing [m]')
    plt.grid()
    plt.legend()
    ylim = 0.05
    ax.set(xlim=(-ylim, ylim), ylim=(-ylim, ylim))
    plt.axis('equal')

plotFileFormat = 'eps'
plotFileName = '.'.join(('test_ppprtk', plotFileFormat))

plt.savefig(plotFileName, format=plotFileFormat, bbox_inches='tight', dpi=300)
# plt.show()
