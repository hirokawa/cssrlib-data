"""
 kinematic test for PPP-RTK (QZSS CLAS)
"""
import matplotlib.pyplot as plt
import numpy as np
import sys

from cssrlib.cssrlib import cssr
import cssrlib.rinex as rn
import cssrlib.gnss as gn
from cssrlib.gnss import rSigRnx, time2str, sys2str
from cssrlib.ppprtk import rtkinit, ppprtkpos
from cssrlib.peph import atxdec, searchpcv

ep = [2021, 9, 22, 6, 30, 0]
time = gn.epoch2time(ep)

atxfile = '../data/igs14.atx'
navfile = '../data/SEPT2650.21P'
obsfile = '../data/SEPT265G.21O'
l6file = '../data/2021265G.l6'
griddef = '../data/clas_grid.def'

xyz_ref = gn.pos2ecef([35.342058098, 139.521986657, 47.5515], True)

# Intial position guess
#
rr0 = xyz_ref  # from pntpos

cs = cssr()
cs.monlevel = 1
cs.week = 2176  # 2021/9/22
cs.read_griddef(griddef)

nep = 360
t = np.zeros(nep)
enu = np.zeros((nep, 3))
smode = np.zeros(nep, dtype=int)
# rr0 = [-3961951.1326752,  3381198.11019757,  3668916.0417232]  # from pntpos
pos_ref = gn.ecef2pos(xyz_ref)

# Load ANTEX data for satellites and stations
#
atx = atxdec()
atx.readpcv(atxfile)

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

# rover
rnx = rn.rnxdec()
rnx.setSignals(sigs)

nav = gn.Nav()
rnx.decode_nav(navfile, nav)

if rnx.decode_obsh(obsfile) >= 0:

    # Auto-substitute signals
    #
    rnx.autoSubstituteSignals()

    # Initialize position
    #
    rtkinit(nav, rnx.pos, 'test_ppprtk2.log')
    nav.armode = 3
    nav.excl_sat = [5, 58]  # [5, 58, 65]

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
    nav.rcv_ant = searchpcv(atx.pcvr, rnx.ant,  rnx.ts)
    if nav.rcv_ant is None:
        nav.fout.write("ERROR: missing antenna type <{}> in ANTEX file!\n"
                       .format(rnx.ant))

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

    pos = gn.ecef2pos(rr0)
    inet = cs.find_grid_index(pos)

    fc = open(l6file, 'rb')
    if not fc:
        nav.fout.write("ERROR: cannot open L6 messsage file {}!"
                       .format(l6file))
        sys.exit(-1)

    # t_obs 06:29:30
    fc.seek(250*(29*60+30))  # seek to 06:29:30
    if True:
        for k in range(30):  # read 30 sec
            cs.decode_l6msg(fc.read(250), 0)
            if cs.fcnt == 5:  # end of sub-frame
                cs.decode_cssr(bytes(cs.buff), 0)

    # Skip epoch until start time
    #
    obs = rnx.decode_obs()
    while time > obs.t and obs.t.time != 0:
        obs = rnx.decode_obs()

    for ne in range(nep):

        week, tow = gn.time2gpst(obs.t)

        cs.decode_l6msg(fc.read(250), 0)
        if cs.fcnt == 5:  # end of sub-frame
            cs.week = week
            cs.decode_cssr(bytes(cs.buff), 0)

        if ne == 0:
            t0 = nav.t = obs.t
            cs.time = obs.t
            nav.time_p = t0

        cstat = cs.chk_stat()

        if cstat:
            ppprtkpos(nav, obs, cs)
        t[ne] = gn.timediff(nav.t, t0)
        sol = nav.xa[0:3] if nav.smode == 4 else nav.x[0:3]
        enu[ne, :] = gn.ecef2enu(pos_ref, sol-xyz_ref)
        smode[ne] = nav.smode

        nav.fout.write("{} {:14.4f} {:14.4f} {:14.4f} {:14.4f} {:14.4f} {:14.4f} {:2d}\n"
                       .format(time2str(obs.t),
                               sol[0], sol[1], sol[2],
                               enu[ne, 0], enu[ne, 1], enu[ne, 2],
                               smode[ne]))

        # Get new epoch, exit after last epoch
        #
        obs = rnx.decode_obs()
        if obs.t.time == 0:
            break

    rnx.fobs.close()

    # Close output file
    #
    if nav.fout is not None:
        nav.fout.close()

# Plot results
#
idx4 = np.where(smode == 4)[0]
idx5 = np.where(smode == 5)[0]
idx1 = np.where(smode == 1)[0]

ms = 8
lbl_t = ['East [m]', 'North [m]', 'Up [m]']

plotFileFormat = 'eps'

fig = plt.figure(figsize=(6, 10))
fig.set_rasterized(True)

for k in range(3):

    plt.subplot(3, 1, k+1)
    plt.plot(t, enu[:, k], '-', color='gray')
    # plt.plot(t[idx1], enu[idx1, k], 'm.', label='stdpos')
    plt.plot(t[idx5], enu[idx5, k], 'y.', markersize=ms, label='float')
    plt.plot(t[idx4], enu[idx4, k], 'g.', markersize=ms, label='fix')
    plt.xticks(np.arange(0, nep+1, step=30))
    plt.ylabel(lbl_t[k])
    plt.xlabel('time[s]')
    if k == 0:
        plt.legend()
    plt.grid()

plotFileName = '.'.join(('test_ppprtk2_1', plotFileFormat))
plt.savefig(plotFileName, format=plotFileFormat, bbox_inches='tight', dpi=300)

# plt.show()

fig = plt.figure(figsize=(6, 10))
fig.set_rasterized(True)

plt.plot(enu[:, 0], enu[:, 1], '-', color='gray')
# plt.plot(enu[idx1, 0], enu[idx1, 1], 'm.', markersize=ms, label='stdpos')
plt.plot(enu[idx5, 0], enu[idx5, 1], 'y.', markersize=ms, label='float')
plt.plot(enu[idx4, 0], enu[idx4, 1], 'g.', markersize=ms, label='fix')

plt.xlabel('Easting [m]')
plt.ylabel('Northing [m]')
plt.grid()
plt.axis('equal')
plt.legend()

plotFileName = '.'.join(('test_ppprtk2_2', plotFileFormat))
plt.savefig(plotFileName, format=plotFileFormat, bbox_inches='tight', dpi=300)

# plt.show()
