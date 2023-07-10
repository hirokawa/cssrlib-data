"""
 kinematic test for PPP-RTK (QZSS CLAS)
"""
import matplotlib.pyplot as plt
import numpy as np
import sys

from cssrlib.cssrlib import cssr
import cssrlib.rinex as rn
import cssrlib.gnss as gn
from cssrlib.gnss import rSigRnx
from cssrlib.ppprtk import rtkinit, ppprtkpos
from cssrlib.peph import atxdec, searchpcv

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

# Load ANTEX data for satellites and stations
#
atx = atxdec()
atx.readpcv(atxfile)

# Define signals to be processed
#
sigs = [rSigRnx("GC1C"), rSigRnx("GC2W"),
        rSigRnx("EC1C"), rSigRnx("EC5Q"),
        rSigRnx("JC1C"), rSigRnx("JC2L"),
        rSigRnx("GL1C"), rSigRnx("GL2W"),
        rSigRnx("EL1C"), rSigRnx("EL5Q"),
        rSigRnx("JL1C"), rSigRnx("JL2L"),
        rSigRnx("GS1C"), rSigRnx("GS2W"),
        rSigRnx("ES1C"), rSigRnx("ES5Q"),
        rSigRnx("JS1C"), rSigRnx("JS2L")]

# rover
dec = rn.rnxdec()
dec.setSignals(sigs)
nav = gn.Nav()
dec.decode_nav(navfile, nav)
dec.decode_obsh(obsfile)


if 'UNKNOWN' in dec.ant:
    dec.ant = "{:16s}{:4s}".format("JAVRINGANT_DM", "SCIS")

# Get equipment information
#
print("Receiver:", dec.rcv)
print("Antenna :", dec.ant)

if 'UNKNOWN' in dec.ant or dec.ant.strip() == "":
    print("ERROR: missing antenna type in RINEX OBS header!")
    sys.exit(-1)

# Set PCO/PCV information
#
nav.rcv_ant = searchpcv(atx.pcvr, dec.ant,  dec.ts)
if nav.rcv_ant is None:
    print("ERROR: missing antenna type <{}> in ANTEX file!".format(dec.ant))
    sys.exit(-1)

nep = 360
t = np.zeros(nep)
enu = np.zeros((nep, 3))
smode = np.zeros(nep, dtype=int)
# rr0 = [-3961951.1326752,  3381198.11019757,  3668916.0417232]  # from pntpos
pos_ref = gn.ecef2pos(xyz_ref)

nav.excl_sat = [5, 58, 65]

if True:

    rtkinit(nav, rr0)
    nav.armode = 3
    pos = gn.ecef2pos(rr0)
    inet = cs.find_grid_index(pos)

    fc = open(l6file, 'rb')
    if not fc:
        print("L6 messsage file cannot open.")
        exit(-1)

    # t_obs 06:29:30
    fc.seek(250*(29*60+30))  # seek to 06:29:30
    if True:
        for k in range(30):  # read 30 sec
            cs.decode_l6msg(fc.read(250), 0)
            if cs.fcnt == 5:  # end of sub-frame
                cs.decode_cssr(cs.buff, 0)

    for ne in range(nep):
        obs = dec.decode_obs()
        week, tow = gn.time2gpst(obs.t)

        cs.decode_l6msg(fc.read(250), 0)
        if cs.fcnt == 5:  # end of sub-frame
            cs.week = week
            cs.decode_cssr(cs.buff, 0)

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

    dec.fobs.close()

if True:
    idx4 = np.where(smode == 4)[0]
    idx5 = np.where(smode == 5)[0]
    idx1 = np.where(smode == 1)[0]

    ms = 8
    lbl_t = ['east [m]', 'north [m]', 'up [m]']
    fig = plt.figure(figsize=(6, 10))

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
    plt.show()

    plt.plot(enu[:, 0], enu[:, 1], '-', color='gray')
    # plt.plot(enu[idx1, 0], enu[idx1, 1], 'm.', markersize=ms, label='stdpos')
    plt.plot(enu[idx5, 0], enu[idx5, 1], 'y.', markersize=ms, label='float')
    plt.plot(enu[idx4, 0], enu[idx4, 1], 'g.', markersize=ms, label='fix')

    plt.xlabel('easting [m]')
    plt.ylabel('northing [m]')
    plt.grid()
    plt.axis('equal')
    plt.legend()
    plt.show()
