"""
 static test for PPP-RTK (QZSS CLAS)
"""
import matplotlib.pyplot as plt
import numpy as np
import sys

import cssrlib.gnss as gn
from cssrlib.cssrlib import cssr
from cssrlib.gnss import ecef2pos, Nav, time2gpst, timediff
from cssrlib.gnss import rSigRnx
from cssrlib.peph import atxdec, searchpcv
from cssrlib.ppprtk import rtkinit, ppprtkpos
from cssrlib.rinex import rnxdec

#from cssrlib.pntpos import stdinit, pntpos

atxfile = '../data/igs14.atx'
l6file = '../data/2021078M.l6'
griddef = '../data/clas_grid.def'
navfile = '../data/SEPT078M.21P'
obsfile = '../data/SEPT078M.21O'

# based on GSI F5 solution
xyz_ref = [-3962108.673,   3381309.574,   3668678.638]
pos_ref = ecef2pos(xyz_ref)

cs = cssr()
cs.monlevel = 1
cs.week = 2149
cs.read_griddef(griddef)

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

dec = rnxdec()
dec.setSignals(sigs)

nav = Nav()
nav = dec.decode_nav(navfile, nav)
nep = 180

# Load ANTEX data for satellites and stations
#
atx = atxdec()
atx.readpcv(atxfile)

t = np.zeros(nep)
tc = np.zeros(nep)
enu = np.ones((nep, 3))*np.nan
sol = np.zeros((nep, 4))
dop = np.zeros((nep, 4))
smode = np.zeros(nep, dtype=int)
if dec.decode_obsh(obsfile) >= 0:

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

    rr = dec.pos
    rtkinit(nav, dec.pos)
    pos = ecef2pos(rr)
    inet = cs.find_grid_index(pos)

    fc = open(l6file, 'rb')
    if not fc:
        print("L6 messsage file cannot open.")
        sys.exit(-1)

    for ne in range(nep):
        obs = dec.decode_obs()
        week, tow = time2gpst(obs.t)

        cs.decode_l6msg(fc.read(250), 0)
        if cs.fcnt == 5:  # end of sub-frame
            cs.week = week
            cs.decode_cssr(cs.buff, 0)

        if ne == 0:
            t0 = nav.t = obs.t
            t0.time = t0.time//30*30
            cs.time = obs.t
            nav.time_p = t0

        cstat = cs.chk_stat()
        if cstat:
            ppprtkpos(nav, obs, cs)

        t[ne] = timediff(nav.t, t0)
        tc[ne] = timediff(cs.time, t0)

        sol = nav.xa[0:3] if nav.smode == 4 else nav.x[0:3]
        enu[ne, :] = gn.ecef2enu(pos_ref, sol-xyz_ref)
        smode[ne] = nav.smode

    fc.close()
    dec.fobs.close()

fig_type = 1
ylim = 0.4

idx4 = np.where(smode == 4)[0]
idx5 = np.where(smode == 5)[0]
idx0 = np.where(smode == 0)[0]

fig = plt.figure(figsize=[7, 9])

if fig_type == 1:

    lbl_t = ['east[m]', 'north[m]', 'up[m]']
    for k in range(3):
        plt.subplot(3, 1, k+1)
        plt.plot(t[idx0], enu[idx0, k], 'r.')
        plt.plot(t[idx5], enu[idx5, k], 'y.')
        plt.plot(t[idx4], enu[idx4, k], 'g.')

        plt.xticks(np.arange(0, nep+1, step=30))
        if k == 2:
            plt.xlabel('time[s]')
        plt.ylabel(lbl_t[k])
        plt.grid()
        #plt.axis([0, ne, -ylim, ylim])

elif fig_type == 2:

    ax = fig.add_subplot(111)

    #plt.plot(enu[idx0, 0], enu[idx0, 1], 'r.', label='stdpos')
    plt.plot(enu[idx5, 0], enu[idx5, 1], 'y.', label='float')
    plt.plot(enu[idx4, 0], enu[idx4, 1], 'g.', label='fix')

    plt.xlabel('easting [m]')
    plt.ylabel('northing [m]')
    plt.grid()
    # plt.legend()
    ylim = 0.05
    ax.set(xlim=(-ylim, ylim), ylim=(-ylim, ylim))
    # plt.axis('equal')

plt.show()

if nav.fout is not None:
    nav.fout.close()
