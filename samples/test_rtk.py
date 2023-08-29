"""
 static test for RTK
"""

import matplotlib.pyplot as plt
import numpy as np
import sys

import cssrlib.rinex as rn
import cssrlib.gnss as gn

from cssrlib.gnss import rSigRnx
from cssrlib.peph import atxdec, searchpcv
from cssrlib.rtk import rtkinit, relpos

bdir = '../data/'
atxfile = bdir+'igs14.atx'
ngsantfile = bdir+'GSI_PCV.TXT'

if False:
    navfile = bdir+'SEPT078M.21P'
    obsfile = bdir+'SEPT078M.21O'
    basefile = bdir+'3034078M.21O'
    xyz_ref = [-3962108.673, 3381309.574, 3668678.638]
    nav.rb = [-3959400.631, 3385704.533, 3667523.111]  # GSI 3034 fujisawa
else:
    navfile = bdir+'SEPT238A.23P'
    obsfile = bdir+'SEPT238A.23O'
    basefile = bdir+'3034238A.23O'
    xyz_ref = [-3962108.7007, 3381309.5532, 3668678.6648]
    nav.rb = [-3959400.6443, 3385704.4948, 3667523.1275]  # GSI 3034 fujisawa

pos_ref = gn.ecef2pos(xyz_ref)

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

# Define signals to be processed
#
sigsb = [rSigRnx("GC1C"), rSigRnx("GC2W"),
         rSigRnx("EC1X"), rSigRnx("EC5X"),
         rSigRnx("JC1X"), rSigRnx("JC2X"),
         rSigRnx("GL1C"), rSigRnx("GL2W"),
         rSigRnx("EL1X"), rSigRnx("EL5X"),
         rSigRnx("JL1X"), rSigRnx("JL2X"),
         rSigRnx("GS1C"), rSigRnx("GS2W"),
         rSigRnx("ES1X"), rSigRnx("ES5X"),
         rSigRnx("JS1X"), rSigRnx("JS2X")]

# rover
dec = rn.rnxdec()
dec.setSignals(sigs)
nav = gn.Nav()
dec.decode_nav(navfile, nav)

# base
decb = rn.rnxdec()
decb.setSignals(sigsb)
decb.decode_obsh(basefile)
dec.decode_obsh(obsfile)

nep = 300


t = np.zeros(nep)
enu = np.zeros((nep, 3))
smode = np.zeros(nep, dtype=int)
# ecef_ = np.zeros((nep, 3))

rtkinit(nav, dec.pos)
rr = dec.pos

# Load ANTEX data for satellites and stations
#
atx = atxdec()
atx.readpcv(atxfile)
atx.readngspcv(ngsantfile)

# Set PCO/PCV information
#
nav.rcv_ant = searchpcv(atx.pcvr, dec.ant,  dec.ts)
if nav.rcv_ant is None:
    print("ERROR: missing antenna type <{}> in ANTEX file!".format(dec.ant))
    sys.exit(-1)
nav.rcv_ant_b = searchpcv(atx.pcvr, decb.ant,  dec.ts)
if nav.rcv_ant_b is None:
    print("ERROR: missing antenna type <{}> in ANTEX file!".format(decb.ant))
    sys.exit(-1)

# Get equipment information
#
print("Rover:")
print("  Receiver:", dec.rcv)
print("  Antenna :", dec.ant)
print()
print("Base:")
print("  Receiver:", decb.rcv)
print("  Antenna :", decb.ant)
print()

for ne in range(nep):
    obs, obsb = rn.sync_obs(dec, decb)
    if ne == 0:
        t0 = nav.t = obs.t
    relpos(nav, obs, obsb)
    t[ne] = gn.timediff(nav.t, t0)
    sol = nav.xa[0:3] if nav.smode == 4 else nav.x[0:3]
    enu[ne, :] = gn.ecef2enu(pos_ref, sol-xyz_ref)
    # ecef_[ne, :] = sol
    smode[ne] = nav.smode

dec.fobs.close()
decb.fobs.close()

# xyz_ref = np.mean(ecef_, axis=0)

fig_type = 1
ylim = 0.2

if fig_type == 1:
    plt.plot(t, enu)
    plt.xticks(np.arange(0, nep+1, step=30))
    plt.ylabel('position error [m]')
    plt.xlabel('time[s]')
    plt.legend(['east', 'north', 'up'])
    plt.grid()
    plt.axis([0, ne, -ylim, ylim])
    plt.show()
else:
    plt.plot(enu[:, 0], enu[:, 1])
    plt.xlabel('easting [m]')
    plt.ylabel('northing [m]')
    plt.grid()
    plt.axis([-ylim, ylim, -ylim, ylim])
    plt.show()
