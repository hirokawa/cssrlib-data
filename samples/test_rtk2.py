"""
 kinematic test for RTK
"""
import matplotlib.pyplot as plt
import numpy as np
from sys import exit as sys_exit

import cssrlib.rinex as rn
import cssrlib.gnss as gn

from cssrlib.gnss import rSigRnx, time2doy, epoch2time
from cssrlib.peph import atxdec, searchpcv
from cssrlib.rtk import rtkpos

ep = [2021, 9, 22, 6, 30, 0]

time = epoch2time(ep)
year = ep[0]
doy = int(time2doy(time))
let = chr(ord('A')+ep[3])

atxfile = '../data/antex/'
if time > epoch2time([2022, 11, 27, 0, 0, 0]):
    atxfile += 'igs20.atx'
else:
    atxfile += 'igs14.atx'

bdir = '../data/doy{:04d}-{:03d}/'.format(year, doy)
navfile = bdir+'SEPT{:03d}{}.{:02d}P'.format(doy, '0', year % 2000)
obsfile = bdir+'SEPT{:03d}{}.{:02d}O'.format(doy, let, year % 2000)
basefile = bdir+'3034{:03d}{}.{:02d}O'.format(doy, let, year % 2000)

xyz_ref = gn.pos2ecef([35.342058098, 139.521986657, 47.5515], True)

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

nep = 360

# GSI 3034 fujisawa
nav.rb = [-3959400.631, 3385704.533, 3667523.111]
t = np.zeros(nep)
enu = np.zeros((nep, 3))
smode = np.zeros(nep, dtype=int)

rr0 = [-3961951.1326752,  3381198.11019757,  3668916.0417232]  # from pntpos
pos_ref = gn.ecef2pos(xyz_ref)

rtk = rtkpos(nav, rr0, 'test_rtk2.log')
rr = rr0

# Load ANTEX data for satellites and stations
#
atx = atxdec()
atx.readpcv(atxfile)

# Set PCO/PCV information
#
nav.rcv_ant = searchpcv(atx.pcvr, dec.ant,  dec.ts)
if nav.rcv_ant is None:
    print("ERROR: missing antenna type <{}> in ANTEX file!".format(dec.ant))
    sys_exit(-1)

nav.rcv_ant_b = searchpcv(atx.pcvr, decb.ant,  dec.ts)
if nav.rcv_ant_b is None:
    print("ERROR: missing antenna type <{}> in ANTEX file!".format(decb.ant))
    sys_exit(-1)

# Print equipment information
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
    rtk.process(obs, obsb=obsb)
    t[ne] = gn.timediff(nav.t, t0)
    sol = nav.xa[0:3] if nav.smode == 4 else nav.x[0:3]
    enu[ne, :] = gn.ecef2enu(pos_ref, sol-xyz_ref)
    smode[ne] = nav.smode

dec.fobs.close()
decb.fobs.close()


idx4 = np.where(smode == 4)[0]
idx5 = np.where(smode == 5)[0]
idx1 = np.where(smode == 1)[0]

# East-north-up position error plot
#

lbl_t = ['east [m]', 'north [m]', 'up [m]']
fig = plt.figure(figsize=(6, 10))

for k in range(3):
    plt.subplot(3, 1, k+1)
    plt.plot(t, enu[:, k], '-', color='gray')
    # plt.plot(t[idx1], enu[idx1, k], 'm.', label='stdpos')
    plt.plot(t[idx5], enu[idx5, k], 'r.', markersize=1, label='float')
    plt.plot(t[idx4], enu[idx4, k], 'b.', markersize=1, label='fix')
    plt.xticks(np.arange(0, nep+1, step=30))
    plt.ylabel(lbl_t[k])
    plt.xlabel('time[s]')
    if k == 0:
        plt.legend()
    plt.grid()

plotFileFormat = 'eps'
plotFileName = '.'.join(('test_rtk2_1', plotFileFormat))
plt.savefig(plotFileName, format=plotFileFormat, bbox_inches='tight', dpi=300)
# plt.show()

# East-north trajectory plot
#
fig = plt.figure(figsize=(6, 10))
plt.plot(enu[:, 0], enu[:, 1], '-', color='gray')
# plt.plot(enu[idx1, 0], enu[idx1, 1], 'm.', markersize=1, label='stdpos')
plt.plot(enu[idx5, 0], enu[idx5, 1], 'r.', markersize=1, label='float')
plt.plot(enu[idx4, 0], enu[idx4, 1], 'b.', markersize=1, label='fix')

plt.xlabel('easting [m]')
plt.ylabel('northing [m]')
plt.grid()
plt.axis('equal')
plt.legend()

plotFileFormat = 'eps'
plotFileName = '.'.join(('test_rtk2_2', plotFileFormat))
plt.savefig(plotFileName, format=plotFileFormat, bbox_inches='tight', dpi=300)
# plt.show()
