"""
 static test for RTK
"""

import matplotlib.pyplot as plt
import numpy as np
from sys import exit as sys_exit

import cssrlib.rinex as rn
import cssrlib.gnss as gn

from cssrlib.gnss import rSigRnx, time2str
from cssrlib.peph import atxdec, searchpcv
from cssrlib.rtk import rtkpos

bdir = '../data/'
ngsantfile = bdir+'GSI_PCV.TXT'

nav = gn.Nav()

if False:
    navfile = bdir+'SEPT078M.21P'
    obsfile = bdir+'SEPT078M.21O'
    basefile = bdir+'3034078M.21O'
    xyz_ref = [-3962108.673, 3381309.574, 3668678.638]
    nav.rb = [-3959400.631, 3385704.533, 3667523.111]  # GSI 3034 fujisawa
    atxfile = bdir+'igs14.atx'
else:
    navfile = bdir+'SEPT238A.23P'
    obsfile = bdir+'SEPT238A.23O'
    basefile = bdir+'3034238A.23O'
    xyz_ref = [-3962108.7007, 3381309.5532, 3668678.6648]
    nav.rb = [-3959400.6443, 3385704.4948, 3667523.1275]  # GSI 3034 fujisawa
    atxfile = bdir+'igs20.atx'

pos_ref = gn.ecef2pos(xyz_ref)

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
dec = rn.rnxdec()
dec.setSignals(sigs)

dec.decode_nav(navfile, nav)

# base
decb = rn.rnxdec()
decb.setSignals(sigs)

decb.decode_obsh(basefile)
dec.decode_obsh(obsfile)

decb.autoSubstituteSignals()
dec.autoSubstituteSignals()

nep = 300

t = np.zeros(nep)
enu = np.zeros((nep, 3))
smode = np.zeros(nep, dtype=int)

rtk = rtkpos(nav, dec.pos, 'test_rtk.log')
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
    sys_exit(-1)
nav.rcv_ant_b = searchpcv(atx.pcvr, decb.ant,  dec.ts)
if nav.rcv_ant_b is None:
    print("ERROR: missing antenna type <{}> in ANTEX file!".format(decb.ant))
    sys_exit(-1)

# nav.excl_sat = [20]
# nav.cnr_min_gpy = 20

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

    rtk.process(obs, obsb=obsb)
    t[ne] = gn.timediff(nav.t, t0)
    sol = nav.xa[0:3] if nav.smode == 4 else nav.x[0:3]
    enu[ne, :] = gn.ecef2enu(pos_ref, sol-xyz_ref)
    smode[ne] = nav.smode

    nav.fout.write("{} {:14.4f} {:14.4f} {:14.4f} "
                   "ENU {:7.4f} {:7.4f} {:7.4f}, 2D {:6.4f}, mode {:1d}\n"
                   .format(time2str(obs.t),
                           sol[0], sol[1], sol[2],
                           enu[ne, 0], enu[ne, 1], enu[ne, 2],
                           np.sqrt(enu[ne, 0]**2+enu[ne, 1]**2),
                           smode[ne]))

    # Log to standard output
    #
    sys.stdout.write('\r {} ENU {:7.4f} {:7.4f} {:7.4f}, 2D {:6.4f}, mode {:1d}'
                     .format(time2str(obs.t),
                             enu[ne, 0], enu[ne, 1], enu[ne, 2],
                             np.sqrt(enu[ne, 0]**2+enu[ne, 1]**2),
                             smode[ne]))

# Send line-break to stdout
#
sys.stdout.write('\n')

dec.fobs.close()
decb.fobs.close()

fig_type = 1
ylim = 0.2

fig = plt.figure(figsize=[7, 9])
fig.set_rasterized(True)

if fig_type == 1:
    plt.plot(t, enu)
    plt.xticks(np.arange(0, nep+1, step=30))
    plt.ylabel('position error [m]')
    plt.xlabel('time[s]')
    plt.legend(['east', 'north', 'up'])
    plt.grid()
    plt.axis([0, ne, -ylim, ylim])
else:
    plt.plot(enu[:, 0], enu[:, 1])
    plt.xlabel('easting [m]')
    plt.ylabel('northing [m]')
    plt.grid()
    plt.axis([-ylim, ylim, -ylim, ylim])

plotFileFormat = 'eps'
plotFileName = '.'.join(('test_rtk', plotFileFormat))

plt.savefig(plotFileName, format=plotFileFormat, bbox_inches='tight', dpi=300)
# plt.show()
