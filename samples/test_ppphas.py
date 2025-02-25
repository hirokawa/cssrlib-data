"""
 static test for PPP (Galileo HAS SIS)
"""

from binascii import unhexlify
import bitstruct as bs
from copy import deepcopy
import matplotlib.pyplot as plt
import matplotlib.dates as md
import numpy as np
from sys import exit as sys_exit
from sys import stdout

import cssrlib.gnss as gn
from cssrlib.gnss import ecef2pos, Nav
from cssrlib.gnss import time2gpst, time2doy, time2str, timediff, epoch2time
from cssrlib.gnss import rSigRnx, prn2sat, uGNSS
from cssrlib.gnss import sys2str
from cssrlib.peph import atxdec, searchpcv
from cssrlib.cssr_has import cssr_has
from cssrlib.pppssr import pppos
from cssrlib.rinex import rnxdec


# Select test case
#
dataset = 2
excl_sat = []
# Start epoch and number of epochs
#
fromSbfConvert = False
if dataset == 0:
    ep = [2023, 7, 8, 4, 0, 0]
    xyz_ref = [-3962108.7007, 3381309.5532, 3668678.6648]
    navfile = '../data/doy2023-189/SEPT1890.23P'
    obsfile = '../data/doy2023-189/SEPT1890.23O'
    file_has = '../data/doy2023-189/gale6_189e.txt'
elif dataset == 1:
    ep = [2023, 8, 11, 21, 0, 0]
    xyz_ref = [-3962108.7007, 3381309.5532, 3668678.6648]
    navfile = '../data/doy2023-223/NAV223.23p'
    # obsfile = '../data/doy2023-223/SEPT223Z.23O'  # MOSAIC-CLAS
    obsfile = '../data/doy2023-223/SEPT223Y.23O'  # PolaRX5
    file_has = '../data/doy2023-223/223v_gale6.txt'
elif dataset == 2:
    ep = [2025, 2, 15, 17, 0, 0]
    xyz_ref = [-3962108.6836, 3381309.5672, 3668678.6720]
    navfile = '../data/doy2025-046/046r_rnx.nav'  # Mosaic-X5
    obsfile = '../data/doy2025-046/046r_rnx.obs'  # Mosaic-X5
    file_has = '../data/doy2025-046/046r_gale6.txt'
    excl_sat = [prn2sat(uGNSS.GAL, 29)]  # E29

# Convert epoch and user reference position
#
pos_ref = ecef2pos(xyz_ref)
time = epoch2time(ep)
year = ep[0]
doy = int(time2doy(time))

nep = 900*4

# Load SSR correction file
#
if fromSbfConvert:

    dtype = [('tow', 'float64'), ('wn', 'int'),  ('prn', 'S3'),
             ('validity', 'object'), ('num1', 'int'), ('signal', 'object'),
             ('num2', 'int'), ('num3', 'int'), ('nav', 'S144')]
    v = np.genfromtxt(file_has, dtype=dtype, delimiter=',')
    v = v[v['validity'] == b'Passed']
    # Eliminate whitespace
    for i, nav in enumerate(v['nav']):
        v[i]['nav'] = (b''.join(nav.split()))

else:

    dtype = [('wn', 'int'), ('tow', 'int'), ('prn', 'int'),
             ('type', 'int'), ('len', 'int'), ('nav', 'S124')]
    v = np.genfromtxt(file_has, dtype=dtype)

# Define signals to be processed
#
gnss = "GE"
sigs = []
if 'G' in gnss:
    sigs.extend([rSigRnx("GC1C"), rSigRnx("GC2W"),
                 rSigRnx("GL1C"), rSigRnx("GL2W"),
                 rSigRnx("GS1C"), rSigRnx("GS2W")])
if 'E' in gnss:
    sigs.extend([rSigRnx("EC1C"), rSigRnx("EC7Q"),
                 rSigRnx("EL1C"), rSigRnx("EL7Q"),
                 rSigRnx("ES1C"), rSigRnx("ES7Q")])

rnx = rnxdec()
rnx.setSignals(sigs)

nav = Nav()

# Positioning mode
# 0:static, 1:kinematic
#
nav.pmode = 0

if len(excl_sat) > 0:
    nav.excl_sat = excl_sat

# Decode RINEX NAV data
#
nav = rnx.decode_nav(navfile, nav)

cs = cssr_has()
cs.monlevel = 0
"""
cs = cssr_has('has.log')
cs.monlevel = 2
"""

file_gm = "Galileo-HAS-SIS-ICD_1.0_Annex_B_Reed_Solomon_Generator_Matrix.txt"
gMat = np.genfromtxt(file_gm, dtype="u1", delimiter=",")

# Load ANTEX data for satellites and stations
#
atxfile = '../data/antex/igs14.atx'
atx = atxdec()
atx.readpcv(atxfile)

# Initialize data structures for results
#
t = np.zeros(nep)
enu = np.ones((nep, 3))*np.nan
sol = np.zeros((nep, 4))
ztd = np.zeros((nep, 1))
smode = np.zeros(nep, dtype=int)

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
    ppp = pppos(nav, rnx.pos, 'test_ppphas.log')
    nav.elmin = np.deg2rad(5.0)

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

    mid_ = -1
    ms_ = -1
    icnt = 0
    rec = []
    mid_decoded = []
    has_pages = np.zeros((255, 53), dtype=int)

    # Skip epochs until start time
    #
    obs = rnx.decode_obs()
    while time > obs.t and obs.t.time != 0:
        obs = rnx.decode_obs()

    # Loop over number of epoch from file start
    #
    for ne in range(nep):

        week, tow = time2gpst(obs.t)
        cs.week = week
        cs.tow0 = tow//3600*3600

        # Set initial epoch
        #
        if ne == 0:
            nav.t = deepcopy(obs.t)
            t0 = deepcopy(obs.t)
            t0.time = t0.time//30*30
            nav.time_p = t0

        vi = v[v['tow'] == tow]
        for vn in vi:
            buff = unhexlify(vn['nav'])
            i = 14
            if bs.unpack_from('u24', buff, i)[0] == 0xaf3bc3:
                continue
            hass, res = bs.unpack_from('u2u2', buff, i)
            i += 4
            if hass >= 2:  # 0:test,1:operational,2:res,3:dnu
                continue
            mt, mid, ms, pid = bs.unpack_from('u2u5u5u8', buff, i)

            cs.msgtype = mt
            ms += 1
            i += 20

            if mid_ == -1 and mid not in mid_decoded:
                mid_ = mid
                ms_ = ms
            if mid == mid_ and pid-1 not in rec:
                page = bs.unpack_from('u8'*53, buff, i)
                rec += [pid-1]
                has_pages[pid-1, :] = page

            # print(f"{mt} {mid} {ms} {pid}")

        if len(rec) >= ms_:
            if cs.monlevel >= 2:
                print("data collected mid={:2d} ms={:2d} tow={:.0f}"
                      .format(mid_, ms_, tow))
            HASmsg = cs.decode_has_page(rec, has_pages, gMat, ms_)
            cs.decode_cssr(HASmsg)
            rec = []

            mid_decoded += [mid_]
            mid_ = -1
            if len(mid_decoded) > 10:
                mid_decoded = mid_decoded[1:]
        else:
            icnt += 1
            if icnt > 2*ms_ and mid_ != -1:
                icnt = 0
                if cs.monlevel >= 2:
                    print(f"reset mid={mid_} ms={ms_} tow={tow}")
                rec = []
                mid_ = -1

        # Call PPP module with HAS corrections
        #
        if (cs.lc[0].cstat & 0xf) == 0xf:
            ppp.process(obs, cs=cs)

        # Save output
        #
        t[ne] = timediff(nav.t, t0)/86400.0

        sol = nav.xa[0:3] if nav.smode == 4 else nav.x[0:3]
        enu[ne, :] = gn.ecef2enu(pos_ref, sol-xyz_ref)

        ztd[ne] = nav.xa[ppp.IT(nav.na)] \
            if nav.smode == 4 else nav.x[ppp.IT(nav.na)]
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

    # Close RINEX observation file
    #
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

fmt = '%H:%M'

if fig_type == 1:

    lbl_t = ['East [m]', 'North [m]', 'Up [m]']

    for k in range(3):
        plt.subplot(4, 1, k+1)
        plt.plot(t[idx0], enu[idx0, k], 'r.')
        plt.plot(t[idx5], enu[idx5, k], 'y.')
        plt.plot(t[idx4], enu[idx4, k], 'g.')

        plt.ylabel(lbl_t[k])
        plt.grid()
        plt.ylim([-ylim, ylim])
        plt.gca().xaxis.set_major_formatter(md.DateFormatter(fmt))

    plt.subplot(4, 1, 4)
    plt.plot(t[idx0], ztd[idx0]*1e2, 'r.', markersize=8, label='none')
    plt.plot(t[idx5], ztd[idx5]*1e2, 'y.', markersize=8, label='float')
    plt.plot(t[idx4], ztd[idx4]*1e2, 'g.', markersize=8, label='fix')
    plt.ylabel('ZTD [cm]')
    plt.grid()
    plt.gca().xaxis.set_major_formatter(md.DateFormatter(fmt))

    plt.xlabel('Time [HH:MM]')
    plt.legend()

elif fig_type == 2:

    ax = fig.add_subplot(111)

    plt.plot(enu[idx0, 0], enu[idx0, 1], 'r.', label='none')
    plt.plot(enu[idx5, 0], enu[idx5, 1], 'y.', label='float')
    plt.plot(enu[idx4, 0], enu[idx4, 1], 'g.', label='fix')

    plt.xlabel('Easting [m]')
    plt.ylabel('Northing [m]')
    plt.grid()
    plt.axis('equal')
    plt.legend()
    # ax.set(xlim=(-ylim, ylim), ylim=(-ylim, ylim))

plotFileFormat = 'eps'
plotFileName = '.'.join(('test_ppphas', plotFileFormat))

plt.savefig(plotFileName, format=plotFileFormat, bbox_inches='tight', dpi=300)
# plt.show()
