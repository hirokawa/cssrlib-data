"""
 static test for PPP (Galileo HAS)
"""
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
from os.path import exists

from cssrlib.ephemeris import eph2pos, findeph
from cssrlib.gnss import Nav, sat2prn, sat2id, vnorm
from cssrlib.gnss import time2gpst, time2doy, timeadd, timediff, epoch2time, time2str
from cssrlib.gnss import uGNSS, rSigRnx, rCST, sys2str
from cssrlib.peph import atxdec
from cssrlib.peph import peph, biasdec, apc2com
from cssrlib.cssr_has import cssr_has
from cssrlib.rinex import rnxdec
from binascii import unhexlify
import bitstruct.c as bs

# Start epoch and number of epochs
#
ep = [2023, 7, 8, 4, 0, 0]

time = epoch2time(ep)
year = ep[0]
doy = int(time2doy(time))

nep = 3600
step = 1

#navfile = '../data/SEPT1890.23P'
navfile = '../data/BRDC00IGS_R_20231890000_01D_MN.rnx'

orbfile = '../data/COD0OPSFIN_{:4d}{:03d}0000_01D_15M_ORB.SP3'\
    .format(year, doy)

clkfile = '../data/COD0OPSFIN_{:4d}{:03d}0000_01D_30S_CLK.CLK'\
    .format(year, doy)

bsxfile = '../data/COD0OPSFIN_{:4d}{:03d}0000_01D_01D_OSB.BIA'\
    .format(year, doy)

if not exists(orbfile):
    orbfile = orbfile.replace('_15M_', '_05M_')

# Read Galile HAS corrections file
#
file_has = '../data/gale6_189e.txt'
dtype = [('wn', 'int'), ('tow', 'int'), ('prn', 'int'),
         ('type', 'int'), ('len', 'int'), ('nav', 'S124')]
v = np.genfromtxt(file_has, dtype=dtype)

# Define signals to be processed
#
sigs = [rSigRnx("GC1C"), rSigRnx("GC2W"),
        rSigRnx("GL1C"), rSigRnx("GL2W"),
        rSigRnx("GS1C"), rSigRnx("GS2W"),
        rSigRnx("EC1C"), rSigRnx("EC7Q"),
        rSigRnx("EL1C"), rSigRnx("EL7Q"),
        rSigRnx("ES1C"), rSigRnx("ES7Q")]

"""
if time > epoch2time([2022, 11, 22, 0, 0, 0]):
    atxfile = '../data/igs20.atx'
else:
    atxfile = '../data/igs14.atx'
"""

# NOTE: igs14 values seem to be yield better consistency with
#       CODE reference orbits
atxfile = '../data/igs14.atx'

rnx = rnxdec()
nav = Nav()
orb = peph()

# Decode RINEX NAV data
#
nav = rnx.decode_nav(navfile, nav)

# Load precise orbits and clock offsets
#
nav = orb.parse_sp3(orbfile, nav)
nav = rnx.decode_clk(clkfile, nav)

# Load code and phase biases from Bias-SINEX
#
bsx = biasdec()
bsx.parse(bsxfile)

# Setup SSR decoder
#
cs = cssr_has()
cs.mon_level = 2

file_gm = "Galileo-HAS-SIS-ICD_1.0_Annex_B_Reed_Solomon_Generator_Matrix.txt"
gMat = np.genfromtxt(file_gm, dtype="u1", delimiter=",")

# Load ANTEX data for satellites and stations
#
atx = atxdec()
atx.readpcv(atxfile)

# Set PCO/PCV information
#
nav.sat_ant = atx.pcvs

# Intialize data structures for results
#
t = np.zeros(nep)
orb_r = np.zeros((nep, uGNSS.MAXSAT))*np.nan
orb_a = np.zeros((nep, uGNSS.MAXSAT))*np.nan
orb_c = np.zeros((nep, uGNSS.MAXSAT))*np.nan
clk = np.zeros((nep, uGNSS.MAXSAT))*np.nan

# Initialize HAS decoding
#
mid_ = -1
ms_ = -1
icnt = 0
rec = []
mid_decoded = []
has_pages = np.zeros((255, 53), dtype=int)

# Loop over number of epochs from start time
#
for ne in range(nep):

    print(time2str(time))

    week, tow = time2gpst(time)
    cs.week = week
    cs.tow0 = tow//3600*3600

    # Set initial epoch
    #
    if ne == 0:
        t0 = deepcopy(time)
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

        #print(f"{mt} {mid} {ms} {pid}")

    if len(rec) >= ms_:
        print("data collected mid={:2d} ms={:2d}".format(mid_, ms_))
        HASmsg = cs.decode_has_page(rec, has_pages, gMat, ms_)
        cs.decode_cssr(HASmsg)
        rec = []

        mid_decoded += [mid_]
        mid_ = -1
        if len(mid_decoded) > 10:
            mid_decoded = mid_decoded[1:]
    else:
        icnt += 1
        if icnt > 10 and mid_ != -1:
            icnt = 0
            print(f"reset mid={mid_} ms={ms_}")
            rec = []
            mid_ = -1

    # Call PPP module with HAS corrections
    #
    if (cs.lc[0].cstat & 0xf) == 0xf:

        ns = len(cs.sat_n)

        rs0 = np.ones((ns, 3))*np.nan
        vs0 = np.ones((ns, 3))*np.nan
        dts0 = np.ones((ns, 1))*np.nan

        rs = np.ones((ns, 3))*np.nan
        vs = np.ones((ns, 3))*np.nan
        dts = np.ones((ns, 1))*np.nan

        d_rs = np.ones((ns, 3))*np.nan
        d_dts = np.ones((ns, 1))*np.nan

        for j, sat in enumerate(cs.sat_n):

            sys, _ = sat2prn(sat)

            # Precise reference orbit and clock
            #
            rs_, dts_, _ = orb.peph2pos(time, sat, nav)
            rs0[j, :] = rs_[0:3]
            vs0[j, :] = rs_[3:6]
            dts0[j] = dts_[0]

            """
            print("{} {} prec xyz [m] {:14.3f} {:14.3f}m {:14.3f} clk [ms] {:12.6f}"
                  .format(time2str(time), sat2id(sat),
                          rs0[j,0], rs0[j,1], rs0[j,2], dts0[j]*1e6))
            """

            # Broadcast ephemeris IODE, orbit and clock corrections for HAS
            #
            idx = cs.sat_n.index(sat)
            iode = cs.lc[0].iode[idx]
            dorb = cs.lc[0].dorb[idx, :]

            if sat not in cs.sat_n_p:
                continue
            idx = cs.sat_n_p.index(sat)
            dclk = cs.lc[0].dclk[idx]

            mode = cs.nav_mode[sys]

            # Get position, velocity and clock offset from broadcast ephemerides
            #
            eph = findeph(nav.eph, time, sat, iode, mode=mode)
            if eph is None:
                print("ERROR: cannot find BRDC for {}".format(sat2id(sat)))
                continue

            rs[j, :], vs[j, :], dts[j] = eph2pos(time, eph, True)

            """
            print("{} {} brdc xyz [m] {:14.3f} {:14.3f}m {:14.3f} clk [ms] {:12.6f}"
                  .format(time2str(time), sat2id(sat),
                          rs[j,0], rs[j,1], rs[j,2], dts[j]*1e6))
            """

            # Convert to CoM using ANTEX PCO corrections
            #

            # Select PCO reference signals
            #
            if sys == uGNSS.GPS:
                sigs = (rSigRnx("GC1W"), rSigRnx("GC2W"))
            elif sys == uGNSS.GAL:
                sigs = (rSigRnx("EC1C"), rSigRnx("EC7Q"))
            else:
                print("ERROR: invalid sytem {}".format(sys2str(sys)))
                continue

            rs[j, :] += apc2com(nav, sat, time, rs[j, :], sigs)

            # Along-track, cross-track and radial conversion
            #
            ea = vnorm(vs[j, :])
            rc = np.cross(rs[j, :], vs[j, :])
            ec = vnorm(rc)
            er = np.cross(ea, ec)
            A = np.array([er, ea, ec])

            # Convert orbit corrections from orbital frame to ECEF
            #
            dorb_e = dorb@A

            # Apply SSR correction
            #
            # NOTE: For Galileo HAS, add the orbit corrections to the BRDC orbit
            #
            rs[j, :] += dorb_e
            dts[j] += dclk/rCST.CLIGHT

            """
            print("{} {} hasc xyz [m] {:14.3f} {:14.3f}m {:14.3f} clk [ms] {:12.6f}"
                  .format(time2str(time), sat2id(sat),
                          rs[j,0], rs[j,1], rs[j,2], dts[j]*1e6))
            """

            # HAS vs. precise
            #
            d_rs[j, :] = (rs[j, :] - rs0[j, :])@A.T
            d_dts[j, 0] = dts[j, 0] - dts0[j, 0]

            print("{} {} diff rac [m] {:14.3f} {:14.3f} {:14.3f} clk [ns] {:12.6f}"
                  .format(time2str(time), sat2id(sat),
                          d_rs[j, 0], d_rs[j, 1], d_rs[j, 2], d_dts[j, 0]*1e9))

            orb_r[ne, sat-1] = d_rs[j, 0]
            orb_a[ne, sat-1] = d_rs[j, 1]
            orb_c[ne, sat-1] = d_rs[j, 2]
            clk[ne, sat-1] = d_dts[j, 0]*1e9

    # Save output
    #
    t[ne] = timediff(time, t0)/60

    # Next time-step
    #
    time = timeadd(time, step)

fig = plt.figure(figsize=[7, 9])
fig.set_rasterized(True)

lbl_t = ['Radial [m]', 'Along [m]', 'Cross [m]', 'Clock [ns]']

plt.subplot(4, 1, 1)
plt.plot(t, orb_r, 'r.')
plt.ylabel(lbl_t[0])
plt.grid()

plt.subplot(4, 1, 2)
plt.plot(t, orb_a, 'r.')
plt.ylabel(lbl_t[1])
plt.grid()

plt.subplot(4, 1, 3)
plt.plot(t, orb_c, 'r.')
plt.ylabel(lbl_t[2])
plt.grid()

plt.subplot(4, 1, 4)
plt.plot(t, clk, 'r.')
plt.ylabel(lbl_t[3])
plt.grid()

plotFileFormat = 'eps'
plotFileName = '.'.join(('test_has_sisre', plotFileFormat))

plt.savefig(plotFileName, format=plotFileFormat, bbox_inches='tight', dpi=300)

plt.show()
