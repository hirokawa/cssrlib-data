"""
 SSR correction conversion to SP3 file format
"""

from binascii import unhexlify
import bitstruct as bs
from copy import deepcopy
import numpy as np
import sys

from cssrlib.ephemeris import satpos
from cssrlib.gnss import Nav, sat2prn, sys2str
from cssrlib.gnss import time2gpst, time2doy, epoch2time, time2str
from cssrlib.gnss import timeadd
from cssrlib.gnss import uGNSS as ug, rSigRnx
from cssrlib.gnss import rCST
from cssrlib.peph import atxdec
from cssrlib.peph import peph, peph_t, apc2com
from cssrlib.cssr_has import cssr_has
from cssrlib.rinex import rnxdec


# Start epoch and number of epochs
#
ep = [2023, 8, 11, 21, 0, 0]
#ep = [2023, 7, 8, 4, 0, 0]

time = epoch2time(ep)
year = ep[0]
doy = int(time2doy(time))

nep = 900*4
step = 1

navfile = '../data{}/BRD400DLR_S_{:4d}{:03d}0000_01D_MN.rnx'\
    .format('/doy223' if doy == 223 else '', year, doy)

orbfile = '{}_{:4d}{:03d}0000_01D_05M_ORB.SP3'\
    .format('ESA0HASOPS', year, doy)

# Read Galile HAS corrections file
#
if doy == 223:
    file_has = '../data/doy{:3d}/{:3d}v_gale6.txt'.format(doy, doy)
else:
    file_has = '../data/gale6_{:3d}e.txt'.format(doy)

dtype = [('wn', 'int'), ('tow', 'int'), ('prn', 'int'),
         ('type', 'int'), ('len', 'int'), ('nav', 'S124')]
v = np.genfromtxt(file_has, dtype=dtype)

# NOTE: igs14 values seem to be yield better consistency with
#       CODE reference orbits
atxfile = '../data/igs14.atx'

rnx = rnxdec()
nav = Nav()
orb = peph()

# Decode RINEX NAV data
#
nav = rnx.decode_nav(navfile, nav)

# Setup SSR decoder
#
cs = cssr_has()
cs.monlevel = 0

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
orb_r = np.zeros((nep, ug.MAXSAT))*np.nan
orb_a = np.zeros((nep, ug.MAXSAT))*np.nan
orb_c = np.zeros((nep, ug.MAXSAT))*np.nan
clk = np.zeros((nep, ug.MAXSAT))*np.nan
sats = set()

# Initialize HAS decoding
#
mid_ = -1
ms_ = -1
icnt = 0
rec = []
mid_decoded = []
has_pages = np.zeros((255, 53), dtype=int)

ns2m = rCST.CLIGHT*1e-9

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

        # print(f"{mt} {mid} {ms} {pid}")

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

    # Convert HAS corrections
    #
    if (cs.lc[0].cstat & 0xf) == 0xf:

        ns = len(cs.sat_n)

        rs = np.ones((ns, 3))*np.nan
        vs = np.ones((ns, 3))*np.nan
        dts = np.ones((ns, 1))*np.nan

        # Store in SP3 dataset
        #
        peph = peph_t(time)

        for j, sat in enumerate(cs.sat_n):

            sys, _ = sat2prn(sat)

            rs, vs, dts, svh = satpos(sat, time, nav, cs)

            # Select PCO reference signals for Galileo HAS
            #
            if sys == ug.GPS:
                sig0 = (rSigRnx("GC1W"), rSigRnx("GC2W"))
            elif sys == ug.GAL:
                sig0 = (rSigRnx("EC1C"), rSigRnx("EC7Q"))
            else:
                print("ERROR: invalid sytem {}".format(sys2str(sys)))
                continue

            # Convert to CoM using ANTEX PCO corrections
            #
            rs[0, :] += apc2com(nav, sat, time, rs[0, :], sig0)

            for i in range(3):
                peph.pos[sat-1, i] = rs[0, i]
            peph.pos[sat-1, 3] = dts[0]

            # Store satellite in set
            #
            if sat not in sats:
                sats.add(sat)

        # Save and store
        #
        nav.peph.append(peph)

    # Next time-step
    #
    time = timeadd(time, step)

# Write results to output file
#
orb.write_sp3(orbfile, nav, sats)
