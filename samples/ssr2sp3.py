"""
SSR correction conversion to SP3 file format
"""

from binascii import unhexlify
import bitstruct as bs
from copy import deepcopy
import numpy as np
import os
import sys

from cssrlib.ephemeris import satpos
from cssrlib.gnss import Nav, sat2prn, sys2str, sat2id
from cssrlib.gnss import time2gpst, time2doy, epoch2time, time2epoch, time2str
from cssrlib.gnss import timeadd, timeget
from cssrlib.gnss import uGNSS as ug, rSigRnx
from cssrlib.gnss import rCST
from cssrlib.peph import atxdec
from cssrlib.peph import peph, peph_t, apc2com
from cssrlib.cssrlib import cssr
from cssrlib.cssr_bds import cssr_bds
from cssrlib.cssr_has import cssr_has
from cssrlib.rinex import rnxdec


def time2bsxstr(t):
    """
    Time conversion to year, day-of-year and seconds-of-day
    """

    year = time2epoch(t)[0]
    doy = time2doy(t)
    sec = (doy-int(doy))*86400

    return "{:04d}:{:03d}:{:5d}".format(year, int(doy), int(sec))


def write_bsx(bsxfile, ac, data):
    """
    Write Bias-SINEX file
    """

    lines = []

    lines.append("%= <to be replaced>")
    lines.append("*"+79*'-')
    lines.append("* Bias Solution INdependent EXchange Format (Bias-SINEX)")
    lines.append("*"+79*'-')

    lines.append("+BIAS/DESCRIPTION")
    lines.append(
        "*KEYWORD________________________________ VALUE (S) _____________________________")
    lines.append(" BIAS_MODE                               ABSOLUTE")
    lines.append("-BIAS/DESCRIPTION")
    lines.append("*"+79*'-')

    lines.append("+BIAS/SOLUTION")
    lines.append("*BIAS SVN_ PRN STATION__ OBS1 OBS2 BIAS_START____ "
                 "BIAS_END______ UNIT __ESTIMATED_VALUE____ _STD_DEV___")

    nVal = 0
    tNow = timeget()
    tFirst = None
    tLast = None

    for sat in data.keys():
        for sig in data[sat].keys():
            for ts, te, osb in data[sat][sig]:
                lines.append(" OSB  {:4s} {:3s} {:9s} {:3s}  {:3s}  {} {} ns   {:21.4f} {:11.4f}"
                             .format(sat2id(sat)[0]+"999", sat2id(sat), "",
                                     sig.str(), "",
                                     time2bsxstr(ts), time2bsxstr(te),
                                     osb*ns2m, 0.0))
                nVal += 1
                tFirst = ts if tFirst is None or ts < tFirst else tFirst
                tLast = te if tLast is None or te < tLast else tLast

    lines.append("-BIAS/SOLUTION")
    lines.append("%=ENDBIA")

    # Replacefirst line
    #
    lines[0] = "%=BIA 1.00 {ac} {tn}   {ac} {ts} {te} R {nVal:08d}"\
        .format(ac=ac, tn=time2bsxstr(tNow), ts=time2bsxstr(tFirst),
                te=time2bsxstr(tLast), nVal=nVal)

    # Write to file
    #
    with open(bsxfile, 'w') as f:
        for line in lines:
            f.write(line+'\n')


baseDirName = os.path.dirname(os.path.abspath(__file__))+"/"

# SSR file for conversion
#
if len(sys.argv) > 1:
    ssrfile = sys.argv[1]
else:
    #ssrfile = '../data/gale6_189e.txt'
    ssrfile = '../data/bdsb2b_189e.txt'
    #ssrfile = '../data/qzsl6_189e.txt'

# Start epoch and number of epochs
#
if "_189e" in ssrfile:
    ep = [2023, 7, 8, 4, 0, 0]
else:
    ep = [2023, 8, 11, 21, 0, 0]

time = epoch2time(ep)
year = ep[0]
doy = int(time2doy(time))

nep = 900*4
step = 1

navfile = baseDirName+'../data{}/BRD400DLR_S_{:4d}{:03d}0000_01D_MN.rnx'\
    .format('/doy223' if doy == 223 else '', year, doy)


if "qzsl6_" in ssrfile:

    name = 'QZS0CLSOPS'

    dtype = [('wn', 'int'), ('tow', 'int'), ('prn', 'int'),
             ('type', 'int'), ('len', 'int'), ('nav', 'S500')]

    prn_ref = 199  # QZSS PRN
    l6_ch = 1  # 0:L6D, 1:L6E
    atxfile = baseDirName+'../data/igs20.atx'


elif "gale6_" in ssrfile:

    name = 'ESA0HASOPS'

    dtype = [('wn', 'int'), ('tow', 'int'), ('prn', 'int'),
             ('type', 'int'), ('len', 'int'), ('nav', 'S124')]

    # NOTE: igs14 values seem to be yield better consistency with
    #       CODE reference orbits
    atxfile = baseDirName+'../data/igs14.atx'

elif "bdsb2b_" in ssrfile:

    name = 'BDS0PPPOPS'

    dtype = [('wn', 'int'), ('tow', 'int'), ('prn', 'int'),
             ('type', 'int'), ('len', 'int'), ('nav', 'S124')]

    prn_ref = 59  # satellite PRN to receive BDS PPP collection

    # NOTE: igs14 values seem to be yield better consistency with
    #       CODE reference orbits
    atxfile = baseDirName+'../data/igs14.atx'

else:

    print("ERROR: unkown SSR format for {}!".format(ssrfile))
    sys.exit(1)

# Load SSR corrections
#
v = np.genfromtxt(ssrfile, dtype=dtype)

# Output files
#
orbfile = '{}_{:4d}{:03d}0000_01D_01S_ORB.SP3'\
    .format(name, year, doy)
bsxfile = '{}_{:4d}{:03d}0000_01D_00U_OSB.BIA'\
    .format(name, year, doy)

# Initialize objects
#
rnx = rnxdec()
nav = Nav()
orb = peph()

# Decode RINEX NAV data
#
nav = rnx.decode_nav(navfile, nav)

# Setup SSR decoder
#
if 'gale6_' in ssrfile:
    cs = cssr_has()
    file_gm = baseDirName+'Galileo-HAS-SIS-ICD_1.0_Annex_B_Reed_Solomon_Generator_Matrix.txt'
    gMat = np.genfromtxt(file_gm, dtype="u1", delimiter=",")
elif 'qzsl6_' in ssrfile:
    cs = cssr()
elif "bdsb2b_" in ssrfile:
    cs = cssr_bds()
else:
    print("ERROR: unkown SSR format for {}!".format(ssrfile))
    sys.exit(1)

cs.monlevel = 0

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

biases = {}

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

    if 'gale6_' in ssrfile:

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

    elif 'qzsl6' in ssrfile:

        vi = v[(v['tow'] == tow) & (v['type'] == l6_ch)
               & (v['prn'] == prn_ref)]
        if len(vi) > 0:
            msg = unhexlify(vi['nav'][0])
            cs.decode_l6msg(msg, 0)
            if cs.fcnt == 5:  # end of sub-frame
                cs.decode_cssr(bytes(cs.buff), 0)

    elif "bdsb2b_" in ssrfile:

        vi = v[(v['tow'] == tow) & (v['prn'] == prn_ref)]
        if len(vi) > 0:
            buff = unhexlify(vi['nav'][0])
            # prn, rev = bs.unpack_from('u6u6', buff, 0)
            cs.decode_cssr(buff, 0)

    else:

        continue

    # Convert SSR corrections
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

            # Select PCO reference signals
            #
            if sys == ug.GPS:
                sig0 = (rSigRnx("GC1W"), rSigRnx("GC2W"))
            elif sys == ug.GAL:
                sig0 = (rSigRnx("EC1C"), rSigRnx("EC7Q"))
            elif sys == ug.QZS:
                sig0 = (rSigRnx("JC1C"), rSigRnx("JC2S"))
            elif sys == ug.GLO:
                sig0 = (rSigRnx("RC1C"), rSigRnx("RC2C"))
            elif sys == ug.BDS:
                sig0 = (rSigRnx("CC6I"),)
            else:
                print("ERROR: invalid sytem {}".format(sys2str(sys)))
                continue

            # Convert to CoM using ANTEX PCO corrections
            #
            if np.linalg.norm(rs[0, :]) > 0:
                rs[0, :] += apc2com(nav, sat, time, rs[0, :], sig0, k=0)

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

        # Get SSR biases
        #
        for sat_, dat_ in cs.lc[0].cbias.items():

            for sig_, val_ in dat_.items():

                # Fix GPS L2 P(Y) signal code for Galileo HAS
                #
                if 'gale6' in ssrfile and rSigRnx('GC2P') == sig_:
                    sig_ = rSigRnx('GC2W')

                if sat_ not in biases.keys():
                    biases.update({sat_: {}})
                if sig_ not in biases[sat_].keys():
                    biases[sat_].update({sig_: []})

                if len(biases[sat_][sig_]) == 0 or \
                        biases[sat_][sig_][-1][2] != val_:
                    biases[sat_][sig_].append([time, None, val_])
                else:
                    biases[sat_][sig_][-1][1] = time

    # Next time-step
    #
    time = timeadd(time, step)

# Write results to output file
#
orb.write_sp3(orbfile, nav, sats)

# Write biases ot Bias-SINEX
#
write_bsx(bsxfile, name[0:3], biases)
