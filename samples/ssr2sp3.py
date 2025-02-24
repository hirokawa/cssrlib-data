"""
SSR correction conversion to SP3 file format
"""

from binascii import unhexlify
import bitstruct as bs
from itertools import chain
import numpy as np
import os
from sys import argv as sys_argv
from sys import exit as sys_exit


from cssrlib.ephemeris import satpos
from cssrlib.gnss import Nav, sat2prn, sys2str, sat2id
from cssrlib.gnss import time2doy, epoch2time, time2epoch, time2str
from cssrlib.gnss import timeadd, timeget, gpst2time
from cssrlib.gnss import uGNSS as ug, rSigRnx
from cssrlib.gnss import rCST
from cssrlib.peph import atxdec
from cssrlib.peph import peph, peph_t, apc2com
from cssrlib.cssrlib import sCSSRTYPE as sc
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

    return "{:04d}:{:03d}:{:05d}".format(year, int(doy), int(sec))


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

    for sat in sorted(data.keys()):
        for sig in data[sat].keys():
            for ts, te, osb in data[sat][sig]:
                lines.append(" OSB  {:4s} {:3s} {:9s} {:3s}  {:3s}  {} {} ns   {:21.4f} {:11.4f}"
                             .format(sat2id(sat)[0]+"999", sat2id(sat), "",
                                     sig.str(), "",
                                     time2bsxstr(ts), time2bsxstr(te),
                                     -osb*m2ns, 0.0))
                nVal += 1
                tFirst = ts if tFirst is None or ts < tFirst else tFirst
                tLast = te if tLast is None or te < tLast else tLast

    lines.append("-BIAS/SOLUTION")
    lines.append("%=ENDBIA")

    # Replace first line
    #
    lines[0] = "%=BIA 1.00 {ac} {tn}   {ac} {ts} {te} R {nVal:08d}"\
        .format(ac=ac, tn=time2bsxstr(tNow), ts=time2bsxstr(tFirst),
                te=time2bsxstr(tLast), nVal=nVal)

    # Write to file
    #
    with open(bsxfile, 'w') as f:
        for line in lines:
            f.write(line+'\n')


def file2time(year, fileName):
    """
    Convert hourly SBF filename to epoch
    """

    doy = os.path.basename(fileName).split('_')[0]
    hour = ord(doy[3])-ord('a')
    doy = int(doy[0:3])
    time = epoch2time([year, 1, 1, hour, 0, 0])

    return timeadd(time, (doy-1)*86400)


# Base directory of this script
#
baseDirName = os.path.dirname(os.path.abspath(__file__))+"/"

# SSR file for conversion
#
ssrfiles = []
if len(sys_argv) > 1:
    ssrfiles = sys_argv[1:]
else:
    ssrfiles = ['../data/doy2023-189/gale6_189e.txt', ]

# Start time
#
if "_189e" in ssrfiles[0]:
    time = epoch2time([2023, 7, 8, 4, 0, 0])
else:
    time = file2time(2023, ssrfiles[0])
    """
    print("ERROR: unknown epoch!")
    sys_exit(1)
    """

ep = time2epoch(time)
year = ep[0]
hour = ep[3]
doy = int(time2doy(time))

if "qzsl6" in ssrfiles[0]:

    name = 'QZS0CLSOPS'
    step = "10S"

    dtype = [('wn', 'int'), ('tow', 'int'), ('prn', 'int'),
             ('type', 'int'), ('len', 'int'), ('nav', 'S500')]

    prn_ref = 199  # QZSS PRN
    l6_ch = 1  # 0:L6D, 1:L6E
    atxfile = baseDirName+'../data/antex/igs20.atx'

elif "gale6" in ssrfiles[0]:

    name = 'ESA0HASOPS'
    step = "10S"

    dtype = [('wn', 'int'), ('tow', 'int'), ('prn', 'int'),
             ('type', 'int'), ('len', 'int'), ('nav', 'S124')]

    # NOTE: igs14 values seem to be yield better consistency with
    #       CODE reference orbits
    atxfile = baseDirName+'../data/antex/igs14.atx'

elif "bdsb2b" in ssrfiles[0]:

    name = 'BDS0PPPOPS'
    step = "10S"

    dtype = [('wn', 'int'), ('tow', 'int'), ('prn', 'int'),
             ('type', 'int'), ('len', 'int'), ('nav', 'S124')]

    prn_ref = 59  # satellite PRN to receive BDS PPP collection
    atxfile = baseDirName+'../data/antex/igs20.atx'

else:

    print("ERROR: unknown SSR format for {}!".format(ssrfiles[0]))
    sys_exit(1)

# Output files
#
orbfile = '{}_{:4d}{:03d}{:02d}00_01D_{}_ORB.SP3'\
    .format(name, year, doy, hour, step)
bsxfile = '{}_{:4d}{:03d}{:02d}00_01D_00U_OSB.BIA'\
    .format(name, year, doy, hour)

# Initialize objects
#
rnx = rnxdec()
nav = Nav()
orb = peph()

# Load RINEX navigation files
#
navfiles = []
for dt in (-1, 0, +1):

    t = timeadd(time, dt*86400)
    ep = time2epoch(t)
    year = ep[0]
    hour = ep[3]
    doy = int(time2doy(t))

    navfile = baseDirName+'../data/brdc/BRD400DLR_S_{:4d}{:03d}0000_01D_MN.rnx'\
        .format(year, doy, year, doy)

    if os.path.exists(navfile):
        navfiles.append(navfile)
    else:
        print("WARNING: cannot find  {}".format(navfile))

# Decode RINEX NAV data
#
for navfile in navfiles:
    nav = rnx.decode_nav(navfile, nav, append=True)

# Setup SSR decoder
#
if 'gale6' in ssrfiles[0]:
    cs = cssr_has()
    file_gm = baseDirName+'Galileo-HAS-SIS-ICD_1.0_Annex_B_Reed_Solomon_Generator_Matrix.txt'
    gMat = np.genfromtxt(file_gm, dtype="u1", delimiter=",")
elif 'qzsl6' in ssrfiles[0]:
    cs = cssr()
    cs.cssrmode = sc.QZS_MADOCA
elif "bdsb2b" in ssrfiles[0]:
    cs = cssr_bds()
else:
    print("ERROR: unknown SSR format for {}!".format(ssrfiles[0]))
    sys_exit(1)

cs.monlevel = 0

# Load SSR corrections
#
v = np.array([], dtype=dtype)
for ssrfile in ssrfiles:
    v = np.append(v, np.genfromtxt(ssrfile, dtype=dtype))

# Load ANTEX data for satellites and stations
#
atx = atxdec()
atx.readpcv(atxfile)

# Set PCO/PCV information
#
nav.sat_ant = atx.pcvs

# Initialize data structures for results
#
t0 = None
biases = {}
sats = set()

# Set flags for adding additional biases
#
extClasBiases = True
extBdsBiases = True

# Initialize HAS decoding
#
mid_ = -1
ms_ = -1
icnt = 0
rec = []
mid_decoded = []
has_pages = np.zeros((255, 53), dtype=int)

# Meter-to-nanosecond conversion
#
m2ns = 1e9/rCST.CLIGHT

# Loop over SSR packages
#
for vi in v:

    week, tow = vi['wn'], vi['tow']
    time = gpst2time(week, tow)
    cs.week = week
    cs.tow0 = tow//3600*3600

    hasNew = False

    if cs.cssrmode == sc.GAL_HAS_SIS:

        buff = unhexlify(vi['nav'])
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

        if len(rec) >= ms_:
            HASmsg = cs.decode_has_page(rec, has_pages, gMat, ms_)
            cs.decode_cssr(HASmsg)
            hasNew = (ms_ == 2)  # only clock messages
            time = cs.time

            rec = []
            mid_decoded += [mid_]
            mid_ = -1
            if len(mid_decoded) > 10:
                mid_decoded = mid_decoded[1:]
        else:
            icnt += 1
            if icnt > 10 and mid_ != -1:
                icnt = 0
                rec = []
                mid_ = -1

    elif cs.cssrmode == sc.QZS_MADOCA:

        if vi['type'] != l6_ch or vi['prn'] != prn_ref:
            continue

        msg = unhexlify(vi['nav'])
        cs.decode_l6msg(msg, 0)

        if cs.fcnt == 5:  # end of sub-frame
            cs.decode_cssr(bytes(cs.buff), 0)
            hasNew = True
            time = cs.time

    elif cs.cssrmode == sc.BDS_PPP:

        if vi['prn'] != prn_ref:
            continue

        buff = unhexlify(vi['nav'])
        cs.decode_cssr(buff, 0)
        hasNew = (tow % 10 == 0)

    else:

        continue

    # Convert SSR corrections
    #
    if (cs.lc[0].cstat & 0xf) == 0xf and hasNew:

        print(time2str(time))

        hasNew = False

        ns = len(cs.sat_n)

        rs = np.ones((ns, 3))*np.nan
        vs = np.ones((ns, 3))*np.nan
        dts = np.ones((ns, 1))*np.nan

        # Store in SP3 data set
        #
        peph = peph_t(time)

        for j, sat in enumerate(cs.sat_n):

            sys, _ = sat2prn(sat)

            rs, vs, dts, svh = satpos(sat, time, nav, cs)

            if cs.cssrmode == sc.QZS_MADOCA:

                if sys == ug.GPS:
                    sig0 = (rSigRnx("GC1C"), rSigRnx("GC2W"))
                elif sys == ug.GLO:
                    sig0 = (rSigRnx("RC1C"), rSigRnx("RC2C"))
                elif sys == ug.GAL:
                    sig0 = (rSigRnx("EC1C"), rSigRnx("EC5Q"))
                elif sys == ug.QZS:
                    sig0 = (rSigRnx("JC1C"), rSigRnx("JC2S"))
                else:
                    print("ERROR: invalid system {}".format(sys2str(sys)))
                    continue

            elif cs.cssrmode == sc.GAL_HAS_SIS:

                if sys == ug.GPS:
                    sig0 = (rSigRnx("GC1C"), rSigRnx("GC2W"))
                elif sys == ug.GAL:
                    sig0 = (rSigRnx("EC1C"), rSigRnx("EC7Q"))
                else:
                    print("ERROR: invalid system {}".format(sys2str(sys)))
                    continue

            elif cs.cssrmode == sc.BDS_PPP:

                if sys == ug.GPS:
                    sig0 = (rSigRnx("GC1C"), rSigRnx("GC2W"))
                elif sys == ug.BDS:
                    sig0 = (rSigRnx("CC6I"),)
                else:
                    print("ERROR: invalid system {}".format(sys2str(sys)))
                    continue

            # Skip invalid positions
            #
            if np.isnan(rs[0, :]).any():
                continue

            # Convert to CoM using ANTEX PCO corrections
            #
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

        # Get SSR code and phase biases
        #
        for sat_, dat_ in chain(cs.lc[0].cbias.items(), cs.lc[0].pbias.items()):

            for sig_, val_ in dat_.items():

                # Skip invalid biases
                #
                if np.isnan(val_):
                    continue

                # Fix GPS L2 P(Y) signal code for Galileo HAS
                #
                if cs.cssrmode == sc.GAL_HAS_SIS and rSigRnx('GC2P') == sig_:
                    sig_ = sig_.toAtt('W')

                if sat_ not in biases.keys():
                    biases.update({sat_: {}})
                if sig_ not in biases[sat_].keys():
                    biases[sat_].update({sig_: []})

                # Add first entry if empty
                #
                if len(biases[sat_][sig_]) == 0:
                    biases[sat_][sig_].append([time, time, val_])

                # Extend previous record with end time of current record
                #
                biases[sat_][sig_][-1][1] = time

                # Add new value if bias has changed
                #
                if biases[sat_][sig_][-1][2] != val_:
                    biases[sat_][sig_].append([time, time, val_])

                # Add additional Galileo biases for QZSS CLAS
                #
                if cs.cssrmode == sc.QZS_MADOCA and extClasBiases:

                    if rSigRnx('EC1X') == sig_ or rSigRnx('EL1X') == sig_:
                        sig_ = sig_.toAtt('C')
                    elif rSigRnx('EC5X') == sig_ or rSigRnx('EL5X') == sig_:
                        sig_ = sig_.toAtt('Q')

                    if sig_ not in biases[sat_].keys():
                        biases[sat_].update({sig_: []})

                    # Add first entry if empty
                    #
                    if len(biases[sat_][sig_]) == 0:
                        biases[sat_][sig_].append([time, time, val_])

                    # Extend previous record with end time of current record
                    #
                    biases[sat_][sig_][-1][1] = time

                    # Add new value if bias has changed
                    #
                    if biases[sat_][sig_][-1][2] != val_:
                        biases[sat_][sig_].append([time, time, val_])

        # Add fake GPS code biases for BeiDou B2b-PPP
        #
        if cs.cssrmode == sc.BDS_PPP and extBdsBiases:

            for sat_ in range(1, 33):

                sigs = [rSigRnx('GC1C'), rSigRnx('GC2W')]
                for sig_ in sigs:

                    val_ = 0.0

                    if sat_ not in biases.keys():
                        biases.update({sat_: {}})
                    if sig_ not in biases[sat_].keys():
                        biases[sat_].update({sig_: []})

                    # Add first entry if empty
                    #
                    if len(biases[sat_][sig_]) == 0:
                        biases[sat_][sig_].append([time, time, val_])

                    # Extend previous record with end time of current record
                    #
                    biases[sat_][sig_][-1][1] = time

                    # Add new value if bias has changed
                    #
                    if biases[sat_][sig_][-1][2] != val_:
                        biases[sat_][sig_].append([time, time, val_])

# Write results to output file
#
orb.write_sp3(orbfile, nav, sats)

# Write biases to Bias-SINEX
#
write_bsx(bsxfile, name[0:3], biases)
