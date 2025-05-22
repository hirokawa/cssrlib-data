"""
SSR correction conversion to SP3 file format
"""

from itertools import chain
import numpy as np
import os
from sys import argv as sys_argv
from sys import exit as sys_exit


from cssrlib.ephemeris import satpos
from cssrlib.gnss import Nav, sat2prn, sys2str, sat2id
from cssrlib.gnss import time2doy, epoch2time, time2epoch
from cssrlib.gnss import timeadd, timeget, gpst2time, time2gpst, time2str, timediff
from cssrlib.gnss import uGNSS as ug, rSigRnx
from cssrlib.gnss import rCST
from cssrlib.peph import atxdec
from cssrlib.peph import peph, peph_t, apc2com
from cssrlib.cssrlib import sCSSRTYPE as sc
from cssrlib.cssr_pvs import cssr_pvs, decode_sinca_line
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


def file2time(fileName):
    """
    Convert hourly filename to epoch
    """
    folder = next((s for s in
                   os.path.dirname(fileName).split('/') if "doy" in s), None)
    year = int(folder[3:7])
    doy = int(folder[8:11])
    hour = ord(os.path.basename(fileName)[-5])-ord('a')
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
    ssrfiles = ['../data/doy2025-110/DAS2025110f.txt', ]

# Start time and interval length
#
time = file2time(ssrfiles[0])
itv = "{:02d}H".format(len(ssrfiles))
if itv == "24H":
    itv = "01D"

ep = time2epoch(time)

year = ep[0]
hour = ep[3]
doy = int(time2doy(time))
dur = len(ssrfiles)

if "DAS" in ssrfiles[0]:

    name = 'PVS0PPPOPS'
    step = "10S"
    atxfile = baseDirName+'../data/antex/igs20.atx'

else:

    print("ERROR: unknown SSR format for {}!".format(ssrfiles[0]))
    sys_exit(1)

# Output files
#
orbfile = '{}_{:4d}{:03d}{:02d}00_{}_{}_ORB.SP3'\
    .format(name, year, doy, hour, itv, step)
bsxfile = '{}_{:4d}{:03d}{:02d}00_{}_00U_OSB.BIA'\
    .format(name, year, doy, hour, itv)

print("Processing SSR corrections")
print("  SSR files: {}".format(" ".join([f for f in ssrfiles])))
print("  Output SP3 file: {}".format(orbfile))
print("  Output BIA file: {}".format(bsxfile))
print()

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
        pass

# Decode RINEX NAV data
#
for navfile in navfiles:
    nav = rnx.decode_nav(navfile, nav, append=True)

# Setup SSR decoder
#
if 'DAS' in ssrfiles[0]:
    cs = cssr_pvs("test.log")
else:
    print("ERROR: unknown SSR format for {}!".format(ssrfiles[0]))
    sys_exit(1)

cs.monlevel = 4

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

# Meter-to-nanosecond conversion
#
m2ns = 1e9/rCST.CLIGHT

# Loop over SSR packages
#
for ssrfile in ssrfiles:

    # Load SSR corrections
    #
    fc = open(ssrfile, 'rt')
    for line in fc:

        time, buff = decode_sinca_line(line)

        week, tow = time2gpst(time)
        cs.week = week
        cs.tow0 = tow//86400*86400
        cs.time0 = time

        if cs.cssrmode == sc.PVS_PPP:

            cs.decode_cssr(buff, 0)
            cs.check_validity(time)
            hasNew = (tow % 30 == 0)

        else:

            continue

        # Convert SSR corrections
        #
        if (cs.lc[0].cstat & 0x6) == 0x6 and hasNew:

            print(time2str(time))

            hasNew = False

            ns = len(cs.sat_n)

            rs = np.ones((ns, 3))*np.nan
            vs = np.ones((ns, 3))*np.nan
            dts = np.ones((ns, 1))*np.nan

            # Get SSR code and phase biases
            #
            for sat_, dat_ in chain(cs.lc[0].cbias.items(), cs.lc[0].pbias.items()):

                for sig_, val_ in dat_.items():

                    # Skip invalid biases
                    #
                    if np.isnan(val_):
                        continue

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

            # Store orbit and clock offset in SP3
            #
            peph = peph_t(time)

            for j, sat in enumerate(cs.sat_n):

                sys, _ = sat2prn(sat)

                rs, vs, dts, svh = satpos(sat, time, nav, cs)

                if cs.cssrmode == sc.PVS_PPP:

                    if sys == ug.GPS:
                        sig0 = (rSigRnx("GC1C"), rSigRnx("GC5Q"))
                    elif sys == ug.GAL:
                        sig0 = (rSigRnx("EC1C"), rSigRnx("EC5Q"))
                    else:
                        print("ERROR: invalid system {}".format(sys2str(sys)))
                        continue

                # Skip invalid positions
                #
                if np.isnan(rs[0, :]).any():
                    continue

                # Convert to CoM using ANTEX PCO corrections
                #
                pco = apc2com(nav, sat, time, rs[0, :], sig0, k=0)
                if pco is None:
                    continue

                # Clock reference signals for reference product
                #
                if sys == ug.GPS:
                    sigClk = (rSigRnx("GC1W"), rSigRnx("GC2W"))
                elif sys == ug.GLO:
                    sigClk = (rSigRnx("RC1C"), rSigRnx("RC2C"))
                elif sys == ug.GAL:
                    sigClk = (rSigRnx("EC1C"), rSigRnx("EC5Q"))
                elif sys == ug.BDS:
                    sigClk = (rSigRnx("CC2I"), rSigRnx("CC6I"))
                elif sys == ug.QZS:
                    sigClk = (rSigRnx("JC1C"), rSigRnx("JC5Q"))
                else:
                    print("ERROR: invalid system {}".format(sys2str(sys)))
                    continue

                freq = [s.frequency() for s in sigClk]
                facs = (+freq[0]**2/(freq[0]**2-freq[1]**2),
                        -freq[1]**2/(freq[0]**2-freq[1]**2))

                # Compute ionosphere-free combination of biases
                #
                bias = 0.0
                """
                if sat in biases.keys() and all(s in biases[sat] for s in sigClk):
                    bias = facs[0]*biases[sat][sigClk[0]][-1][2] + \
                        facs[1]*biases[sat][sigClk[1]][-1][2]
                else:
                    #print("ERROR: missing bias for {} {}".format(
                        sat2id(sat), sigClk))
                    continue

                # Adjust sign of biases
                if cs.cssrmode in [sc.GAL_HAS_SIS, sc.QZS_MADOCA]:
                    bias = -bias
                """

                # Adjust orbit and clock offset and store in SP3
                #
                for i in range(3):
                    peph.pos[sat-1, i] = rs[0, i] + pco[i]
                peph.pos[sat-1, 3] = dts[0] - bias*1e-9

                # Store satellite in set
                #
                if sat not in sats:
                    sats.add(sat)

            # Save and store
            #
            nav.peph.append(peph)

# Write results to output file
#
orb.write_sp3(orbfile, nav, sats)

# Write biases to Bias-SINEX
#
if len(biases) > 0:
    write_bsx(bsxfile, name[0:3], biases)
