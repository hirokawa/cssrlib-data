"""
 Signal-In-Space Range Error for Galileo HAS
"""
from binascii import unhexlify
import bitstruct as bs
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
from os.path import exists
import sys

from cssrlib.ephemeris import eph2pos, findeph
from cssrlib.gnss import Nav, sat2prn, sat2id, vnorm
from cssrlib.gnss import time2gpst, time2doy, epoch2time, time2str
from cssrlib.gnss import timeadd, timediff
from cssrlib.gnss import uGNSS as ug
from cssrlib.gnss import rSigRnx, rCST, sys2str
from cssrlib.peph import atxdec
from cssrlib.peph import peph, biasdec, apc2com
from cssrlib.cssr_has import cssr_has
from cssrlib.cssrlib import sCType
from cssrlib.cssrlib import sCSSRTYPE as sc
from cssrlib.pppssr import find_corr_idx
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

ac = 'COD0MGXFIN'

orbfile = '../data/{}_{:4d}{:03d}0000_01D_05M_ORB.SP3'\
    .format(ac, year, doy)

clkfile = '../data/{}_{:4d}{:03d}0000_01D_30S_CLK.CLK'\
    .format(ac, year, doy)

bsxfile = '../data/{}_{:4d}{:03d}0000_01D_01D_OSB.BIA'\
    .format(ac, year, doy)

if not exists(clkfile):
    orbfile = orbfile.replace('COD0OPSFIN', 'COD0OPSRAP')
    clkfile = clkfile.replace('COD0OPSFIN', 'COD0OPSRAP')
    bsxfile = bsxfile.replace('COD0OPSFIN', 'COD0OPSRAP')
if not exists(orbfile):
    orbfile = orbfile.replace('_05M_', '_15M_')

# Read Galile HAS corrections file
#
if doy == 223:
    file_has = '../data/doy223/223v_gale6.txt'
else:
    file_has = '../data/gale6_189e.txt'
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
cs = cssr_has("test_sisre_has_ssr.log")
cs.monlevel = 2

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
            if rs_ is None or dts is None:
                continue

            rs0[j, :] = rs_[0:3]
            vs0[j, :] = rs_[3:6]
            dts0[j] = dts_[0]

            """
            print("{} {} prec xyz [m] {:14.3f} {:14.3f}m {:14.3f} clk [ms] {:12.6f}"
                  .format(time2str(time), sat2id(sat),
                          rs0[j,0], rs0[j,1], rs0[j,2], dts0[j]*1e6))
            """

            # Broadcast ephemeris IODE, orbit and clock corrections
            #
            if cs.iodssr_c[sCType.ORBIT] == cs.iodssr:
                if sat not in cs.sat_n:
                    continue
                idx = cs.sat_n.index(sat)
            else:
                if cs.iodssr_c[sCType.ORBIT] == cs.iodssr_p:
                    if sat not in cs.sat_n_p:
                        continue
                    idx = cs.sat_n_p.index(sat)
                else:
                    continue

            iode = cs.lc[0].iode[idx]
            dorb = cs.lc[0].dorb[idx, :]  # radial,along-track,cross-track

            if cs.cssrmode == sc.BDS_PPP:  # consitency check for IOD corr
                if cs.lc[0].iodc[idx] == cs.lc[0].iodc_c[idx]:
                    dclk = cs.lc[0].dclk[idx]
                else:
                    if cs.lc[0].iodc[idx] == cs.lc[0].iodc_c_p[idx]:
                        dclk = cs.lc[0].dclk_p[idx]
                    else:
                        continue

            else:

                if cs.cssrmode == sc.GAL_HAS_SIS:  # HAS only

                    if cs.mask_id != cs.mask_id_clk:  # mask has changed
                        if sat not in cs.sat_n_p:
                            continue
                        idx = cs.sat_n_p.index(sat)

                else:

                    if cs.iodssr_c[sCType.CLOCK] == cs.iodssr:
                        if sat not in cs.sat_n:
                            continue
                        idx = cs.sat_n.index(sat)
                    else:
                        if cs.iodssr_c[sCType.CLOCK] == cs.iodssr_p:
                            if sat not in cs.sat_n_p:
                                continue
                            idx = cs.sat_n_p.index(sat)
                        else:
                            continue

                dclk = cs.lc[0].dclk[idx]

            if np.isnan(dclk) or np.isnan(dorb@dorb):
                continue

            mode = cs.nav_mode[sys]

            # Get position, velocity and clock offset from broadcast ephemerides
            #
            eph = findeph(nav.eph, time, sat, iode=iode, mode=mode)
            if eph is None:
                """
                print("ERROR: cannot find BRDC for {} mode {} iode {} at {}"
                      .format(sat2id(sat), mode, iode, time2str(time)))
                """
                continue

            rs[j, :], vs[j, :], dts[j] = eph2pos(time, eph, True)
            """
            print("{} {} brdc xyz [m] {:14.3f} {:14.3f} {:14.3f} clk [ms] {:12.6f} iode {:3d}"
                  .format(time2str(time), sat2id(sat),
                          rs[j, 0], rs[j, 1], rs[j, 2], dts[j, 0]*1e6, eph.iode))
            """

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
            rs[j, :] += apc2com(nav, sat, time, rs[j, :], sig0)

            # Select user reference signals
            #
            if sys == ug.GPS:
                sigs = (rSigRnx("GC1C"), rSigRnx("GC2W"))
            elif sys == ug.GAL:
                sigs = (rSigRnx("EC1C"), rSigRnx("EC5Q"))
            else:
                print("ERROR: invalid sytem {}".format(sys2str(sys)))
                continue

            freq = [s.frequency() for s in sigs]
            facs = (+freq[0]**2/(freq[0]**2-freq[1]**2),
                    -freq[1]**2/(freq[0]**2-freq[1]**2))

            # Get HAS biases
            #
            cbias = np.ones(len(sigs))*np.nan

            if cs.lc[0].cstat & (1 << sCType.CBIAS) == (1 << sCType.CBIAS):
                nsig, idx_n, kidx = find_corr_idx(cs, nav.nf, sCType.CBIAS,
                                                  sigs, sat)

                if nsig >= nav.nf:
                    cbias = cs.lc[0].cbias[idx_n][kidx]
                elif nav.monlevel > 1:
                    print("skip cbias for sat={:d}".format(sat))

                # - IS-QZSS-MDC-001 sec 5.5.3.3
                # - HAS SIS ICD sec 7.4, 7.5
                # - HAS IDD ICD sec 3.3.4
                if cs.cssrmode in [sc.GAL_HAS_IDD, sc.GAL_HAS_SIS, sc.QZS_MADOCA]:
                    cbias = -cbias

            # Get CODE biases
            #
            cbias_ = np.ones(len(sigs))*np.nan
            for i, sig in enumerate(sigs):
                cbias_[i] = bsx.getosb(sat, time, sig)*ns2m

            dcb = cbias[0]-cbias[1]
            dcb_ = cbias_[0]-cbias_[1]

            osbIF = facs[0]*cbias[0]+facs[1]*cbias[1]
            osbIF_ = facs[0]*cbias_[0]+facs[1]*cbias_[1]

            # Apply ionosphere-free bias correction to clock offsets
            #
            dts[j] -= osbIF/rCST.CLIGHT
            dts0[j] -= osbIF_/rCST.CLIGHT

            # Along-track, cross-track and radial conversion
            #
            ea = vnorm(vs[j, :])
            rc = np.cross(rs[j, :], vs[j, :])
            ec = vnorm(rc)
            er = np.cross(ea, ec)
            A = np.array([er, ea, ec])

            """
            print("{} {} dorb rac [m] {:14.3f} {:14.3f} {:14.3f} clk [ms] {:12.6f} iode {:3d}"
                  .format(time2str(time), sat2id(sat),
                          dorb[0], dorb[1], dorb[2], dclk/rCST.CLIGHT*1e6, iode))
            """

            # Convert orbit corrections from orbital frame to ECEF
            #
            dorb_e = dorb@A

            # Apply SSR correction
            #
            rs[j, :] -= dorb_e
            dts[j] += dclk/rCST.CLIGHT  # [m] -> [s]

            """
            print("{} {} hasc xyz [m] {:14.3f} {:14.3f} {:14.3f} clk [ms] {:12.6f}"
                  .format(time2str(time), sat2id(sat),
                          rs[j, 0], rs[j, 1], rs[j, 2], dclk))
            """

            # HAS vs. precise
            #
            d_rs[j, :] = (rs[j, :] - rs0[j, :])@A.T
            d_dts[j, 0] = dts[j, 0] - dts0[j, 0]

            print("{} {} diff rac [m] {:8.3f} {:8.3f} {:8.3f} "
                  "clk [m] {:12.6f} iode {:3d} "
                  "bias SSR [m] {} {:7.3f} {} {:7.3f} DCB {:7.3f} "
                  "bias COD [m] {} {:7.3f} {} {:7.3f} DCB {:7.3f} "
                  .format(time2str(time), sat2id(sat),
                          d_rs[j, 0], d_rs[j, 1], d_rs[j, 2],
                          d_dts[j, 0]*rCST.CLIGHT, iode,
                          sigs[0], cbias[0], sigs[1], cbias[1], dcb,
                          sigs[0], cbias_[0], sigs[1], cbias_[1], dcb_))

            orb_r[ne, sat-1] = d_rs[j, 0]
            orb_a[ne, sat-1] = d_rs[j, 1]
            orb_c[ne, sat-1] = d_rs[j, 2]
            clk[ne, sat-1] = d_dts[j, 0]*rCST.CLIGHT

    print()

    # Save output
    #
    t[ne] = timediff(time, t0)/60

    # Next time-step
    #
    time = timeadd(time, step)

fig = plt.figure(figsize=[7, 9])
fig.set_rasterized(True)

lbl_t = ['Radial [m]', 'Along [m]', 'Cross [m]', 'Clock [m]']

idx_G = np.arange(ug.GPSMIN, ug.GPSMIN+ug.GPSMAX)
idx_E = np.arange(ug.GALMIN, ug.GALMIN+ug.GALMAX)

plt.subplot(4, 1, 1)
plt.plot(t, orb_r[:, idx_G], 'k.', label='GPS')
plt.plot(t, orb_r[:, idx_E], 'b.', label='GAL')
plt.ylabel(lbl_t[0])
plt.grid()

plt.subplot(4, 1, 2)
plt.plot(t, orb_a[:, idx_G], 'k.', label='GPS')
plt.plot(t, orb_a[:, idx_E], 'b.', label='GAL')
plt.ylabel(lbl_t[1])
plt.grid()

plt.subplot(4, 1, 3)
plt.plot(t, orb_c[:, idx_G], 'k.', label='GPS')
plt.plot(t, orb_c[:, idx_E], 'b.', label='GAL')
plt.ylabel(lbl_t[2])
plt.grid()

plt.subplot(4, 1, 4)
plt.plot(t, clk[:, idx_G], 'k.', label='GPS')
plt.plot(t, clk[:, idx_E], 'b.', label='GAL')
plt.ylabel(lbl_t[3])
plt.grid()

plotFileFormat = 'eps'
plotFileName = '.'.join(('test_sisre_has', plotFileFormat))

plt.savefig(plotFileName, format=plotFileFormat, bbox_inches='tight', dpi=300)
# plt.show()
