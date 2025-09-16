import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from cssrlib.rinex import rnxdec
from cssrlib.gnss import Nav, epoch2time, prn2sat, uGNSS, sat2prn, \
    timeadd, ecef2pos
from cssrlib.ephemeris import findeph, eph2pos, geph2pos

rnx_ver = 4

if rnx_ver == 3:  # RINEX 3
    navfile = '../data/brdc/BRDC00IGS_R_20231890000_01D_MN.rnx'
    t0 = epoch2time([2023, 7, 8, 0, 0, 0])
    sys_ref = uGNSS.QZS
    prn_ref = 194
    mode_ref = 0  # 0: LNAV, 1: CNAV, 2: CNAV2
elif rnx_ver == 4:  # RINEX 4
    navfile = '../data/brdc/BRD400DLR_S_20231890000_01D_MN.rnx'
    t0 = epoch2time([2023, 7, 8, 0, 0, 0])
    sys_ref = uGNSS.GLO
    prn_ref = 7
    mode_ref = 0  # 0: LNAV, 1: CNAV, 2: CNAV2

dec = rnxdec()
nav = dec.decode_nav(navfile, Nav())

step = 300
n = 24*3600//step
flg_plot = True


if True:
    t = t0
    sat = prn2sat(sys_ref, prn_ref)
    if sys_ref == uGNSS.GLO:
        eph = findeph(nav.geph, t, sat, mode=mode_ref)
    else:
        eph = findeph(nav.eph, t, sat, mode=mode_ref)
    if eph is not None:
        if sys_ref == uGNSS.GLO:
            rs, vs, dts = geph2pos(t, eph, True)
        else:
            rs, vs, dts = eph2pos(t, eph, True)

if flg_plot:

    lon0 = 135
    plt.figure(figsize=(6, 6))
    ax = plt.axes(projection=ccrs.Orthographic(central_longitude=lon0,
                                               central_latitude=0))
    # ax.coastlines(resolution='50m')
    ax.stock_img()
    ax.gridlines()

    pos = np.zeros((n, 3))

    for k in range(uGNSS.MAXSAT):
        sat = k+1
        sys, prn = sat2prn(sat)
        if sys != sys_ref:
            continue
        """
        if prn != prn_ref:
            continue
        """
        for i in range(n):
            t = timeadd(t0, i*step)
            if sys_ref == uGNSS.GLO:
                eph = findeph(nav.geph, t, sat, mode=mode_ref)
            else:
                eph = findeph(nav.eph, t, sat, mode=mode_ref)
            if eph is None:
                continue
            if sys_ref == uGNSS.GLO:
                rs, dts = geph2pos(t, eph)
            else:
                rs, dts = eph2pos(t, eph)
            pos[i, :] = ecef2pos(rs)
            pos[i, 0] = np.rad2deg(pos[i, 0])
            pos[i, 1] = np.rad2deg(pos[i, 1])

        plt.plot(pos[:, 1], pos[:, 0], 'm-', transform=ccrs.Geodetic())
    plt.show()
