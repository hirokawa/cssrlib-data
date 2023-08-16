import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from cssrlib.rinex import rnxdec
from cssrlib.gnss import Nav, epoch2time, prn2sat, uGNSS, sat2prn,\
    timeadd, ecef2pos
from cssrlib.ephemeris import findeph, eph2pos

rnx_ver = 4

if rnx_ver == 3:  # RINEX 3
    navfile = '../data/30340780.21q'
    t0 = epoch2time([2021, 3, 19, 0, 0, 0])
    sys_ref = uGNSS.QZS
    prn_ref = 194
    mode_ref = 0  # 0: LNAV, 1: CNAV, 2: CNAV2
elif rnx_ver == 4:  # RINEX 4
    navfile = '../data/BRD400DLR_S_20231890000_01D_MN.rnx'
    t0 = epoch2time([2023, 7, 8, 4, 0, 0])
    sys_ref = uGNSS.BDS
    prn_ref = 35
    mode_ref = 1  # 0: LNAV, 1: CNAV, 2: CNAV2

dec = rnxdec()
nav = dec.decode_nav(navfile, Nav())

n = 24*3600//300
flg_plot = True


if True:
    t = t0
    sat = prn2sat(sys_ref, prn_ref)
    eph = findeph(nav.eph, t, sat, mode=mode_ref)
    if eph is not None:
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
        for i in range(n):
            t = timeadd(t0, i*300)
            if eph is None:
                continue
            rs, dts = eph2pos(t, eph)
            pos[i, :] = ecef2pos(rs)
            pos[i, 0] = np.rad2deg(pos[i, 0])
            pos[i, 1] = np.rad2deg(pos[i, 1])

        plt.plot(pos[:, 1], pos[:, 0], 'm-', transform=ccrs.Geodetic())
    plt.show()
