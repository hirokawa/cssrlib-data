# -*- coding: utf-8 -*-
"""
Created on Sun Aug 22 21:01:49 2021

@author: ruihi
"""

import numpy as np
import matplotlib.pyplot as plt
from cssrlib.gnss import Nav
from cssrlib.gnss import rSigRnx
from cssrlib.gnss import time2gpst, epoch2time, time2epoch, timeadd
from cssrlib.gnss import tropmodel, tropmapf
from cssrlib.gnss import pos2ecef
from cssrlib.peph import atxdec, searchpcv, antModelRx

bdir = '../data/'
atxfile = bdir+"igs14.atx"

igfont = {'family': 'Meiryo'}

t = epoch2time([2021, 3, 19, 12, 0, 0])
ep = time2epoch(t)
week, tow = time2gpst(t)
t1 = timeadd(t, 300)

# Receiver position
#
lat = np.deg2rad(45)
lon = np.deg2rad(11)
hgt = 0

rr = pos2ecef(np.array([lat, lon, hgt]))

atx = atxdec()
atx.readpcv(atxfile)

# Retrieve station antennas
#
antr = "{:16s}{:4s}".format("JAVRINGANT_DM", "SCIS")
antb = "{:16s}{:4s}".format("TRM59800.80", "NONE")

nav = Nav()
nav.rcv_ant = searchpcv(atx.pcvr, antr, t1)
nav.rcv_ant_b = searchpcv(atx.pcvr, antb, t1)

sigs = [rSigRnx("GC1C"), rSigRnx("GC2W")]

el_t = np.arange(10, 90)
ofst_r = np.zeros((len(el_t), len(sigs)))
ofst_b = np.zeros((len(el_t), len(sigs)))
for k, el in enumerate(el_t):
    el = np.deg2rad(el)
    e = np.array([0, np.cos(el), np.sin(el)])
    ofst_r[k, :] = antModelRx(nav, rr, e, sigs, 0)
    ofst_b[k, :] = antModelRx(nav, rr, e, sigs, 1)

flg_ant = True  # False
flg_trop = True  # True

plt.figure()
if flg_ant:
    plt.plot(el_t, ofst_b[:, 0]*100, label='Trimble TRM59800.80')
    plt.plot(el_t, ofst_r[:, 0]*100, label='JAVAD RINGANT')
    plt.grid()
    plt.legend()
    plt.xlabel('elevation[deg]')
    plt.ylabel('range correction for antenna offset [cm]')
    plt.show()

if flg_trop:
    ep = [2021, 4, 1, 0, 0, 0]
    t = epoch2time(ep)
    el_t = np.arange(0.01, np.pi/2, 0.01)
    n = len(el_t)
    trop = np.zeros(n)
    trop_hs = np.zeros(n)
    trop_wet = np.zeros(n)
    lat_t = [45]
    for lat in lat_t:
        pos = [np.deg2rad(lat), 0, 0]
        for k, el in enumerate(el_t):
            ths, twet, z = tropmodel(t, pos, el)
            mapfh, mapfw = tropmapf(t, pos, el)
            trop_hs[k] = mapfh*ths
            trop_wet[k] = mapfw*twet

        trop = trop_hs+trop_wet
        plt.plot(np.rad2deg(el_t), trop)
        plt.plot(np.rad2deg(el_t), trop_hs)
        plt.plot(np.rad2deg(el_t), trop_wet)

    plt.grid()
    plt.axis([0, 90, 0, 10])
    plt.legend(['total', 'dry', 'wet'], prop=igfont)
    plt.xlabel('Elevation angle [deg]', **igfont)
    plt.ylabel('Tropospheric delay [m]', **igfont)
    plt.show()
