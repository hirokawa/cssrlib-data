"""
 kinematic test for stand-alone positioning
"""
import matplotlib.pyplot as plt
import numpy as np
from cssrlib.rinex import rnxdec
from cssrlib.gnss import ecef2pos, timediff, ecef2enu, pos2ecef, xyz2enu
from cssrlib.gnss import rSigRnx, Nav
from cssrlib.pntpos import stdpos


rr0 = xyz_ref = pos2ecef([35.342058098, 139.521986657, 47.5515], True)
pos_ref = ecef2pos(xyz_ref)
nep = 360

E = xyz2enu(pos_ref)

navfile = '../data/doy2021-265/SEPT2650.21P'
obsfile = '../data/doy2021-265/SEPT265G.21O'

# Define signals to be processed
#
sigs = [rSigRnx("GC1C"), rSigRnx("EC1C"),
        rSigRnx("GL1C"), rSigRnx("EL1C"),
        rSigRnx("GS1C"), rSigRnx("ES1C")]

rnx = rnxdec()
rnx.setSignals(sigs)

nav = Nav(nf=1)
nav.pmode = 1  # 0: static, 1: kinematic
nav = rnx.decode_nav(navfile, nav)

t = np.zeros(nep)
enu = np.zeros((nep, 3))
dop = np.zeros((nep, 4))
nsat = np.zeros(nep, dtype=int)

if rnx.decode_obsh(obsfile) >= 0:

    # Auto-substitute signals
    #
    rnx.autoSubstituteSignals()

    # Initialize position
    #
    std = stdpos(nav, rr0, 'test_stdpos2.log')
    nav.elmin = np.deg2rad(5.0)

    sol = np.zeros((nep, nav.nx))

    for ne in range(nep):

        obs = rnx.decode_obs()
        if ne == 0:
            t0 = nav.t = obs.t
        t[ne] = timediff(obs.t, t0)

        std.process(obs, cs=None)

        sol[ne, :] = nav.x
        enu[ne, :] = ecef2enu(pos_ref, sol[ne, 0:3]-xyz_ref)
        dop[ne, :] = std.dop

    rnx.fobs.close()


if True:
    lbl_t = ['east [m]', 'north [m]', 'up [m]']
    fig = plt.figure(figsize=(6, 10))

    for k in range(3):
        plt.subplot(3, 1, k+1)
        plt.plot(t, enu[:, k])
        plt.ylabel(lbl_t[k])
        plt.xlabel('time[s]')
        plt.grid()
    plt.show()

    venu = sol[:, 3:6]@E.T

    plt.figure()
    plt.plot(t, venu)
    plt.ylabel('velocity [m/s]')
    plt.xlabel('time[s]')
    plt.legend(['east', 'north', 'up'])
    plt.grid()
    plt.axis([0, nep, -10, 10])
    plt.show()

    sol[0, 7] = np.nan
    plt.figure()
    plt.subplot(211)
    plt.plot(t, sol[:, 6]-sol[0, 6])
    plt.ylabel('clock bias [m]')
    plt.grid()
    plt.subplot(212)
    plt.plot(t, sol[:, 7])
    plt.ylabel('clock drift [m/s]')
    plt.xlabel('time[s]')
    plt.grid()
    plt.show()

if True:
    plt.figure()
    plt.plot(enu[:, 0], enu[:, 1], '-', color='gray')
    plt.plot(enu[:, 0], enu[:, 1], 'm.', markersize=1)
    plt.xlabel('easting[m]')
    plt.ylabel('northing[m]')
    plt.grid()
    plt.axis('equal')
    plt.show()

    plt.figure()
    plt.plot(t, dop[:, 1:])
    plt.legend(['pdop', 'hdop', 'vdop'])
    plt.grid()
    plt.axis([0, nep, 0, 3])
    plt.xlabel('time[s]')
    plt.show()
