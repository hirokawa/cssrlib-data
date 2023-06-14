import numpy as np
from cssrlib.gnss import Nav, ecef2pos, geodist, satazel, timediff, uGNSS
from cssrlib.gnss import rSigRnx
from cssrlib.ephemeris import findeph, eph2pos
from cssrlib.plot import skyplot, plot_elv
from cssrlib.rinex import rnxdec

navfile = '../data/SEPT078M.21P'
obsfile = '../data/SEPT078M.21O'

# Define signals to be processed
#
sigs = [rSigRnx("GC1C"), rSigRnx("GC2W"),
        rSigRnx("GL1C"), rSigRnx("GL2W"),
        rSigRnx("GS1C"), rSigRnx("GS2W"),
        rSigRnx("EC1C"), rSigRnx("EC5Q"),
        rSigRnx("EL1C"), rSigRnx("EL5Q"),
        rSigRnx("ES1C"), rSigRnx("ES5Q"),
        rSigRnx("JC1C"), rSigRnx("JC5Q"),
        rSigRnx("JL1C"), rSigRnx("JL5Q"),
        rSigRnx("JS1C"), rSigRnx("JS5Q")]


dec = rnxdec()
dec.setSignals(sigs)
nav = Nav()
nav = dec.decode_nav(navfile, nav)

nep = 60*15
elv = np.ones((nep, uGNSS.MAXSAT))*np.nan
azm = np.ones((nep, uGNSS.MAXSAT))*np.nan

t = np.zeros(nep)*np.nan
if dec.decode_obsh(obsfile) >= 0:
    rr = dec.pos
    pos = ecef2pos(rr)
    for ne in range(nep):
        obs = dec.decode_obs()
        if ne == 0:
            t0 = obs.t
        t[ne] = timediff(obs.t, t0)
        for k, sat in enumerate(obs.sat):
            eph = findeph(nav.eph, obs.t, sat)
            rs, dts = eph2pos(obs.t, eph)
            r, e = geodist(rs, rr)
            azm[ne, sat-1], elv[ne, sat-1] = satazel(pos, e)

    dec.fobs.close()

plot_elv(t, elv)
skyplot(azm, elv)
