"""
 test of RINEX decoder
"""

from cssrlib.rinex import rnxdec
from cssrlib.gnss import uTYP,  rSigRnx
from cssrlib.gnss import sat2id, sat2prn
from cssrlib.gnss import time2str

from os.path import expanduser

#obsfile = '../data/SEPT078M.21O'
obsfile = '../data/SEPT1890.23O'
#obsfile = '../data/3034078M.21O'

sigs = [rSigRnx("GC1C"), rSigRnx("EC1C"), rSigRnx("CC2I"), rSigRnx("JC1C"),
        rSigRnx("GC2W"), rSigRnx("EC5Q"), rSigRnx("CC6I"), rSigRnx("JC2L"),
        rSigRnx("GC5Q"), rSigRnx("EC6C"), rSigRnx("CC5P"), rSigRnx("JC5Q"),
        rSigRnx("GL1C"), rSigRnx("EL1C"), rSigRnx("CL2I"), rSigRnx("JL1C"),
        rSigRnx("GL2W"), rSigRnx("EL5Q"), rSigRnx("CL6I"), rSigRnx("JL2L"),
        rSigRnx("GL5Q"), rSigRnx("EL6C"), rSigRnx("CL5P"), rSigRnx("JL5Q"), ]

dec = rnxdec()
dec.setSignals(sigs)

nep = 1
if dec.decode_obsh(expanduser(obsfile)) >= 0:

    dec.autoSubstituteSignals()

    for ne in range(nep):

        obs = dec.decode_obs()

        print(time2str(obs.t))
        for i, sat in enumerate(obs.sat):

            txt = "{} ".format(sat2id(sat))

            sys, _ = sat2prn(sat)
            sigs = dec.getSignals(sys, uTYP.C)
            for j, sig in enumerate(sigs):
                txt += "{} {:13.3f}  ".format(sig.str(), obs.P[i, j])
            sigs = dec.getSignals(sys, uTYP.L)
            for j, sig in enumerate(sigs):
                txt += "{} {:13.3f}  ".format(sig.str(), obs.L[i, j])

            print(txt)

        print()

    dec.fobs.close()
