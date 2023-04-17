from cssrlib.gnss import rSigRnx, uGNSS, uTYP, uSIG
from cssrlib.gnss import sat2id, sat2prn
from cssrlib.rinex import rnxdec

obsfile = '../data/SEPT078M.21O'

sigs = [rSigRnx(uGNSS.GPS, uTYP.C, uSIG.L1C),
        rSigRnx(uGNSS.GPS, uTYP.C, uSIG.L2W),
        rSigRnx(uGNSS.GPS, uTYP.L, uSIG.L1C),
        rSigRnx(uGNSS.GPS, uTYP.L, uSIG.L2W)]

dec = rnxdec()
dec.setSignals(sigs)

nep = 1
if dec.decode_obsh(obsfile) >= 0:

    for ne in range(nep):

        obs = dec.decode_obs()

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

    dec.fobs.close()
