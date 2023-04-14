"""
 test of RINEX decoder
"""

from cssrlib.rinex import rnxdec
from cssrlib.gnss import uGNSS, uTYP, uSIG, rSigRnx

obsfile = '/home/andre/GNSS_OBS/IGS/DAILY/2021/078/CHOF00JPN_S_20210780000_01D_30S_MO.rnx'

dec = rnxdec()

nep = 1

if dec.decode_obsh(obsfile) >= 0:

    sigs = [rSigRnx(uGNSS.GPS, uTYP.C, uSIG.L1C),
            rSigRnx(uGNSS.GPS, uTYP.C, uSIG.L2W),
            rSigRnx(uGNSS.GPS, uTYP.L, uSIG.L1C),
            rSigRnx(uGNSS.GPS, uTYP.L, uSIG.L2W)]
    dec.setSignals(sigs)

    for ne in range(nep):
        obs = dec.decode_obs()
        print(obs)

    dec.fobs.close()
