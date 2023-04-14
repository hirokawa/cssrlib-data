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

    for sig in sigs:
        if sig.gns not in dec.sig_tab:
            dec.sig_tab.update({sig.gns: {}})
        if sig.typ not in dec.sig_tab[sig.gns]:
            dec.sig_tab[sig.gns].update({sig.typ: []})
        if sig.sig not in dec.sig_tab[sig.gns][sig.typ]:
            dec.sig_tab[sig.gns][sig.typ].append(sig.sig)

    for ne in range(nep):
        obs = dec.decode_obs()
        print(obs)

    dec.fobs.close()
