"""
RTCM 3 messages decoder

[1] RTCM Standard 10403.4 Differential GNSS Services -
    Version 3 with Amendment 1, November, 2024

@author: Rui Hirokawa
"""

import argparse
from glob import glob
import multiprocessing as mp
import numpy as np
import os
from pathlib import Path

from cssrlib.gnss import uGNSS, uTYP, rSigRnx, Obs, gtime_t, timediff
from cssrlib.rawnav import rcvDec, rcvOpt

from cssrlib.rtcm import rtcm


class rtcmDec(rcvDec):
    tow = -1
    week = -1
    sat = []
    sys = []
    prn = []
    nsat = 0
    nmax = 96
    nsigmax = 7

    pr = []
    cp = []
    dp = []

    def __init__(self, opt=None, prefix='', gnss_t='GECJ'):
        super().__init__(opt, prefix, gnss_t)

        self.sig_tab = {}

        if 'G' in gnss_t:
            self.sig_tab[uGNSS.GPS] = {
                uTYP.C: [rSigRnx('GC1C'), rSigRnx('GC1W'), rSigRnx('GC2W'),
                         rSigRnx('GC2L'), rSigRnx('GC5Q'), rSigRnx('GC1L')],
                uTYP.L: [rSigRnx('GL1C'), rSigRnx('GL1W'), rSigRnx('GL2W'),
                         rSigRnx('GL2L'), rSigRnx('GL5Q'), rSigRnx('GL1L')],
                uTYP.S: [rSigRnx('GS1C'), rSigRnx('GS1W'), rSigRnx('GS2W'),
                         rSigRnx('GS2L'), rSigRnx('GS5Q'), rSigRnx('GS1L')],
            }
        if 'R' in gnss_t:
            self.sig_tab[uGNSS.GLO] = {
                uTYP.C: [rSigRnx('RC1C'), rSigRnx('RC1P'), rSigRnx('RC2C'),
                         rSigRnx('RC2P'), rSigRnx('RC3X')],
                uTYP.L: [rSigRnx('RL1C'), rSigRnx('RL1P'), rSigRnx('RL2C'),
                         rSigRnx('RL2P'), rSigRnx('RL3X')],
                uTYP.S: [rSigRnx('RS1C'), rSigRnx('RS1P'), rSigRnx('RS2C'),
                         rSigRnx('RS2P'), rSigRnx('RS3X')],
            }
        if 'E' in gnss_t:
            self.sig_tab[uGNSS.GAL] = {
                uTYP.C: [rSigRnx('EC1C'), rSigRnx('EC5Q'), rSigRnx('EC7Q'),
                         rSigRnx('EC8Q'), rSigRnx('EC6C')],
                uTYP.L: [rSigRnx('EL1C'), rSigRnx('EL5Q'), rSigRnx('EL7Q'),
                         rSigRnx('EL8Q'), rSigRnx('EL6C')],
                uTYP.S: [rSigRnx('ES1C'), rSigRnx('ES5Q'), rSigRnx('GS7Q'),
                         rSigRnx('ES8Q'), rSigRnx('ES6C')],
            }
        if 'C' in gnss_t:
            self.sig_tab[uGNSS.BDS] = {
                uTYP.C: [rSigRnx('CC1P'), rSigRnx('CC2I'), rSigRnx('CC5P'),
                         rSigRnx('CC6I'), rSigRnx('CC7D'), rSigRnx('CC7I')],
                uTYP.L: [rSigRnx('CL1P'), rSigRnx('CL2I'), rSigRnx('CL5P'),
                         rSigRnx('CL6I'), rSigRnx('CL7D'), rSigRnx('CL7I')],
                uTYP.S: [rSigRnx('CS1P'), rSigRnx('CS2I'), rSigRnx('CS5P'),
                         rSigRnx('CS6I'), rSigRnx('CS7D'), rSigRnx('CS7I')],
            }
        if 'J' in gnss_t:
            self.sig_tab[uGNSS.QZS] = {
                uTYP.C: [rSigRnx('JC1C'), rSigRnx('JC1L'), rSigRnx('JC2L'),
                         rSigRnx('JC5Q'), rSigRnx('JC6X'), rSigRnx('JC1E')],
                uTYP.L: [rSigRnx('JL1C'), rSigRnx('JL1L'), rSigRnx('JL2L'),
                         rSigRnx('JL5Q'), rSigRnx('JL6X'), rSigRnx('JL1E')],
                uTYP.S: [rSigRnx('JS1C'), rSigRnx('JS1L'), rSigRnx('JS2L'),
                         rSigRnx('JS5Q'), rSigRnx('JS6X'), rSigRnx('JS1E')],
            }
        if 'S' in gnss_t:
            self.sig_tab[uGNSS.SBS] = {
                uTYP.C: [rSigRnx('SC1C'), rSigRnx('SC5X')],
                uTYP.L: [rSigRnx('SL1C'), rSigRnx('SL5X')],
                uTYP.S: [rSigRnx('SS1C'), rSigRnx('SS5X')],
            }
        if 'I' in gnss_t:
            self.sig_tab[uGNSS.IRN] = {
                uTYP.C: [rSigRnx('IC5A'), rSigRnx('IC1X')],
                uTYP.L: [rSigRnx('IL5A'), rSigRnx('IL1X')],
                uTYP.S: [rSigRnx('IS5A'), rSigRnx('IS1X')],
            }

        self.rtcm = rtcm()
        self.time_p = gtime_t()
        self.obs = None

        if opt is not None:
            self.init_param(opt=opt, prefix=prefix)

    def add_obs(self, obs):
        self.obs.sat = np.hstack((self.obs.sat, obs.sat))
        nsat = len(obs.sat)

        sys = list(obs.sig)[0]
        if sys not in self.obs.sig:
            self.obs.sig[sys] = obs.sig[sys]

        obs_ = Obs()
        obs_.P = np.zeros((nsat, self.nsig[uTYP.C]))
        obs_.L = np.zeros((nsat, self.nsig[uTYP.L]))
        obs_.S = np.zeros((nsat, self.nsig[uTYP.S]))
        obs_.D = np.zeros((nsat, self.nsig[uTYP.D]))
        obs_.lli = np.zeros((nsat, self.nsig[uTYP.C]), dtype=int)

        obs.P[np.isnan(obs.L)] = 0.0
        obs.L[np.isnan(obs.L)] = 0.0

        for k, sig in enumerate(obs.sig[sys][uTYP.L]):
            if sig in self.sig_tab[sys][uTYP.L]:
                idx = self.sig_tab[sys][uTYP.L].index(sig)
                obs_.P[:, idx] = obs.P[:, k]
                obs_.L[:, idx] = obs.L[:, k]
                obs_.S[:, idx] = obs.S[:, k]
                obs_.D[:, idx] = obs.D[:, k]
                obs_.lli[:, idx] = obs.lli[:, k]

        self.obs.P = np.vstack((self.obs.P, obs_.P))
        self.obs.L = np.vstack((self.obs.L, obs_.L))
        self.obs.S = np.vstack((self.obs.S, obs_.S))
        self.obs.D = np.vstack((self.obs.D, obs_.D))
        self.obs.lli = np.vstack((self.obs.lli, obs_.lli))

    def init_obs(self, time):
        nsig_max = 0
        for s in self.sig_tab:
            if len(self.sig_tab[s][uTYP.L]) > nsig_max:
                nsig_max = len(self.sig_tab[s][uTYP.L])

        self.nsig[uTYP.C] = nsig_max
        self.nsig[uTYP.L] = nsig_max
        self.nsig[uTYP.D] = nsig_max
        self.nsig[uTYP.S] = nsig_max

        self.obs = Obs()
        self.obs.time = time
        self.obs.P = np.empty((0, nsig_max))
        self.obs.L = np.empty((0, nsig_max))
        self.obs.D = np.empty((0, nsig_max))
        self.obs.S = np.empty((0, nsig_max))
        self.obs.lli = np.empty((0, nsig_max), dtype=int)
        self.obs.sat = np.empty(0, dtype=int)
        self.obs.sig = {}

    def decode(self, buff, len_, sys=[], prn=[]):

        _, obs, eph, geph, seph = self.rtcm.decode(buff, len_)

        if self.flg_rnxobs and obs is not None:
            self.time = obs.time
            self.re.rnx_obs_header(obs.time, self.fh_rnxobs)

            if timediff(self.time, self.time_p) != 0.0:
                if self.obs is not None:
                    self.re.rnx_obs_body(self.obs, self.fh_rnxobs)
                self.init_obs(obs.time)

            self.add_obs(obs)

            self.time_p = self.time

        if eph is not None:
            self.re.rnx_nav_body(eph, self.fh_rnxnav)

        if geph is not None:
            self.re.rnx_gnav_body(geph, self.fh_rnxnav)

        if seph is not None:
            self.re.rnx_snav_body(seph, self.fh_rnxnav)


def decode(f, opt, args):

    print("Decoding {}".format(f))

    bdir, fname = os.path.split(f)

    prefix = fname[4:].removesuffix('.rtcm3')+'_'
    prefix = str(Path(bdir) / prefix) if bdir else prefix
    rtcmdec = rtcmDec(opt=opt, prefix=prefix, gnss_t=args.gnss)
    rtcmdec.monlevel = 1

    rtcmdec.rtcm.week = args.weekref

    path = str(Path(bdir) / fname) if bdir else fname
    blen = os.path.getsize(path)
    with open(path, 'rb') as f:
        msg = f.read(blen)
        maxlen = len(msg)-5
        # maxlen = 400000
        k = 0
        for _ in range(maxlen):
            if k > maxlen:
                break
            stat = rtcmdec.rtcm.sync(msg, k)
            if stat is False:
                k += 1
                continue
            if not rtcmdec.rtcm.checksum(msg, k, maxlen):
                k += 1
                continue

            len_ = rtcmdec.rtcm.len+3
            rtcmdec.decode(msg[k:k+len_], len_)
            k += rtcmdec.rtcm.dlen

    rtcmdec.file_close()


def main():

    # Parse command line arguments
    #
    parser = argparse.ArgumentParser(description="RTCM3 converter")

    # Input file and folder
    #
    parser.add_argument("inpFileName",
                        help="Input RTCM3 file(s) (wildcards allowed)")

    parser.add_argument("--receiver", default='unknown',
                        help="Receiver type [unknown]")
    parser.add_argument("--antenna", default='unknown',
                        help="Antenna type [unknown]")

    parser.add_argument("-g", "--gnss", default='GRECJ',
                        help="GNSS [GRECJ]")

    parser.add_argument("--weekref", type=int, required=True,
                        help="GPS week number reference [unknown]")

    parser.add_argument("-j", "--jobs", default=int(mp.cpu_count() / 2),
                        type=int, help='Max. number of parallel processes')

    # Retrieve all command line arguments
    #
    args = parser.parse_args()

    opt = rcvOpt()

    opt.flg_rnxobs = True
    opt.flg_rnxnav = True

    # Start processing pool
    #
    with mp.Pool(processes=args.jobs) as pool:
        pool.starmap(decode, [(f, opt, args) for f in glob(args.inpFileName)])


# Call main function
#
if __name__ == "__main__":
    main()
