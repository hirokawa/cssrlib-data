"""
 interopetabity test for SSR Integrity messages MT11,12,13,MT54.*

 @author Rui Hirokawa
"""

import os
import copy
from cssrlib.gnss import uGNSS, prn2sat
from cssrlib.rtcm import rtcm, rtcme, Integrity
from random import randint, seed, sample
from binascii import unhexlify


def read_asc(file):
    b = bytearray()
    with open(file) as fh:
        for line in fh:
            b += unhexlify(''.join(line.split()))

    return b


def read_bin(file):
    with open(file, 'rb') as fh:
        return fh.read()
    return None


def gen_data(mt, sys_t, svid_t):
    """ generate random message data """
    intr = Integrity()

    if mt in (11, 12, 13):  # SSR integrity messages

        for sys in sys_t:
            intr.iod_sys[sys] = randint(0, 3)
            intr.nid[sys] = randint(0, 255)  # network id DFi071 0-255
            intr.flag[sys] = {}
            for svid in svid_t[sys]:
                sat = prn2sat(sys, svid)
                intr.flag[sys][sat] = randint(0, 2)
                # integrity flag DFi068 0:use,1:DNU,2:not monitored,3:reserved

        # issue of GNSS satellite mask DFi010
        intr.pid = randint(0, 4095)  # provider id DFi027 (0-4095)
        intr.tow = randint(0, 604799)*1e-3  # tow
        intr.vp = intr.vp_tbl[randint(0, 15)]  # validity period DFi065 (0-15)
        intr.uri = randint(0, 65535)*0.1  # update rate interval DFi067

        intr.pidssr = randint(0, 65535)  # SSR Provider ID DFi078 (0-65535)
        intr.sidssr = randint(0, 15)  # SSR Solution Type DFi076 (0-15)
        intr.iodssr = randint(0, 15)  # SSR IOD DFi077 (0-15)

    return intr


def gen_sat_list(sys_t, prn_rng_t):
    prn_t = {}
    for sys in sys_t:
        nsat = randint(1, nsatmax)
        prn_t[sys] = []
        smin = prn_rng_t[sys][0]
        smax = prn_rng_t[sys][1]
        prn_t[sys] = sample(range(smin, smax), nsat)
        prn_t[sys].sort()
    return prn_t


def write_rtcm(file_rtcm, msg_t, intr, nep=1):

    cs = rtcme()
    cs.integ = copy.deepcopy(intr)
    cs.iodssr = intr.iodssr
    cs.pid = intr.pidssr
    cs.sid = intr.sidssr

    fc = open(file_rtcm, 'wb')
    if not fc:
        print("RTCM message file cannot open.")
    k = 0
    msg = bytearray(maxlen)
    buff = bytearray(1024)

    for ne in range(nep):
        cs.msgtype = msg_t[ne % len(msg_t)]
        i = cs.encode(buff)

        len_ = (i+7)//8
        if k+len_+6 >= maxlen:
            break

        cs.set_body(msg, buff, k, len_)
        cs.set_sync(msg, k)
        cs.set_len(msg, k, len_)
        cs.set_checksum(msg, k)
        k += cs.dlen

    fc.write(msg[:k])
    fc.close()

    return msg[:k]


def decode_rtcm(msg, intr=None, nep=1, logfile=None, maxlen=1024):
    cs = rtcm(foutname=logfile)
    cs.monlevel = 2

    k = 0
    for ne in range(nep):

        # loop for RTCM
        while k < maxlen:
            stat = cs.sync(msg, k)
            if stat is False:
                k += 1
                continue
            if not cs.checksum(msg, k, maxlen):
                print("checksum failed.")
                # k += 1
                # continue

            cs.decode(msg[k:k+cs.len+3])
            k += cs.dlen

            if intr is not None:
                print(f"pid: {cs.integ.pid == intr.pid}")
                print(f"tow: {cs.integ.tow == intr.tow}")
                print(f"flag: {cs.integ.flag == intr.flag}")
                # print(f"iodsys: {cs.integ.iod_sys == intr.iod_sys}")

                print(f"vp: {cs.integ.vp == intr.vp}")
                print(f"uri: {cs.integ.uri == intr.uri}")

                print(f"pidssr: {cs.pid == intr.pidssr}")
                print(f"sidssr: {cs.sid == intr.sidssr}")
                print(f"iodssr: {cs.iodssr == intr.iodssr}")
    return cs


def read_rtcm(file_rtcm, intr, nep=1, logfile=None):
    """ read test script for SC-134 messages """

    fc = open(file_rtcm, 'rb')
    if not fc:
        print("RTCM message file cannot open.")

    blen = os.path.getsize(file_rtcm)
    msg = fc.read(blen)
    maxlen = len(msg)-5
    fc.close()

    return decode_rtcm(msg, intr, nep, logfile, maxlen)


if __name__ == "__main__":
    bdir = '../data/sc134/msg/'
    flg_sim = False

    if flg_sim:  # generate test data
        file_rtcm = bdir+'test.rtcm'
        file_log = bdir+'test.log'
        nep = 1
        maxlen = 1024
        nsatmax = 10

        seed_ = 1
        # parameters
        # msg_t = [11, 12, 13]
        msg_t = [11]
        # msg_t = [12]
        # msg_t = [13]

        # 0:GPS,1:GLO,2:GAL,3:BDS,4:QZS,5:IRN
        sys_t = [uGNSS.GPS, uGNSS.GLO, uGNSS.GAL, uGNSS.QZS]

        # GNSS satellite mask DFi009
        prn_rng_t = {uGNSS.GPS: [1, 32],  # Table 8.5-1
                     uGNSS.GLO: [1, 27],  # Table 8.5-3
                     uGNSS.GAL: [1, 36],  # Table 8.5-5
                     uGNSS.BDS: [1, 63],
                     uGNSS.QZS: [193, 209],
                     uGNSS.IRN: [1, 14]}

        mt = msg_t[0]

        seed(seed_)
        prn_t = gen_sat_list(sys_t, prn_rng_t)  # generate random sat list
        intr = gen_data(mt, sys_t, prn_t)  # generate random message data
        msg = write_rtcm(file_rtcm, msg_t, intr, nep)
        cs = read_rtcm(file_rtcm, intr, nep, logfile=file_log)

    else:  # decode using sample dataset (*.bin)

        file_rtcm = ['MT54_9', 'MT54_10_DFi209=0',
                     'MT54_10_DFi209=1', 'MT54_10_DFi209=2']

        # file_rtcm = ['sampledataMT03-04-05-06-07']

        # msg = read_asc(file_asc)
        for f in file_rtcm:
            file_log = bdir+f+'.dlg'
            msg = read_bin(bdir+f+'.bin')
            decode_rtcm(msg, logfile=file_log, maxlen=len(msg))
