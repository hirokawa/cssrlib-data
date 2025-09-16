"""
 interopetabity test for SSR Integrity messages MT11,12,13

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
        intr.tow = randint(0, 604799)  # tow
        intr.vp = randint(0, 15)  # validity period DFi065 (0-15)
        intr.uri = randint(0, 65535)  # update rate interval DFi067

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
                print(cs.integ.pid == intr.pid)
                print(cs.integ.tow == intr.tow)
                print(cs.integ.flag == intr.flag)
                print(cs.integ.nid == intr.nid)
                print(cs.integ.iod_sys == intr.iod_sys)

                if cs.msgtype == 11:
                    print(cs.integ.vp == intr.vp)
                    print(cs.integ.uri == intr.uri)

    return cs


def read_rtcm(file_rtcm, intr, nep=1, logfile=None):

    fc = open(file_rtcm, 'rb')
    if not fc:
        print("RTCM message file cannot open.")

    blen = os.path.getsize(file_rtcm)
    msg = fc.read(blen)
    maxlen = len(msg)-5
    fc.close()

    return decode_rtcm(msg, intr, nep, logfile, maxlen)


if __name__ == "__main__":
    file_rtcm = '../data/sample.rtcm'
    file_log = '../data/sample.log'

    file_asc = '../data/sc134/MT05_DFi56=00.txt'

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
                 uGNSS.GLO: [1, 24],  # Table 8.5-3
                 uGNSS.GAL: [1, 36],  # Table 8.5-5
                 uGNSS.BDS: [1, 63],
                 uGNSS.QZS: [193, 202],
                 uGNSS.IRN: [1, 14]}

    mt = msg_t[0]

    seed(seed_)
    prn_t = gen_sat_list(sys_t, prn_rng_t)  # generate random sat list
    intr = gen_data(mt, sys_t, prn_t)  # generate random message data
    msg = write_rtcm(file_rtcm, msg_t, intr, nep)
    cs = read_rtcm(file_rtcm, intr, nep, logfile=file_log)

    # msg = read_asc(file_asc)
    # cs = rtcm(foutname=file_log)
    # decode_rtcm(msg)
