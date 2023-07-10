# -*- coding: utf-8 -*-
"""
Created on Fri May  7 20:14:35 2021

@author: ruihi
"""

import numpy as np
import matplotlib.pyplot as plt

CA_LEN = 1023


def rotate(r):
    sz = len(r)
    idx_r = np.arange(sz-1, 0, -1)
    c = r[sz-1]
    r[idx_r] = r[idx_r-1]
    r[0] = c
    return r


def codegen(prn=1):
    DW_SF = 10
    g2_delay_gps = [
        5,   6,   7,   8,  17,  18, 139, 140, 141, 251,
        252, 254, 255, 256, 257, 258, 469, 470, 471, 472,
        473, 474, 509, 512, 513, 514, 515, 516, 859, 860,
        861, 862]
    g2_delay_qzs = [
        339, 208, 711, 189, 263, 537, 663, 942, 173, 900, 30, 500, 935, 556]
    ca = np.zeros(CA_LEN, dtype=int)

    if prn >= 1 and prn <= 32:
        g2_0 = g2_delay_gps[prn-1]
    elif prn >= 193 and prn <= 206:
        g2_0 = g2_delay_qzs[prn-193]
    else:
        return ca

    r1 = -1*np.ones(DW_SF, dtype=int)
    r2 = -1*np.ones(DW_SF, dtype=int)

    g1 = np.zeros(CA_LEN, dtype=int)
    g2 = np.zeros(CA_LEN, dtype=int)

    for i in range(CA_LEN):
        # g1/g2 tap
        g1[i] = r1[9]
        g2[i] = r2[9]
        c1 = r1[2]*r1[9]
        c2 = r2[1]*r2[2]*r2[5]*r2[7]*r2[8]*r2[9]
        # shift register
        r1 = rotate(r1)
        r1[0] = c1
        r2 = rotate(r2)
        r2[0] = c2

    idx_ca = np.arange(CA_LEN)
    ca = (1-g1*g2[(CA_LEN-g2_0+idx_ca) % CA_LEN])//2
    ca[ca == 0] = -1
    return ca


if __name__ == '__main__':

    CLIGHT = 299792458.0
    crate = 1.023e6    # chip rate [chip/sec]

    prn = 193
    ca = codegen(prn)  # C/A code generation

    P = 36952979.472  # Suspected distance [m]
    # P=36953010.338    # Pseudorange (rounded) [m]
    dt = P/CLIGHT  # propagation time [s]
    chips = crate*dt
    delay = round(chips) % CA_LEN  # 268
    n = chips//CA_LEN    # 123
    for i in range(delay):  # code delay shift
        ca = rotate(ca)

    # Processing on the receiver side
    prn = 193
    ca_ = codegen(prn)  # C/A code replication of PRN193
    c = []
    for i in range(CA_LEN):  # Correlate with Delayed Replicas
        c.append(ca@ca_)
        ca_ = rotate(ca_)

    delay_max = np.argmax(c)  # Find correlation peaks: 268
    P_ = (n*CA_LEN+delay_max)*CLIGHT/crate

    plt.figure()
    plt.plot(c)
    plt.grid()
    plt.xlabel('code delay [chips]')
    plt.show()

    plt.figure()
    plt.step(np.arange(CA_LEN), ca)
    plt.axis([0, 50, -1.2, 1.2])
    plt.xlabel('time [chips]')
    plt.grid()
    plt.show
