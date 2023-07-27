# -*- coding: utf-8 -*-
"""
example of Galileo HAS correction data decoder

@author: ruihi

[1] Galileo High Accuracy Service Signal-in-Space
  Interface Control Document (HAS SIS ICD), Issue 1.0, May 2022

"""

import numpy as np
from cssrlib.cssr_has import cssr_has


def parse_has_data_sample(file, ex=1):
    """ load HAS pages from file attached with Galileo HAS ICD """
    rec = []
    valid = False
    nd = 256
    nc = 53
    has_pages = np.zeros((nd-1, nc), dtype=int)

    with open(file, "r") as fh:
        for line in fh:
            if "* HAS MESSAGE DECODING EXAMPLE {:d} *".format(ex) in line:
                valid = True
            if "* END HAS MESSAGE DECODING EXAMPLE {:d} *".format(ex) in line:
                valid = False

            if valid:
                if '// PID' in line:
                    pid = int(fh.readline())
                    rec += [pid-1]
                elif '// HAS encoded page' in line:
                    line = fh.readline()
                    has_pages[pid-1, :] = \
                        np.genfromtxt(line.split()[1:-1], dtype="u1")
    return rec, has_pages


# available from Galileo HAS ICD
file = 'Galileo-HAS-SIS-ICD-1.0_Annex_D_HAS_Message_Decoding_Example.txt'
file_gm = "Galileo-HAS-SIS-ICD_1.0_Annex_B_Reed_Solomon_Generator_Matrix.txt"

gMat = np.genfromtxt(file_gm, dtype="u1", delimiter=",")

dec = cssr_has()
dec.mon_level = 2

# example 1 MS=15 mask/orbit/cbias/pbias
rec, has_pages = parse_has_data_sample(file, 1)
HASmsg = dec.decode_has_page(rec, has_pages, gMat, 15)
dec.decode_cssr(HASmsg)

# example 2 MS=2 clock
rec, has_pages = parse_has_data_sample(file, 2)
HASmsg = dec.decode_has_page(rec, has_pages, gMat, 2)
dec.decode_cssr(HASmsg)
