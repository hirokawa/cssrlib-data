"""
 static test for PPP (MADOCA PPP)
"""
from binascii import unhexlify
import numpy as np
from sys import exit as sys_exit
from cssrlib.gnss import ecef2pos, Nav, load_config
from cssrlib.gnss import time2doy, epoch2time
from cssrlib.cssrlib import sCType as sc
from cssrlib.cssr_mdc import cssr_mdc
from cssrlib.pppssr import pppos
from cssrlib.rinex import rnxdec
from cssrlib.utils import process

config = load_config('config_ppp.yml')


def decode_msg(v, tow, prn_ref, l6_ch=0, prn_ref_ext=0, l6_ch_ext=0):
    """ find valid correction message """

    msg, msg_e = None, None

    vi_ = v[v['tow'] == tow]
    vi = vi_[(vi_['type'] == l6_ch) & (vi_['prn'] == prn_ref)]
    if len(vi) > 0:
        msg = unhexlify(vi['nav'][0])

    # load regional STEC info (experimental)
    if prn_ref_ext > 0:
        vi = vi_[(vi_['type'] == l6_ch_ext) &
                 (vi_['prn'] == prn_ref_ext)]
        if len(vi) > 0:
            msg_e = unhexlify(vi['nav'][0])

    return msg, msg_e


# Select test case
#
ttl = 'test_pppmdc'
l6_mode = 0  # 0: from receiver log, 1: from archive on QZSS
dataset = 1
navfile = None
file_l6 = None
file_stec = None

prn_ref = 199  # QZSS PRN
l6_ch = 1  # 0:L6D, 1:L6E

prn_ref_ext = -1
l6_ch_ext = 0  # 0:L6D,1:L6E

# Start epoch and number of epochs
#
if dataset == 0:
    ep = [2025, 2, 15, 17, 0, 0]
    xyz_ref = [-3962108.6836, 3381309.5672, 3668678.6720]
    if l6_mode == 1:
        file_l6 = '../data/doy2025-046/2025046R.l6'
elif dataset == 1:
    ep = [2025, 8, 21, 7, 0, 0]
    xyz_ref = [-3962108.6836, 3381309.5672, 3668678.6720]  # Kamakura
elif dataset == 2:  # MADOCA-PPP with iono correction
    ep = [2025, 8, 21, 7, 0, 0]
    xyz_ref = [-3962108.6836, 3381309.5672, 3668678.6720]  # Kamakura
    file_stec = '../data/qzsl6/2025233H.201.l6'  # STEC correction
    # prn_ref_ext = 201  # QZSS PRN for iono-correction


time = epoch2time(ep)
year = ep[0]
doy = int(time2doy(time))
ses = chr(ord('a')+ep[3])

if navfile is None:
    navfile = f'../data/doy{year}-{doy:03d}/{doy:03d}{ses}_rnx.nav'
    obsfile = f'../data/doy{year}-{doy:03d}/{doy:03d}{ses}_rnx.obs'
if file_l6 is None:
    file_l6 = f'../data/doy{year}-{doy:03d}/{doy:03d}{ses}_qzsl6.txt'

nep = 3600-5
# nep = 120

if l6_mode == 1:
    fc = open(file_l6, 'rb')
    if not fc:
        print("ERROR: cannot open L6 message file {}!".format(file_l6))
        sys_exit(-1)
else:
    dtype = [('wn', 'int'), ('tow', 'int'), ('prn', 'int'),
             ('type', 'int'), ('len', 'int'), ('nav', 'S500')]
    v = np.genfromtxt(file_l6, dtype=dtype)

if file_stec is not None:
    fc = open(file_stec, 'rb')
    if not fc:
        print("ERROR: cannot open L6 message file {}!".format(file_stec))
        sys_exit(-1)

iono_opt = 2 if prn_ref_ext > 0 or (file_stec is not None) else 1

rnx = rnxdec()
nav = Nav()
proc = process(nav, rnx, config, nep=nep, xyz_ref=xyz_ref)

# Define signals to be processed

# sig_t = {'G': ['1C', '2W'], 'E': ['1C', '5Q']} # GE
sig_t = {'G': ['1C', '2W'], 'E': ['1C', '5Q'], 'J': ['1C', '5Q'],
         'R': ['1C', '2C'], 'C': ['2I', '5P']}  # GEJRC

sigs, nav.nf = proc.init_sig(sig_t)
rnx.setSignals(sigs)

# Positioning mode
# 0:static, 1:kinematic
#
nav.pmode = 0

# Decode RINEX NAV data
#
nav = rnx.decode_nav(navfile, nav)

cs = cssr_mdc()
cs.monlevel = 0

if prn_ref_ext > 0 or file_stec is not None:
    cs_ = cssr_mdc()
    cs_.monlevel = 0
else:
    cs_ = None

"""
cs = cssr_mdc('../data/madoca_cssr.log')
cs.monlevel = 2
cs_ = cssr_mdc('../data/madoca_cssr_ex.log')
cs_.monlevel = 2
"""

# Load RINEX OBS file header
#
rnx.decode_obsh(obsfile)

# Initialize position
#
ppp = pppos(nav, rnx.pos, f'{ttl}.log', iono_opt=iono_opt)
# nav.armode = 3
# nav.thresar = 2.0

proc.prepare_signal(obsfile)

# Skip epochs until start time
#
obs = rnx.decode_obs()
while time > obs.t and obs.t.time != 0:
    obs = rnx.decode_obs()

# Loop over number of epoch from file start
#
for ne in range(nep):
    week, tow = cs.set_time(obs.t)  # set time for reference

    # Set initial epoch
    #
    if ne == 0:
        proc.init_time(obs.t)

    if l6_mode == 1:  # from log file
        cs.decode_l6msg(fc.read(250), 0)
        if cs.fcnt == 5:  # end of sub-frame
            cs.week = week
            cs.decode_cssr(cs.buff, 0)
    else:  # from L6 log
        msg, msg_e = decode_msg(v, tow, prn_ref, l6_ch, prn_ref_ext, l6_ch_ext)
        if msg is not None:
            cs.decode_l6msg(msg, 0)
            if cs.fcnt == 5:  # end of sub-frame
                cs.decode_cssr(bytes(cs.buff), 0)

        # load regional STEC info (experimental)
        if prn_ref_ext > 0 and msg_e is not None:
            cs_.decode_l6msg(msg_e, 0)
            if cs_.sid == 1:  # end of sub-frame
                cs.decode_cssr(bytes(cs_.buff_p), 0)

        if file_stec is not None:  # STEC read from file
            cs_.decode_l6msg(fc.read(250), 0)
            if cs_.sid == 1:  # end of sub-frame
                cs.decode_cssr(bytes(cs_.buff_p), 0)

        cs.inet = cs.find_grid_index(ecef2pos(nav.x[:3]))

    # Call PPP module
    #
    if (cs.lc[0].cstat & 0xf) == 0xf:
        if iono_opt == 2:  # STEC is available
            mask_s = 1 << sc.STEC
            if cs.inet > 0 and \
                    cs.lc[cs.inet].cstat & mask_s == mask_s:
                ppp.process(obs, cs=cs)
        else:
            ppp.process(obs, cs=cs)

    proc.save_output(obs.t, ne, ppp)  # save output

    # Get new epoch, exit after last epoch
    #
    obs = rnx.decode_obs()
    if obs.t.time == 0:
        break

proc.close()
proc.plot(ttl, fig_type=1)
