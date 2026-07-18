"""
 static test for SBAS (L1 or DFMC)
"""
from binascii import unhexlify
import numpy as np
from sys import exit as sys_exit
from cssrlib.gnss import Nav, load_config, uIonoModel
from cssrlib.gnss import time2doy, timediff, epoch2time
from cssrlib.pntpos import stdpos
from cssrlib.sbas import sbasDec
from cssrlib.rinex import rnxdec
from cssrlib.cssr_pvs import decode_sinca_line
from cssrlib.utils import process


def decode_msg(v, tow, prn_ref, sbas_type=0):
    """ find valid correction message """

    if len(prn_ref) == 1:
        vi = v[(v['tow'] == tow) & (v['prn'] == prn_ref)]
    else:
        vi = v[(v['tow'] == tow) & (v['prn'] >= prn_ref[0]) &
               (v['prn'] <= prn_ref[1])]
    if sbas_type == 0:  # L1
        vi = vi[vi['type'] <= 28]
    else:  # DFMC L5
        vi = vi[(vi['type'] == 31) | (vi['type'] == 32) |
                ((vi['type'] >= 34) & (vi['type'] <= 37))]
    if len(vi) > 0:
        msg = {}
        for vi_ in vi:
            msg[vi_['prn']] = unhexlify(vi_['nav'])
    else:
        msg = None

    return msg


# Select test case
#
nep = 3600
ttl = 'test_sbas'  # title for case
dataset = 1
navfile = None

config = load_config('config.yml')

# Start epoch and number of epochs
#
if dataset == 1:  # MSAS, L1 SBAS
    ep = [2025, 8, 21, 7, 0, 0]
    xyz_ref = [-3962108.6836, 3381309.5672, 3668678.6720]
    prn_ref = [137]  # satellite PRN for SBAS correction
    sbas_type = 0  # L1: 0, L5: 1
    nf = 1

elif dataset == 2:  # QZSS, L5 DFMC
    ep = [2025, 8, 21, 7, 0, 0]
    xyz_ref = [-3962108.6836, 3381309.5672, 3668678.6720]
    # prn_ref = [193, 202]  # satellite PRN for SBAS correction
    prn_ref = [199]
    sbas_type = 1  # L1: 0, L5: 1
    nf = 2

elif dataset == 3:  # SouthPAN L5 DFMC (DAS)
    ep = [2025, 4, 20, 5, 0, 0]
    # navfile = '../data/doy2025-110/ALIC00AUS0110f.nav'
    navfile = '../data/doy2025-110/BRD400DLR_S_20251100000_01D_MN.rnx'
    obsfile = '../data/doy2025-110/ALIC00AUS0110f.obs'
    file_sbas = '../data/doy2025-110/DAS2025110f.txt'
    xyz_ref = [-4052052.9320,  4212835.9496, -2545104.3074]
    sbas_type = 1
    nf = 2

elif dataset == 4:  # SouthPAN L5 DFMC (SIS)
    ep = [2025, 8, 21, 7, 0, 0]
    xyz_ref = [-3962108.6836, 3381309.5672, 3668678.6720]  # Kamakura
    prn_ref = [122]
    sbas_type = 1  # L1: 0, L5: 1
    nf = 2

time = epoch2time(ep)
year = ep[0]
doy = int(time2doy(time))
ses = chr(ord('a')+ep[3])

if navfile is None:
    navfile = f'../data/doy{year}-{doy:03d}/{doy:03d}{ses}_rnx.nav'
    obsfile = f'../data/doy{year}-{doy:03d}/{doy:03d}{ses}_rnx.obs'
    file_sbas = f'../data/doy{year}-{doy:03d}/{doy:03d}{ses}_sbas.txt'


if sbas_type == 0:  # L1 SBAS
    sig_t = {'G': ['1C']}

else:  # DFMC
    sig_t = {'G': ['1C', '5Q'], 'E': ['1C', '5Q']}

rnx = rnxdec()
nav = Nav()
proc = process(nav, rnx, config, nep=nep, xyz_ref=xyz_ref)

# Define signals to be processed
#
sigs, nav.nf = proc.init_sig(sig_t)
rnx.setSignals(sigs)

# Decode RINEX NAV data
#
nav = rnx.decode_nav(navfile, nav)

cs = sbasDec('test_sbas_cs.log', nf=nf)
cs.monlevel = config['cs']['monlevel']

rnx.decode_obsh(obsfile)  # Load RINEX OBS file header

# Initialize position
#
std = stdpos(nav, rnx.pos, f'{ttl}.log')
std.monlevel = 1
std.ionoModel = uIonoModel.SBAS

proc.prepare_signal(obsfile)

nav.rmode = 2 if nav.nf == 2 else 0  # L1/L5 iono-free combination

# Skip epochs until start time
#
obs = rnx.decode_obs()
while time > obs.t and obs.t.time != 0:
    obs = rnx.decode_obs()

if 'sbas' in file_sbas:  # SIS
    dtype = [('wn', 'int'), ('tow', 'float'), ('prn', 'int'),
             ('type', 'int'), ('marker', 'S2'), ('nav', 'S124')]
    v = np.genfromtxt(file_sbas, dtype=dtype)
elif 'DAS' in file_sbas:  # DAS
    fc = open(file_sbas, 'rt')
else:
    print("ERROR: unknown file format for correction data")
    sys_exit(1)

# Loop over number of epoch from file start
#
for ne in range(nep):
    _, tow = cs.set_time(obs.t)  # set time for reference

    # Set initial epoch
    #
    if ne == 0:
        proc.init_time(obs.t)

    if 'sbas' in file_sbas:  # SIS
        msgs = decode_msg(v, tow, prn_ref, sbas_type)
        if msgs is not None:
            for prn, msg in msgs.items():
                cs.decode_cssr(msg, 0, src=sbas_type, prn=prn)

    else:  # DAS
        for line in fc:
            tc, buff = decode_sinca_line(line)
            cs.decode_cssr(buff, 0, src=sbas_type)
            if timediff(obs.t, tc) >= 0.0:
                break

    # Call standard positioning module with SBAS corrections
    #
    if (cs.lc[0].cstat & 0x6) == 0x6:  # wait for orbit/clock correction
        std.process(obs, cs=cs)

    proc.save_output(obs.t, ne)  # save output

    # Get new epoch, exit after last epoch
    #
    obs = rnx.decode_obs()
    if obs.t.time == 0:
        break

proc.close()
proc.plot(ttl, fig_type=1, ylim=2, ylim_v=4)
