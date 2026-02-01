"""
 static test for PPP (PVS PPP)
"""
from binascii import unhexlify
import numpy as np
from sys import exit as sys_exit

from cssrlib.gnss import Nav, time2doy, timediff, epoch2time, load_config
from cssrlib.cssr_pvs import cssr_pvs
from cssrlib.pppssr import pppos
from cssrlib.rinex import rnxdec
from cssrlib.cssr_pvs import decode_sinca_line
from cssrlib.utils import process

config = load_config('config_ppp.yml')

# Select test case
#
dataset = 2
navfile = None
file_pvs = None
ttl = 'test_ppppvs'


def decode_msg(v, tow, prn_ref, sbas_type=1):
    """ find valid correction message """
    msg = None
    vi = v[(v['tow'] == tow) & (v['prn'] == prn_ref)]
    if sbas_type == 0:  # L1
        vi = vi[vi['type'] <= 30]
    else:  # DFMC L5
        vi = vi[vi['type'] > 30]

    if len(vi) > 0:
        msg = unhexlify(vi['nav'][0])

    return msg


# Start epoch and input files
#
if dataset == 0:  # DAS

    ep = [2025, 4, 20, 5, 0, 0]
    navfile = '../data/doy2025-110/BRD400DLR_S_20251100000_01D_MN.rnx'
    obsfile = '../data/doy2025-110/ALIC00AUS0110f.obs'
    file_pvs = '../data/doy2025-110/DAS2025110f.txt'
    xyz_ref = [-4052052.9320,  4212835.9496, -2545104.3074]

elif dataset == 1:  # SIS

    ep = [2025, 2, 15, 17, 0, 0]
    xyz_ref = [-3962108.6836, 3381309.5672, 3668678.6720]

elif dataset == 2:  # SIS

    ep = [2025, 8, 21, 7, 0, 0]
    xyz_ref = [-3962108.6836, 3381309.5672, 3668678.6720]  # Kamakura

elif dataset == 3:  # SIS (Australia)

    ep = [2025, 8, 21, 7, 0, 0]
    navfile = '../data/doy2025-233/alby233h_rnx.nav'
    obsfile = '../data/doy2025-233/alby233h_rnx.obs'  # SEPT POLARX5
    xyz_ref = [-2441715.2741, 4629128.6896, -3633362.5218]  # Albany, AUSTRALIA

time = epoch2time(ep)
year = ep[0]
doy = int(time2doy(time))
ses = chr(ord('a')+ep[3])

if navfile is None:
    navfile = f'../data/doy{year}-{doy:03d}/{doy:03d}{ses}_rnx.nav'
    obsfile = f'../data/doy{year}-{doy:03d}/{doy:03d}{ses}_rnx.obs'
if file_pvs is None:
    file_pvs = f'../data/doy{year}-{doy:03d}/{doy:03d}{ses}_sbas.txt'

nep = 900*4
prn_ref = 122  # satellite PRN for PRN122
sbas_type = 1  # L1: 0, L5: 1

rnx = rnxdec()
nav = Nav()
proc = process(nav, rnx, config, nep=nep, xyz_ref=xyz_ref)

sig_t = {'G': ['1C', '5Q'], 'E': ['1C', '5Q']}  # GE

sigs, nav.nf = proc.init_sig(sig_t)
rnx.setSignals(sigs)


# Positioning mode
# 0:static, 1:kinematic
#
nav.pmode = 0

# Decode RINEX NAV data
#
nav = rnx.decode_nav(navfile, nav)

# cs = cssr_pvs()
# cs.monlevel = 0
cs = cssr_pvs(f'{ttl}_ssr.log')
cs.monlevel = 2

# Load RINEX OBS file header
#
rnx.decode_obsh(obsfile)

# Initialize position
#
ppp = pppos(nav, rnx.pos, f'{ttl}.log')

proc.prepare_signal(obsfile)

# Skip epochs until start time
#
obs = rnx.decode_obs()
while time > obs.t and obs.t.time != 0:
    obs = rnx.decode_obs()

if 'sbas' in file_pvs:  # SIS
    dtype = [('wn', 'int'), ('tow', 'float'), ('prn', 'int'),
             ('type', 'int'), ('marker', 'S2'), ('nav', 'S124')]
    v = np.genfromtxt(file_pvs, dtype=dtype)
elif 'DAS' in file_pvs:  # DAS
    fc = open(file_pvs, 'rt')
else:
    print("ERROR: unknown file format for correction data")
    sys_exit(1)

# Loop over number of epoch from file start
#
for ne in range(nep):

    week, tow = cs.set_time(obs.t)  # set time for reference

    # Set initial epoch
    #
    if ne == 0:
        proc.init_time(obs.t)

    if 'sbas' in file_pvs:  # SIS
        msg = decode_msg(v, tow, prn_ref, sbas_type)
        if msg is not None:
            cs.decode_cssr(msg, 0)

    else:  # DAS
        for line in fc:
            tc, buff = decode_sinca_line(line)
            cs.decode_cssr(buff, 0)
            if timediff(obs.t, tc) >= 0.0:
                break

    cs.check_validity(obs.t)

    # Call PPP module with PVS corrections
    #
    if (cs.lc[0].cstat & 0x6) == 0x6:
        ppp.process(obs, cs=cs)

    proc.save_output(obs.t, ne, ppp)  # save output

    # Get new epoch, exit after last epoch
    #
    obs = rnx.decode_obs()
    if obs.t.time == 0:
        break

proc.close()
proc.plot(ttl, fig_type=1)
