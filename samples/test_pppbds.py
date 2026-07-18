"""
 static test for PPP (BeiDou PPP)
"""
from binascii import unhexlify
import numpy as np

from cssrlib.gnss import Nav, load_config, time2doy, epoch2time
from cssrlib.cssr_bds import cssr_bds
from cssrlib.pppssr import pppos
from cssrlib.rinex import rnxdec
from cssrlib.utils import process


def decode_msg(v, tow, prn_ref):
    """ find valid correction message """
    msg = None
    vi = v[(v['tow'] == tow) & (v['prn'] == prn_ref)]
    if len(vi) > 0:
        msg = unhexlify(vi['nav'][0])
    return msg


config = load_config('config_ppp.yml')

# Select test case
#
nep = 900*4
prn_ref = 59  # satellite PRN to receive BDS PPP collection
dataset = 1
navfile = None
file_bds = None
ttl = 'test_pppbds'

# Start epoch and number of epochs
#
if dataset == 0:
    ep = [2025, 2, 15, 17, 0, 0]
    xyz_ref = [-3962108.6836, 3381309.5672, 3668678.6720]

elif dataset == 1:
    ep = [2025, 8, 21, 7, 0, 0]
    xyz_ref = [-3962108.6836, 3381309.5672, 3668678.6720]  # Kamakura

time = epoch2time(ep)
year = ep[0]
doy = int(time2doy(time))
ses = chr(ord('a')+ep[3])

if navfile is None:
    navfile = f'../data/doy{year}-{doy:03d}/{doy:03d}{ses}_rnx.nav'
    obsfile = f'../data/doy{year}-{doy:03d}/{doy:03d}{ses}_rnx.obs'
if file_bds is None:
    file_bds = f'../data/doy{year}-{doy:03d}/{doy:03d}{ses}_bdsb2b.txt'


dtype = [('wn', 'int'), ('tow', 'int'), ('prn', 'int'),
         ('type', 'int'), ('len', 'int'), ('nav', 'S124')]
v = np.genfromtxt(file_bds, dtype=dtype)

rnx = rnxdec()
nav = Nav()
proc = process(nav, rnx, config, nep=nep, xyz_ref=xyz_ref)

sig_t = {'G': ['1C', '2W'], 'C': ['1P', '5P']}  # GC

sigs, nav.nf = proc.init_sig(sig_t)
rnx.setSignals(sigs)

# Positioning mode
# 0:static, 1:kinematic
#
nav.pmode = 0

# Decode RINEX NAV data
#
nav = rnx.decode_nav(navfile, nav)

cs = cssr_bds()
cs.monlevel = 0
"""
cs = cssr_bds(f'{ttl}_ssr.log')
cs.monlevel = 2
"""

# Load RINEX OBS file header
#
rnx.decode_obsh(obsfile)

# Initialize position
#
ppp = pppos(nav, rnx.pos, f'{ttl}.log')

if time < epoch2time([2022, 11, 27, 0, 0, 0]):
    config['atxfile'] = '../data/antex/igs14.atx'

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

    msg = decode_msg(v, tow, prn_ref)
    if msg is not None:
        cs.decode_cssr(msg, 0)

    # Call PPP module with BDS-PPP corrections
    #
    if (cs.lc[0].cstat & 0xf) == 0xf:  # wait for mask/orb/clk/cbias
        ppp.process(obs, cs=cs)

    proc.save_output(obs.t, ne, ppp)  # save output

    # Get new epoch, exit after last epoch
    #
    obs = rnx.decode_obs()
    if obs.t.time == 0:
        break

proc.close()
proc.plot(ttl, fig_type=1)
