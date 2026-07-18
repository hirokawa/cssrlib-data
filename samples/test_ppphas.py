"""
 static test for PPP (Galileo HAS SIS)
"""

import numpy as np
from cssrlib.gnss import Nav, time2doy, epoch2time
from cssrlib.gnss import prn2sat, uGNSS, load_config
from cssrlib.cssr_has import cssr_has, cnav_msg
from cssrlib.pppssr import pppos
from cssrlib.rinex import rnxdec
from cssrlib.utils import process

config = load_config('config_ppp.yml')

# Select test case
#
ttl = 'test_ppphas'
dataset = 3
nep = 900*4

excl_sat = []
navfile = None
file_has = None
# Start epoch and number of epochs
#
fromSbfConvert = False
if dataset == 0:
    ep = [2023, 7, 8, 4, 0, 0]
    xyz_ref = [-3962108.7007, 3381309.5532, 3668678.6648]
    navfile = '../data/doy2023-189/SEPT1890.23P'
    obsfile = '../data/doy2023-189/SEPT1890.23O'
    file_has = '../data/doy2023-189/gale6_189e.txt'
elif dataset == 1:
    ep = [2023, 8, 11, 21, 0, 0]
    xyz_ref = [-3962108.7007, 3381309.5532, 3668678.6648]
    navfile = '../data/doy2023-223/NAV223.23p'
    # obsfile = '../data/doy2023-223/SEPT223Z.23O'  # MOSAIC-CLAS
    obsfile = '../data/doy2023-223/SEPT223Y.23O'  # PolaRX5
    file_has = '../data/doy2023-223/223v_gale6.txt'
elif dataset == 2:
    ep = [2025, 2, 15, 17, 0, 0]
    xyz_ref = [-3962108.6836, 3381309.5672, 3668678.6720]
    excl_sat = [prn2sat(uGNSS.GAL, 29)]  # E29
elif dataset == 3:
    ep = [2025, 8, 21, 7, 0, 0]
    # obsfile = '../data/doy2025-233/sept233h_rnx.obs'  # SEPT POLARX5
    # obsfile = '../data/doy2025-233/ux2233h_rnx.obs'  # u-blox X20P
    # obsfile = '../data/doy2025-233/jav3233h_rnx.obs'  # javad DELTA-3S
    xyz_ref = [-3962108.6836, 3381309.5672, 3668678.6720]  # Kamakura


# Convert epoch and user reference position
#
time = epoch2time(ep)
year = ep[0]
doy = int(time2doy(time))
ses = chr(ord('a')+ep[3])

if navfile is None:
    navfile = f'../data/doy{year}-{doy:03d}/{doy:03d}{ses}_rnx.nav'
    obsfile = f'../data/doy{year}-{doy:03d}/{doy:03d}{ses}_rnx.obs'
if file_has is None:
    file_has = f'../data/doy{year}-{doy:03d}/{doy:03d}{ses}_gale6.txt'

# Load SSR correction file
#
if fromSbfConvert:

    dtype = [('tow', 'float64'), ('wn', 'int'),  ('prn', 'S3'),
             ('validity', 'object'), ('num1', 'int'), ('signal', 'object'),
             ('num2', 'int'), ('num3', 'int'), ('nav', 'S144')]
    v = np.genfromtxt(file_has, dtype=dtype, delimiter=',')
    v = v[v['validity'] == b'Passed']
    # Eliminate whitespace
    for i, nav in enumerate(v['nav']):
        v[i]['nav'] = (b''.join(nav.split()))

else:

    dtype = [('wn', 'int'), ('tow', 'int'), ('prn', 'int'),
             ('type', 'int'), ('len', 'int'), ('nav', 'S124')]
    v = np.genfromtxt(file_has, dtype=dtype)

# Define signals to be processed
#
if dataset == 3:
    sig_t = {'G': ['1C', '2L'], 'E': ['1C', '5Q']}  # GE
else:
    sig_t = {'G': ['1C', '2W'], 'E': ['1C', '5Q']}  # GE

rnx = rnxdec()
nav = Nav()
proc = process(nav, rnx, config, nep=nep, xyz_ref=xyz_ref)

sigs, nav.nf = proc.init_sig(sig_t)
rnx.setSignals(sigs)

# Positioning mode
# 0:static, 1:kinematic
#
nav.pmode = 0

if len(excl_sat) > 0:
    nav.excl_sat = excl_sat

# Decode RINEX NAV data
#
nav = rnx.decode_nav(navfile, nav)

cs = cssr_has()
cs.monlevel = 0

"""
cs = cssr_has('has.log')
cs.monlevel = 2
"""

file_gm = "Galileo-HAS-SIS-ICD_1.0_Annex_B_Reed_Solomon_Generator_Matrix.txt"

# Galileo CNAV message parser
cnav = cnav_msg()
cnav.load_gmat(file_gm)

# Load RINEX OBS file header
#
rnx.decode_obsh(obsfile)

# Initialize position
#
ppp = pppos(nav, rnx.pos, f'{ttl}.log')

if time < epoch2time([2025, 5, 15, 17, 18, 0]):
    config['atxfile'] = '../data/antex/has14_2345.atx'

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

    vi = v[v['tow'] == tow]

    HASmsg = cnav.decode_cnav(tow, vi)  # decode CNAV pages
    if HASmsg is not None:
        cs.msgtype = cnav.msgtype
        cs.decode_cssr(HASmsg)  # decode HAS messages

    # Call PPP module with HAS corrections
    #
    if (cs.lc[0].cstat & 0xf) == 0xf:
        ppp.process(obs, cs=cs)

    proc.save_output(obs.t, ne, ppp)  # save output

    # Get new epoch, exit after last epoch
    #
    obs = rnx.decode_obs()
    if obs.t.time == 0:
        break

proc.close()
proc.plot(ttl, fig_type=1)
