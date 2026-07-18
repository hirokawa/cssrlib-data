"""
 static test for PPP (IGS)
"""
from copy import deepcopy
from cssrlib.gnss import ecef2pos, Nav, load_config
from cssrlib.gnss import time2doy, epoch2time
from cssrlib.peph import peph, biasdec
from cssrlib.pppssr import pppos
from cssrlib.rinex import rnxdec
from cssrlib.utils import process

config = load_config('config_ppp.yml')

# Start epoch and number of epochs
#
dataset = 4
ttl = 'test_pppigs'

if dataset == 0:  # SETP078M.21O
    ep = [2021, 3, 19, 12, 0, 0]
    xyz_ref = [-3962108.6617, 3381309.5232, 3668678.6410]
elif dataset == 1:  # SETP1890.23O
    ep = [2023, 7, 8, 4, 0, 0]
    xyz_ref = [-3962108.7063, 3381309.5703, 3668678.6690]
elif dataset == 2:  # SETP223Z.23O
    ep = [2023, 8, 11, 21, 0, 0]
    xyz_ref = [-3962108.7063, 3381309.5703, 3668678.6690]
elif dataset == 3:  # 046r_rnx.obs
    ep = [2025, 2, 15, 17, 0, 0]
    xyz_ref = [-3962108.6819, 3381309.5707, 3668678.6750]
elif dataset == 4:  # 046r_rnx.obs
    ep = [2025, 8, 21, 7, 0, 0]
    xyz_ref = [-3962108.6836, 3381309.5672, 3668678.6720]  # Kamakura
else:
    print("ERROR: no RINEX data set selected!")
    exit(1)

time = epoch2time(ep)
year = ep[0]
doy = int(time2doy(time))

nep = 900*4

pos_ref = ecef2pos(xyz_ref)

bdir = '../data/doy{:04d}-{:03d}/'.format(year, doy)

if dataset == 0:
    let = chr(ord('A')+ep[3])
    navfile = bdir+'SEPT{:03d}{}.{:02d}P'.format(doy, let, year % 2000)
    obsfile = bdir+'SEPT{:03d}{}.{:02d}O'.format(doy, let, year % 2000)
if dataset == 1:
    let = '0'
    navfile = bdir+'SEPT{:03d}{}.{:02d}P'.format(doy, let, year % 2000)
    obsfile = bdir+'SEPT{:03d}{}.{:02d}O'.format(doy, let, year % 2000)
elif dataset == 2:
    let = 'Z'
    navfile = '../data/brdc/' + \
        'BRD400DLR_S_{:04d}{:03d}0000_01D_MN.rnx'.format(year, doy)
    obsfile = bdir+'SEPT{:03d}{}.{:02d}O'.format(doy, let, year % 2000)
else:
    let = chr(ord('a')+ep[3])
    navfile = bdir+'{:03d}{}_rnx.nav'.format(doy, let)
    obsfile = bdir+'{:03d}{}_rnx.obs'.format(doy, let)

ac = 'COD0OPSFIN'

orbfile = '../data/igs/{}_{:4d}{:03d}0000_01D_05M_ORB.SP3'\
    .format(ac, year, doy)

clkfile = '../data/igs/{}_{:4d}{:03d}0000_01D_30S_CLK.CLK'\
    .format(ac, year, doy)

bsxfile = '../data/igs/{}_{:4d}{:03d}0000_01D_01D_OSB.BIA'\
    .format(ac, year, doy)

# Define signals to be processed
#
sig_t = {'G': ['1C', '2W'], 'E': ['1C', '5Q'], 'J': ['1C', '5Q']}  # GEJ

rnx = rnxdec()
nav = Nav()
orb = peph()

proc = process(nav, rnx, config, nep=nep, xyz_ref=xyz_ref)

sigs, nav.nf = proc.init_sig(sig_t)
rnx.setSignals(sigs)

# Set site code for GLONASS code biases
#
site = "SUWN"

# Positioning mode
# 0:static, 1:kinematic
#
nav.pmode = 0

# Decode RINEX NAV data
#
nav = rnx.decode_nav(navfile, nav)

# Load precise orbits and clock offsets
#
nav = orb.parse_sp3(orbfile, nav)
nav = rnx.decode_clk(clkfile, nav)

# Load code and phase biases from Bias-SINEX
#
bsx = biasdec()
bsx.parse(bsxfile, site)

# Load ANTEX data for satellites and stations
#
atxfile = '../data/antex/'
if time > epoch2time([2022, 11, 27, 0, 0, 0]):
    atxfile += 'I20.ATX' if 'COD0MGXFIN' in ac else 'igs20_2353.atx'
elif time > epoch2time([2021, 5, 2, 0, 0, 0]):
    atxfile += 'M20.ATX' if 'COD0MGXFIN' in ac else 'igs14.atx'
else:
    atxfile += 'M14.ATX' if 'COD0MGXFIN' in ac else 'igs14.atx'

config['atxfile'] = atxfile


# Load RINEX OBS file header
#
rnx.decode_obsh(obsfile)

# Initialize position
#
ppp = pppos(nav, rnx.pos, f'{ttl}.log')
nav.ephopt = 4  # IGS
nav.armode = 3  # 1: continuous, 3: fix-and-hold
nav.parmode = 1  # 1: normal, 2: partial ambiguity resolution
nav.thresar = 2.0

proc.prepare_signal(obsfile)

# Skip epochs until start time
#
obs = rnx.decode_obs()
while time > obs.t and obs.t.time != 0:
    obs = rnx.decode_obs()

# Loop over number of epoch from file start
#
for ne in range(nep):

    # Set initial epoch
    #
    if ne == 0:
        nav.t = deepcopy(obs.t)
        t0 = deepcopy(obs.t)

    # Call PPP module with IGS products
    #
    ppp.process(obs, orb=orb, bsx=bsx)

    proc.save_output(obs.t, ne, ppp)  # save output

    # Get new epoch, exit after last epoch
    #
    obs = rnx.decode_obs()
    if obs.t.time == 0:
        break

proc.close()
proc.plot(ttl, fig_type=1)
