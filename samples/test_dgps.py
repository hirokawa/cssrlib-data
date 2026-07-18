"""
 static test for DGPS (QZSS SLAS)
"""
from binascii import unhexlify
import numpy as np
from cssrlib.gnss import Nav, time2gpst, time2doy, epoch2time
from cssrlib.gnss import rSigRnx, uIonoModel, load_config
from cssrlib.pntpos import stdpos
from cssrlib.dgps import dgpsDec
from cssrlib.rinex import rnxdec
from cssrlib.utils import process


def decode_msg(v, tow, prn_ref):
    """ find valid correction message """

    vi = v[(v['tow'] == tow) & (v['prn'] == prn_ref)
           & ((v['type'] >= 47) & (v['type'] <= 51))]
    if len(vi) > 0:
        msg = unhexlify(vi['nav'][0])
    else:
        msg = None

    return msg


config = load_config('config.yml')

# Select test case
#
dataset = 0
nep = 3600
ttl = 'test_dgps'  # title for case

navfile = None

# Start epoch and number of epochs
#
if dataset == 0:  # QZSS SLAS
    ep = [2025, 8, 21, 7, 0, 0]
    xyz_ref = [-3962108.6836, 3381309.5672, 3668678.6720]  # Kamakura

    prn_ref = 199  # satellite PRN for SBAS correction
    sbas_type = 0  # L1: 0, L5: 1
    nf = 1

time = epoch2time(ep)
year = ep[0]
doy = int(time2doy(time))
ses = chr(ord('a')+ep[3])

if navfile is None:
    navfile = f'../data/doy{year}-{doy:03d}/{doy:03d}{ses}_rnx.nav'
    obsfile = f'../data/doy{year}-{doy:03d}/{doy:03d}{ses}_rnx.obs'
    file_sbas = f'../data/doy{year}-{doy:03d}/{doy:03d}{ses}_sbas.txt'

if sbas_type == 0:  # DGNSS
    sig_t = {'G': ['1C'], 'J': ['1C']}

else:  # dual frequency
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

cs = dgpsDec(f'{ttl}_cs.log')
cs.monlevel = config['cs']['monlevel']

# Load RINEX OBS file header
#
rnx.decode_obsh(obsfile)

std = stdpos(nav, rnx.pos, f'{ttl}.log', trop_opt=2, iono_opt=2)
std.ionoModel = uIonoModel.SBAS

proc.prepare_signal(obsfile)

# Skip epochs until start time
#
obs = rnx.decode_obs()
while time > obs.t and obs.t.time != 0:
    obs = rnx.decode_obs()

dtype = [('wn', 'int'), ('tow', 'float'), ('prn', 'int'),
         ('type', 'int'), ('marker', 'S2'), ('nav', 'S124')]
v = np.genfromtxt(file_sbas, dtype=dtype)

# Loop over number of epoch from file start
#
for ne in range(nep):
    _, tow = cs.set_time(obs.t)  # set time for reference

    # Set initial epoch
    #
    if ne == 0:
        proc.init_time(obs.t)

    msg = decode_msg(v, tow, prn_ref)
    if msg is not None:
        cs.decode_cssr(msg, 0)  # decode DGPS correction

    if (cs.lc[0].cstat & 0x7) == 0x7:  # wait for mask/clock/orbit
        std.process(obs, cs=cs)  # standalone positioning

    proc.save_output(obs.t, ne)  # save output

    obs = rnx.decode_obs()  # get new epoch
    if obs.t.time == 0:  # exit after last epoch
        break

proc.close()
proc.plot(ttl, fig_type=1, ylim=2, ylim_v=4)
