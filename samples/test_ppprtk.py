"""
 static test for PPP-RTK (QZSS CLAS)
"""
import numpy as np
from sys import exit as sys_exit

from cssrlib.cssrlib import cssr
from cssrlib.gnss import ecef2pos, Nav, time2gpst, time2doy
from cssrlib.gnss import epoch2time, load_config
from cssrlib.ppprtk import ppprtkpos
from cssrlib.rinex import rnxdec
from binascii import unhexlify
from cssrlib.utils import process


def decode_msg(v, tow, l6_ch, prn_ref):
    """ find valid correction message """

    msg = None
    vi = v[(v['tow'] == tow) & (v['type'] == l6_ch) & (v['prn'] == prn_ref)]
    if len(vi) > 0:
        msg = unhexlify(vi['nav'][0])

    return msg


config = load_config('config_ppprtk.yml')

l6_mode = 0  # 0: from receiver log, 1: from archive on QZSS
dataset = 3
nep = 900*4
# nep = 60

navfile = None
file_l6 = None
ttl = 'test_ppprtk'
l6_ch = 0  # 0:L6D, 1:L6E
prn_p1 = 199
prn_p2 = -1

sig_t = {'G': ['1C', '2W'], 'E': ['1C', '5Q'], 'J': ['1C', '2L']}  # GEJ

if l6_mode == 1:  # from archive

    if dataset == 1:

        ep = [2021, 3, 19, 12, 0, 0]
        xyz_ref = [-3962108.673, 3381309.574, 3668678.638]
        navfile = '../data/doy2021-078/SEPT078M.21P'
        obsfile = '../data/doy2021-078/SEPT078M.21O'
        l6file = '../data/doy2021-078/2021078M.l6'

    elif dataset == 2:

        ep = [2025, 2, 15, 17, 0, 0]
        xyz_ref = [-3962108.6836, 3381309.5672, 3668678.6720]

else:  # from receiver log

    if dataset == 0:

        ep = [2023, 8, 11, 21, 0, 0]
        xyz_ref = [-3962108.7007, 3381309.5532, 3668678.6648]
        navfile = '../data/doy2023-223/NAV223.23p'
        obsfile = '../data/doy2023-223/SEPT223Y.23O'  # PolaRX5

    elif dataset == 1:  # single channel

        ep = [2025, 8, 21, 7, 0, 0]
        xyz_ref = [-3962108.6836, 3381309.5672, 3668678.6720]

    elif dataset == 2:  # two channel

        ep = [2025, 8, 21, 7, 0, 0]
        xyz_ref = [-3962108.6836, 3381309.5672, 3668678.6720]

        # from Tab 4.1.1-1 of IS-QZSS-L6
        prn_p1 = 199  # QZSS PRN pattern 1 (195, 197, 199)
        prn_p2 = 194  # QZSS PRN pattern 2 (194, 196)

    elif dataset == 3:  # single channel

        ep = [2025, 8, 21, 7, 0, 0]
        xyz_ref = [-3962108.6836, 3381309.5672, 3668678.6720]

time = epoch2time(ep)
year = ep[0]
doy = int(time2doy(time))
ses = chr(ord('a')+ep[3])

if navfile is None:
    navfile = f'../data/doy{year}-{doy:03d}/{doy:03d}{ses}_rnx.nav'
    obsfile = f'../data/doy{year}-{doy:03d}/{doy:03d}{ses}_rnx.obs'
if file_l6 is None:
    file_l6 = f'../data/doy{year}-{doy:03d}/{doy:03d}{ses}_qzsl6.txt'
    if l6_mode == 1:
        file_l6 = f'../data/doy{year}-{doy:03d}/{year}{doy:03d}{ses.upper()}.txt'


time = epoch2time(ep)

if time < epoch2time([2022, 11, 27, 0, 0, 0]):
    config['atxfile'] = '../data/antex/igs14.atx'

griddef = config['griddef']

pos_ref = ecef2pos(xyz_ref)

cs = cssr()
cs.monlevel = config['cs']['monlevel']
cs.week = time2gpst(time)[0]
cs.read_griddef(griddef)

cs_ = cssr()
cs_.monlevel = config['cs']['monlevel']
cs_.week = cs.week
cs_.read_griddef(griddef)

rnx = rnxdec()
nav = Nav()
proc = process(nav, rnx, config, nep=nep, xyz_ref=xyz_ref)

sigs, nav.nf = proc.init_sig(sig_t)
rnx.setSignals(sigs)

nav = rnx.decode_nav(navfile, nav)
rnx.decode_obsh(obsfile)

# Initialize position
#
ppprtk = ppprtkpos(nav, rnx.pos, logfile=f'{ttl}.log', config=config)

proc.prepare_signal(obsfile)

# Get grid location
#
pos = ecef2pos(rnx.pos)
inet = cs.find_grid_index(pos)

if l6_mode == 1:
    fc = open(file_l6, 'rb')
    if not fc:
        nav.fout.write(f"ERROR: cannot open L6 message file {file_l6}!")
        sys_exit(-1)
else:
    dtype = [('wn', 'int'), ('tow', 'int'), ('prn', 'int'),
             ('type', 'int'), ('len', 'int'), ('nav', 'S500')]
    v = np.genfromtxt(file_l6, dtype=dtype)

# Skip epoch until start time
#
obs = rnx.decode_obs()
while time > obs.t and obs.t.time != 0:
    obs = rnx.decode_obs()

msg, msg2 = None, None

for ne in range(nep):
    week, tow = cs.set_time(obs.t)  # set time for reference

    if ne == 0:
        proc.init_time(obs.t)

    if l6_mode == 1:  # from log file
        cs.decode_l6msg(fc.read(250), 0)
        if cs.fcnt == 5:  # end of sub-frame
            cs.week = week
            cs.decode_cssr(cs.buff, 0)
    else:
        msg = decode_msg(v, tow, l6_ch, prn_p1)

        if msg is not None:
            cs.decode_l6msg(msg, 0)
            if cs.fcnt == 5:  # end of sub-frame
                cs.decode_cssr(bytes(cs.buff), 0)

        if prn_p2 > 0:
            msg2 = decode_msg(v, tow, l6_ch, prn_p2)

        if msg2 is not None:
            cs_.decode_l6msg(msg2, 0)
            if cs_.fcnt == 5:  # end of sub-frame
                cs_.decode_cssr(bytes(cs_.buff), 0)
                cs.merge_cssr(cs_)

    cstat = cs.chk_stat()
    if cstat:
        ppprtk.process(obs, cs=cs)

    proc.save_output(obs.t, ne)  # save output

    # Get new epoch, exit after last epoch
    #
    obs = rnx.decode_obs()
    if obs.t.time == 0:
        break

if l6_mode == 1:
    fc.close()

proc.close()
proc.plot(ttl, fig_type=1)
