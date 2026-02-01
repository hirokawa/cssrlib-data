"""
 static test for PPP (Galileo HAS IDD, JPL GDGPS HAS)
"""

import os
from cssrlib.gnss import Nav, load_config, time2doy, timediff, epoch2time
from cssrlib.cssrlib import sCSSRTYPE, sCType
from cssrlib.rtcm import rtcm
from cssrlib.pppssr import pppos
from cssrlib.rinex import rnxdec
from cssrlib.utils import process

config = load_config('config.yml')

# Select test case
#
icase = 3
nep = 900*4-5

ttl = 'test_ppprtcm'
navfile = None
file_rtcm = None
# Start epoch and number of epochs
#
if icase == 1:  # Galileo HAS IDD

    ep = [2023, 8, 17, 2, 0, 0]
    navfile = '../data/doy2023-229/OBE42023229c.nav'
    # navfile = '../data/brdc/BRD400DLR_S_20232290000_01D_MN.rnx'
    obsfile = '../data/doy2023-229/OBE42023229c.obs'
    xyz_ref = [4186704.2262, 834903.7677, 4723664.9337]
    file_rtcm = '../data/doy2023-229/idd2023229c.rtc'
    file_rtcm_log = '../data/doy2023-229/idd2023229c.log'
    cs_mask = 1 << sCType.CLOCK | 1 << sCType.ORBIT | 1 << sCType.CBIAS

elif icase == 2:  # JPL GDGPS  Mosaic-X5

    ep = [2024, 2, 12, 7, 0, 0]
    xyz_ref = [-3962108.7007, 3381309.5532, 3668678.6648]
    file_rtcm = '../data/doy2024-043/JPL32T2043h.rtcm3'
    file_rtcm_log = '../data/doy2024-043/JPL32T2043h.log'
    cs_mask = 1 << sCType.CLOCK | 1 << sCType.ORBIT | 1 << sCType.CBIAS

elif icase == 3:  # JPL GDGPS (w/o code bias) JAVAD DELTA-3S

    ep = [2025, 8, 21, 7, 0, 0]
    xyz_ref = [-3962108.6836, 3381309.5672, 3668678.6720]
    # SSRA11JPL0 GPS+GAL orbit+clock corrs
    file_rtcm = '../data/doy2025-233/jpl233h.rtcm3'
    file_rtcm_log = '../data/doy2025-233/jpl233h.log'
    cs_mask = 1 << sCType.CLOCK | 1 << sCType.ORBIT

time = epoch2time(ep)
year = ep[0]
doy = int(time2doy(time))
ses = chr(ord('a')+ep[3])

if navfile is None:
    navfile = f'../data/doy{year}-{doy:03d}/{doy:03d}{ses}_rnx.nav'
    obsfile = f'../data/doy{year}-{doy:03d}/{doy:03d}{ses}_rnx.obs'

# Define signals to be processed
#
if icase in [1, 2]:
    sig_t = {'G': ['1C', '2W'], 'E': ['1C', '7Q']}  # GE


elif icase in [3, 4]:
    sig_t = {'G': ['1W', '2W'], 'E': ['1C', '7Q']}  # GE

rnx = rnxdec()
nav = Nav()
proc = process(nav, rnx, config, nep=nep, xyz_ref=xyz_ref)

sigs, nav.nf = proc.init_sig(sig_t)
rnx.setSignals(sigs)

# Positioning mode
# 0:static, 1:kinematic
#
nav.pmode = 0

# Decode RINEX NAV data
#
nav = rnx.decode_nav(navfile, nav)

cs = rtcm(file_rtcm_log)
cs.monlevel = 1
cs.cssrmode = sCSSRTYPE.RTCM3_SSR
cs.inet = 0

if icase in [2, 3]:  # mask phase-bias for JPL GDGPS
    cs.mask_pbias = True

if True:
    fc = open(file_rtcm, 'rb')
    if not fc:
        print("RTCM message file cannot open.")

    blen = os.path.getsize(file_rtcm)
    msg = fc.read(blen)
    maxlen = len(msg)-5
    fc.close()

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

k = 0
# Loop over number of epoch from file start
#
for ne in range(nep):
    week, tow = cs.set_time(obs.t)  # set time for reference

    # Set initial epoch
    #
    if ne == 0:
        # nav.t = deepcopy(obs.t)
        # t0 = deepcopy(obs.t)
        # t0.time = t0.time//30*30
        # nav.time_p = t0
        proc.init_time(obs.t)

    while True:
        stat = cs.sync(msg, k)
        if stat is False:
            k += 1
            continue
        if not cs.checksum(msg, k, maxlen):
            k += 1
            continue

        tc = cs.decode_time(msg[k:k+cs.len+3])
        if (tc is not False) and timediff(tc, obs.t) > 0:
            break

        _, _, eph, geph, seph = cs.decode(msg[k:k+cs.len+3])
        k += cs.dlen

        if cs.msgtype in cs.eph_t.values():
            nav.eph.append(eph)

    # Call PPP module with HAS corrections
    #
    if (cs.lc[0].cstat & cs_mask) == cs_mask:
        ppp.process(obs, cs=cs)

    proc.save_output(obs.t, ne, ppp)  # save output

    # Get new epoch, exit after last epoch
    #
    obs = rnx.decode_obs()
    if obs.t.time == 0:
        break

proc.close()
proc.plot(ttl, fig_type=1)
