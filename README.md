# cssrlib-data
Sample scripts and dataset for CSSRlib
Pre-installation of CSSRlib is required.

## Ephemeris: RINEX/TLE
- test_eph.py reading/plotting ephemeris from RINEX 3
- test_tlesim.py reading/plotting TLE orbit

## Model
- test_trop.py tropospheric delay model
- cacode.py GPS/QZSS C/A code simulation

## Positioning: Static
- test_pntpos.py Standalone positioning
- test_rtk.py RTK positioning
- test_ppprtk.py PPP-RTK positioning

## Positioning: Kinematic
- test_pntpos2.py Standalone positioning
- test_rtk2.py RTK positioning
- test_ppprtk2.py PPP-RTK positioning

## ToDo-List

- [ ] Implement pole tide displacements
- [ ] Check and improve observation noise settings
- [ ] Check slant iono process noise settings
- [ ] Add residual output
- [ ] Add check for observations, must be on different frequencies