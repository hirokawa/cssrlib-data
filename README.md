# cssrlib-data
Sample scripts and dataset for CSSRlib
Pre-installation of CSSRlib is required.

## Prerequisites:

Additional python packages are required as prerequisites and can be installed via the following command

```
pip install -r requirements.txt
```

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
- [ ] Add residual output
- [ ] Add check for observations, first two observations must be on different frequencies
- [ ] Number of frequencies `nav.nf` should be set automatically depending on specified signals
- [x] Check triple-frequency processing for PPP
- [x] Code/phase error ratio `nav.erratio` should also work for > 2 signals
