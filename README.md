# CSSRlib-data

Sample scripts and datasets for **CSSRlib**. The scripts require an installation of the **CSSRlib** package.

## Prerequisites:

Pre-installation of **CSSRlib** is required. Install the [officially released package](https://pypi.org/project/cssrlib/) with `pip` using

```bash
pip install cssrlib
```

Additional python packages are required as prerequisites and can be installed via the following command

```bash
pip install -r requirements.txt
```

It is recommended to use a virtual environment to avoid version conflicts with the system-wide python installation.

## Installation in Virtual Environment

Create a virtual environment `cssrlib-venv` in the current directory with

```bash
python3 -m venv cssrlib-venv

```
Activate the environment with

```bash
source cssrlib-venv/bin/activate
```

To make the base **CSSRlib** module availabe in the virtual environment, there are several options:

***1. Install the official package***

Install the [officially released package](https://pypi.org/project/cssrlib/) with `pip`

```bash
pip install cssrlib
```

This will be sufficient for most cases, but typically the version from the github repository is more up-to-date.

***2. Install from github repository***

Install `cssrlib` from the `main` branch of the [github repository](https://github.com/hirokawa/cssrlib.git) with this command

```bash
pip install -U git+https://github.com/hirokawa/cssrlib.git@main
```

If you like to use the development branch with the most recent changes, substitute `@main`with `@devel`.

*IMPORTANT*: make sure to always use consistent branches for `cssrlib` and `cssrlib-data`!

***3. Install a local copy of the github repository***

For this installation guide, it is assumed that the **CSSRlib** base repository has been cloned to a folder `cssrlib` in the same directory as the samples repository `cssrlib-data`.  The command `ls | grep cssrlib` should then return

```bash
cssrlib
cssrlib-data
cssrlib-venv
```

Install the dependencies from the `requirements.txt` file of `cssrlib`:

```bash
pip install -r cssrlib/requirements.txt
```

Install the local module in editing mode using the `-e` flag:


```bash
pip install -e ./cssrlib
```

***4. Make a local copy available without installation***

*NOTE:* this is not necessry if the `cssrlib` package has been installed with `pip` through one of the steps above!

Make sure the path to the `src` folder of the `cssrlib` base repository appears in the python path.

```bash
echo $PYTHONPATH
```

If not, add it with the following export command, where `<path-to-cssrlib>` must be replaced with the full path to **CSSRlib** base repository.

```bash
export PYTHONPATH="$PYTHONPATH:<path-to-cssrlib>/src"
```

## Check **CSSRlib** base package installation

Check if the package inclusion works with

```bash
python -c "import cssrlib; print(cssrlib.__version__); print(cssrlib.__file__)"
```

If the `cssrlib` package was installed with `pip`, you can display the installation details with

```bash
pip show cssrlib
```

## Install dependencies for data repository

Install the dependencies from the `requirements.txt` file of `cssrlib-data`:

```bash
pip install -r cssrlib-data/requirements.txt
```

## Handling Warnings from the PySolid module

The module `pysolid` is used for the computation of solid Earth tides. It contains a hard-coded leap second table with an expiration date, which is set to the next possible injection date of a leap second at the time of the last update. The table is frequently updated by the package maintainers. The following warning is issued when the expiration date is exceeded:

> Mild Warning -- time crossed leap second table boundaries.  Boundary edge value used instead

If you encounter this warning when executing **CSSRlib** scripts, it can most likely be removed by updating `pysolid` to the most recent version using

```bash
pip install --upgrade pysolid
```

If the warnings persist, installing the latest version from the repository can help

```bash
pip install -U git+https://github.com/insarlab/PySolid.git@main
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
