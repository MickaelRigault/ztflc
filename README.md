# ztflc
LightCurve Estimation from ZTF data


[![PyPI](https://img.shields.io/pypi/v/ztflc.svg?style=flat-square)](https://pypi.python.org/pypi/ztflc)


### Credit

M. Rigault (corresponding author, m.rigault@ipnl.in2p3.fr, CNRS/IN2P3).
A similar code exists [here](https://github.com/yaoyuhan/ForcePhotZTF), used for the ZTF high cadence SNeIa paper [Yao et al. 2019](http://cdsads.u-strasbg.fr/abs/2019arXiv191002967Y).

`ztflc` is however solely based on [ztfquery](https://github.com/MickaelRigault/ztfquery).

### Acknowledgment

If you have used `ztflc` for a research you are publishing, please **include the following in your acknowledgments**:
_"The ztflc code was funded by the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (grant agreement n°759194 - USNAC, PI: Rigault)."_

# Usage
### First time using `ztflc`

1) You may want to make sure `$ZTFDATA` is defined (see [ztfquery](https://github.com/MickaelRigault/ztfquery).)
2) You should first run *Storing or updating marshal data.* for your marshal's science program.


### Example 1, force photometry from coordinates:

For a target with coordinates `ra`, `dec`, you want to do force photometry on ZTF data from julian data `jdmin` up to julian date `jdmax`:

```python
from ztflc import forcephotometry
# Setup the target
fp = forcephotometry.ForcePhotometry.from_coords(ra, dec, jdmin=jdmin,jdmax=jdmax)
# Load the ztf metadata
fp.load_metadata()
# Load the ztf file you are going to need (diffimg and psfimg)
fp.load_filepathes()
# and run the forcephotometry
fp.run_forcefit(verbose=True)
```
Once the force photometry has run, data are stored as:
```python
fp.data_forcefit
```
and use `fp.show_lc()` to plot the lightcurve.

To store the `date_forcefit`, simply do:
```python
fp.store()
```
the dataframe will be store in `$ZTFDATA/forcephotometry/name.csv` (c.f. `ztfquery` for `$ZTFDATA`). you can also directly provide the filename if you want to dataframe stored elsewhere `fp.store(filename)`.


### Example 2, from ztf target name:
Say you want to get the force photometry from a ztf target named `ZTF18aaaaaaa`. Use the same code as example 1, except that, this time, you load `fp` using:
```python
from ztflc import forcephotometry
# Setup the target
fp = forcephotometry.ForcePhotometry.from_name(ZTF18aaaaaaa)
```

the code will use `ztfquery.marshal` to recover the information you need (ra, dec jdmin, jdmax).

# Storing or updating marshal data.

You can store ztf mashal target information locally using `ztfquery`. For instance, if you want to store the "Cosmology" program, simply do:

```python
from ztfquery import marshal
m = marshal.MarshalAccess()
m.load_target_sources("Cosmology")
m.store()
```
Once you did that, the code will use the locally stored data when you use `forcephotometry.ForcePhotometry.from_name(ZTF18aaaaaaa)`

# Downloading the files you need.

The first time you are going to run the forcephotometry for a given target, you most likely need to download the associated data. For instance for `ZTF18aaaaaaa`

```python
from ztflc import forcephotometry
# Setup the target
fp = forcephotometry.ForcePhotometry.from_name(ZTF18aaaaaaa)
# Download the data:
fp.io.download_data()
```
