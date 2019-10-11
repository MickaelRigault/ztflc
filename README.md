# ztflc
LightCurve Estimation from ZTF data


### Credit

M. Rigault (corresponding author, m.rigault@ipnl.in2p3.fr, CNRS/IN2P3). 
A similar code exists [here](https://github.com/yaoyuhan/ForcePhotZTF), used for the ZTF high cadence SNeIa paper [Yao et al. 2019](http://cdsads.u-strasbg.fr/abs/2019arXiv191002967Y). 

`ztflc` is however solely based on [ztfquery](https://github.com/MickaelRigault/ztfquery).

### Acknowledgment

If you have used `ztflc` for a research you are publishing, please **include the following in your acknowledgments**:
_"The ztflc code was funded by the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (grant agreement n°759194 - USNAC, PI: Rigault)."_

# Usage

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
fp.fp.data_forcefit
```
and use `fp.show_lc()` to plot the lightcurve.


### Example 2, from ztf target name:
Say you want to get the force photometry from a ztf target named `ZTF18aaaaaaa`. Use the same code as example 1, except that, this time, you load `fp` using:
```python 
from ztflc import forcephotometry
# Setup the target
fp = forcephotometry.ForcePhotometry.from_name(ZTF18aaaaaaa)
```

the code will use `ztfquery.marshal` to recover the information you need (ra, dec jdmin, jdmax).

