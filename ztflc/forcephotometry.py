#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas
import numpy as np
from .diffdata import DiffData


def run_forcephotometry(target, download_files=True, nprocess=4, verbose=True,
                        filecheck=True, ignore_warnings=False,  **kwargs):
    """ """
    fp = ForcePhotometry.from_name(target)
    fp.load_metadata()
    if download_files:
        fp.io.download_data(nprocess=nprocess)
    fp.load_filepathes(filecheck=filecheck, ignore_warnings=ignore_warnings)
    fp.run_forcefit(update_diffdata=False,verbose=verbose, **kwargs)
    
    

class ForcePhotometry():
    """ """
    # =============== #
    #    Init         #
    # =============== #
    def __init__(self, ztftarget):
        """ """
        self.io = ztftarget
        
    @classmethod
    def from_name(cls, name):
        """ Load the ForcePhotometry from Target Name """
        from .io import ZTFTarget
        return cls(ZTFTarget.from_name(name))

    @classmethod
    def from_coords(cls, ra, dec, jdmin=None, jdmax=None,  name="unknown"):
        """ Load the ForcePhotometry from Target Name """
        from .io import ZTFTarget
        return cls(ZTFTarget.from_coords(ra, dec, jdmin=jdmin, jdmax=jdmax, name=name))
        
    # =============== #
    #   Method        #
    # =============== #
    # -------- #
    # LOADDER  #
    # -------- #
    def load_metadata(self, **kwargs):
        """ """
        if not hasattr(self.io,"marshal"):
            self.io.load_marshal(**kwargs)
        self.io.load_metadata()

    def load_filepathes(self, **kwargs):
        """ 
        **kwargs eventually goes to ztfquery.Query.get_local_data()
                 -> exists=True, indexes=None, filecheck=True, ignore_warnings=False, etc.

        """
        self._filepathes = self.io.get_diffimg_forcepsf_filepath(**kwargs)
        self._diffdata = None

    def store(self, filename=None):
        """ """
        if filename is None:
            from .io import LOCALDATA
            filename = LOCALDATA+"/%s.csv"%self.io.name
        self._data_forcefit.to_csv(filename, index=False)
        
    # -------- #
    # FITTER   #
    # -------- #
    def run_forcefit(self, indexes=None, update_diffdata=False,
                     store=False, verbose=True, no_badsub=False):
        """ """
        import gc
        if indexes is None:
            indexes = np.arange(self.nsources)

        dataout = {}
        if verbose:
            print("Starting run_forcefit() for %d image differences"%len(indexes))
            
        for i in indexes:            
            if verbose:
                print("running %d "%i)
                print(self.filepathes[i][0].split("/")[-1])
            diffdata = self.get_ith_diffdata(i, update=update_diffdata)
            has_nan = np.any(np.isnan(diffdata.diffimg))
            if has_nan and no_badsub:
                print("NaNs in the image, skipped")
            else:
                fitresults = diffdata.fit_flux()
                datainfo   = diffdata.get_main_info()
                dataout[i] = {**fitresults,**datainfo}
                dataout[i]["data_hasnan"] = has_nan
            del diffdata
            gc.collect()
            
        self._data_forcefit = pandas.DataFrame(dataout).T
        if store:
            self.store()
            
    def get_ith_diffdata(self, index, update=False, rebuild=False, **kwargs):
        """ loads and returns a DiffData object corresponding 
        at the ith-entry of the self.filepathes 

        Parameters
        ----------
        index: [int]
            which entry iof self.filepathes/self.diffdata do you want ?

        update: [bool] -optional-
            if the DiffData is derived (because unknown or rebuild), should this
            be updated to self.diffdata

        rebuild: [bool] -optional-
            Even if self.DiffData[index] exists, shall this remeasure it ?

        Returns
        -------
        DiffData
        """
        if self.diffdata is None or self.diffdata[index] is None or rebuild:
            diffdata = DiffData(*self.filepathes[index], self.io.get_coordinate(), **kwargs)
            if update:
                self.diffdata[index] = diffdata
            return diffdata

        return self.diffdata[index]
            

    # --------- #
    #  PLOTTER  #
    # --------- #
    def show_lc(self, ax=None, scalezp=25, reference=False, **kwargs):
        """ """
        # - Figure
        if ax is None:
            import matplotlib.pyplot as mpl
            fig = mpl.figure()
            ax = fig.add_subplot(111)
        else:
            fig = ax.figure

        x, y, err = self.data_forcefit[["obsmjd","ampl","ampl.err"]].values.T
        f0coef = 10**(-(self.data_forcefit["magzp"]-scalezp)/2.5) if scalezp is not None else 1
    
        for i,band_ in [[1,"C0"],[2,"C2"],[3,"C1"]]:
            if i not in self.data_forcefit["filterid"]:
                continue
            flag = np.asarray(self.data_forcefit["filterid"]==i, dtype="bool")
            
            prop = {**dict(marker="o", color=band_, zorder=5), **kwargs}
            ax.scatter(x[flag],(y*f0coef)[flag], **prop)
            prop["zorder"] = prop.get("zorder", 5)-1
            ax.errorbar(x[flag],(y*f0coef)[flag], yerr=err[flag], 
                                ls="None",ecolor="0.6",**prop)

        ax.axhline(0, ls="--", color="0.5")
        return ax
    # =============== #
    #  Properties     #
    # =============== #
    @property
    def data_forcefit(self):
        """ DataFrame containing the forcefit results and data info """
        if not hasattr(self,"_data_forcefit"):
            raise AttributeError("data_forcefit not loaded. run self.run_forcefit()")
        return self._data_forcefit
    
    @property
    def nsources(self):
        """ """
        return len(self.filepathes)
    
    @property
    def metadata(self):
        """ """
        return self.io.zquery.metatable

    @property
    def filepathes(self):
        """ List of DiffImage and PSFDiffImage file path """
        if not hasattr(self,"_filepathes"):
            self.load_filepathes()
        return self._filepathes

    def had_filepathes(self):
        """ """
        return hasattr(self,"_filepathes") and self._filepathes is not None
    
    @property
    def diffdata(self):
        """ DiffData object associated with the filepathes"""
        if not hasattr(self, "_diffdata") or self._diffdata is None:
            if self.had_filepathes():
                self._diffdata = [None for i in range(self.nsources)]
            else:
                self._diffdata = None
                
        return self._diffdata
