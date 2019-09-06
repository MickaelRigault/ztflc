#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas
import numpy as np
from .diffdata import DiffData


def run_forcephotometry(target, download_files=True, nprocess=4, verbose=True, **kwargs):
    """ """
    fp = ForcePhotometry.from_name(target)
    fp.load_metadata()
    if download_files:
        fp.io.download_data(nprocess=nprocess)
    fp.load_filepathes()
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
    def from_coords(cls, ra, dec, jdmin=None, jdmax=None):
        """ Load the ForcePhotometry from Target Name """
        from .io import ZTFTarget
        return cls(ZTFTarget.from_coords(ra, dec, jdmin=jdmin, jdmax=jdmax))
        
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

    def load_filepathes(self,  **kwargs):
        """ """
        self._filepathes = self.io.get_diffimg_forcepsf_filepath(**kwargs)
        self._diffdata = None
        
    # -------- #
    # FITTER   #
    # -------- #
    def run_forcefit(self, indexes=None, update_diffdata=False, store=True, verbose=True):
        """ """
        import gc
        if indexes is None:
            indexes = np.arange(self.nsources)

        dataout = {}
        if verbose:
            print("Starting run_forcefit() for %d image differences"%len(indexes))
            
        for i in indexes:
            if i>244:
                print("skip for no")
                break
            
            if verbose:
                print("running %d "%i)
                print(self.filepathes[i][0].split("/")[-1])
            diffdata = self.get_ith_diffdata(i, update=update_diffdata)
            fitresults = diffdata.fit_flux()
            datainfo   = diffdata.get_main_info()
            dataout[i] = {**fitresults,**datainfo}
            
            del diffdata
            gc.collect()
            
        self._data_forcefit = pandas.DataFrame(dataout).T
        if store:
            from .io import LOCALDATA
            self._data_forcefit.to_csv(LOCALDATA+"/%s.csv"%self.io.target, index=False)
            
    def get_ith_diffdata(self, index, update=False, rebuild=False):
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
    da            Even if self.DiffData[index] exists, shall this remeasure it ?

        Returns
        -------
        DiffData
        """
        if self.diffdata is None or self.diffdata[index] is None or rebuild:
            diffdata = DiffData(*self.filepathes[index], self.io.get_coordinate())
            if update:
                self.diffdata[index] = diffdata
            return diffdata

        return self.diffdata[index]
            
        
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
