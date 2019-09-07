#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
from astropy.io import fits

from ztfquery import query, marshal
from ztfquery.io import LOCALSOURCE
LOCALDATA = LOCALSOURCE+"forcephotometry/"

class ZTFTarget( object ):
    """ """

    # ================ #
    #   Init           #
    # ================ #
    def __init__(self, target=None):
        """ """
        self._zquery = query.ZTFQuery()
        if target is not None:
            self.set_target(target)

    @classmethod
    def from_name(cls, name):
        """ """
        return cls(name)

    @classmethod
    def from_coords(cls, ra, dec, jdmin=None, jdmax=None):
        """ """
        this = cls()
        this.set_coordinate(ra, dec)
        this.set_jdrange(jdmin, jdmax)
        return this
    
    # ================ #
    #  Methods         #
    # ================ #

    # ------- #
    # SETTER  #
    # ------- #
    def set_target(self, targetname, load_metadata=True):
        """ """
        self._target = targetname

    def set_coordinate(self, ra, dec):
        """ provide the coordinate.
        Get them using self.get_coordinate().
        """
        self._radec = [ra,dec]

    def set_jdrange(self, jdmin, jdmax):
        """ """
        self._jdrange = [jdmin, jdmax]
        
    # ------- #
    # LOADDER #
    # ------- #    
    def load_marshal(self, source="local", program="Cosmology"):
        """ """
        if source == "local":
            self._marshal = marshal.MarshalAccess.load_local(program)
        else:
            self.update_marshal()

    def update_marshal(self, program="Cosmology", store=True):
        """ """
        self._marshal = marshal.MarshalAccess()
        self._marshal.load_target_sources(program)
        if store:
            self._marshal.store()
        
    def load_metadata(self, fromtarget=True):
        """ """
        if fromtarget and self.has_target():
            dictmetadata = self.marshal.get_target_metadataquery(self.target)
        else:
            dictmetadata = self.marshal.get_metadataquery(*self._radec, *self._jdrange, size=0.01)
            
        self.zquery.load_metadata(**dictmetadata)

    def download_data(self, which=["scimrefdiffimg.fits.fz","diffimgpsf.fits"],
                          nprocess=4, **kwargs ):
        """ 

        **kwargs goes to zquery.download_data()
              - source='IRSA',
              - indexes=None,
              - download_dir=None,
              - show_progress=False,
              - notebook=False,
              - verbose=True,
              - nodl=False,
              - overwrite=False,
              - auth=None,
              **kwargs
        """
        for suffix in np.atleast_1d(which):
            self.zquery.download_data(suffix, nprocess=nprocess, **{**{"show_progress":False},
                                                                    **kwargs})
        
    # ------- #
    # GETTER  #
    # ------- #
    def get_coordinate(self, forcetarget=False):
        """ """
        if hasattr(self,"_radec") and not forcetarget:
            return self._radec
        
        if len( np.atleast_1d(self.target) )==1:
            return self.marshal.get_target_coordinates(self.target).values[0]
        return self.marshal.get_target_coordinates(self.target).values

    def get_jdrange(self, forcetarget=False):
        """ """
        if hasattr(self,"_jdrange") and not forcetarget:
            return self._jdrange
        
        if len( np.atleast_1d(self.target) )==1:
            return self.marshal.get_target_jdrange(self.target)
        return self.marshal.get_target_jdrange(self.target)
    
    # Others 
    def get_diffimg_forcepsf_filepath(self, exists=True, indexes=None):
        """ """
        diffimg = self.get_diffimg_filepath(exists=exists, indexes=indexes)
        diffpsf = self.get_diffpsf_filepath(exists=exists, indexes=indexes)
        out = []
        for diff in diffimg:
            expected_psf = diff.replace("scimrefdiffimg.fits.fz", "diffimgpsf.fits")
            if expected_psf not in diffpsf:
                out.append([diff, None])
            else:
                out.append([diff, expected_psf])
        return out

    def get_diffimg_filepath(self, exists=True, indexes=None):
        """ """
        return self.zquery.get_local_data("scimrefdiffimg.fits.fz", exists=exists, indexes=indexes)

    def get_diffpsf_filepath(self, exists=True, indexes=None):
        """ """
        return self.zquery.get_local_data("diffimgpsf.fits", exists=exists, indexes=indexes)
    
    # ================ #
    #  Properties      #
    # ================ #
    @property
    def zquery(self):
        """ """
        return self._zquery

    @property
    def marshal(self):
        """ """
        if not hasattr(self,"_marshal"):
            raise AttributeError("Marshal not loaded. See load_marshal()")
        return self._marshal
    
    @property
    def target(self):
        """ """
        if not self.has_target():
            raise AttributeError("Target not loaded. See set_target()")
        return self._target
    
    def has_target(self):
        """ Test if the current instance has a target set """
        return hasattr(self,"_target") and self._target is not None
