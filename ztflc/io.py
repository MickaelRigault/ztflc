#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import pandas
from astropy.io import fits

from ztfquery import query, marshal
from ztfquery.io import LOCALSOURCE

LOCALDATA = LOCALSOURCE + "forcephotometry/"

import os

if os.path.isfile(LOCALSOURCE + "/externaldata/target_position.csv"):
    TARGET_IO = pandas.read_csv(LOCALSOURCE + "/externaldata/target_position.csv")
else:
    TARGET_IO = pandas.DataFrame(
        columns=["name", "mean_ra", "mean_dec", "median_ra", "median_dec"]
    )


class ZTFTarget(object):
    """ """

    # ================ #
    #   Init           #
    # ================ #
    def __init__(self, target=None):
        """ """
        self._zquery = query.ZTFQuery()
        if target is not None:
            self.set_name(target)

    @classmethod
    def from_name(cls, name):
        """ """
        return cls(name)

    @classmethod
    def from_coords(cls, ra, dec, jdmin=None, jdmax=None, name="unknown"):
        """ """
        this = cls()
        this.set_coordinate(ra, dec)
        this.set_jdrange(jdmin, jdmax)
        this.set_name(name)
        return this

    # ================ #
    #  Methods         #
    # ================ #

    # ------- #
    # SETTER  #
    # ------- #
    def set_name(self, targetname):
        """ """
        self._name = targetname

    def set_coordinate(self, ra, dec):
        """ provide the coordinate.
        Get them using self.get_coordinate().
        """
        self._radec = [ra, dec]

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
            self.update_marshal(program=program)

    def update_marshal(self, program="Cosmology", store=True):
        """ """
        self._marshal = marshal.MarshalAccess()
        self._marshal.load_target_sources(program)
        if store:
            self._marshal.store()

    def load_metadata(self, fromname=False):
        """ """
        if self.has_radec() and not fromname:
            # Test if one radec value if given. If not, use 
            # the first one (it does not matter here which
            # one is used)
            if isinstance(self.radec[0], list):
                radec = [self.radec[0][0], self.radec[1][0]]
            else:
                radec = self.radec

            dictmetadata = self.marshal.get_metadataquery(
                *radec, *self.jdrange, size=0.01
            )

        elif self.has_name():
            dictmetadata = self.marshal.get_target_metadataquery(self.name)
        else:
            raise AttributeError("Cannot load the metadata, no self.name and no _radec")

        self.zquery.load_metadata(**dictmetadata)

    def download_data(
        self, which=["scimrefdiffimg.fits.fz", "diffimgpsf.fits"], nprocess=4, **kwargs
    ):
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
            self.zquery.download_data(
                suffix, nprocess=nprocess, **{**{"show_progress": False}, **kwargs}
            )

    # ------- #
    # GETTER  #
    # ------- #
    def get_coordinate(self, forcename=False, use="mean"):
        """ """
        if self.has_radec() and not forcename:
            return self.radec

        if len(np.atleast_1d(self.name)) == 1:
            if self.name in TARGET_IO["name"].values:
                return TARGET_IO[TARGET_IO["name"] == self.name][
                    ["%s_ra" % use, "%s_dec" % use]
                ].values[0]

            return self.marshal.get_target_coordinates(self.name).values[0]
        return self.marshal.get_target_coordinates(self.name).values

    def get_jdrange(self, forcetarget=False):
        """ """
        if hasattr(self, "_jdrange") and not forcetarget:
            return self._jdrange

        if len(np.atleast_1d(self.name)) == 1:
            return self.marshal.get_target_jdrange(self.name)
        return self.marshal.get_target_jdrange(self.name)

    # Others
    def get_diffimg_forcepsf_filepath(self, exists=True, indexes=None, **kwargs):
        """
        **kwargs eventually goes to ztfquery.Query.get_local_data()
                 -> filecheck=True, ignore_warnings=False, etc.
        """
        diffimg = self.get_diffimg_filepath(exists=exists, indexes=indexes, **kwargs)
        diffpsf = self.get_diffpsf_filepath(exists=exists, indexes=indexes, **kwargs)
        out = []
        for diff in diffimg:
            expected_psf = diff.replace("scimrefdiffimg.fits.fz", "diffimgpsf.fits")
            if expected_psf not in diffpsf:
                out.append([diff, None])
            else:
                out.append([diff, expected_psf])
        return out

    def get_diffimg_filepath(self, exists=True, indexes=None, **kwargs):
        """
        **kwargs goes to ztfquery.Query.get_local_data()
                 -> filecheck=True, ignore_warnings=False, etc.
        """
        return self.zquery.get_local_data(
            "scimrefdiffimg.fits.fz", exists=exists, indexes=indexes, **kwargs
        )

    def get_diffpsf_filepath(self, exists=True, indexes=None, **kwargs):
        """
        **kwargs goes to ztfquery.Query.get_local_data()
                 -> filecheck=True, ignore_warnings=False, etc.
        """
        return self.zquery.get_local_data(
            "diffimgpsf.fits", exists=exists, indexes=indexes, **kwargs
        )

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
        if not hasattr(self, "_marshal"):
            raise AttributeError("Marshal not loaded. See load_marshal()")
        return self._marshal

    @property
    def name(self):
        """ """
        if not self.has_name():
            raise AttributeError("Target not loaded. See set_name()")
        return self._name

    @property
    def jdrange(self):
        """ """
        if not self.has_jdrange():
            return None
        return self._jdrange

    @property
    def radec(self):
        """ """
        if not self.has_radec():
            return None
        return self._radec

    def has_name(self):
        """ Test if the current instance has a target set """
        return hasattr(self, "_name") and self._name is not None

    def has_radec(self):
        """ """
        return hasattr(self, "_radec") and self._radec is not None

    def has_jdrange(self):
        """ """
        return hasattr(self, "_jdrange")
