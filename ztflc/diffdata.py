#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import warnings
import numpy as np
from astropy.io import fits
from scipy import stats


class DiffData(object):
    """ """

    def __init__(self, diffimgpath, psfimgpath, coords, inpixels=False, clean=True):
        """ """
        self.set_data(diffimgpath, psfimgpath, coords, inpixels=inpixels, clean=clean)

    # ================ #
    #    Methods       #
    # ================ #
    # ------- #
    # FITTER  #
    # ------- #
    def fit_flux(self, **kwargs):
        """ """
        from .fitter import DiffImgFitter
        from astropy.stats import mad_std

        # Load the fitter
        self.fitter = DiffImgFitter(self.diffimg, self.psfimg, 0, shape=self.psfshape)
        robust_nmad = mad_std(self.fitter.data[~np.isnan(self.fitter.data)])
        fit_prop = {
            "ampl_guess": np.nanmax(self.fitter.data)
            * (np.sqrt(2 * np.pi * 2**2)),  # 2 sigma guess
            "sigma_guess": robust_nmad,
        }
        fit_prop["ampl_boundaires"] = [
            -2 * np.nanmin(self.fitter.data) * (np.sqrt(2 * np.pi * 2**2)),
            5 * np.nanmax(self.fitter.data) * (np.sqrt(2 * np.pi * 2**2)),
        ]
        fit_prop["sigma_boundaries"] = [
            robust_nmad / 10.0,
            np.nanstd(self.fitter.data) * 2,
        ]

        # Return fitter output
        return self.fitter.fit(**{**fit_prop, **kwargs})

    # ------- #
    # SETTER  #
    # ------- #
    def set_data(
        self, diffimgpath, psfimgpath, coords, inpixels=False, clean=True, **kwargs
    ):
        """ """
        from astropy import wcs
        from astropy.io import fits
        from scipy import ndimage
        from astropy.stats import mad_std

        self._xy = coords if inpixels else None
        #
        # Handeling too many open files
        #
        # PSF
        with fits.open(psfimgpath) as psf:
            self._psfimg_raw = psf[0].data.copy()

        with fits.open(diffimgpath) as fdiff:
            # If coords is a list, then assume different coordinates for each
            # filter
            filterid = fdiff[1].header["FILTERID"]
            imgshape = fdiff[1].data.shape

            if isinstance(coords[0], list):
                ra = coords[0][filterid - 1]
                dec = coords[1][filterid - 1]
                coords = [ra, dec]
            # x,y position
            if self._xy is None:
                wcs_ = wcs.WCS(fdiff[1].header)
                self._xy = wcs_.all_world2pix([coords], 0)[0]
                del wcs_

            # data position
            x, y = np.asarray(self._xy, dtype="int")
            buffer = (np.asarray(self.psfshape)) / 2 - 0.5
            xmin, xmax = int(x - buffer[0]), int(x + buffer[0] + 1)
            ymin, ymax = int(y - buffer[1]), int(y + buffer[1] + 1)
            self._diffimg = fdiff[1].data[ymin:ymax, xmin:xmax].copy()
            self._datatmp = fdiff[1].data.copy()
            self._diffimg_targetpos = self._xy - [xmin, ymin]
            self._header = fdiff[1].header.copy()
            self._istarget_in = (
                (x - buffer[0] > 0)
                and (x + buffer[0] < imgshape[1])
                and (y - buffer[1] > 0)
                and (y + buffer[1] < imgshape[0])
            )

        if clean:
            self._iscleaned = True
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                flagout = self._diffimg < -mad_std(self._diffimg) * 10
            self._diffimg[flagout] = np.NaN
        else:
            self._iscleaned = False

        self._psf_shift = self._diffimg_targetpos - np.asarray(self.psfshape) / 2 + 0.5
        # ::-1 because numpy arrays are not images but matrices
        self._psfimg = ndimage.interpolation.shift(
            self._psfimg_raw, self._psf_shift[::-1], order=5
        )

    # ------- #
    # GETTER  #
    # ------- #
    def get_main_info(self, backup_value=None):
        """ """
        header = {
            k.lower(): self.header.get(k) if k in self.header else backup_value
            for k in self._main_columns
        }
        header["target_x"], header["target_y"] = self.target_position
        return header

    def estimate_noise(self, buffer=10, use="nmad"):
        """ """
        raise NotImplementedError("estimate_noise is not ready yet")
        """
        from shapely import geometry
        from astropy.stats import mad_std

        circle = geometry.Point(*self.target_position).buffer(buffer)
        pixelsout, flagbuffer = self.grid_diffimg.get_pixels_in(circle, invert=True)

        dataout = self.grid_diffimg.geodataframe[flagbuffer]["data"]
        self.noise_data = {"flagnoise":flagbuffer,
                            "mean": np.mean(dataout),
                            "std": np.std(dataout),
                            "median": np.median(dataout),
                            "nmad":mad_std(dataout),
                            "buffer":buffer,
                            }

        self.noise_data["mean.err"] = self.noise_data["nmad"]/np.sqrt(len(dataout)-1)
        self.noise_data["p(zero)"]  = stats.norm.pdf(self.noise_data["mean"], loc=0, scale=self.noise_data["mean.err"])
        self.noise_data["differr"]  = self.noise_data[use]
        return self.noise_data["differr"]
        """

    # ------- #
    # PLOTTER #
    # ------- #
    def show(self):
        """ """
        import matplotlib.pyplot as mpl

        fig = mpl.figure(figsize=[8, 3])
        prop = dict(origin="lower")  # vmin=-20,vmax=60)

        ax = fig.add_subplot(131)
        ax.imshow(self.diffimg, **prop)

        ax = fig.add_subplot(132)
        ax.imshow(self.psfimg_raw, **prop)

        ax = fig.add_subplot(133)
        ax.imshow(self.psfimg, **prop)

    def show_alignment(self, ax=None):
        """ """
        # - Data
        ax = self.grid_diffimg.show("data", ax=ax, lw=0)
        ax.scatter(*self.target_position, marker="x", color="C1", s=80)
        # - Grid out
        _ = self.grid_psf.show(None, ax=ax, lw=0.5)

    def show_noise(self, ax=None):
        """ """
        import matplotlib.pyplot as mpl
        from matplotlib.patches import Circle

        fig = mpl.figure(figsize=[9, 3])
        axleft = fig.add_axes([0.05, 0.12, 0.3, 0.8], frameon=False)
        ax = fig.add_axes([0.4, 0.12, 0.45, 0.78])

        dataout = self.grid_diffimg.geodataframe[self.noise_data["flagnoise"]]["data"]
        datain = self.grid_diffimg.geodataframe[~self.noise_data["flagnoise"]]["data"]

        _ = self.grid_diffimg.show("data", ax=axleft, lw=0)
        axleft.add_patch(
            Circle(
                self.target_position,
                self.noise_data["buffer"],
                facecolor="None",
                edgecolor="C1",
                lw=2,
            )
        )
        axleft.scatter(*self.target_position, marker="x", color="C1", s=80)

        range = [
            int(self.noise_data["median"] - 5 * self.noise_data["nmad"]),
            int(self.noise_data["median"] + 5 * self.noise_data["nmad"]),
        ]
        ax.hist(
            dataout,
            density=True,
            range=range,
            bins=15,
            facecolor="C0",
            histtype="step",
            fill=True,
            alpha=0.4,
            label="out data",
        )
        ax.hist(
            dataout,
            density=True,
            range=range,
            bins=15,
            edgecolor="C0",
            histtype="step",
            lw=2,
        )
        ax.hist(
            datain,
            histtype="step",
            density=True,
            range=[range[0], range[1] * 2],
            bins=20,
            color="C1",
            lw=2,
            label="in data",
        )

        mean, std = np.mean(dataout), np.std(dataout)
        xx = np.linspace(*range, 100)
        ax.plot(xx, stats.norm.pdf(xx, mean, std), ls="--", color="C0")
        ax.axvline(
            mean,
            label=r"$%.2f \pm %.2f$"
            % (self.noise_data["mean"], self.noise_data["mean.err"]),
            color="C0",
        )

        fig.text(
            0.01,
            0.98,
            "Estimated error: %.2f" % self.estimated_noise,
            va="top",
            ha="left",
        )

        ax.legend()

    # ================ #
    #  Properties      #
    # ================ #
    #
    # Data Image
    #
    @property
    def diffimg(self):
        """Difference data patch around the target."""
        if not hasattr(self, "_diffimg"):
            self._diffimg = None
        return self._diffimg

    @property
    def psfimg_raw(self):
        """PSF data as given by the IRSA pipeline"""
        if not hasattr(self, "_psfimg_raw"):
            self._psfimg_raw = None
        return self._psfimg_raw

    @property
    def psfimg(self):
        """PSF data re-aligned with the target"""
        if not hasattr(self, "_psfimg_raw"):
            self._psfimg = None
        return self._psfimg

    @property
    def psfshape(self):
        """shape of the psf image"""
        return np.shape(self.psfimg_raw)

    #
    # target
    #
    @property
    def target_position(self):
        """ """
        if not hasattr(self, "_xy"):
            self._xy = None
        return self._xy

    @property
    def target_in(self):
        """Test if the target is enough within the image"""
        if not hasattr(self, "_istarget_in"):
            return None
        return self._istarget_in

    #
    # Additional information
    #
    @property
    def header(self):
        """Data header"""
        if not hasattr(self, "_header"):
            self._header = None
        return self._header

    @property
    def _main_columns(self):
        """ """
        return [
            "FILENAME",
            "EXPTIME",
            "HUMIDITY",
            "FILTER",
            "OBSMJD",
            "CCDID",
            "AMP_ID",
            "GAIN",
            "READNOI",
            "DARKCUR",
            "MAGZP",
            "MAGZPUNC",
            "MAGZPRMS",
            "CLRCOEFF",
            "CLRCOUNC",
            "ZPCLRCOV",
            "ZPMED",
            "ZPAVG",
            "ZPRMSALL",
            "CLRMED",
            "CLRAVG",
            "CLRRMS",
            "QID",
            "RCID",
            "SEEING",
            "AIRMASS",
            "NMATCHES",
            "MAGLIM",
            "STATUS",
            "INFOBITS",
            "FILTERID",
            "FILTER",
            "FIELDID",
            "MOONALT",
            "MOONILLF",
        ]
