#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import astrobject
from astropy.io import fits
from scipy import stats
from pixelproject import grid

class DiffData():
    """ """
    def __init__(self, diffimgpath, psfimgpath, radec=None, inpixels=False):
        """ """
        self.set_data(diffimgpath, psfimgpath)
        self.set_target_coords(*radec, inpixels=inpixels)

    # ================ #
    #    Methods       #
    # ================ #
    # ------- #
    # FITTER  #
    # ------- #
    def fit_flux(self, use_datanoise=False):
        """ """
        from .fitter import DiffImgFitter
        if not hasattr(self,"grid_diffaligned"):
            _ = self.get_aligned_grid()
            
        # Load the fitter
        self.fitter = DiffImgFitter(self.grid_diffaligned.geodataframe["data"].values, 
                                    np.ravel(self.psfimg),
                                    self.estimated_noise if use_datanoise else 0,
                                    shape= self.psfshape)
        
        # Load the fitter
        amplitude_guess = np.max(self.fitter.data)*(np.sqrt(2*np.pi*2**2)) # 2 sigma guess
        sigma_int_guess = 1 if use_datanoise else np.std(self.fitter.data)
        # Return fitter output
        return self.fitter.fit(ampl_guess=amplitude_guess, sigma_guess=sigma_int_guess)
    
    # ------- #
    # SETTER  #
    # ------- #
    def set_target_coords(self, ra,dec, inpixels=False):
        """ set the target coordinates.
        
        Parameters
        ----------
        ra,dec: [float, float]
            Coordinates (in deg or in pixels)

        inpixels: [bool] -optional-
            Are the coordinates in pixels [True] or in RaDEC [False] ? 

        Returns
        -------
        None
        """
        if not inpixels:
            self.radec = ra,dec
        else:
            self._xy = ra,dec

    def set_data(self, diffimgpath, psfimgpath, **kwargs):
        """ """
        self.diffimg = astrobject.get_image(diffimgpath, index=1, background=0, ascopy=True, **kwargs)
        self.psfimg  = fits.getdata(psfimgpath)

    # ------- #
    # GETTER  #
    # ------- #
    def get_main_info(self):
        """ """
        header = {k.lower():self.header[k] for k in self._main_columns}
        x,y = self.target_position
        header["target_x"] = x
        header["target_y"] = y
        return header
    
    def get_aligned_grid(self, update=True):
        """ """
        aligned_diffgrid = self.grid_diffimg.project_to(self.grid_psf)
        if update:
            self._grid_diffaligned = aligned_diffgrid
        return aligned_diffgrid

    def estimate_noise(self, buffer=10, use="nmad"):
        """ """
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
        
    # ------- #
    # PLOTTER #
    # ------- #
    def show_alignment(self, ax=None):
        """ """
        # - Data
        ax = self.grid_diffimg.show("data",ax=ax, lw=0)
        ax.scatter(*self.target_position, marker="x", color="C1", s=80)
        # - Grid out
        _ = self.grid_psf.show(None, ax=ax, lw=0.5)

    def show_noise(self, ax=None):
        """ """
        import matplotlib.pyplot as mpl
        from matplotlib.patches import Circle
        fig    = mpl.figure(figsize=[9,3])
        axleft = fig.add_axes([0.05,0.12,0.3, 0.8], frameon=False)
        ax     = fig.add_axes([0.4,0.12,0.45, 0.78])
        
        dataout = self.grid_diffimg.geodataframe[ self.noise_data["flagnoise"]]["data"]
        datain  = self.grid_diffimg.geodataframe[~self.noise_data["flagnoise"]]["data"]

        _ = self.grid_diffimg.show("data", ax=axleft, lw=0)
        axleft.add_patch(Circle(self.target_position,self.noise_data["buffer"], facecolor="None", edgecolor="C1", lw=2))
        axleft.scatter(*self.target_position, marker="x", color="C1", s=80)

        range = [int(self.noise_data["median"]-5*self.noise_data["nmad"]), int(self.noise_data["median"]+5*self.noise_data["nmad"])]
        ax.hist(dataout, density=True, range=range, bins=15, facecolor="C0", histtype="step", fill=True, alpha=0.4,
                    label="out data")
        ax.hist(dataout, density=True, range=range, bins=15, edgecolor="C0", histtype="step", lw=2)
        ax.hist(datain, histtype="step", density=True, range=[range[0],range[1]*2], bins=20, color="C1", lw=2,
                    label="in data")

        mean,std = np.mean(dataout),np.std(dataout)
        xx = np.linspace(*range, 100)
        ax.plot(xx, stats.norm.pdf(xx, mean, std), ls="--", color="C0")
        ax.axvline(mean, label=r"$%.2f \pm %.2f$"%(self.noise_data["mean"],self.noise_data["mean.err"]), color="C0")

        fig.text(0.01,0.98, "Estimated error: %.2f"%self.estimated_noise, va="top", ha="left")
        
        ax.legend()
    # ================ #
    #  Properties      #
    # ================ #
    # Grids
    @property
    def grid_diffaligned(self):
        """ Re-aligned diffimg (aligned with PSF image, that is the target is in the center of a pixel)  """
        if not hasattr(self,"_grid_diffaligned"):
            raise AttributeError("No aligned grid derived yet. Run self.get_aligned_grid()")
        return self._grid_diffaligned
        
    @property
    def grid_diffimg(self):
        """ PixelProjet Grid of the input difference image"""
        if not hasattr(self, "_grid_diffimg"):
            x,y = self.target_position
            buffer = (np.asarray(self.psfshape))/2+5
            xmin,xmax = int(x-buffer[0]),int(x+buffer[0])
            ymin,ymax = int(y-buffer[1]),int(y+buffer[1])
            # Data Grid
            self._grid_diffimg = grid.Grid.from_stamps(self.diffimg.data[ymin:ymax,xmin:xmax ], origin=[xmin,ymin])
            
        # - Defined, let's go
        return self._grid_diffimg

    @property
    def grid_psf(self):
        """ PixelProjet Grid of the input PSF image """
        if not hasattr(self, "_grid_psf"):
            x,y = self.target_position
            # Data Grid Centered at the target location
            self._grid_psf = grid.get_simple_grid(*self.psfshape, shift_origin=(np.asarray([x,y]) - np.asarray(self.psfshape)/2+0.5))
            
        # - Defined, let's go
        return self._grid_psf
            
    # Others
    @property
    def estimated_noise(self):
        """ """
        if not hasattr(self,"noise_data"):
            self.estimate_noise(10)
        return self.noise_data["differr"]
    
    @property
    def target_position(self):
        """ """
        if hasattr(self,"_xy"):
            return self._xy
        
        if not hasattr(self,"radec"):
            raise AttributeError("No coordinate set. see self.set_target_radec")

        if not hasattr(self,"diffimg"):
            raise AttributeError("No diff image (that has the wcs solution)")

        self._xy = self.diffimg.coords_to_pixel(*self.radec)
        return self._xy

    @property
    def psfshape(self):
        """ shape of the psf image """
        return np.shape(self.psfimg)
    
    @property
    def header(self):
        """ Data header"""
        return self.diffimg.header

    @property
    def _main_columns(self):
        """ """
        return ["FILENAME","HUMIDITY","FILTER","OBSMJD","CCDID","AMP_ID","GAIN","READNOI","DARKCUR",
                "MAGZP","MAGZPUNC","MAGZPRMS", "CLRCOEFF","CLRCOUNC","ZPCLRCOV","ZPMED", "ZPAVG",
                "ZPRMSALL","CLRMED","CLRAVG","CLRRMS","QID","RCID","SEEING","MAGLIM","STATUS",
                "FILTERID","FILTER", "FIELDID"]
