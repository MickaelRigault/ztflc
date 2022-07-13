#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
from scipy import stats
import warnings
from iminuit import Minuit

from .utils import make_method


class DiffImgFitter(object):

    FREEPARAMETERS = ["sigma", "ampl"]

    def __new__(cls, *arg, **kwargs):
        """ """
        obj = super(DiffImgFitter, cls).__new__(cls)
        exec(
            "@make_method(DiffImgFitter)\n"
            + "def _minuit_chi2_(self,%s): \n" % (", ".join(obj.FREEPARAMETERS))
            + "    parameters = %s \n" % (", ".join(obj.FREEPARAMETERS))
            + "    return self.get_chi2(parameters)\n"
        )
        return obj

    def __init__(self, diffimg, psfimg, diffvar=None, shape=None):
        """ """
        self.set_data(diffimg, psfimg, diffvar=diffvar, shape=shape)

    # ============== #
    #   FITTER       #
    # ============== #
    def fit(self, verbose=False, **kwargs):
        """ """
        self._hfitverbose = verbose
        # Setting guesses
        self._setup_minuit_(**kwargs)
        # fitting
        self.minuit.migrad()

        # checking fit
        if not self.minuit.valid:
            warnings.warn("migrad is not valid")
            self.fit_ok = False
        else:
            self.fit_ok = True

        # saving results
        self.best_parameters = np.asarray(
            [self.minuit.values[k] for k in self.FREEPARAMETERS]
        )
        # - readout -
        self.fitvalues = {}
        covmat = self.get_fit_covmatrix()

        for i, name in enumerate(self.FREEPARAMETERS):
            self.fitvalues[name] = self.best_parameters[i]
            self.fitvalues[name + ".err"] = np.sqrt(covmat[i, i])

        # -- Additional data -- #
        self.fitvalues["fval"] = self.minuit.fval
        self.fitvalues["chi2"] = self.get_chi2(simple_chi2=True)
        self.fitvalues["chi2dof"] = self.fitvalues["chi2"] / self.dof

        return self.fitvalues

    def _setup_minuit_(self, step=1, print_level=0, **kwargs):
        """ """
        # == Minuit Keys == #
        minuit_kwargs = {}
        for param in self.FREEPARAMETERS:
            minuit_kwargs[param] = kwargs.get(
                "%s_guess" % param, 0 if not "sigma" in param else 1
            )

        self._minuit_input = minuit_kwargs

        self.minuit = Minuit(self._minuit_chi2_, **minuit_kwargs)
        self.minuit.errordef = step
        self.minuit.print_level = print_level

        for param in self.FREEPARAMETERS:
            limits = kwargs.get(
                "%s_boundaries" % param, None if not "sigma" in param else [0, None]
            )
            if limits:
                self.minuit.limits[param] = (limits[0], limits[1])
            fixed = kwargs.get("%s_fixed" % param, False)
            self.minuit.fixed[param] = fixed

    # -------- #
    # SETTER   #
    # -------- #
    def set_data(self, diffimg, psfimg, diffvar=None, shape=None):
        """ """
        if np.shape(diffimg) != np.shape(psfimg):
            raise ValueError(
                "diffimg and psfimg do not have the same shape. They must."
            )

        if diffvar is not None:
            if len(np.atleast_1d(diffvar)) == 1:
                self.diffvar = float(np.atleast_1d(diffvar)[0])
            elif np.shape(diffimg) == np.shape(diffvar):
                self.diffvar = np.asarray(diffvar)
            else:
                raise ValueError(
                    "diffvar is not a single value and diffimg and diffvar do not have the same shape. "
                )
        else:
            self.diffvar = 0

        self._inputshape = np.shape(diffimg)
        self.data = np.ravel(diffimg)
        self.flagout = np.isnan(self.data)
        self.profile = np.ravel(psfimg)
        self._shape = shape

    def set_parameters(self, param):
        """ """
        self.parameters = {k: v for k, v in zip(self.FREEPARAMETERS, param)}

    # ============== #
    #   GETTER       #
    # ============== #
    def get_chi2(self, param=None, simple_chi2=False, use_prior=True):
        """ """
        if param is not None:
            self.set_parameters(param)

        if self._fitverbose:
            print(self.parameters)
        if simple_chi2:
            return np.nansum((self.data - self.scaled_model) ** 2 / self.variance_tot)

        if not use_prior:
            return -2 * np.nansum(self.get_loglikelihood())
        return -2 * np.nansum(self.get_logproba())

    # // Posterior
    def get_logproba(self):
        """ """
        val = self.get_loglikelihood() + self.get_logprior()
        if self._fitverbose:
            print("get_logproba:", val)
        return val

    # // Likelihood
    def get_loglikelihood(self):
        """ """
        val = stats.norm.logpdf(
            self.data, loc=self.scaled_model, scale=np.sqrt(self.variance_tot)
        )
        if self._fitverbose:
            print("get_loglikelihood:", val)
        return val

    # // Priors
    def get_logprior(self):
        """ """
        return 0

    def get_fit_covmatrix(self):
        """coveriance matrix after the fit"""

        def _read_hess_(hess):
            """ """
            if len(hess) == len(self.FREEPARAMETERS):
                return hess

            indexFixed = [
                i
                for i, name in enumerate(self.FREEPARAMETERS)
                if "%s_fixed" % name in self._minuit_input
            ]
            for i in indexFixed:
                newhess = np.insert(hess, i, 0, axis=0)
                newhess = np.insert(newhess, i, 0, axis=1)
                hess = newhess

            return hess

        if self.minuit.valid:
            return _read_hess_(np.asarray(self.minuit.covariance))
        else:
            fakeMatrix = np.zeros(
                (len(self.best_parameters), len(self.best_parameters))
            )
            for i, k in enumerate(self.FREEPARAMETERS):
                fakeMatrix[i, i] = self.minuit.errors[k] ** 2
            warnings.warn("Inaccurate covariance Matrix. Only trace defined")
            return _read_hess_(fakeMatrix)

    # --------- #
    #  PLOTTER  #
    # --------- #
    def show(self, axes=None, res_aspull=True):
        """ """
        if axes is not None:
            axd, axm, axr = axes
            fig = axd.figure
        else:
            import matplotlib.pyplot as mpl

            fig = mpl.figure(figsize=[9, 3])
            axd = fig.add_axes([0.08, 0.1, 0.25, 0.8])
            axm = fig.add_axes([0.38, 0.1, 0.25, 0.8])
            axr = fig.add_axes([0.68, 0.1, 0.25, 0.8])

        data_ = self.data
        axd.imshow(data_.reshape(self._shape), origin="lower")

        model_ = self.scaled_model
        axm.imshow(model_.reshape(self._shape), origin="lower")
        # Residual
        if res_aspull:
            norm = np.sqrt(self.variance_tot)
            prop = dict(vmin=-3, vmax=3)
        else:
            norm = 1
            prop = {}
        axr.imshow(
            ((data_ - model_) / norm).reshape(self._shape), origin="lower", **prop
        )
        if res_aspull:
            axr.text(
                0.05,
                0.95,
                "Pull [-3,3]",
                color="w",
                va="top",
                ha="left",
                transform=axr.transAxes,
                fontdict=dict(weight="bold"),
            )
        return fig

    # ============== #
    #  Properties    #
    # ============== #
    @property
    def variance_tot(self):
        """ """
        return self.diffvar + self.parameters["sigma"] ** 2

    @property
    def scaled_model(self):
        """ """
        return self.profile * self.parameters["ampl"]

    @property
    def ndata(self):
        """ """
        return len(self.data)

    @property
    def nfreeparameters(self):
        """ """
        fixed = (
            0
            if not hasattr(self, "_minuit_input")
            else len([k for k, v in self._minuit_input.items() if "fixed" in k and v])
        )
        return len(self.FREEPARAMETERS) - fixed

    @property
    def dof(self):
        """ """
        return self.ndata - self.nfreeparameters

    @property
    def _fitverbose(self):
        """ """
        if not hasattr(self, "_hfitverbose"):
            self._hfitverbose = False
        return self._hfitverbose
