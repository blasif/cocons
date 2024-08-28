\name{cocons-package}
\alias{cocons-package}
\alias{cocons}
\docType{package}
\title{Covariate-based Covariance Functions for Nonstationary Gaussian Processes}
\description{
Provides routines and methods for estimation and prediction of nonstationary Gaussian process models with modular covariate-based covariance functions. Several sources of nonstationarity can be modeled based on spatial information, including a trend, marginal standard deviation, local geometric anisotropy, local nugget, and spatially varying smoothness. Each of these components is modeled independently. For large datasets, a sparse-induced version of the nonstationary covariance function is provided to speed up computations. Models are estimated via maximum likelihood. A variety of functions are also included to compute prediction metrics and to visualize, simulate, and summarize these types of models.
}
\section{Disclaimer}{
This package is provided "as is" without warranty of any kind, either express or implied.
Backwards compatibility will not be offered until version 0.1.5.
}
\author{
Federico Blasi [aut, cre], \email{federico.blasi@uzh.ch}
}
\examples{
\dontrun{
    vignette("cocons", package = "cocons")
    methods(class = "coco")
  }
}