#' Holes Data Set
#'
#' The synthetic "holes" provides a set of training and test data.frame of a Gaussian process realization with a (inherently dense) nonstationary covariance function. 
#' Four holes are present in the training dataset, and the task is to predict them. 
#'
#' @format A list with training and test data.frame with rows and variables:
#' \describe{
#'   \item{x}{first spatial coordinate}
#'   \item{y}{second spatial coordinate}
#'   \item{cox_x}{first spatial characteristic}
#'   \item{cov_y}{second spatial characteristic}
#'   \item{z}{response variable}
#' }
#' @source Source of the data
#' @examples
#' data(holes)
"holes"

#' Holes with trend + multiple realizations Data Set
#'
#' The synthetic "holes_bm" provides a set of training and test data.frame of a Gaussian process realization with a (inherently dense) nonstationary covariance function. 
#' Four holes are present in the training dataset, and the task is to predict them. This version provides ten independent realizaitons of the process, as well as considers
#' a spatial mean effect.
#'
#' @format A list with training and test data.frame with rows and variables:
#' \describe{
#'   \item{x}{first spatial coordinate}
#'   \item{y}{second spatial coordinate}
#'   \item{cox_x}{first spatial characteristic}
#'   \item{cov_y}{second spatial characteristic}
#'   \item{cov_z}{second spatial characteristic}
#'   \item{z}{response variable}
#' }
#' @source Source of the data
#' @examples
#' data(holes_bm)
"holes_bm"

#' Stripes Data Set
#'
#' The synthetic "stripes" provides a set of training and test data.frame of a Gaussian process realization with a (inherently sparse) nonstationary covariance function. 
#' Several stripes are present in the training dataset, and the task is to predict them. 
#'
#' @format A list with training and test data.frame with rows and variables:
#' \describe{
#'   \item{x}{first spatial coordinate}
#'   \item{y}{second spatial coordinate}
#'   \item{cox_x}{first spatial characteristic}
#'   \item{cov_y}{second spatial characteristic}
#'   \item{cov_xy}{third spatial characteristic}
#'   \item{z}{response variable}
#' }
#' @source Source of the data
#' @examples
#' data(stripes)
"stripes"