#' Creates a coco S4 object
#' @description Creates an S4 object of class \code{coco}, which is the centerpiece of the \pkg{cocons} package. The function provides a set of consistency checks for ensuring the suitability of the different objects involved.
#' @details 
#' Two types of \code{coco} objects are available, each assuming a different type of covariance matrix for the Gaussian process. 
#' Type \code{"dense"} builds dense covariance matrices (non zero elements), while type \code{"sparse"} builds sparse covariance 
#' matrices by tapering the dense covariance matrix with a compact isotropic compact-supported correlation matrix \[1\]. 
#' Type \code{"sparse"} allows a set of efficient algorithms, thus making it more suitable for large sample sizes.
#' 
#' An important component of the \code{coco} S4 class is the \code{model.list} specification, involving individual formulas provided as a list, where each of them specifies a covariate-based parametric model for a specific source of nonstationarity. 
#' It involves \code{"mean"} for the spatial mean, the \code{"std.dev"} for the marginal standard deviation, 
#' \code{"scale"}, \code{"aniso"} and \code{"tilt"}, each of them shaping specific aspects of the local spatial geometrically anisotropy structure,
#' \code{"smooth"} handling local smoothness, and \code{"nugget"} handling the local nugget effect. The models are defined as:
#' 
#' 
#' | Source       | Related to    | Description                             | Model                                                                    |
#' |--------------|---------------|-----------------------------------------|--------------------------------------------------------------------------|
#' | \emph{mean}      | \eqn{\mu}       | Spatial mean function             | \eqn{\boldsymbol{X}_1\boldsymbol{\beta}}                                               |
#' | \emph{std.dev}   | \eqn{\sigma^{X}}| Marginal standard deviation       | \eqn{\text{exp}(0.5 \boldsymbol{X}_2 \boldsymbol{\alpha})}                                    |
#' | \emph{scale}     | \eqn{\boldsymbol{\Sigma}^{X}}| Local scale          | \eqn{\text{exp}(\boldsymbol{X}_3 \boldsymbol{\theta}_1)}                                      |
#' | \emph{aniso}     | \eqn{\boldsymbol{\Sigma}^{X}}| Local geometric anisotropy        | \eqn{\text{exp}(\boldsymbol{X}_4 \boldsymbol{\theta}_2)}                         |
#' | \emph{tilt}      | \eqn{\boldsymbol{\Sigma}^{X}}| (Restricted) local tilt           | \eqn{\cos(\text{logit}^{-1}(\boldsymbol{X}_5 \boldsymbol{\theta}_3))}          |
#' | \emph{smooth}    | \eqn{\nu^{X}}    | Local smoothness                  | \eqn{(\nu_{u} - \nu_{l})/(1+\text{exp}(-\boldsymbol{X}_6 \boldsymbol{\phi})) + \nu_{l}}      |
#' | \emph{nugget}    | \eqn{\sigma^{X}_{\epsilon}}| Local micro-scale variability   | \eqn{\text{exp}(\boldsymbol{X}_7 \boldsymbol{\zeta})}                                |
#' 
#' where \eqn{\boldsymbol{\beta}}, \eqn{\boldsymbol{\alpha}}, \eqn{\boldsymbol{\theta}_1}, \eqn{\boldsymbol{\theta}_2}, \eqn{\boldsymbol{\theta}_3}, \eqn{\boldsymbol{\phi}}, and \eqn{\boldsymbol{\zeta}} are the parameter vectors of each model,
#' \eqn{\nu_{l}}, and \eqn{\nu_{u}} are the lower and upper bounds limiting the range of variation of the spatially-varying smoothness, and where \eqn{\boldsymbol{X}_{\ell}} relates to the design matrix defined by the specific models for each of the source of nonstationarity.
#' 
#' Lastly, arguments for the \code{"info"} list argument involve: \itemize{
#' \item \code{"lambda.reg"}: (\code{numeric}) a positive scalar specifying the regularization parameter. Larger values regularize highly-smoothed long-tailed covariance functions. 
#' \item \code{"lambda.Sigma"}: (\code{numeric}) a positive scalar specifying the penalization parameter for the covariate-driven covariance parameters. 
#' \item \code{"lambda.betas"}: (\code{numeric}) a positive scalar specifying the penalization parameter for the covariate-driven spatial mean parameters. 
#' \item \code{"sparse.point"}: (\code{numeric}) a cutting point for which smaller coefficients in absolute value will be set to zero after the smoothed L1 penalization optimization. Used in combination with \code{lambda.Sigma} and \code{lambda.betas}. By default, it is set to \code{1e-4}.
#' \item \code{"smooth.limits"}: (\code{numeric vector}) specifying the range of variation for the spatially varying smoothness (e.g. c(0.5, 2.5)).
#' \item \code{"taper"}: (\code{numeric}) specifying the desired taper function from the spam package (only for "sparse" coco objects).
#' \item \code{"delta"}: (\code{numeric}) specifying the taper range/scale (only for "sparse" coco objects).
#' \item \code{"skip.scale"}: (\code{integer vector}) By default, all covariates are scaled. \code{skip.scale} allows to specify the index of those variables in \code{data} that should not be scaled during the optimization.
#' }
#' 
#' @usage coco(type, data, locs, z, model.list, info, output = list())
#' @param type (\code{character}) One of two available types \code{"dense"} or \code{"sparse"}. See description.
#' @param data (\code{data.frame}) A \code{data.frame} with covariates information, where \code{colnames(data)} matches model.list specification.
#' @param locs (\code{matrix}) A \code{matrix} with spatial locations. 
#' @param z (\code{vector} or \code{matrix}) A matrix of \eqn{n \times r} response realizations, one realization per column. When considering only one realization, a vector can also be provided.
#' @param model.list (\code{list}) A \code{list} specifying a model for each source of nonstationarity.
#' @param info (\code{list} or \code{NULL}) A \code{list} specifying characteristics of the coco object.
#' @param output (\code{list} or \code{NULL}) Empty or the resulting object from running \link[optimParallel]{optimParallel}, adding to this a list with boundary information (see \link{getBoundaries} to check the expected structure).
#' 
#' @returns (\code{S4}) An S4 object of class \code{coco}.
#' 
#' @examples
#' \dontrun{
#' locs <- expand.grid(seq(0,1,length.out = 10),
#' seq(0,1,length.out = 10))
#' 
#' toydata <- data.frame('x' = locs[,1])
#' 
#' set.seed(1)
#' z <- rnorm(100)
#' 
#' model.list <- list('mean' = 0,
#'                    'std.dev' = formula( ~ 1),
#'                    'scale' = formula( ~ 1 + x),
#'                    'aniso' = 0,
#'                    'tilt' = 0,
#'                    'smooth' = 3/2,
#'                    'nugget' = -Inf)
#'                    
#' coco_object <- coco(type = 'dense',
#'                     data = toydata,
#'                     locs = as.matrix(locs),
#'                     z = z,
#'                     model.list = model.list)
#' 
#' coco_object
#' }
#' 
#' @author Federico Blasi
#' @seealso [spam::cov.wend1()]
#' @references 
#' \[1\] Furrer, Reinhard, Marc G. Genton, and Douglas Nychka. 
#' \emph{"Covariance tapering for interpolation of large spatial datasets."} 
#' Journal of Computational and Graphical Statistics 15.3 (2006): 502-523.
#' 
coco <- function(type,
                  data,
                  locs,
                  z,
                  model.list,
                  info = list(),
                  output = list()) {
  
  .cocons.check.type(type)
  .cocons.check.data(data)
  
  if(is.data.frame(locs)){
    locs <- as.matrix(locs)
  }
  
  .cocons.check.locs(locs)
  
  if(is.vector(z)){
    z <- matrix(z, ncol = 1)
  }
  
  .cocons.check.z(z, data)
  
  .cocons.check.model.list(model.list,
                          data)
  
  #
  
  if(is.null(model.list$mean)){
    model.list$mean <- 0
  }
  
  if(is.null(model.list$aniso)){
    model.list$aniso <- 0
  }
  
  if(is.null(model.list$tilt)){
    model.list$tilt <- 0
  }
  
  if(is.null(model.list$smooth)){
    model.list$smooth <- 0.5
  }
  
  if(is.null(model.list$nugget)){
    model.list$nugget <- -Inf
  }
  
  model.list <- model.list[getOption("cocons.Dictionary")]
  
  .cocons.check.info(type = type, 
                    info = info,
                    model.list = model.list,
                    data = data)
  
  .cocons.check.output(output)
  
  if(is.null(info$lambda.reg)){
    info$lambda.reg <- 0
  }
  
  if(is.null(info$lambda.betas)){
    info$lambda.betas <- 0
  }
  
  if(is.null(info$lambda.Sigma)){
    info$lambda.Sigma <- 0
  }
  
  if(any(c(info$lambda.betas,info$lambda.Sigma) > 0) && (is.null(info$sparse.point))){
    info$sparse.point <- 1e-4
  }
  
  if(T){
    if(!is.formula(model.list$smooth)){
      info$smooth.limits <- c(model.list$smooth[1],
                              model.list$smooth[1])
    }
  }
  
  methods::new(
    Class = "coco",
    type = type,
    data = data,
    locs = locs,
    z = z,
    model.list = model.list,
    info = info,
    output = output
  )
  
}
