
# is.formula() -----------------------------------------------------------------

#' check whether an R object is a formula
#' @description check whether an R object is a formula
#' @usage is.formula(x)
#' @param x (ANY) an R object.
#' @returns TRUE/FALSE
#' @author Federico Blasi
is.formula <- function(x) {
  methods::is(x,"formula")
}
