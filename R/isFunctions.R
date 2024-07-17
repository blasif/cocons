
# is.formula() -----------------------------------------------------------------

#' check whether an object belongs to a formula class
#' @description check whether an object belongs to a formula class
#' @usage is.formula(x)
#' @param x an R object
#' @returns TRUE/FALSE
#' @author Federico Blasi
is.formula <- function(x) {
  ifelse(any("formula" %in% class(x)), yes = TRUE, no = FALSE)
}
