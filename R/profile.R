.onLoad <- function(libname, pkgname) {
  
  default.options <- list(
  
    cocons.Dictionary = c("mean", "std.dev", "scale", 
                           "aniso", "tilt",
                           "smooth", "nugget"),
    
    cocons.Optim.Control = list("parallel" = list(forward = FALSE,
                                                   loginfo = TRUE),
                                 "control" = list(factr = 1e-8/.Machine$double.eps, # tolerance of 1e-8
                                                  trace = 1, 
                                                  maxit = 200,
                                                  lmm = 100), # we save more information than the standard 5 for lmm
                                 "hessian" = F)
    
  )
  options(default.options)
}