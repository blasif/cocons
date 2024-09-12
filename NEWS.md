cocons 0.1.2
----------------------------------------------------------------

* Improve overall help files
* Add examples for coco, cocoOptim, cocoPredict, and cocoSim
* update Vignette
* some bug fixes and overall code polishing
* Some important changes in R files:
* \code{smooth_limits} from \code{info} from the coco function is now called \code{smooth.limits} to match the style of other arguments (no backward compatibility)
* new names for C++ sparse covariance functions
* add warnings of non-convergence of the LBFGSB for \code{cocoOptim}
* bug fixes for coco "methods"
* less redundant code