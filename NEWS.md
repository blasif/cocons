## cocons 0.1.5

- Vignette updated
- changed default max.iters from 200 to 500.
- now `coco` objects store information about optim.control arguments for full reproducibility.
- introduce a two-step optimization approach for automatic model selection. See `help(coco)`.

## cocons 0.1.4

- Documentation polished

### Enhances

- `GetNeg2loglikelihood`, `GetNeg2loglikelihoodTaper`, and `GetNeg2loglikelihoodTaperProfile` now much faster (~35\% faster)
- `cocoOptim`:
  - new `REML` estimation for dense `coco` type.
- now cpp functions `cov_rns` and `cov_rns_taper` include as special cases the well-known shapes when nu = 0.5, 1.5, and 2.5, yielding computational speed-ups.
- now `getCIs` display names of covariate instead of index of Design Matrix.
- now `.cocons.check.info` also checks for smoothness model and smooth.limits.
  
### Changes

- `coco` now accepts `data.frame` for locs argument, which is then converted to matrix.
- `getTrend` now is called `getSpatMean`.
- `cocoPredict` renamed output spatial mean vectors: `trend` is now called `systematic`, while `mean` is called `stochastic`.
- `cocoOptim` : 
  - reordering of arguments.
  - `mle` estimation method now called `ml`.
  - `pmle` estimation method now called `pml`.
- `getHessian` :
  - `mle` to `ml` and `pmle` to `pml`.
- `getCIs` now `alpha` argument reflects confidence level instead of 1-confidence level.

### Fixes

- fix a bug for plotOptimInfo when handling `pml` or `reml` objects.

## cocons 0.1.3

- Vignette updates, watermark removed
- Documentation polished
- New `holes_bm` dataset with independent realizations and spatial trend

### Enhances

- added a NEWS.md file with version updates / modifications / enhances / etc
- automatized `delta` for method `plot` for coco class.
- Better visualization for `plot(cocoOptim object, type = "ellipse")`
- `coco`:
  - now it is not necessary to provide all models for each source of nonstationarity. Those not specified will be set to those referenced
  with a stationary 0-mean model (i.e. tilt = 0 , aniso = 0, nugget = -Inf (because of log-parameterization)). If 'smooth' is not specified, then it is set to 0.5.
- `cocoOptim`:
  - "auto" option for `ncores` argument for `cocoOptim`, providing a convenient number of threads based on the number of parameters to estimate, available threads, and settings of the LBFGSB routine
  - "safe" argument, which prevents crashes due to ill-posed covariance matrices (Choelsky factorization error)
  - `.cocons.check.convergence` now checks and reports at which iteration ill-posed covariance matrices have been found during the optimization. 
  - now "pmle" works with multiple independent realizations for coco types `dense` and `sparse`
  - safer parallel handling
- `getHessian` more memory efficient
- `cocoSim`: 
  - if provided a fitted coco object, then `pars` argument can be `NULL` (default), and `coco.object@output$par` is used instead (and also `type` is set to `diff`).
  - more memory efficient
- safer parallel handling for `getHessian`
- polishing of neg2loglikelihood functions, leading to more efficient code
- new and more polished internal functions to assess the validity of arguments (`stopifnot()` instead of `if() stop()`)
- `getCondNumber` optimized
- small improvement over cpp functions

### Changes

- switched `getPen` as an internal function
- method `plot` for `coco` objects shows rotation angle of the kernel w.r.t x-axis
- `GetSpateffects` now provides angle w.r.t to x-axis
- method "summary" for `coco` objecets (former "print" method)
- more proper naming of objects inside functions
- renaming of "cat.vars" to "skip.scale" + associated checks and optimization
- `getCondNumber` removed, which can be replaced with `kappa` function from base R (i.e. `kappa(getCovMatrix(coco.object),exact = TRUE)`)

## cocons 0.1.2

-   improve overall help files
-   add examples for `coco`, `cocoOptim`, `cocoPredict`, and `cocoSim`
-   update Vignette
-   some bug fixes and overall code polishing
-   `smooth_limits` from `info` from the coco function is now called `smooth.limits` to match the style of other arguments (no backward compatibility)
-   new names for C++ sparse covariance functions
-   add warnings of non-convergence of the LBFGSB for `cocoOptim`
-   bug fixes for coco "methods"
-   less redundant code
