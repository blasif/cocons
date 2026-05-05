# tests/coco_test.R

## --- Load package ----------------------------------------------------------

suppressPackageStartupMessages({
  library(cocons)
})

## Optional: check version at least 0.1.5
if (packageVersion("cocons") < "0.1.5") {
  stop("cocons >= 0.1.5 required for these checks.", call. = FALSE)
}

## --------

data_test <- cocons::holes$training[1:50,]

model.list <- list(
  'mean'    = 0,
  'std.dev' = ~ 1,
  'scale'   = ~ 1,
  'aniso'   = 0,
  'tilt'    = 0,
  'smooth'  = 1.5,
  'nugget'  = -Inf
)

# Run function
obj <- coco(type = "dense", 
            data = data_test[,1:4], 
            locs = data_test[,1:2], 
            z = data_test$z,
            model.list = model.list)

# Start manual checking
stopifnot(inherits(obj, "coco"))

test_optim <- cocoOptim(obj, ncores = 1)

cmat <- getCovMatrix(test_optim)
stopifnot(is.matrix(cmat), nrow(cmat) == 50, ncol(cmat) == 50)
stopifnot(all(eigen(cmat)$values > 0))

predss <- cocoPredict(test_optim, newdataset = cocons::holes$test[1:50,1:4],newlocs = as.matrix(holes$test[1:50,1:2]))
stopifnot(all(predss$systematic == 0))
stopifnot(all(!is.na(predss$stochastic)))

model.list <- list(
  'mean'    = formula(~ 1 + cov_x + cov_y),
  'std.dev' = ~ 1,
  'scale'   = ~ 1,
  'aniso'   = 0,
  'tilt'    = 0,
  'smooth'  = 1.5,
  'nugget'  = -Inf
)

obj <- coco(type = "dense", 
            data = data_test[,1:4], 
            locs = data_test[,1:2], 
            z = data_test$z,
            model.list = model.list)

test_optim <- cocoOptim(obj, ncores = 1, optim.type = 'pml')

cmat <- getCovMatrix(test_optim)
stopifnot(is.matrix(cmat), nrow(cmat) == 50, ncol(cmat) == 50)
stopifnot(all(eigen(cmat)$values > 0))

predss <- cocoPredict(test_optim, newdataset = cocons::holes$test[1:50,1:4], newlocs = as.matrix(holes$test[1:50,1:2]))
stopifnot(all(!is.na(predss$systematic)))
stopifnot(all(!is.na(predss$stochastic)))

# Penalization

model.list <- list(
  'mean'    = formula(~ 1 + cov_x + cov_y),
  'std.dev' = formula(~ 1 + cov_x + cov_y),
  'scale'   = formula(~ 1 + cov_x + cov_y),
  'aniso'   = formula(~ 1 + cov_x + cov_y),
  'tilt'    = formula(~ 1 + cov_x + cov_y),
  'smooth'  = formula(~ 1 + cov_x + cov_y),
  'nugget'  = -Inf
)

set.seed(181222)
sample_index <- sample(5570, 100)

data_test <- cocons::holes$training[sample_index,]

# Run function
obj <- coco(type = "dense", 
            data = data_test[,1:4], 
            locs = data_test[,1:2], 
            z = data_test$z,
            model.list = model.list,
            info = list('lambda.Sigma' = 2,
                        'lambda.betas' = 2,
                        'lambda.reg' = 0.3,
                        'smooth.limits' = c(0.5,2.5)))

test_optim <- cocoOptim(obj, ncores = 1, optim.type = 'ml')

obj <- coco(type = "dense", 
            data = data_test[,1:4], 
            locs = data_test[,1:2], 
            z = data_test$z,
            model.list = model.list,
            info = list('lambda.Sigma' = 0.25,
                        'lambda.betas' = 0.25,
                        'lambda.reg' = 0.3,
                        'smooth.limits' = c(0.5,2.5)))

test_optim <- cocoOptim(obj, ncores = 1, 
                        boundaries = getBoundariesV4(coco.object = obj,lower.bound = -0.5,upper.bound = 0.5), optim.type = 'ml')

getEstims(test_optim)

plot(test_optim)

test_Hess <- getHessian(test_optim)

getCIs(test_optim,solve(test_Hess))

cocoSim(test_optim)

# Spam and sparse coco type

# Run function

data_test <- cocons::holes$training[sample_index,]

model.list <- list(
  'mean'    = 0,
  'std.dev' = ~ 1,
  'scale'   = ~ 1,
  'aniso'   = 0,
  'tilt'    = 0,
  'smooth'  = 1.5,
  'nugget'  = -Inf
)

obj <- coco(type = "sparse", 
            data = data_test[,1:4], 
            locs = data_test[,1:2], 
            z = data_test$z,
            model.list = model.list,
            info = list('taper'= spam::cov.wend1,
                        'delta'= 1))

# Start manual checking
stopifnot(inherits(obj, "coco"))

test_optim <- cocoOptim(obj, ncores = 1)

test_optim <- cocoOptim(coco.object = obj, ncores = 1,optim.type = 'ml')

model.list <- list(
  'mean'    = formula(~ 1 + cov_x + cov_y),
  'std.dev' = formula(~ 1 + cov_x + cov_y),
  'scale'   = ~ 1,
  'aniso'   = 0,
  'tilt'    = 0,
  'smooth'  = 1.5,
  'nugget'  = -Inf
)

obj <- coco(type = "sparse", 
            data = data_test[,1:4], 
            locs = data_test[,1:2], 
            z = data_test$z,
            model.list = model.list,
            info = list('taper'= spam::cov.wend1,
                        'delta'= 1.5,
                        'lambda.Sigma' = 0.1,
                        'lambda.betas' = 0.1,
                        'lambda.reg' = 0.1))

test_optim <- cocoOptim(obj,boundaries = getBoundariesV4(obj, lower.bound = -0.25, upper.bound = 0.25) , ncores = 1)

getEstims(test_optim)

plot(test_optim)

getHessian(test_optim)

cmat <- as.matrix(getCovMatrix(test_optim))
stopifnot(is.matrix(cmat), nrow(cmat) == 100, ncol(cmat) == 100)
stopifnot(all(eigen(cmat)$values > 0))

predss <- cocoPredict(test_optim, newdataset = cocons::holes$test[1:200,1:4],newlocs = as.matrix(holes$test[1:200,1:2]))
#stopifnot(all(predss$systematic == 0))
stopifnot(all(!is.na(predss$stochastic)))

test_optim <- cocoOptim(obj, ncores = 1, optim.type = "pml")

cmat <- as.matrix(getCovMatrix(test_optim))
stopifnot(is.matrix(cmat), nrow(cmat) == 100, ncol(cmat) == 100)
stopifnot(all(eigen(cmat)$values > 0))

predss <- cocoPredict(test_optim, newdataset = cocons::holes$test[1:200,1:4],newlocs = as.matrix(holes$test[1:200,1:2]))
stopifnot(all(!is.na(predss$stochastic)))

aeq <- function(a, b, tol = 1e-8) {
  if (!isTRUE(all.equal(a, b, tolerance = tol))) {
    stop(sprintf("Numeric mismatch:\nexpected: %s\ngot:      %s",
                 paste0(head(b,10), collapse=", "),
                 paste0(head(a,10), collapse=", ")),
         call. = FALSE)
  }
  invisible(TRUE)
}

expect_error <- function(expr, pattern = NULL) {
  err <- NULL
  tryCatch(eval.parent(substitute(expr)), error = function(e) err <<- e)
  if (is.null(err)) stop("Expected an error, but none was thrown.", call. = FALSE)
  if (!is.null(pattern) && !grepl(pattern, conditionMessage(err))) {
    stop(sprintf("Error did not match /%s/; got: %s",
                 pattern, conditionMessage(err)), call. = FALSE)
  }
  invisible(TRUE)
}

set.seed(1)
locs <- expand.grid(seq(0, 1, length.out = 5),
                    seq(0, 1, length.out = 5))
toydata <- data.frame(x = locs[, 1])
z <- rnorm(nrow(locs))

model.list <- list(
  mean    = 0,
  std.dev = formula(~ 1),
  scale   = formula(~ 1 + x),
  aniso   = 0,
  tilt    = 0,
  smooth  = 3/2,
  nugget  = -Inf
)

coco_dense <- coco(
  type = "dense",
  data = toydata,
  locs = as.matrix(locs),
  z = z,
  model.list = model.list
)

## S4 class + slots
stopifnot(isS4(coco_dense))
stopifnot(methods::is(coco_dense, "coco"))
needed_slots <- c("type","data","locs","z","model.list","info","output")
stopifnot(all(needed_slots %in% slotNames(coco_dense)))
stopifnot(identical(coco_dense@type, "dense"))
stopifnot(nrow(coco_dense@locs) == nrow(toydata))
stopifnot(nrow(coco_dense@locs) == length(z))

## --- 2) invalid 'type' should error ---------------------------------------

expect_error(
  coco(type = "weird",
       data = toydata,
       locs = as.matrix(locs),
       z = z,
       model.list = model.list),
  pattern = "type|dense|sparse"
)

## --- 3) is.formula() helper ------------------------------------------------

stopifnot(isTRUE( cocons::is.formula(~ 1) ))
stopifnot(isTRUE(!cocons::is.formula(1)))
stopifnot(isTRUE(!cocons::is.formula("~ 1 + x")))

## --- 4) Scoring rules: getLogScore / getCRPS -------------------------------

z.pred    <- c(0.0,  1.0, -0.5)
mean.pred <- c(0.1,  0.8, -0.4)
sd.pred   <- c(0.5,  1.2,  2.0)

## reference formulas (match package code)
logscore_ref <- (log(2*pi) + ((z.pred - mean.pred)/sd.pred)^2)/2 + log(sd.pred)
v <- (mean.pred - z.pred)/sd.pred
crps_ref <- sd.pred * ( v * (2 * stats::pnorm(v) - 1) +
                          2 * stats::dnorm(v) - 1 / sqrt(pi) )

aeq(cocons::getLogScore(z.pred, mean.pred, sd.pred), logscore_ref, tol = 1e-10)
aeq(cocons::getCRPS(z.pred, mean.pred, sd.pred),    crps_ref,    tol = 1e-8)

## --- 5) getDensityFromDelta() behavior ------------------------------------

## Dense objects should error
expect_error(cocons::getDensityFromDelta(coco_dense, delta = 0.2),
             pattern = "only for sparse coco objects")

## Optional: sparse path (skipped if 'spam' not available)
if (requireNamespace("spam", quietly = TRUE)) {
  coco_sparse <- coco(
    type = "sparse",
    data = toydata,
    locs = as.matrix(locs),
    z = z,
    model.list = model.list,
    info = list(
      taper = spam::cov.wend1,
      delta = 0.15,
      smooth.limits = c(0.5, 2.5)
    )
  )
  
  dens <- cocons::getDensityFromDelta(coco_sparse, delta = 0.20)
  stopifnot(is.numeric(dens), length(dens) == 1L, is.finite(dens))
  ## For a 25x25 grid with short taper, density should be in (0, 1)
  stopifnot(dens > 0, dens < 1)
} else {
  message("Skipping sparse/taper density check (package 'spam' not available).")
}