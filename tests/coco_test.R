# tests/coco_test.R

# Load your package
library(cocons)

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

# Spam and sparse coco type

# Run function
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

cmat <- getCovMatrix(test_optim)
stopifnot(is.matrix(cmat), nrow(cmat) == 50, ncol(cmat) == 50)
stopifnot(all(eigen(cmat)$values > 0))

predss <- cocoPredict(test_optim, newdataset = cocons::holes$test[1:50,1:4],newlocs = as.matrix(holes$test[1:50,1:2]))
stopifnot(all(predss$systematic == 0))
stopifnot(all(!is.na(predss$stochastic)))

test_optim <- cocoOptim(obj, ncores = 1,optim.type = "pml")

cmat <- getCovMatrix(test_optim)
stopifnot(is.matrix(cmat), nrow(cmat) == 50, ncol(cmat) == 50)
stopifnot(all(eigen(cmat)$values > 0))

predss <- cocoPredict(test_optim, newdataset = cocons::holes$test[1:50,1:4],newlocs = as.matrix(holes$test[1:50,1:2]))
stopifnot(all(predss$systematic == 0))
stopifnot(all(!is.na(predss$stochastic)))