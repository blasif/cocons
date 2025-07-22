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
