# tests/coco_test.R

# Load your package
library(cocons)

data_test <- cocons::holes$training[1:50,]

# Create a minimal reproducible example
locs <- matrix(c(0, 1, 0, 1), ncol = 2)
data <- data.frame(x = locs[, 1])
z <- rnorm(2)

# Define a trivial model.list
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

# Basic sanity check: covariance matrix has correct dimensions
cmat <- getCovMatrix(test_optim)

stopifnot(all(eigen(cmat)$values > 0))
stopifnot(is.matrix(cmat), nrow(cmat) == 50, ncol(cmat) == 50)