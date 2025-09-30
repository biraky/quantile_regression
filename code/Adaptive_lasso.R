library(glmnet)
library(dplyr)

x <- model.matrix(dependent_variable ~ ., data = dataset)
y <- dataset$dependent_variable

set.seed(2025)

cv_ridge <- cv.glmnet(
  x, y,
  family = "gaussian",
  alpha = 0,        # ridge
  nfolds = 10,
  standardize = TRUE
)

# Coefficients from the initial estimator at lambda.min
# coef() returns a sparse vector including the intercept in the first entry.
beta_init <- as.numeric(coef(cv_ridge, s = "lambda.min"))
beta_init_no_int <- beta_init[-1]  # drop intercept

# Step 2b: build adaptive weights w_j = 1 / (|beta_init_j| + eps)^gamma
gamma <- 1.0
eps   <- 1e-6
w <- 1 / (abs(beta_init_no_int) + eps)^gamma

# Safety check: length must match number of columns in x
stopifnot(length(w) == ncol(x))

# Step 2c: weighted LASSO via penalty.factor
cv_adalasso <- cv.glmnet(
  x, y,
  family = "gaussian",
  alpha = 1,                 # LASSO
  nfolds = 10,
  standardize = TRUE,
  penalty.factor = w         # <- adaptive part
)

# Coefficients at the two usual lambdas
coef_ada_min <- coef(cv_adalasso, s = "lambda.min")
coef_ada_1se <- coef(cv_adalasso, s = "lambda.1se")
