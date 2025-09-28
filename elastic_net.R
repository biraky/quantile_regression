# --- Packages
library(glmnet)
library(dplyr)

# --- Design matrix (drop explicit intercept; glmnet adds it by default)
x <- model.matrix(dependent_variable ~ . - 1, data = dataset)
y <- dataset$dependent_variable

set.seed(2025)

# ---- Elastic net with K-fold CV at a chosen alpha
cv_en <- cv.glmnet(
  x, y,
  family      = "gaussian",
  alpha       = 0.5,      # 0<alpha<1 = elastic net
  nfolds      = 10,
  standardize = TRUE
)

# Key lambdas
cv_en$lambda.min   # lambda with minimal CV error
cv_en$lambda.1se   # largest lambda within 1 SE of the minimum

# Coefficients at common choices
coef_min <- coef(cv_en, s = "lambda.min")
coef_1se <- coef(cv_en, s = "lambda.1se")

