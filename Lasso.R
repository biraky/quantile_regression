library(glmnet)
library(dplyr)

x <- model.matrix(dependent_variable ~ ., data = dataset)
y <- dataset$dependent_variable

set.seed(2025)

# ---- LASSO with K-fold cross-validation ----
cv_lasso <- cv.glmnet(
  x, y,
  family = "gaussian",  # outcome is continuous
  alpha = 1,            # 1 = LASSO (0 = ridge; in-between = elastic net)
  nfolds = 10,          # K in K-fold CV; with n=25, 5 or 10 is typical
  standardize = TRUE    # standardise predictors internally
)

# Inspect the two key lambdas:
cv_lasso$lambda.min   # lambda giving minimal CV error
cv_lasso$lambda.1se   # largest lambda within 1 SE of the minimum

# ---- Coefficients at the two lambdas ----
coef_min  <- coef(cv_lasso, s = "lambda.min")
coef_1se  <- coef(cv_lasso, s = "lambda.1se")

# Non-zero variables selected
vars_min <- setdiff(rownames(coef_min)[as.vector(coef_min[,1] != 0)], "(Intercept)")
vars_1se <- setdiff(rownames(coef_1se)[as.vector(coef_1se[,1] != 0)], "(Intercept)")

vars_min
vars_1se

# ---- Quick in-sample performance (just for a feel) ----
pred_min <- as.numeric(predict(cv_lasso, newx = x, s = "lambda.min"))
pred_1se <- as.numeric(predict(cv_lasso, newx = x, s = "lambda.1se"))

rmse <- function(a, b) sqrt(mean((a - b)^2))
rmse_min <- rmse(y, pred_min)
rmse_1se <- rmse(y, pred_1se)

list(
  lambda_min = cv_lasso$lambda.min,
  lambda_1se = cv_lasso$lambda.1se,
  cv_mse_min = min(cv_lasso$cvm),  # the minimal average CV error (MSE for gaussian)
  rmse_in_sample_min = rmse_min,
  rmse_in_sample_1se = rmse_1se,
  selected_vars_lambda_min = vars_min,
  selected_vars_lambda_1se = vars_1se
)

# ---- Optional: plot CV curve and coefficient paths ----
# CV error vs log(lambda):
 plot(cv_lasso)

# Coefficient paths (as lambda decreases, more vars enter):
 fit_lasso <- glmnet(x, y, alpha = 1, standardize = TRUE, family = "gaussian")
 plot(fit_lasso, xvar = "lambda")
