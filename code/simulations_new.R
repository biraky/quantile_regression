library(Matrix)
library(MASS)
library(gWQS)

library(corrplot)
library(ggplot2)
library(dplyr)

generate_correlation_matrix <- \(n_correlated = 10,   correlation_min = 0.6,
                        correlation_max = 0.9, seed = 12345) {
  set.seed(seed)
  # Build correlation matrix for correlated mixtures
  corr_vals <- runif(n_correlated * (n_correlated - 1) / 2,
                     correlation_min, correlation_max)
  sigma_corr <- diag(1, n_correlated)
  sigma_corr[upper.tri(sigma_corr)] <- corr_vals
  sigma_corr[lower.tri(sigma_corr)] <- t(sigma_corr)[lower.tri(sigma_corr)]

  # Ensure the matrix is positive definite
  sigma_corr <- as.matrix(nearPD(sigma_corr)$mat)
  sigma_corr
}


vcor <- outer(c(seq(0.6,0.8, length.out = 5), seq(0.92,0.98, length.out = 5)),
              c(seq(0.6,0.8, length.out = 5), seq(0.92,0.98, length.out = 5)))
diag(vcor) <- rep(1, 10)

corrplot(vcor, method = "color",type = "full",
         col = rev(rep(viridis(200), 3)),
         col.lim = c(0.3,1),
         tl.col = "black", addCoef.col = "white",
         number.digits = 2)

generate_simulated_dataset <- \(
  n_obs = 5,
  n_total = 25,
  n_correlated = 10,
  correlation_matrix,
  weights_feature = rep(1, n_correlated),
  betas_features =  rep(1, n_correlated),
  betas_biomarkers = rep(1, n_total - n_correlated)
)
{
  n_independent <- n_total - n_correlated

  # Generate the mixture components
  correlated_mix <- mvrnorm(n_obs, mu = rep(0, n_correlated), Sigma = sigma)
  independent_mix <- matrix(rnorm(n_obs * n_independent, 0, 1), ncol = n_independent)

  mixtures <- cbind(correlated_mix, independent_mix)
  colnames(mixtures) <- c(
    paste0("feature_",   seq_len(n_correlated)),
    paste0("biomarker_", seq.int(n_correlated + 1, n_total))
  )

  betas <- c(betas_features, betas_biomarkers)

  weights <- c(weights_feature, rep(1, n_independent))

  dependent_variable <- as.matrix(mixtures) %*% (betas * weights) + rnorm(n_obs)

  mixtures <- as.data.frame(mixtures)
  mixtures$dependent_variable <- dependent_variable
  attr(mixtures, "true_weights") <- weights
  mixtures
}

weights <- c(5:1, rep(0, 5))
weights <- weights/sum(weights)

dataset <- generate_simulated_dataset(n_obs = 1000, n_total = 15,
                                      correlation_matrix = sigma,
                                      weights_feature = weights)
names(dataset)
feature_names <- names(dataset)[1:10]



dataset |> dplyr::select(starts_with("feature_")) |> cor(use = "pairwise.complete.obs") |>
  corrplot(
    method = "number",
    type = "upper",
    tl.col = "black"
  )

summary(lm(dependent_variable ~ . , data = dataset))

form <- paste0("dependent_variable ~ wqs + ", paste0("biomarker_", 11:15, collapse = "+"))

fit_gwqs <- gwqs(as.formula(form),
                  mix_name = feature_names, data = as.data.frame(dataset),
                  q = 4, validation = 0.6, b = 50, rh = 10,
                  family = "gaussian", seed = 2016)

summary(fit_gwqs)


gwqs_barplot(fit_gwqs)

names(dataset)

dataset <- generate_simulated_dataset(n_samples = 1000,
                                      correlation_min = 0.8,
                                      correlation_max = 0.95)

feature_names <- names(dataset)[1:25]

fit_gwqs <- gwqs(dependent_variable ~ pwqs + nwqs, mix_name =  feature_names, data = dataset,
                 q = 4, validation = 0.6, b = 50, rh = 20,
                 family = "gaussian", seed = 2016)
summary(fit_gwqs)

cor(dataset$feature_5, dataset$feature_6)
cor(dataset$feature_4, dataset$feature_3)

features <- dataset[paste0("feature_", 1:25)]
cor_matrix <- cor(features, use = "pairwise.complete.obs")
View(cor_matrix)
gwqs_barplot(fit_gwqs)


### =======================================================
library(glmnet)

# Build design matrix (all predictors except dependent_variable)
X <- model.matrix(dependent_variable ~ ., data = dataset)[, -1]  # drop intercept
y <- dataset$dependent_variable

# Fit Lasso path (alpha = 1 = Lasso)
fit_lasso <- glmnet(X, y, alpha = 1, standardize = TRUE)

# Cross-validation to choose lambda
set.seed(123)
cvfit <- cv.glmnet(X, y, alpha = 1, standardize = TRUE)

# Best lambda
best_lambda <- cvfit$lambda.min

# Coefficients at best lambda
lasso_coefs <- coef(cvfit, s = "lambda.min")

print(best_lambda)
print(lasso_coefs)

fit_lasso_best_lambda <- glmnet(X, y, alpha = 1, lambda = best_lambda,
                                standardize = TRUE)
set.seed(123)
cvfit <- cv.glmnet(X, y, alpha = 1, lambda = best_lambda)

lasso_coefs <- coef(cvfit, s = "lambda.min")
