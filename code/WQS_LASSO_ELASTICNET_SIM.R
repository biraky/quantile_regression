# ============================================================
# Comparative Sensitivity Analysis: gWQS vs Regularized Regressions
# Using Simulated Dataset (25 mixtures, 10 correlated)
# Monte Carlo Simulation with 1000 Replications
# RMSE-focused comparison with comprehensive visualizations
# ============================================================

# Load required libraries
library(gWQS)
library(glmnet)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(MASS)
library(Matrix)
library(parallel)

# -------------------------------
# 1. Define Core Functions
# -------------------------------
generate_simulated_dataset <- function(n_samples = 5,
                                       correlation_min = 0.4,
                                       correlation_max = 0.7)
                                       {
  n_total <- 25
  n_correlated <- 10
  n_independent <- n_total - n_correlated

  # Build correlation matrix for correlated mixtures
  corr_vals <- runif(n_correlated * (n_correlated - 1) / 2,
                     correlation_min, correlation_max)
  sigma_corr <- diag(1, n_correlated)
  sigma_corr[upper.tri(sigma_corr)] <- corr_vals
  sigma_corr[lower.tri(sigma_corr)] <- t(sigma_corr)[lower.tri(sigma_corr)]

  # Ensure the matrix is positive definite
  sigma_corr <- as.matrix(nearPD(sigma_corr)$mat)

  # Generate the mixture components
  correlated_mix <- MASS::mvrnorm(n_samples, mu = rep(0, n_correlated), Sigma = sigma_corr)
  independent_mix <- matrix(rnorm(n_samples * n_independent, 0, 1), ncol = n_independent)

  mixtures <- cbind(correlated_mix, independent_mix)
  colnames(mixtures) <- paste0("feature_", 1:n_total)
  mixtures <- as.data.frame(mixtures)

  # Generate true weights
  true_weights <- rep(0, n_total)
  active_idx <- sample(1:n_total, 10)
  raw_weights <- runif(10, 0.05, 0.15)
  normalized_weights <- raw_weights / sum(raw_weights)
  true_weights[active_idx] <- normalized_weights

  # Generate outcome
  wqs_true <- as.matrix(mixtures) %*% true_weights
  y <- 2 + 0.5 * wqs_true + rnorm(n_samples, 0, 0.5)

  mixtures$yLBX <- y
  return(list(data = mixtures, true_weights = true_weights))
}

run_gwqs <- function(data, n, seed) {
  set.seed(seed)
  sub_data <- data[sample(1:nrow(data), n, replace = FALSE), ]
  PCBs <- names(sub_data)[1:25]

  result <- tryCatch({
    gwqs(
      formula = yLBX ~ wqs,
      mix_name = PCBs,
      data = sub_data,
      q = 4,
      validation = 0.6,
      b = 100,
      b1_pos = TRUE,
      b_constr = FALSE,
      family = "gaussian",
      seed = seed
    )
  }, error = function(e) {
    return(NULL)
  })
  return(result)
}

run_glmnet <- function(X, y, n, alpha, seed) {
  set.seed(seed)
  sample_idx <- sample(1:nrow(X), n, replace = FALSE)
  X_sub <- X[sample_idx, ]
  y_sub <- y[sample_idx]

  train_idx <- sample(1:n, size = round(0.6 * n), replace = FALSE)
  X_train <- X_sub[train_idx, ]
  X_test <- X_sub[-train_idx, ]
  y_train <- y_sub[train_idx]
  y_test <- y_sub[-train_idx]

  cv_fit <- cv.glmnet(X_train, y_train, alpha = alpha, nfolds = 5)
  coefs <- as.matrix(coef(cv_fit, s = "lambda.min"))

  predictions <- predict(cv_fit, newx = X_test, s = "lambda.min")
  test_mse <- mean((y_test - predictions)^2)
  test_rmse <- sqrt(test_mse)
  test_r2 <- 1 - (sum((y_test - predictions)^2) / sum((y_test - mean(y_test))^2))

  return(list(
    coefficients = coefs[-1, 1],
    test_rmse = test_rmse,
    test_r2 = test_r2
  ))
}

run_adaptive_lasso <- function(X, y, n, seed, gamma = 1) {
  set.seed(seed)
  sample_idx <- sample(1:nrow(X), n, replace = FALSE)
  X_sub <- X[sample_idx, ]
  y_sub <- y[sample_idx]

  train_idx <- sample(1:n, size = round(0.6 * n), replace = FALSE)
  X_train <- X_sub[train_idx, ]
  X_test <- X_sub[-train_idx, ]
  y_train <- y_sub[train_idx]
  y_test <- y_sub[-train_idx]

  ols_fit <- lm(y_train ~ X_train)
  init_coefs <- coef(ols_fit)[-1]
  adaptive_weights <- 1 / (abs(init_coefs)^gamma + 1e-6)

  cv_fit <- cv.glmnet(X_train, y_train, alpha = 1, penalty.factor = adaptive_weights, nfolds = 5)
  coefs <- as.matrix(coef(cv_fit, s = "lambda.min"))
  predictions <- predict(cv_fit, newx = X_test, s = "lambda.min")

  test_mse <- mean((y_test - predictions)^2)
  test_rmse <- sqrt(test_mse)
  test_r2 <- 1 - (sum((y_test - predictions)^2) / sum((y_test - mean(y_test))^2))

  return(list(
    coefficients = coefs[-1, 1],
    test_rmse = test_rmse,
    test_r2 = test_r2
  ))
}

run_one_simulation <- function(rep_id, n) {
  sim_data <- generate_simulated_dataset(n_samples = 5000)
  data_pool <- sim_data$data
  X_pool <- as.matrix(data_pool[, 1:25])
  y_pool <- data_pool$yLBX
  current_true_weights <- sim_data$true_weights

  results <- list()
  seed <- 123 + rep_id

  gwqs_result <- run_gwqs(data_pool, n, seed)
  if (!is.null(gwqs_result)) {
    sumry <- summary(gwqs_result)
    fitted_vals <- gwqs_result$fit$fitted.values
    test_rmse <- sqrt(mean((gwqs_result$data$yLBX - fitted_vals)^2))
    gwqs_weights <- gwqs_result$final_weights$mean_weight
    names(gwqs_weights) <- gwqs_result$final_weights$mix_name
    gwqs_weights_ordered <- gwqs_weights[paste0("PCB", 1:25)]

    spearman_rho <- cor(gwqs_weights_ordered, current_true_weights, method = "spearman", use = "complete.obs")

    results$gwqs <- list(
      beta = sumry$coefficients["wqs", "Estimate"],
      pval = sumry$coefficients["wqs", "Pr(>|t|)"],
      test_rmse = test_rmse,
      spearman_rho = spearman_rho,
      weights = gwqs_weights_ordered,
      coefficients = gwqs_weights_ordered
    )
  } else {
    results$gwqs <- list(
      beta = NA, pval = NA, test_rmse = NA, spearman_rho = NA,
      weights = setNames(rep(NA, 25), paste0("PCB", 1:25)),
      coefficients = setNames(rep(NA, 25), paste0("PCB", 1:25))
    )
  }

  results$lasso <- run_glmnet(X_pool, y_pool, n, alpha = 1, seed)
  results$elastic_net <- run_glmnet(X_pool, y_pool, n, alpha = 0.5, seed)
  results$adaptive_lasso <- run_adaptive_lasso(X_pool, y_pool, n, seed)

  return(results)
}

# -------------------------------
# 2. Set Simulation Parameters
# -------------------------------
n_replications <- 1000
sample_sizes <- c(100, 500, 1000, 2500, 5000)

# -------------------------------
# 3. Run the Monte Carlo Simulation and Return DataFrames
# -------------------------------
run_simulation_and_return_dfs <- function() {
  performance_all_df <- data.frame()
  coefficients_all_df <- data.frame()

  for (n in sample_sizes) {
    cat("Starting simulations for sample size n =", n, "\n")

    num_cores <- detectCores() - 1
    cl <- makeCluster(num_cores)

    clusterExport(cl, varlist = c("n", "generate_simulated_dataset", "run_gwqs",
                                  "run_glmnet", "run_adaptive_lasso", "run_one_simulation",
                                  "nearPD", "gwqs", "cv.glmnet", "lm"))

    clusterEvalQ(cl, {
      library(gWQS)
      library(glmnet)
      library(MASS)
      library(Matrix)
      library(mgcv)
    })

    results_n <- parLapply(cl, 1:n_replications, function(rep_id) {
      run_one_simulation(rep_id, n)
    })

    stopCluster(cl)
    cat("Completed", n_replications, "replications for n =", n, "\n")

    performance_temp_df <- data.frame()
    coefficients_temp_df <- data.frame()

    for (r in 1:n_replications) {
      rep_result <- results_n[[r]]

      perf_row <- data.frame(
        replication = r,
        sample_size = n,
        model = c("gWQS", "LASSO", "ElasticNet", "AdaptiveLASSO"),
        rmse = c(rep_result$gwqs$test_rmse,
                 rep_result$lasso$test_rmse,
                 rep_result$elastic_net$test_rmse,
                 rep_result$adaptive_lasso$test_rmse)
      )
      performance_temp_df <- rbind(performance_temp_df, perf_row)

      for (model_name in c("gwqs", "lasso", "elastic_net", "adaptive_lasso")) {
        model_label <- case_when(
          model_name == "gwqs" ~ "gWQS",
          model_name == "lasso" ~ "LASSO",
          model_name == "elastic_net" ~ "ElasticNet",
          model_name == "adaptive_lasso" ~ "AdaptiveLASSO"
        )

        coefs <- rep_result[[model_name]]$coefficients
        if (length(coefs) == 25) {
          coef_rows <- data.frame(
            replication = r,
            sample_size = n,
            model = model_label,
            variable = paste0("PCB", 1:25),
            coefficient = coefs
          )
          coefficients_temp_df <- rbind(coefficients_temp_df, coef_rows)
        }
      }
    }

    mean_performance <- performance_temp_df %>%
      group_by(sample_size, model) %>%
      summarise(
        mean_rmse = mean(rmse, na.rm = TRUE),
        rmse_sd = sd(rmse, na.rm = TRUE),
        rmse_se = sd(rmse, na.rm = TRUE) / sqrt(n()),
        .groups = 'drop'
      )
    performance_all_df <- rbind(performance_all_df, mean_performance)

    mean_coefficients <- coefficients_temp_df %>%
      group_by(sample_size, model, variable) %>%
      summarise(
        mean_coefficient = mean(coefficient, na.rm = TRUE),
        coefficient_sd = sd(coefficient, na.rm = TRUE),
        .groups = 'drop'
      )
    coefficients_all_df <- rbind(coefficients_all_df, mean_coefficients)

    saveRDS(list(performance = performance_all_df, coefficients = coefficients_all_df),
            paste0("intermediate_results_n", n, ".rds"))
    cat("Processed results for n =", n, "\n")
  }

  return(list(performance_df = performance_all_df, coefficients_df = coefficients_all_df))
}

# -------------------------------
# 4. Run Simulation and Get DataFrames
# -------------------------------
cat("Starting Monte Carlo simulation...\n")
simulation_start_time <- Sys.time()
results <- run_simulation_and_return_dfs()
simulation_end_time <- Sys.time()
cat("Simulation completed in", round(difftime(simulation_end_time, simulation_start_time, units = "mins"), 1), "minutes\n")

performance_df <- results$performance_df
coefficients_df <- results$coefficients_df

# -------------------------------
# 5. Save Results as DataFrames
# -------------------------------
write.csv(performance_df, "simulation_performance_results.csv", row.names = FALSE)
write.csv(coefficients_df, "simulation_coefficient_results.csv", row.names = FALSE)

# -------------------------------
# 6. Display Results
# -------------------------------
cat("\n=== PERFORMANCE RESULTS (RMSE comparison) ===\n")
print(performance_df)

cat("\n=== MODEL RANKING BY SAMPLE SIZE ===\n")
performance_ranking <- performance_df %>%
  dplyr::group_by(sample_size) %>%
  dplyr::arrange(mean_rmse) %>%
  dplyr::mutate(rank = 1:dplyr::n()) %>%
  dplyr::select(sample_size, model, mean_rmse, rank)
print(performance_ranking)

cat("\n=== SUMMARY STATISTICS ===\n")
cat("Performance dataframe dimensions:", dim(performance_df), "\n")
cat("Coefficients dataframe dimensions:", dim(coefficients_df), "\n")
cat("Unique models in performance_df:", toString(unique(performance_df$model)), "\n")
cat("Unique models in coefficients_df:", toString(unique(coefficients_df$model)), "\n")
cat("Sample sizes analyzed:", toString(unique(performance_df$sample_size)), "\n")
