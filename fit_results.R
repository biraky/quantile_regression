#! /nfs/sw/eb/software/R/4.4.3-gfbf-2023b/bin/Rscript
#SBATCH -c 20
#SBATCH -t 01:00:00
#SBATCH --mem 30000

library(dplyr)
library(parallel)
library(Matrix)
library(MASS)
library(gWQS)
library(glmnet)

## Functions ===================================================================
generate_simulated_dataset <- \(
  n_obs = 5,
  n_covariates = 25,
  n_correlated = 10,
  correlation_matrix,
  weights_feature = rep(1, n_correlated),
  betas_features =  rep(1, n_correlated),
  betas_biomarkers = rep(1, n_covariates - n_correlated)
)
{
  n_independent <- n_covariates - n_correlated

  # Generate the mixture components
  correlated_mix <- mvrnorm(n_obs, mu = rep(0, n_correlated), Sigma = correlation_matrix)
  independent_covs <- matrix(rnorm(n_obs * n_independent, 0, 1), ncol = n_independent)

  mixtures <- cbind(correlated_mix, independent_covs)
  colnames(mixtures) <- c(paste0("feature_", seq_len(n_correlated)),
                          paste0("biomarker_", seq.int(n_correlated + 1, n_covariates)))

  betas <- c(betas_features, betas_biomarkers)

  weights <- c(weights_feature, rep(1, n_independent)) # all 1 = no weights by default

  dependent_variable <- as.matrix(mixtures) %*% (betas * weights) + rnorm(n_obs)

  mixtures <- as.data.frame(mixtures)
  mixtures$dependent_variable <- dependent_variable
  attr(mixtures, "true_weights") <- weights
  mixtures
}

fit_wqs_once <- \(
  dataset,
  feature_names,
  biomarker_names = paste0("biomarker_", 11:15),
  beta_intercept_true = 0,
  beta_biomarkers_true = rep(1, length(biomarker_names))
) {
  form <- as.formula(paste0(
    "dependent_variable ~ wqs + ",
    paste(biomarker_names, collapse = " + ")
  ))


  fit <- gwqs(
    as.formula(form),
    mix_name = feature_names,
    data = as.data.frame(dataset),
    q = 4,
    validation = 0.6,
    b = 50,
    rh = 10,
    family = "gaussian",
    seed = 1
  )

  # --- Coefficients and weights
  coefs <- coef(fit)
  coef_intercept <- unname(coefs[["(Intercept)"]])
  coef_biomk     <- unname(coefs[biomarker_names])

  # --- Biases (estimate - truth)
  bias_intercept <- coef_intercept - beta_intercept_true
  # bias_features       <- coef_wqs       - beta_wqs_true
  bias_biomk     <- coef_biomk - beta_biomarkers_true
  names(bias_biomk) <- biomarker_names

  rmse <- NA_real_

  y     <- dataset$dependent_variable
  rmse  <- sqrt(mean((y - unname(predict(
    fit
  )))^2))


  list(
    fit                 = fit,
    estimates           = list(
      intercept = coef_intercept,
      final_weights_table = fit$final_weights,
      # weights and CIs
      correctly_removed_proportion = sum(fit$final_weights[, 2][6:10] <
                                           1 / 10) / 5,
      biomarkers = setNames(coef_biomk, biomarker_names)
    ),
    bias                = list(
      intercept = bias_intercept,
      # wqs       = bias_wqs,
      biomarkers = bias_biomk
    ),
    rmse                = rmse
  )
}

fit_lasso_once <- \(dataset,
                    beta_intercept_true = 0,
                    beta_true = c(rep(1, 5), rep(0, 5), rep(1, 5))) {
  x <- model.matrix(dependent_variable ~ . - 1, data = dataset)
  y <- dataset$dependent_variable

  # --- LASSO with CV
  cv_lasso <- cv.glmnet(
    x,
    y,
    family      = "gaussian",
    alpha       = 1,
    # LASSO
    nfolds      = 10,
    standardize = TRUE
  )

  get_results <- \(lambda_choice) {
    coefs <- coef(cv_lasso, s = lambda_choice)
    coef_intercept <- unname(coefs[1, ])
    coef_pred      <- as.numeric(coefs[-1, ])
    names(coef_pred) <- colnames(x)

    correctly_removed_prop <- sum(coef_pred[6:10] == 0) / 5

    bias_intercept <- coef_intercept - beta_intercept_true
    bias_pred      <- coef_pred - beta_true
    names(bias_pred) <- colnames(x)

    preds <- predict(cv_lasso, newx = x, s = lambda_choice)
    rmse  <- sqrt(mean((y - preds)^2))

    list(
      intercept  = coef_intercept,
      predictors = coef_pred,
      correctly_removed_proportion = correctly_removed_prop,
      bias       = list(intercept  = bias_intercept, predictors = bias_pred),
      rmse = rmse
    )
  }

  list(
    fit   = cv_lasso,
    min   = get_results("lambda.min"),
    oneSE = get_results("lambda.1se")
  )
}

fit_adapt_lasso_once <- \(dataset,
                          beta_intercept_true = 0,
                          beta_true = c(rep(1, 5), rep(0, 5), rep(1, 5))) {
  x <- model.matrix(dependent_variable ~ . - 1, data = dataset)
  y <- dataset$dependent_variable

  # init fit
  cv_ridge <- cv.glmnet(
    x,
    y,
    family = "gaussian",
    alpha = 0,
    # ridge
    nfolds = 10,
    standardize = TRUE
  )
  # Coefficients from the initial estimator at lambda.min
  # coef() vector including the intercept in the first entry.
  beta_init <- as.numeric(coef(cv_ridge, s = "lambda.min"))
  beta_init_no_int <- beta_init[-1]  # drop intercept

  # build adaptive weights w_j = 1 / (|beta_init_j| + eps)^gamma
  gamma <- 1.0
  eps   <- 1e-6
  w <- 1 / (abs(beta_init_no_int) + eps)^gamma

  #weighted LASSO via penalty.factor
  cv_adalasso <- cv.glmnet(
    x,
    y,
    family = "gaussian",
    alpha = 1,
    # LASSO
    nfolds = 10,
    standardize = TRUE,
    penalty.factor = w         # <- adaptive part
  )

  get_results <- \(lambda_choice) {
    coefs <- coef(cv_adalasso, s = lambda_choice)
    coef_intercept <- unname(coefs[1, ])
    coef_pred      <- as.numeric(coefs[-1, ])
    names(coef_pred) <- colnames(x)

    correctly_removed_prop <- sum(coef_pred[6:10] == 0) / 5

    bias_intercept <- coef_intercept - beta_intercept_true
    bias_pred      <- coef_pred - beta_true
    names(bias_pred) <- colnames(x)

    preds <- predict(cv_adalasso, newx = x, s = lambda_choice)
    rmse  <- sqrt(mean((y - preds)^2))

    list(
      intercept  = coef_intercept,
      predictors = coef_pred,
      correctly_removed_proportion = correctly_removed_prop,
      bias       = list(intercept  = bias_intercept, predictors = bias_pred),
      rmse = rmse
    )
  }

  list(
    fit   = cv_adalasso,
    min   = get_results("lambda.min"),
    oneSE = get_results("lambda.1se")
  )
}
fit_elastic_net_once <- \(
  dataset,
  alpha        = 0.5,
  beta_intercept_true = 0,
  beta_true    = c(rep(1, 5), rep(0, 5), rep(1, 5))
) {

  x <- model.matrix(dependent_variable ~ . - 1, data = dataset)
  y <- dataset$dependent_variable

  # ---- Elastic net with K-fold CV at chosen alpha
  cv_en <- cv.glmnet(
    x, y,
    family      = "gaussian",
    alpha       = alpha,
    nfolds      = 10,
    standardize = TRUE
  )

  get_results <- \(lambda_choice) {
    coefs <- coef(cv_en, s = lambda_choice)

    coef_intercept <- unname(coefs[1, ])
    coef_pred      <- as.numeric(coefs[-1, ])
    names(coef_pred) <- colnames(x)

    # Proportion of the truly null block (positions 6:10) set exactly to zero
    correctly_removed_prop <- mean(coef_pred[6:10] == 0)

    bias_intercept <- coef_intercept - beta_intercept_true
    bias_pred      <- coef_pred - beta_true
    names(bias_pred) <- colnames(x)

    preds <- as.numeric(predict(cv_en, newx = x, s = lambda_choice))
    rmse  <- sqrt(mean((y - preds)^2))

    list(
      intercept  = coef_intercept,
      predictors = coef_pred,
      correctly_removed_proportion = correctly_removed_prop,
      bias       = list(
        intercept  = bias_intercept,
        predictors = bias_pred
      ),
      rmse = rmse
    )
  }

  list(
    fit    = cv_en,
    alpha  = alpha,
    min    = get_results("lambda.min"),
    oneSE  = get_results("lambda.1se")
  )
}
# --- helper to collect stats from existing fit_* wrappers
collect_results <- \(dataset, feature_names) {
  lasso <- fit_lasso_once(dataset)
  adapt <- fit_adapt_lasso_once(dataset)
  wqs   <- fit_wqs_once(dataset, feature_names)
  enet  <- fit_elastic_net_once(dataset)


  grab <- \(obj, submodel, model_name) {
    tibble(
      model = model_name,
      bias_intercept   = 100 * obj[[submodel]]$bias$intercept,
      bias_biomarkers  = mean(abs(100 * unlist(
        obj[[submodel]]$bias$predictors[grepl("^biomarker_",
                                              names(obj[[submodel]]$bias$predictors))]
      ))),
      correctly_removed = obj[[submodel]]$correctly_removed_proportion,
      rmse             = obj[[submodel]]$rmse
    )
  }

  grab_wqs <- \(obj, model_name) {
    tibble(
      model = model_name,
      bias_intercept   = 100 * obj$bias$intercept,
      bias_biomarkers  = mean(abs(100 * unlist(obj$bias$biomarkers))),
      correctly_removed = obj$estimates$correctly_removed_proportion,
      rmse             = obj$rmse
    )
  }

  bind_rows(
    grab(lasso, "min", "Lasso_min"),
    grab(lasso, "oneSE", "Lasso_1se"),
    grab(adapt, "min", "AdaptLasso_min"),
    grab(adapt, "oneSE", "AdaptLasso_1se"),
    grab(enet,  "min",   "Elastic_Net_min"),
    grab(enet,  "oneSE", "Elastic_Net_1se"),
    grab_wqs(wqs, "WQS")
  )
}

run_sims <- \(
  sample_sizes,
  reps        = 10,
  cores       = max(1, detectCores() - 2),
  seed        = 1,
  save_dir    = "/nfs/home/biraba/RUNNING/gqws/"
) {
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)

  jobs <- expand.grid(sample_size = sample_sizes, rep = seq_len(reps), KEEP.OUT.ATTRS = FALSE)

  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)

  chunks <- mclapply(
    seq_len(nrow(jobs)),
    \(i) {
      ss  <- jobs$sample_size[i]
      rep <- jobs$rep[i]

      dataset <- generate_simulated_dataset(
        n_obs              = ss,
        n_covariates       = 15,
        correlation_matrix = vcor,
        weights_feature    = weights
      )

      feature_names <- paste0("feature_", 1:10)
      res_tbl <- collect_results(dataset, feature_names) |>
        dplyr::mutate(sample_size = nrow(dataset), rep = rep, .before = 1)

      key <- sprintf("n%d_rep%03d", ss, rep)   # <- list key

      list(name = key, data = dataset, res = res_tbl)
    },
    mc.cores    = cores,
    mc.set.seed = TRUE
  )

  all_res <- dplyr::bind_rows(lapply(chunks, `[[`, "res"))
  datasets <- setNames(lapply(chunks, `[[`, "data"),
                       vapply(chunks, `[[`, character(1), "name"))

  # Save both results and datasets in one .RData
  save(all_res, datasets, file = file.path(save_dir, "model_results_and_datasets.RData"))

 # write.csv(all_res, file.path(save_dir, "model_results.csv"), row.names = FALSE)

  write.csv(dataframe,"~/Downloads/filename.csv", row.names = FALSE)

  list(results = all_res, datasets = datasets)
}


## =============================================================================
vcor <- outer(c(seq(0.6, 0.8, length.out = 5), seq(0.92, 0.98, length.out = 5)),
              c(seq(0.6, 0.8, length.out = 5), seq(0.92, 0.98, length.out = 5)))
diag(vcor) <- rep(1, 10)

weights <- c(rep(1, 5), rep(0, 5))
weights <- weights / sum(weights)

results <- run_sims(
  sample_sizes = c(100, 500, 1000),
  reps = 10,
  cores = 19,
  seed = 1
)

results
