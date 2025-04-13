DR_TMLE_ML <- function(dat, outcome = "continuous", ml_lib = "ranger", glmnet_alpha = NULL) {
  # Check and install packages
  required_pkgs <- c("tidyverse")
  if (ml_lib == "svm") required_pkgs <- c(required_pkgs, "e1071")
  else required_pkgs <- c(required_pkgs, ml_lib)
  
  missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
  if (length(missing_pkgs) > 0) {
    install.packages(missing_pkgs, quiet = TRUE)
  }
  
  # Load required packages
  suppressPackageStartupMessages({
    library(tidyverse)
    if (ml_lib == "svm") {
      if (!require("e1071", quietly = TRUE)) stop("Package e1071 is not available")
    } else if (!require(ml_lib, character.only = TRUE, quietly = TRUE)) {
      stop(paste("Package", ml_lib, "is not available"))
    }
  })
  
  
  # Prepare data
  Y <- dat[, 1]             # First column as the dependent variable
  D <- dat[, 2]             # Second column as treatment indicator
  X <- dat[, -c(1, 2)]      # Remaining columns as controls
  df <- data.frame(Y = Y, D = D, X)  # Rename Y and D
  
  # Validate outcome type
  if (!outcome %in% c("continuous", "binary")) {
    stop("Outcome must be either 'binary' or 'continuous'")
  }
  
  # Normalize Y if continuous
  if (outcome == "continuous") {
    Y_min <- min(Y)
    Y_max <- max(Y)
    df$Y_prob <- (Y - Y_min) / (Y_max - Y_min)
  } else {
    df$Y_prob <- Y
  }
  
  # Outcome Model - using specified ML method
  if (ml_lib == "ranger") {
    mod_Y <- ranger::ranger(Y_prob ~ ., data = df[, c("Y_prob", "D", colnames(X))], 
                            num.trees = 2000, probability = FALSE)
    
    predict_Y <- function(new_data) {
      predict(mod_Y, data = new_data)$predictions
    }
    
  } else if (ml_lib == "glmnet") {
    #X_mat <- model.matrix(~ . , data = df[, c("D", colnames(X))])
    X_mat <- model.matrix(~ . -1, data = data.frame(D=df$D, X))
    
    # If alpha is not specified, perform cross-validation to find best alpha
    if (is.null(glmnet_alpha)) {
      alpha_grid <- seq(0, 1, by = 0.1)
      alpha_table <- data.frame(alpha = numeric(), lambda_min = numeric(),
                                cvm_min = numeric(), stringsAsFactors = FALSE)
      
      for (k in alpha_grid) {
        cv_fit <- glmnet::cv.glmnet(x = X_mat, y = Y, 
                                    family = "gaussian", alpha = k)
        alpha_table <- rbind(alpha_table, 
                             data.frame(alpha = k,
                                        lambda_min = cv_fit$lambda.min,
                                        cvm_min = min(cv_fit$cvm),
                                        stringsAsFactors = FALSE))
      }
      
      # Best alpha for Regularization model 
      best_alpha <- alpha_table$alpha[which.min(alpha_table$cvm_min)]
      cat(paste0("\nBest alpha selected: ", best_alpha, 
                 " (0 = Ridge, 1 = Lasso, (0-1) = Elastic Net)\n"))
      
      mod_Y <- glmnet::cv.glmnet(x = X_mat, y = df$Y_prob, family = "gaussian", alpha = best_alpha)
    
    } else {
      # With alpha specified in model via "glmnet_alpha"
      mod_Y <- glmnet::cv.glmnet(x = X_mat, y = df$Y_prob, family = "gaussian", 
                                 alpha = glmnet_alpha)
      best_alpha <- glmnet_alpha 
    }
    predict_Y <- function(new_data) {
      new_mat <- model.matrix(~ . - 1, data = new_data)
      as.numeric(predict(mod_Y, newx = new_mat, s = "lambda.min"))
    }
    
  } else if (ml_lib == "xgboost") {
    xg_mat <- model.matrix(~ . - 1, data = data.frame(D = df$D, X))
    mod_Y <- xgboost::xgboost(data = xg_mat, label = df$Y_prob, nrounds = 100, 
                              objective = "reg:squarederror", verbose = 0)
    
    predict_Y <- function(new_data) {
      new_mat <- model.matrix(~ . - 1, data = new_data)
      predict(mod_Y, newdata = new_mat)
    }
    
  } else if (ml_lib == "svm") {
    mod_Y <- e1071::svm(Y_prob ~ ., data = df[, c("Y_prob", "D", colnames(X))], 
                        probability = FALSE)
    
    predict_Y <- function(new_data) {
      predict(mod_Y, newdata = new_data)
    }
    
  } else if (ml_lib == "nnet") {
    mod_Y <- nnet::nnet(Y_prob ~ ., data = df[, c("Y_prob", "D", colnames(X))], 
                        size = 5, linout = TRUE, trace = FALSE)
    
    predict_Y <- function(new_data) {
      predict(mod_Y, newdata = new_data)
    }
    
  } else {
    stop("Unsupported ML method specified")
  }
  
  # Get predictions Y_hat, Y1_hat and Y0_hat 
  df$Y_hat <- predict_Y(df[, c("D", colnames(X))])
  df$Y1_hat <- predict_Y(data.frame(D = 1, X))
  df$Y0_hat <- predict_Y(data.frame(D = 0, X))
  
  
  # Propensity score model with different ML method
  if (ml_lib == "ranger") {
    mod_ps <- ranger::ranger(D ~ ., data = df[, c("D", colnames(X))], 
                             num.trees = 2000, probability = TRUE)
    ps_pred <- predict(mod_ps, data = X)$predictions
    df$e_ps <- if (ncol(ps_pred) == 2) ps_pred[, 2] else ps_pred[, 1]
    
  } else if (ml_lib == "glmnet") {
    X_mat_ps <- model.matrix(~ . - 1, data = X)
    mod_ps <- glmnet::cv.glmnet(x = X_mat_ps, y = D, family = "binomial", 
                                alpha = ifelse(is.null(glmnet_alpha), 1, glmnet_alpha))
    ps_pred <- predict(mod_ps, newx = X_mat_ps, type = "response", s = "lambda.min")
    df$e_ps <- as.numeric(ps_pred)
    
  } else if (ml_lib == "xgboost") {
    xg_mat_ps <- model.matrix(~ . - 1, data = X)
    mod_ps <- xgboost::xgboost(data = xg_mat_ps, label = D, nrounds = 100,
                               objective = "binary:logistic", verbose = 0)
    df$e_ps <- predict(mod_ps, newdata = xg_mat_ps)
    
  } else if (ml_lib == "svm") {
    mod_ps <- e1071::svm(D ~ ., data = df[, colnames(X)], probability = TRUE)
    ps_pred <- predict(mod_ps, newdata = X, probability = TRUE)
    df$e_ps <- attr(ps_pred, "probabilities")[, "1"]
    
  } else if (ml_lib == "nnet") {
    mod_ps <- nnet::nnet(D ~ ., data = df[, c("D", colnames(X))], 
                         size = 5,
                         linout = FALSE,   #for binary D
                         trace = FALSE,
                         decay = 0.1, 
                         maxit = 500)
    df$e_ps <- predict(mod_ps, newdata = X)
  }
  
  # Avoid extreme propensity scores 
  #df$e_ps <- pmin(pmax(df$e_ps, 0.01), 0.99)
  df$e_ps <- pmin(pmax(df$e_ps, 0.05), 0.95)
  
  # Define clever covariates
  df$H1_D <- df$D / df$e_ps
  df$H0_D <- (1 - df$D) / (1 - df$e_ps)
  df$H_D <- df$H1_D - df$H0_D
  
 
  # Avoid extreme values of Y_hat 
  safe_qlogis <- function(p) {
    p <- pmin(pmax(p, 0.001), 0.999)
    qlogis(p)
  }
  
  # Fluctuation parameters: two-covariate approach 
  #epsilon_fit <- glm(Y_prob ~ -1 + H0_D + H1_D + offset(qlogis(Y_hat)), 
  #                   data = df, family = 'binomial')
  
  #Ensure Y_hat between 0.001 and 0.999
  epsilon_fit <- glm(Y_prob ~ -1 + H0_D + H1_D + offset(safe_qlogis(Y_hat)), 
                     data = df, family = 'binomial')
  eps0 <- coef(epsilon_fit)["H0_D"]
  eps1 <- coef(epsilon_fit)["H1_D"]
  
  
  # Update predictions on logit scale
  #df$Y1_update <- plogis(qlogis(df$Y1_hat) + eps1 * df$H1_D)
  #df$Y0_update <- plogis(qlogis(df$Y0_hat) + eps0 * df$H0_D)
  df$Y1_update <- plogis(safe_qlogis(df$Y1_hat) + eps1 * df$H1_D)
  df$Y0_update <- plogis(safe_qlogis(df$Y0_hat) + eps0 * df$H0_D)
  
  
  # Rescale if continuous outcome
  if (outcome == "continuous") {
    df$Y1_update_scaled <- (Y_max - Y_min) * df$Y1_update + Y_min
    df$Y0_update_scaled <- (Y_max - Y_min) * df$Y0_update + Y_min
  } else {
    df$Y1_update_scaled <- df$Y1_update
    df$Y0_update_scaled <- df$Y0_update
  }
  
  # TMLE estimates
  ate_tmle <- mean(df$Y1_update_scaled - df$Y0_update_scaled, na.rm = TRUE)
  att_tmle <- mean(df$Y1_update_scaled[df$D == 1] - df$Y0_update_scaled[df$D == 1], 
                   na.rm = TRUE)
  
  # EIC influence curve
  ey1_tmle <- mean(df$Y1_update_scaled)
  ey0_tmle <- mean(df$Y0_update_scaled)
  ic_D1 <- df$H1_D * (df$Y - df$Y1_update_scaled) + df$Y1_update_scaled - ey1_tmle
  ic_D0 <- df$H0_D * (df$Y - df$Y0_update_scaled) + df$Y0_update_scaled - ey0_tmle
  
  eic_ate <- ic_D1 - ic_D0
  ate_se <- sqrt(var(eic_ate, na.rm = TRUE) / length(eic_ate))
  
  eic_att <- eic_ate[df$D == 1]
  att_se <- sqrt(var(eic_att, na.rm = TRUE) / sum(df$D == 1))
  
  return(list(
    ate_tmle = ate_tmle,
    ate_se = ate_se,
    att_tmle = att_tmle,
    att_se = att_se,
    ml_method = ml_lib,
    glmnet_alpha = if (ml_lib == "glmnet") best_alpha else NULL
  ))
}
