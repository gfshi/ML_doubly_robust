DR_AIPW_ML <- function(dat, ml_lib = "ranger", glmnet_alpha = NULL) {
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
  Y <- dat[, 1]          # First column as outcome
  D <- dat[, 2]          # Second column as treatment
  X <- dat[, -c(1, 2)]   # Controls
  df <- data.frame(Y = Y, D = D, X)   # Rename Y and D
  
  # Validate ML method
  valid_methods <- c("ranger", "glmnet", "xgboost", "svm", "nnet")
  if (!ml_lib %in% valid_methods) {
    stop(paste("Invalid ML method. Supported methods are:", paste(valid_methods, collapse = ", ")))
  }
  
  # Outcome Model (mu1 and mu0)
  if (ml_lib == "ranger") {
    mod_Y <- ranger::ranger(Y ~ ., data = df, 
                            num.trees = 2000, probability = FALSE)
  } else if (ml_lib == "glmnet") {
    X_mat <- model.matrix(~ . - 1, data = df[, -1])
    #X_mat <- model.matrix(~ . , data = df[, -1])  #include intercept
    
    # identify best alpha, Lasso, Elastic net, Ridge
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
                 " (0 = Ridge, 1 = Lasso, intermediate = Elastic Net)\n"))
      
      mod_Y <- glmnet::cv.glmnet(x = X_mat, y = Y, family = "gaussian", alpha = best_alpha)
    } else {
      # With alpha specified as "glmnet_alpha"
      mod_Y <- glmnet::cv.glmnet(x = X_mat, y = Y, family = "gaussian", 
                                 alpha = glmnet_alpha)
      best_alpha <- glmnet_alpha
    }
  } else if (ml_lib == "xgboost") {
    xg_mat <- model.matrix(~ . - 1, data = df[, -1])
    mod_Y <- xgboost::xgboost(data = xg_mat, label = Y, nrounds = 100, 
                              objective = "reg:squarederror", verbose = 0)
  } else if (ml_lib == "svm") {
    mod_Y <- e1071::svm(Y ~ ., data = df, probability = FALSE)
  } else if (ml_lib == "nnet") {
    mod_Y <- nnet::nnet(Y ~ ., data = df, size = 5, linout = TRUE, trace = FALSE)
  }
  
  
  # Predict potential outcomes
  predict_outcomes <- function(d_val) {
    new_data <- df[, -1]
    new_data$D <- d_val
    
    if (ml_lib == "ranger") {
      predict(mod_Y, data = new_data)$predictions
    } else if (ml_lib == "glmnet") {
      new_mat <- model.matrix(~ . - 1, data = new_data)
      #new_mat <- model.matrix(~ . , data = new_data)   #include intercept
      as.numeric(predict(mod_Y, newx = new_mat, s = "lambda.min"))
    } else if (ml_lib == "xgboost") {
      new_mat <- model.matrix(~ . - 1, data = new_data)
      predict(mod_Y, newdata = new_mat)
    } else if (ml_lib %in% c("svm", "nnet")) {
      predict(mod_Y, newdata = new_data)
    }
  }
  
  # Get predicted potential outcomes
  mu1_hat <- predict_outcomes(1)
  mu0_hat <- predict_outcomes(0)
  

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
  df$e_ps <- pmin(pmax(df$e_ps, 0.05), 0.95)
  
  # Compute AIPW components
  G_comp <- mu1_hat - mu0_hat
  adj <- (D * (Y - mu1_hat) / df$e_ps) - ((1 - D) * (Y - mu0_hat) / (1 - df$e_ps))
  
  # ATE calculation
  ate_aipw <- mean(G_comp + adj)
  IC <- G_comp + adj - ate_aipw
  ate_se <- sqrt(var(IC) / length(IC))
  
  # ATT calculation
  G_comp_att <- mu1_hat[D == 1] - mu0_hat[D == 1]
  adj_att <- (D[D == 1] * (Y[D == 1] - mu1_hat[D == 1]) / df$e_ps[D == 1]) - 
    ((1 - D[D == 1]) * (Y[D == 1] - mu0_hat[D == 1]) / (1 - df$e_ps[D == 1]))
  
  att_aipw <- mean(G_comp_att + adj_att)
  IC_att <- (G_comp_att + adj_att) - att_aipw
  att_se <- sqrt(var(IC_att) / sum(D == 1))
  
  return(list(
    ate_aipw = ate_aipw,
    ate_se = ate_se,
    att_aipw = att_aipw,
    att_se = att_se,
    ml_method = ml_lib,
    glmnet_alpha = if (ml_lib == "glmnet") best_alpha else NULL
    ))
 
}
