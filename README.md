# ML_doubly_robust

This repository provides R code for conducting **Doubly Robust Estimation** of the **Average Treatment Effect (ATE)** using modern **machine learning techniques**. Developed for a working paper, the package includes two core components:

- `DR_AIPW.R`: Implements Augmented Inverse Probability Weighting (AIPW)
- `DR_TMLE.R`: Implements Targeted Maximum Likelihood Estimation (TMLE)

## Overview

The code enables users to estimate causal effects with **high-dimensional confounders** by combining rigorous statistical methods with flexible machine learning models. This approach provides **consistency** and **efficiency** under less restrictive modeling assumptions, making it especially suitable for **observational studies**.

### Key Features

- 🧠 Implements **two doubly robust estimators**: AIPW and TMLE
- 🔍 Supports **five machine learning models** for nuisance estimation:
  - `glmnet`: Regularized regression (LASSO / Elastic Net)
  - `ranger`: Random Forest
  - `xgboost`: Gradient Boosted Trees
  - `e1071`: Support Vector Machines (SVM)
  - `nnet`: Neural Networks
- 📊 Automatically compares ATE estimates across ML methods
- 📎 Based on theory from an accompanying working paper

## Applications

- Causal inference in observational data
- Policy evaluation
- High-dimensional data analysis

## File Structure

- `DR_AIPW.R`: AIPW estimation with multiple ML learners
- `DR_TMLE.R`: TMLE estimation with multiple ML learners

## Getting Started

To use this package, simply run either file after ensuring the necessary R packages are installed:

```R
source("DR_AIPW.R")  # for AIPW estimation
source("DR_TMLE.R")  # for TMLE estimation
