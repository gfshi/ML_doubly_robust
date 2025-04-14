This R package provides tools for estimating the Average Treatment Effect (ATE) using doubly robust machine learning methods. 
It implements both Augmented Inverse Probability Weighting (AIPW) and Targeted Maximum Likelihood Estimation (TMLE) with difference 
machine learning techniques, allowing users to obtain consistent and efficient causal effect estimates under the presence of high-dimensional confounders.
By combining machine learning algorithms with robust statistical theory, the package offers a flexible and reliable framework for 
causal inference in observational studies.

The code includes five machine learning tools: regularized regression using glmnet, random forest implemented via ranger, 
boosting with xgboost, support vector machines using e1071, and neural networks with nnet.

At the end of the analysis, users can compare the results across these five machine learning methods.
