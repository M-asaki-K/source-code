library(rBayesianOptimization)
Gauss_holdout <- function(a, b){
  Pred <- -a**2 - b**2 + a + b - 1
  list(Score=Pred, Pred=Pred)}

opt_svm <- BayesianOptimization(Gauss_holdout,bounds=list(a=c(0,1),b=c(0,1)),init_points=30, n_iter=5, acq='ei', kappa=2.576,eps=0.0, verbose=TRUE)
