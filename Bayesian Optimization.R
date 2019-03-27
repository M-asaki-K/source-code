library(rBayesianOptimization)
Gauss_holdout <- function(deg, c, q, m){
  model <- GPmodel
  t.Pred <- predict(model, data.frame(deg, c, q, m))
  Pred <- t.Pred$Y_hat + 2*(t.Pred$MSE)**0.5
  list(Score=Pred, Pred=Pred)}

opt_svm <- BayesianOptimization(Gauss_holdout,bounds=list(deg=c(0,1),c=c(0,1),q =c(0,1), m = c(0,1)),init_points=30, n_iter=5, acq='ei', kappa=2.576,eps=0.0, verbose=TRUE)

#---------------------in the case of experimental condition optimization-------------
fn <- function(a){
  a1 <- a[1]
  a2 <- a[2]
  a3 <- a[3]
  a4 <- a[4]
  
  model <- GPmodel
  t.Pred <- predict(model, data.frame(a1, a2, a3, a4))
  (t.Pred$Y_hat + 2*(t.Pred$MSE)**0.5)}

optim(c(0.3,0.3,0.3,0.3), fn, method = "L-BFGS-B", lower = c(0,0,0,0), upper = c(1,1,1,1),control=list(fnscale=-1))
