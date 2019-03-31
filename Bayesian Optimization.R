#--------------------------Bayesian optimization for hyperparameter tuning---------------------------
library(rBayesianOptimization)
Gauss_holdout <- function(deg, c, q, m){
  model <- GPmodel
  t.Pred <- predict(model, data.frame(deg, c, q, m))
  Pred <- t.Pred$Y_hat + 2*(t.Pred$MSE)**0.5
  list(Score=Pred, Pred=Pred)}

opt_svm <- BayesianOptimization(Gauss_holdout,bounds=list(deg=c(0,1),c=c(0,1),q =c(0,1), m = c(0,1)),init_points=30, n_iter=5, acq='ei', kappa=2.576,eps=0.0, verbose=TRUE)

#----------------------------Bayesian optimization_grid search--------------------------------------
parameter_list <- list(a1 = seq(0,1,length.out = 10),
                       a2 = seq(0,1,length.out = 10),
                       a3 = seq(0,1,length.out = 10),
                       a4 = seq(0,1,length.out = 10))

grid <- unname(expand.grid(parameter_list, KEEP.OUT.ATTRS = FALSE)) #form grid
View(grid)
pred <- predict.GP(GPmodel,grid)
acq_func <- function(mu, var, k) {
  acqs <- mu + k * sqrt(var)
  which.max(acqs)
}

ind <-grid[acq_func(pred$Y_hat,pred$MSE,2),]
ind #print optimum parameter
max(pred$Y_hat+2*pred$MSE**0.5) #print maximum

#---------------------Bayesian optimization_find optim-------------
fn <- function(a){
  a1 <- a[1]
  a2 <- a[2]
  a3 <- a[3]
  a4 <- a[4]
  
  model <- GPmodel
  t.Pred <- predict(model, data.frame(a1, a2, a3, a4))
  (t.Pred$Y_hat + 2*(t.Pred$MSE)**0.5)}

optim(c(1,0.7777778,0,0), fn, method = "L-BFGS-B", lower = c(0,0,0,0), upper = c(1,1,1,1),control=list(fnscale=-1))

