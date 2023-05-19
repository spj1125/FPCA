# parallel computing

sim_function <- function(r){ # r = 1
  set.seed(r)
  inverse.weigh <- function(beta, z)  #inverse rate weighting
    exp(-(beta*z))
  
  
  #==============================
  
  #eigenfunction
  phi <- function(t,k) {
    val <- sqrt(2/3)*cos(k*pi*t)
    val
  }
  #mean function
  mu <- function(t)
  {
    val <- sin(t+1/2)
    val
  }
  
  
  
  Z.fun <- function(Xi, t, K=50, nu=2)
  {
    #Xi <- runif(K,-sqrt(3),sqrt(3))
    zeta <- sapply(1:K,function(k)(-1)^(k+1)*(k+1)^(-nu/2))
    val <- 0 
    for(k in 1:K)
      val <- val + zeta[k]*Xi[k]*phi(t,k)
    val <- val + mu(t)
  }
  
  dense.grid <- seq(0.2,2.8,length.out=300)
  M<-length(dense.grid)
  N <- 500
  Y <- matrix(0, nrow=N, ncol=M)
  #generate (latent) response process
  for(i in 1:N)
  {
    # temp <- Z.fun(t=dense.grid, K=50, nu=2)
    score <- runif(50,-sqrt(3),sqrt(3))
    Y[i,] <- Z.fun(Xi=score, t=dense.grid, K=50, nu=2)
  }
  
  true.cov <- cov(Y)
  true.eigenl<- EigenC(true.cov,dense.grid)
  true.eigenfl<-true.eigenl[[1]] #eigenfunction
  true.eigenvl <- true.eigenl[[2]]
  true.eigenfl[,1] <- phi(t=dense.grid, k=1);

  
  
  x.points<-matrix(0,nrow=2,ncol=M*(M+1)/2)
  # cov.dense <- upperTriangle(true.cov, diag = T, byrow=T)
  count<-0
  for (i in 1:M){
    for (j in i:M){
      count<-count+1
      x.points[1,count]<-dense.grid[i]
      x.points[2,count]<-dense.grid[j]
      # cov.dense[count] <- c(true.cov)[count]
    }
  }
  
  h <- dense.grid[2] - dense.grid[1]
  
  
  beta <- 2
  n = 100  #generate 100 curves each time
  ini.Z <- numeric(n)
  Lt <- Lx <- list()
  # generate observed response
  for(i in 1:n)
  { # i = 1
    score <- runif(50,-sqrt(3),sqrt(3))
    ini.Z[i] <- Z.fun(Xi=score, t=0, K=50, nu=2)
    obs_time <- obs_res <-  NULL
    k <- 0; current_time <- 0; current.obs <- ini.Z[i]
    while(current_time < 3)
    { # j = 1
      gap_time <- rexp(1, rate=1*exp(current.obs*beta))
      current_time <- current_time+gap_time
      current.obs <- Z.fun(Xi=score, t=current_time, K=50, nu=2) + rnorm(1, mean=0, sd=0.1)
      obs_time <- c(obs_time, current_time)
      obs_res <- c(obs_res, current.obs)
      k <- k + 1
    }
    
    #observation times
    temp <- sort(obs_time)
    # temp <- 
    Lt[[i]] <- temp[temp < 3]
    #longitudinal observations for each curve
    Lx[[i]] <- obs_res[temp < 3]
  }
  
 # mean(sapply(Lt,length))
  
  Lx.r <- Lx[lapply(Lx,length)>0]; 
  index <- which(lapply(Lx,length)>0)
  Lt.r <- Lt[lapply(Lt,length)>0]
  n.r <- length(Lx.r)
  ini.value <- ini.Z[index]
  
  #id <- rep(1:n.r, times= sapply(Lt.r,length))
  tstart <- tend <- status <- id <- x_cov <- NULL
  for(i in 1:n.r)
  { # i = 1
    temp <- c(0, Lt.r[[i]])
    tstart <- c(tstart, temp)
    temp <- c(Lt.r[[i]], 3)
    tend <- c(tend, temp)
    temp <- c(rep(1, length(Lt.r[[i]])), 0)
    status <- c(status, temp)
    temp <- rep(i, length(Lt.r[[i]])+1)
    id <- c(id, temp)
    temp <- c(ini.value[1], Lx.r[[i]])
    x_cov <- c(x_cov, temp)
  }
  
  surv.dat <- list(start=tstart,stop =tend, event=status, id=id, x_cov=x_cov)
  ss <- coxph(Surv(start, stop, event) ~ cluster(id) + x_cov, data=surv.dat)
  #ss$coefficients
  
  t.max<- max(sapply(Lt.r,length)) 
  Obs.r<-matrix(0,n.r,t.max)
  T.r<-matrix(0,n.r,t.max)
  RWmean.r <- EWmean.r <- matrix(0, n.r, t.max)
  N.r<-sapply(Lt.r,length)
  
  
  for(i in 1:n.r)
  { #i=2
    Obs.r[i,1:N.r[i]]<-Lx.r[[i]]
    T.r[i,1:N.r[i]]<-Lt.r[[i]]
    index <- which(id==i)
    index <- index[-length(index)]
    EWmean.r[i,1:N.r[i]] <- inverse.weigh(beta=as.numeric(ss$coefficients), z=x_cov[index])
    RWmean.r[i,1:N.r[i]] <- inverse.weigh(beta=2,z=x_cov[index])
  }
  
  
  RW_data.r<-TranMtoV(Obs.r,T.r,RWmean.r,N.r)
  y.r<-RW_data.r[,2]  
  t.r<-RW_data.r[,3]  
  RW_mean.r <- RW_data.r[,4] #Raw weight for each observation
  
  EW_data.r<-TranMtoV(Obs.r,T.r,EWmean.r,N.r)
  y.r<-EW_data.r[,2]  #aggregated observations from all curves
  t.r<-EW_data.r[,3]  #thg corresponding aggregated observation times
  EW_mean.r <- EW_data.r[,4] #estimated weight for each observation
  
 
  
  meanf.err <- rep(NA, 3)   #unweighted, raw weighted, estimated weighted
  eigf.err <- rep(NA, 3)
  cov.err <- rep(NA, 3)
  eigv.err <- rep(NA, 3)
  
  
  Ew_est.mean  <- Rw_est.mean <- Uw_est.mean<-rep(NA, M)
  Ew_est.eigf1 <- Rw_est.eigf1 <- Uw_est.eigf1 <- rep(NA, M)
  
  #=============================================
  #unweighted mean estimate 
  #select tuning parameter for estimating the mean function 
  hmu.cv <- h.select(t.r,y.r,method="df")
  uw_mu1 <- sm.regression(t.r,y.r,h=hmu.cv, poly.index=1,eval.points=dense.grid)$estimate 
  uw_mu2 <- sm.regression(t.r,y.r,h=hmu.cv, poly.index=1, eval.points=t.r)$estimate 
  #=============================================
  #raw weighted mean estimate
  hmu.cv<-h.select(t.r,y.r,weights=RW_mean.r,method="df", nbins=0)
  #fit local linear regression 
  Rw_mu1<-sm.regression(t.r,y.r,h=hmu.cv,weights=RW_mean.r, nbins=0, poly.index=1,
                        eval.points=dense.grid)$estimate 
  Rw_mu2 <- sm.regression(t.r,y.r,h=hmu.cv,weights=RW_mean.r, nbins=0, poly.index=1,
                          eval.points=t.r)$estimate 
  
  
  #============================================
  #estimated weighted mean estimate
  hmu.cv<-h.select(t.r,y.r,weights=EW_mean.r,method="df", nbins=0)
  #fit local linear regression 
  Ew_mu1<-sm.regression(t.r,y.r,h=hmu.cv,weights=EW_mean.r, nbins=0, poly.index=1,
                        eval.points=dense.grid)$estimate 
  Ew_mu2 <- sm.regression(t.r,y.r,h=hmu.cv,weights=EW_mean.r, nbins=0, poly.index=1,
                          eval.points=t.r)$estimate 
  
  #=======================
  #calculate MISE of mean estimate
  meanf.err[1] <- trapz(x=dense.grid,y=(uw_mu1-mu(dense.grid))^2)
  meanf.err[2] <- trapz(x=dense.grid, y=(Rw_mu1-mu(dense.grid))^2)
  meanf.err[3] <- trapz(x=dense.grid, y=(Ew_mu1-mu(dense.grid))^2)
  
  
  Uw_est.mean <- uw_mu1
  Rw_est.mean <- Rw_mu1
  Ew_est.mean <- Ew_mu1
  
  
  #============================================================
  #covariance surface estimate smoothing spline
  
  conver<-HatG(uw_mu2,RW_data.r,N.r)  
  Cova<-conver[[2]]  #raw covariance estimate for each pair of time points
  CovT<-conver[[3]] #all pairs of time points 
  
  
  cov_dat <- data.frame(y=Cova, x1=CovT[1,], x2=CovT[2,])
  model_sm <- ssanova(y ~ x1*x2, data=cov_dat, alpha=1.4)
  est <- predict(model_sm, newdata=data.frame(x1=x.points[1,], x2=x.points[2,]))
  
  covmatrix  <- matrix(0,M,M)  
  
  count<-0
  for (i in 1:M){
    for (j in i:M){
      count<-count+1
      covmatrix[i,j]<-est[count]
      covmatrix[j,i]<-covmatrix[i,j] 
    }    
  }
  
  #MISE for Cov estimate
  cov.err[1] <- sum((true.cov-covmatrix)^2)*h^2 
  
  #eigenfunction estimation
  eigenl<-EigenC(covmatrix,dense.grid)
  eigenfl<-eigenl[[1]] #eigenfunction
 
  eigenvl<-eigenl[[2]] #eigenvalue
  a1 <- trapz(x=dense.grid,y=(eigenfl[,1]-true.eigenfl[,1])^2)
  a2 <- trapz(x=dense.grid,y=(eigenfl[,1]+true.eigenfl[,1])^2)
  eigf.err[1] <- min(a1,a2) #MISE for the estimated first eigenfunction
  eigv.err[1] <-  eigenvl[1] - true.eigenvl[1] #estimation error of the estimated first eigenvalue
  Uw_est.eigf1 <- (a1>a2)*(-eigenfl[,1]) + (a1 < a2)*eigenfl[,1]
  
  #======================================
  #raw weigh method 
  conver<-HatG(Rw_mu2,RW_data.r,N.r)  
  Cova<-conver[[2]]  #raw covariance estimate for each pair of time points
  CovT<-conver[[3]] #all pairs of time points 
  CovW<-conver[[4]] #weight for each raw estimate of covariance, given by the product of inverse weights
  
  #directly smoothing
  cov_dat <- data.frame(y=Cova, x1=CovT[1,], x2=CovT[2,])
  model_sm <- ssanova(y ~ x1*x2, data=cov_dat, weights=CovW,alpha=1.4)
  est <- predict(model_sm, newdata=data.frame(x1=x.points[1,], x2=x.points[2,]))
  
  
  
  covmatrix  <- matrix(0,M,M)  
  
  count<-0
  for (i in 1:M){
    for (j in i:M){
      count<-count+1
      covmatrix[i,j]<-est[count]
      covmatrix[j,i]<-covmatrix[i,j] 
    }    
  }
  
  #covariance surface estimate
  cov.err[2] <- sum((true.cov-covmatrix)^2)*h^2
  
  #first eigenfunction estimate
  eigenl<-EigenC(covmatrix,dense.grid)
  eigenfl<-eigenl[[1]] #eigenfunction
  eigenvl<-eigenl[[2]] #eigenvalue
  #eigenfunction estimate error
  b1 <- trapz(x=dense.grid,y=(eigenfl[,1]-true.eigenfl[,1])^2)
  b2 <- trapz(x=dense.grid,y=(eigenfl[,1]+true.eigenfl[,1])^2)
  eigf.err[2] <- min(b1,b2)
  eigv.err[2] <-  eigenvl[1] - true.eigenvl[1]
  Rw_est.eigf1 <- (b1>b2)*(-eigenfl[,1])+ (b1 < b2)*eigenfl[,1]
  
  
  #=====================================
  #estimated weighted method
  
  conver<-HatG(Ew_mu2,EW_data.r,N.r)  
  Cova<-conver[[2]]  #raw covariance estimate for each pair of time points
  CovT<-conver[[3]] #all pairs of time points 
  CovW<-conver[[4]] #weight for each raw estimate of covariance, given by the product of inverse weights
  
  #directly smoothing
  cov_dat <- data.frame(y=Cova, x1=CovT[1,], x2=CovT[2,])
  model_sm <- ssanova(y ~ x1*x2, data=cov_dat, weights=CovW,alpha=1.4)
  est <- predict(model_sm, newdata=data.frame(x1=x.points[1,], x2=x.points[2,]))
  
  
  
  covmatrix  <- matrix(0,M,M)  
  
  count<-0
  for (i in 1:M){
    for (j in i:M){
      count<-count+1
      covmatrix[i,j]<-est[count]
      covmatrix[j,i]<-covmatrix[i,j] 
    }    
  }
  
  #covariance surface estimate
  cov.err[3] <- sum((true.cov-covmatrix)^2)*h^2
  
  #first eigenfunction estimate
  eigenl<-EigenC(covmatrix,dense.grid)
  eigenfl<-eigenl[[1]] #eigenfunction
  eigenvl<-eigenl[[2]] #eigenvalue
  #eigenfunction estimate error
  b1 <- trapz(x=dense.grid,y=(eigenfl[,1]-true.eigenfl[,1])^2)
  b2 <- trapz(x=dense.grid,y=(eigenfl[,1]+true.eigenfl[,1])^2)
  eigf.err[3] <- min(b1,b2)
  eigv.err[3] <-  eigenvl[1] - true.eigenvl[1]
  Ew_est.eigf1 <- (b1>b2)*(-eigenfl[,1])+ (b1 < b2)*eigenfl[,1]
  
  res <- list(Ew_est.mean = Ew_est.mean, Rw_est.mean = Rw_est.mean,
              Uw_est.mean = Uw_est.mean, Ew_est.eigf1 = Ew_est.eigf1,
              Rw_est.eigf1 = Rw_est.eigf1, Uw_est.eigf1 = Uw_est.eigf1,
              meanf.err = meanf.err, cov.err = cov.err, 
              eigf.err = eigf.err, eigv.err = eigv.err)
  return(res)
}



### Setup parallel computing
library(dplyr)
library(pbapply)
library(doParallel)
cl <- makeCluster(min(detectCores() - 1, 100))
registerDoParallel(cl)
invisible(clusterEvalQ(cl, expr = {
  # load packages
  library(NHPoisson)
  library(fdapace)
  library(sm) 
  library(caTools)
  library(gss)
  library(rpart)
  library(survival)
  #source("gen_data.R")
  source("w_loclin.R")
  
}))


R <- 200
res <- pblapply(1:R, cl = cl, FUN = sim_function)
stopCluster(cl)
saveRDS(res, file = "S1_100.rds")

#repead the above procedure with n = 200

#meanf.err, cov.err,  eigf.err,  eigv.err to generate Table 1 of the main text

# take average of  Uw_est.mean, Rw_est.mean,  Ew_est.mean
#                  Uw_est.eigf1, Rw_est.eigf1, Ew_est.eigf1 to generate Figure 1 of the main text.
                   
 
