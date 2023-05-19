#### 04-13-08
#### This file contains  R functions used for FPCA and conditional distance analysis
#### These functions are used in eBay_codes.txt
###########PART I: local linear/quadratic estimation of mu(t)-mean,G(s,t)-covariance,sigma^2--error variance, ect.
###
LocLin.mean<-function(y,t,teval,hmu){
  ##local linear estimation for mu(t)(mean): scatterplot smoother on pooled data
  ##name:LocLin.mean
  ##para:(y,t)--in the TranMtoV form, data<-TranMtoV(Obs,T,N),  y<-data[,2] obs; t<-data[,3] measurement time points;
  ##     teval--time points for evaluation; hmu--bandwidth
  ##result:muhat(teval)
  ssm<-sm.regression(t,y,h=hmu,poly.index=1,eval.points=teval,eval.grid=FALSE)       #sm.regression result:evaluate at teval
  result<-ssm$estimate
  return(result)
}

###
LocLin.G<-function(Cova,CovT,CovTeval,hcov){
  ##local linear estimation for G(s,t)-covariance; use hat(G(tij,til)) and ignore diagonal observations;
  ## scatterplot smoother on pooled data
  ##name:LocLin.G
  ##para: Cova, CovT from formatting function HatGOff; hcov--bandwidth (2*1); CovTeval--evaluation points
  ##result:two dim local linear estimation
  
  ssmcov<-sm.regression(t(CovT),Cova,h=hcov,poly.index=1,eval.points=t(CovTeval),eval.grid=FALSE)
  return(ssmcov$estimate) 
}

###
LocLin<-function(t,y,dime,eval,h,ker="Epan",epsi=1e-8){
  ##local linear regression for 1-dim and 2-dim predictors; Epanechnikov kernel is used
  ##name:LocLin
  ##para:t--predictors (times), y--response,dime--dimension, eval--evaluation grids, h--bandwidth
  ##result:local linear estimate at "eval"
  if(dime==1){
    n<-length(t)
    if(ker=="Gau"){
      w<-sqrt(apply(matrix((t-eval)/h,nrow=1),2,OneDGauKer))  
    }else{ 
      w<-sqrt(apply(matrix((t-eval)/h,nrow=1),2,OneDEpKernel))
    }
    xt<-(eval-t)*w
    xt<-cbind(w,xt)
    yt<-y*w
    beta<-solve(t(xt)%*%xt+epsi*diag(rep(1,2)))%*%t(xt)%*%yt
    return(beta[1])
  }
  if(dime==2){
    n<-nrow(t)
    s1<-eval[1]
    s2<-eval[2]
    kert<-(t-cbind(rep(s1,n),rep(s2,n)))/cbind(rep(h[1],n),rep(h[2],n))
    if(ker=="Gau"){
      w<-sqrt(apply(kert,1,TwoDGauKer))
    }else{
      w<-sqrt(apply(kert,1,TwoDEpKernel))
    }
    xt<-cbind(w,(s1-t[,1])*w,(s2-t[,2])*w)
    yt<-y*w
    beta<-solve(t(xt)%*%xt+epsi*diag(rep(1,3)))%*%t(xt)%*%yt
    return(beta[1]) 
  }
}


###
LocLinErr.diag<-function(y,t,fitmu,teval,hsig){
  ##estimate diagonal with measurement error: G(t,t)+\sigma^2; use local linear estimation
  ##name: locLinErr.diag 
  ##para: (y,t)--in the TranMtoV form,data<-TranMtoV(Obs,T,N),y<-data[,2] obs;t<-data[,3] measurement points;
  ##      fitmu--estimated mu(t) at t; teval--evaluation points; hsig--bandwidth
  ##result: G(t,t)+\sigma^2 at teval 
  ssmd<-sm.regression(t,(y-fitmu)^2,h=hsig,poly.index=1,eval.points=teval,eval.grid=FALSE)
  return(ssmd$estimate) 
}


###
LocLinEst.new<-function(Obs,T,N,indext,hmu,hcov,hsig,epsi=1e-8){
  ##local linear estimation of mu(t), G(s,t) for given bandwidth, and evaluated at fine grids
  ##name:LocLinEst
  ##para:Obs--observed values, T--measurement time points, N--number of measurements;
  ##     indext--1 dim evaluation grid;
  ##     hmu--bandwidth for mu,hcov--bandwidth for G(s,t)(2*1),hd--bandwidth for diagonal,
  ##     hsig--bandwidth for G(t,t)+sigma^2;
  ##result:ssmt$estimate--estimation of mu; covmatrix--estimation of G(s,t)I(s!=t);
  ##       diagcovt--estimation of G(t,t) with quadratic kernel by CovDiag; 
  ##       ssmdt$estimate--estimation of G(t,t)+\sigma^2 by LocLinErr.diag             
  ##       sighat2--estimation of sigma^2*(L2-L1); 
  
  data<-TranMtoV(Obs,T,N)  #transform data to vectors  
  trdata<-TransData(Obs,T,N)
  y<-data[,2]              #obs
  t<-data[,3]              #measurements 
  #(0) two dim evaluation grid
  t.points<-indext    #2-dim
  s.points<-indext
  M<-length(t.points)
  x.points<-matrix(0,2,M*(M+1)/2)
  count<-0
  for (i in 1:M){
    for (j in i:M){
      count<-count+1
      x.points[1,count]<-t.points[i]
      x.points[2,count]<-s.points[j]
    }
  }
  
  #(i) estimate mean by local linear smoother:
  ssm<-sm.regression(t,y,h=hmu,poly.index=1,eval.points=t,eval.grid=FALSE)       #evaluate at t
  ssmt<-sm.regression(t,y,h=hmu,poly.index=1,eval.points=indext,eval.grid=FALSE) #evaluate at equally spaced grid
  ssmtr<-sm.regression(t,y,h=hmu,poly.index=1,eval.points=as.vector(trdata[[3]]),eval.grid=FALSE) 
  
  #(iii)estimate covariance: G(s,t)I(s!=t)
  fitmu<-ssm$estimate   
  conver<-HatGOff(fitmu,data,N)
  Cova<-conver[[2]]
  CovT<-conver[[3]]  
  ssmcovt<-sm.regression(t(CovT),Cova,h=hcov,poly.index=1,eval.points=t(x.points),eval.grid=FALSE) 
  
  #(iv)covert to matrix
  covmatrix<-matrix(0,M,M)  #G(s,t)I(s!=t) (linear diagonal) 
  count<-0
  for (i in 1:M){
    for (j in i:M){
      count<-count+1
      covmatrix[i,j]<-ssmcovt$estimate[count]
      covmatrix[j,i]<-covmatrix[i,j]           
    }    
  }
  
  #(v) estimate G(t,t)+\sigma^2 and \sigma^2   
  ssmdt<-sm.regression(t,(y-fitmu)^2,h=hsig,poly.index=1,eval.points=indext,eval.grid=FALSE) 
  diagsig<-ssmdt$estimate
  
  sig2hat<-sum(diagsig[round(M/4):round(3*M/4)]-diag(covmatrix)[round(M/4):round(3*M/4)])*(indext[2]-indext[1])*2
  
  #(vi) return 
  result<-list(ssmt$estimate,covmatrix,diagsig,sig2hat)
  return(result)
}

###
EstCV<-function(Obs,T,N,indext,hmu,hcov,hd,hsig,epsi=1e-8){
  ##cross-validation and estimation all together
  ##name:EstCV
  ##para:Obs,T,N;indext--evaluation grid;hmu(vec),hcov(L*2),hd(vec),hsig(vec)--bandwidth for selection;
  ##result:LocLinEst result; cross-validation scores
  
  data<-TranMtoV(Obs,T,N)
  y<-data[,2]
  t<-data[,3]
  ##(i)leave one out cross validation to choose bandwidth
  #mu(t)
  cv.mean<-CV.mean(Obs,T,N,hmu)
  hmu.cv<-hmu[order(cv.mean)[1]]
  fitmu<-LocLin.mean(y,t,t,hmu.cv)
  #G(s,t)
  cv.G<-CV.G(Obs,T,N,fitmu,hcov)
  hcov.cv<-hcov[order(cv.G)[1],]
  #G(t,t)
  cv.diag<-CV.diag(Obs,T,N,hmu.cv,hd,epsi=epsi)
  hd.cv<-hd[order(cv.diag)[1]]
  #G(t,t)+\sigma^2
  cv.diagerr<-CV.diagerr(Obs,T,N,hmu.cv,hsig)
  hsig.cv<-hsig[order(cv.diagerr)[1]]
  ##(ii) estimation by using the selected bandwidth
  lin.result<-LocLinEst(Obs,T,N,indext,hmu.cv,hcov.cv,hd.cv,hsig.cv,epsi)
  result<-list(lin.result,cv.mean,cv.G,cv.diag,cv.diagerr)
  return(result)
}


####################### PART II: format data
###
TransData<-function(Obs,L,N){ 
  ##transform data to 2 column matrix (Y_{ij},Y_{il},j!=l)and(L_{ij},L_{il},j<l);
  ##then apply diagonal transformation on Ls;
  ##result used by CovDiag; called by CV.diag
  ##Name:TransData
  ##para: Obs--observed values,L--locations of measurement, N--# of measurements of each curve,
  ##return: TObs--observation(2 columns), TL--tranformed times;
  ##        TLO--original times, TId--subject id for time each pair
  tran<-sqrt(2)/2*matrix(c(1,1,-1,1),2,2,byrow=T)
  num<-sum(N*(N-1)/2) 
  TObs<-matrix(0,num,2)
  TL<-matrix(0,num,2)
  TLO<-matrix(0,num,2)
  TId<-numeric(num)
  count<-0
  for (i in 1:nrow(Obs)){
    for(j in 1:N[i]){
      for(l in j:N[i]){
        if(l !=j){
          count<-count+1   
          TLO[count,]<-c(L[i,j],L[i,l])     
          TL[count,]<-tran%*%matrix(TLO[count,],2,1)                      
          TObs[count,]<-c(Obs[i,j],Obs[i,l])
          TId[count]<-i     
        }     
      }
    }  
  }
  return(list(TObs,TL, TLO,TId))
}

###
TranMtoV<-function(Obs,L,W, N){
  ##transform data into a vector
  ##name:TranMtoV 
  ##para: Obs--observed values,L--locations of measurements, N--# of measurements of each curve, W-- weight of each observation
  ##return: data--3col matrix;1--ID,2--data,3--locations of meaturements
  n<-sum(N)
  result<-matrix(0,n,4)
  count<-0
  for(i in 1:nrow(Obs)){
    for(j in 1:N[i]){
      count<-count+1
      result[count,1]<-i
      result[count,2]<-Obs[i,j]
      result[count,3]<-L[i,j]
      result[count,4] <- W[i,j]
    }
  }
  return(result)
}


###
HatGOff<-function(fitmu,data,N){
  ##standardize data by subtracting mean and then return \hat{G(s,t)} at off diagonal pair (s,t)
  ##result used by LocLin.G; called by CV.G
  ##name:HatGOff
  ##para:fitmu--fitted mean curve at observed points; 
  ##data<-TranMtoV(Obs,T,N) format;N--number of observations for each subject
  ##return:Cova--a vector of \hat{G(tij,til)}; and CovT: time pairs of observation (tij,til),j<l; 
  ##        and  CovW: weight of observation (t_{ij}, t_{il}), j<l;
  ##       index--subject ID for the pair (tij,til), j<l 
  data.mean<-cbind(data,fitmu)
  Cova<-numeric(sum(N*(N-1)/2))     
  CovT<-matrix(0,nrow=2,ncol=sum(N*(N-1)/2))
  CovW<-numeric(sum(N*(N-1)/2))   
  index<-numeric(sum(N*(N-1)/2)) 
  count<-0  
  for (i in 1:length(N)){
    if(N[i]>1){
      med<-data.mean[data[,1]==i,]
      yi<-med[,2]
      ti<-med[,3] 
      wi<-med[,4]
      mi<-med[,5]
      for (j in 1: (N[i]-1)){
        for (l in (j+1): N[i]){
          count<-count+1
          Cova[count]<-(yi[j]-mi[j])*(yi[l]-mi[l])
          CovT[1,count]<-ti[j] 
          CovT[2,count]<-ti[l] 
          CovW[count]<-wi[j]*wi[l]
          index[count]<-i
        }
      }
    }
  }  
  return(list(index,Cova,CovT,CovW))
}




HatG <-function(fitmu,data,N){
  ##standardize data by subtracting mean and then return \hat{G(s,t)} at all pairs (s,t)
  ##result used by LocLin.G; called by CV.G
  ##name:HatGOff
  ##para:fitmu--fitted mean curve at observed points; 
  ##data<-TranMtoV(Obs,T,N) format;N--number of observations for each subject
  ##return:Cova--a vector of \hat{G(tij,til)}; and CovT: time pairs of observation (tij,til),j<l; 
  ##        and  CovW: weight of observation (t_{ij}, t_{il}), j<=l;
  ##        index--subject ID for the pair (tij,til), j<=l 
  
  #data=data.r; N<-N.r
  data.mean<-cbind(data,fitmu)
  Cova<-numeric(sum(N*(N+1)/2))     
  CovT<-matrix(0,nrow=2,ncol=sum(N*(N+1)/2))
  CovW<-numeric(sum(N*(N+1)/2))   
  index<-numeric(sum(N*(N+1)/2)) 
  count<-0  
  for (i in 1:length(N)){ #i=1
      if(N[i] > 1){
        med<-data.mean[data[,1]==i,]
        yi<-med[,2]
        ti<-med[,3] 
        wi<-med[,4]
        mi<-med[,5]
        for (j in 1:N[i]){
          for (l in j: N[i]){
            count<-count+1
            Cova[count]<-(yi[j]-mi[j])*(yi[l]-mi[l])
            CovT[1,count]<-ti[j] 
            CovT[2,count]<-ti[l] 
            if(l!=j){
              CovW[count]<-wi[j]*wi[l]
            } else {
              CovW[count]<-wi[j]
            }
            
            index[count]<-i
          }
        }
        
      } else {
        med<-data.mean[data[,1]==i,]
        yi<-med[2]
        ti<-med[3] 
        wi<-med[4]
        mi<-med[5]
        count<-count+1
        Cova[count]<-(yi-mi)^2
        CovT[1,count]<-ti
        CovT[2,count]<-ti
        CovW[count]<-wi
        index[count]<-i
      }
      
    }
  
  return(list(index,Cova,CovT,CovW))
}

###################PART III: other functions
###
StanResidu<-function(y,fitmu){
  ##caculate residuals for each observation after fitting the mean
  ##then devided by pooled sample standard deviation of the residuals
  ##name:StanResidu
  ##para:y--observations in forms of TranMtoV;
  ##     fitmu--fitted mean curve at oberved points
  ##result: (y-fitmu)/sd(y-fitmu)
  res<-y-fitmu
  ressd<-sd(res)
  result<-res/ressd
  return(result)          
}

###
RemOut<-function(Obs,T,N,fitmu,thre=3){
  ##remove subjects with outliers which are defined as 
  ##absolute values of standardized residuals (from StanResidu) bigger than thre
  ##name:RemOut
  ##para:Obs,T,N;fitmu;thre=3--threshold
  data<-TranMtoV(Obs,T,N)
  index<-data[,1]
  y<-data[,2]
  stanres<-StanResidu(y,fitmu)
  index<-index[abs(stanres)>thre] 
  index.sub<-as.numeric(attr(factor(index),"levels"))
  Obs.new<-Obs[-index.sub,]
  T.new<-T[-index.sub,]
  N.new<-N[-index.sub]
  index.new<-(1:nrow(Obs))[-index.sub]
  return(list(Obs.new,T.new,N.new,index.new))
}



###
IndiKernel<-function(s,t,h){
  ##indicator of I(|s-t|>h)
  ##Name:IndiKernel 
  ##para: s,t,h
  if(abs(s-t)>h){
    return(1)
  }else{
    return(0)
  }
}

###
ConvoKernel<-function(s,t,h,hnor=0.01){
  ##convolution of Indikernel with Gaussian kernel
  ##Name:ConvoKernel 
  ##para:s,t;h--bandwidth for Indikernel, hnor--bandwidth for Gaussian kernel
  result<-pnorm((s-t-h)/hnor)+1-pnorm((s-t+h)/hnor)
  return(result)
}

###
GauDiag<-function(covmatrix,diag, indext,hcov,hnor){
  ##Gausssian kernel convolution for covariance surface G(s,t)
  ##name: GauDiag
  ##para:covmatrix--G(s,t)I(s!=t);diag--G(t,t);indext--evaluation grids;
  ##     hcov--bandwidth for getting G(s,t);hnor--bandwidth for Gaussian convolution
  ##result:covmatrixd--G(s,t)
  #two dim evaluation grid
  t.points<-indext    
  s.points<-indext
  M<-length(t.points)
  x.points<-matrix(0,2,M*(M+1)/2)
  count<-0
  for (i in 1:M){
    for (j in i:M){
      count<-count+1
      x.points[1,count]<-t.points[i]
      x.points[2,count]<-s.points[j]
    }
  }
  #use the ConvoKernel
  covmatrixd<-matrix(0,M,M)  #G(s,t)+diag by ConvoKernel
  count<-0
  for (i in 1:M){
    for (j in i:M){
      count<-count+1    
      ker<-ConvoKernel(x.points[1,count],x.points[2,count],hcov,hnor) 
      covmatrixd[i,j]<-covmatrix[i,j]*ker+diag[round((i+j)/2)]*(1-ker) 
      covmatrixd[j,i]<-covmatrixd[i,j]
    }    
  }
  return(covmatrixd)
}


###
TtoGrid<-function(T,N,indext){
  ##transform observed measurement points to nearest grid points and also return index of the gird points
  ##called by FPcScore, CDist.mtr, ProjFirPc
  ##name:TtoGrid
  ##para:T,N--generated by GenData; indext--grid evaluation points
  ##result:Ts--nearest grid points of T; Tindex--corresponding index
  Ts<-T
  Tindex<-T
  for (i in 1:nrow(T)){
    for (j in 1:N[i]){
      Tindex[i,j]<-order(abs(T[i,j]-indext))[1]
      Ts[i,j]<-indext[Tindex[i,j]]
    }
  }
  return(list(Ts,Tindex))
}



############## PART III:conditional pc scores and conditional distances

###
EigenC<-function(covmatrixd,indext){
  ##get eigenfunctions and eigenvalues given the covariance surface G(s,t)
  ##name:EigenC
  ##para:covmatrixd--estimation of G(s,t);indext--evaluation grid
  ##result:eigenfunctions and eigenvalues
  eigend<-svd(covmatrixd)            ##svd decomposition on the grid of covariance:X=U%*%D%*%t(U)
  eigenvd<-eigend$d                  
  eigenvd<-eigenvd*(indext[2]-indext[1])  ##eigenvalues
  eigenfd<-eigend$u                  
  eigenfd<-eigenfd/sqrt(indext[2]-indext[1]) ##eigenfd[,i] is the ith eigenfunction
  return(list(eigenfd,eigenvd))
}

###
FPcScore<-function(Obs,T,N,indext,muhat,covmatrix,sig2hat,K){
  ##estimated conditional principal component scores: \hat{E(\xi_k|Y)}, k=1,...,K
  ##Name:FPcScore
  ##para:Obs,T,N;indext--evaluation grid;
  ##     muhat,covmatrix,sig2hat--estimation result from LocLinEst/EstCV;
  ##     K--number of eigenfunctions 
  ##result:first K conditional PC scores 
  n<-nrow(Obs)
  result<-matrix(0,n,K)  ##FPC scores for n subjects 
  
  #(i)eigen functions and eigen values based on covmatrix
  eigen<-svd(covmatrix)           
  eigenv<-eigen$d                  
  eigenv<-eigenv*(indext[2]-indext[1])  ##eigenvalues
  eigenv<-eigenv[1:K]
  eigenf<-eigen$u                  
  eigenf<-eigenf/sqrt(indext[2]-indext[1])  ##eigenf[,i] is the ith eigenfunction
  eigenf<-eigenf[,1:K]
  
  #(ii) convert observed time points to nearest grid points
  Tgrid<-TtoGrid(T,N,indext)
  Tindex<-Tgrid[[2]]  
  
  #(iii) first K conditional PC scores for each subject
  for (i in 1:nrow(Obs)){
    Ti<-Tindex[i,1:N[i]]
    X<-Obs[i,1:N[i]]
    mux<-muhat[Ti]
    Phix<-matrix(eigenf[Ti,],ncol=K)       
    result[i,]<-solve(t(Phix)%*%Phix+sig2hat*diag(1/eigenv))%*%t(Phix)%*%(X-mux)    
  }
  return(result)
}


###
CDist.mtr<-function(Obs,T,N,indext,muhat,covmatrix,sig2hat,K){
  ##estimated pairwise conditional L2 distance for n subjects; (d(i,i)=0)
  ##Name:CDist.mtr
  ##para:Obs,T,N;indext--evaluation grid;
  ##     muhat,covmatrix,sig2hat--estimation result from LocLinEst; K--number of eigenfunctions 
  ##result: n*n symmetric matrix of estimated pairwise conditional L2 distance 
  n<-nrow(Obs)
  result<-matrix(0,n,n)
  
  #(i)eigen functions and eigen values
  eigen<-svd(covmatrix)           
  eigenv<-eigen$d                  
  eigenv<-eigenv*(indext[2]-indext[1]) ##eigenvalues
  eigenv<-eigenv[1:K]
  eigenf<-eigen$u                  
  eigenf<-eigenf/sqrt(indext[2]-indext[1]) ##eigenf[,i] is the ith eigenfunction
  eigenf<-eigenf[,1:K]
  
  #(ii) convert observed time points to nearest grid points
  Tgrid<-TtoGrid(T,N,indext)
  Tindex<-Tgrid[[2]]  
  
  #(iii) pairwise conditional L2 distance 
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      Ti<-Tindex[i,1:N[i]]
      X<-Obs[i,1:N[i]]
      mux<-muhat[Ti]
      Phix<-matrix(eigenf[Ti,],ncol=K)    
      Tj<-Tindex[j,1:N[j]]
      Y<-Obs[j,1:N[j]]
      muy<-muhat[Tj]
      Phiy<-matrix(eigenf[Tj,],ncol=K)
      
      result[i,j]<-CDist(eigenv,sig2hat,Phix,Phiy,mux,muy,X,Y)[[1]]
      result[j,i]<-result[i,j]
      print(c(i,j)) 
    }
  }
  return(result)
}




###
CDist<-function(lab,sig2,Phix,Phiy,mux,muy,X,Y){          
  ##conditional L2 distance for two subjects given mu,phi evaluated at observed points;
  ##distance based on conditional second moments and conditional PC scores of subjects X and Y
  ##Name: CDist
  ##para:lab--eigenvalues;sig2--error variance;Phix,Phiy--eigenfunctions;
  ##mux,muy--means;X,Y--observations
  ##return:cm2.distance--conditional L2 distance;cmeanx,cmeany--conditional pc scores
  
  #(i)conditional mean and covariance of pc scores given the observed values
  D<-diag(lab)
  Nx<-length(X)
  Ny<-length(Y)
  cmeanx<-solve(t(Phix)%*%Phix+sig2*diag(1/lab))%*%t(Phix)%*%(X-mux)          ##conditional mean 
  cmeany<-solve(t(Phiy)%*%Phiy+sig2*diag(1/lab))%*%t(Phiy)%*%(Y-muy)
  ccovx<-D-D%*%t(Phix)%*%solve(Phix%*%D%*%t(Phix)+sig2*diag(1,Nx))%*%Phix%*%D ##conditional covariance 
  ccovy<-D-D%*%t(Phiy)%*%solve(Phiy%*%D%*%t(Phiy)+sig2*diag(1,Ny))%*%Phiy%*%D  
  cm2x<-diag(ccovx)+cmeanx^2                                             ##conditional second moments 
  cm2y<-diag(ccovy)+cmeany^2
  
  #(ii)conditional L2 distance based on conditional second moments:\sum_i (cm2x[i]+cm2y[i])
  cm2.distance<-sum(cm2x+cm2y-2*cmeanx*cmeany)
  return(list(cm2.distance,cmeanx,cmeany))
}

###
PreCurve<-function(pcscr,muhat,phihat){
  ##prediction of curves (subjects) giving conditional PC scores
  ##Name:PreCurve
  ##para: pcscores--n*K; muhat(L*1)--means,phihat(K*L)--eigenfunctions; 
  ##result: prediction of curves given conditional pcscores (true or estimated)
  result<-matrix(1,nrow=nrow(pcscr),ncol=1)%*%matrix(muhat,nrow=1)+pcscr%*%phihat
  return(result) 
}
