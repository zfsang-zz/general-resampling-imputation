#'Finding the indicator matrix for missing values
#'
#'This function develops a simple way to generate the indicator matrix for matrix with missing values
#'@param x a n x p matrix of n observation and p predictors, where NAs present as missing values
#'@return a n x p binary indicator matrix where 1 indicates value is missing, 0 indicates value is not missing. 
#'@author Zhifan Sang
#'@details 
#'This function runs NA search and replace on the input matrix, return a binary matrix of the same size.  If the (i,j) element 
#'in the input matrix is NA, then the (i,j) element of output matrix is 1; Otherwise, the (i,j) element of the output matrix is 0   
#'@seealso \code{is.na}
#'@export
detect.missing<-function(x){
  x[!is.na(x)] <- 0 
  x[is.na(x)] <- 1 
  return(x)
}



#'Bootstrap Imputation and Variable Selection
#'
#'This function deploys bootstrap imputation and variable selection on (X,y) data with missing values
#'@param X a n x p matrix of n observation and p predictors, where NAs present as missing values
#'@param Y a n x 1 vector represent the outcome parameter, esponse variable. Quantitative for family="gaussian", or family="poisson" (non-negative counts).
#' For family="binomial" should be either a factor with two levels, or a two-column matrix of counts or proportions (the second column is treated as the target class; 
#' for a factor, the last level in alphabetical order is the target class). 
#' For  "binomial" case, if y is presented as a vector, it will be coerced into a factor.
#'@param missing_col index of the missing column in X (only one column is missing in this case)
#'@param reg.type  the type of regression used in the second step, could be lasso, lasso.adp
#'@param MI.method the type of imputation method.  Can be "mice", or  "balasso".
#'@param BI.size number of bootstrap imputated datasets
#'@param pi threshold value for adaptive lasso variable selection
#'@param nsteps number of steps used in stability selection
#'@param family of regression, could be one of the following values: "gaussian","binomial","poisson")
#'@return a list of output values including Bootstrap imputated dataset(Bimp), Final variable selected indicator(S.fin), 
#'         final parameter estimates for the regression(beta.fin), final selection probability for lasso (S1.prob), 
#'         final selection for variable selection indicator for stability selection with randomized lasso(S2.prob)
#'@author Zhifan Sang
#'@details 
#'This function first performs the bootstrap imputation based on mice method(mice package) and balasso method(MIHD) method on data consists of  
#'outcome parameter and predictors with missing values to obtain bootstrap imputed dataset.  Then this function use lasso and stability selection with randomized lasso to conduct
#'variable selection on the bootstrap imputed datasets, and obtain selection indicators. Lastly, this function calculate selected predictor estimates.    
#'@examples
#'  output_lr_1=BI_SS(lr_1[,-1],lr_1$Y,family="gaussian",link=NULL,missing_col=1,MI.method="mice",BI.size=100,pi=NULL,nsteps=12)

           
#'@seealso \code{MIdurr}
#'@export
BI_SS<-function(X,Y,family,link=NULL,missing_col,MI.method,BI.size,pi,nsteps){
  #indicator<-detect.missing(X)
  if(is.null(link) ){
    family0<-switch (family,
                     gaussian = gaussian(),binomial=binomial(),poisson=poisson())
  }else{
    family0<-switch (family,
                     gaussian = gaussian(link),binomial=binomial(link),poisson=poisson(link))
  }
  
  
  N=dim(X)[1];
  
  K=dim(X)[2];
  M=BI.size; ###generate M bootstrap data sets;
  
  #nsteps=12;
  RLasso.alpha=0.5;
  pct=c(0.6,0.7,0.8,0.9,1.0); 
  q=3;
  nEst = q + 2 * 2 * length(pct);
  
  S.fin = array(0,dim=c(nEst,K)) # final indicator of variable selection
  beta.fin = var.fin = array(NA,dim=c(nEst,1+K))
  
  
  ### naive estiamte using only complete observations and full model ###
  if(missing_col==0){
    all=complete.cases(Y);
  }else{
    all=complete.cases(X);
  }
  
  Yt=as.matrix(Y[all==1]);
  Xt=as.matrix(X[all==1,]);
  
  if(N>K){
    tmp = summary(glm(Yt~Xt,family = family0))$coefficients;
    # print(length(beta.fin[1,]))
    # print(length(tmp[,1]))
    beta.fin[1,] = tmp[,1];
    S.fin[1,] = (tmp[,1] != 0)[-1]*1
  }
  # tmp = summary(lm(Yt~Xt))$coefficients;
  
  
  # naive.betahat[,1]=paste(round(tmp[,1],digits=3)," (",
  #                         round(tmp[,2],digits=3),")",sep="");
  
  ### naive estiamte using only complete observations and lasso ###
  all=complete.cases(X);
  Yt=as.matrix(Y[all==1]);
  Xt=as.matrix(X[all==1,]);
  
  cv=dim(Xt)[1];
  
  
  # 
  # model<- glmnet(Xt,Yt,standardize=FALSE)
  # cvres<-cv.glmnet(Xt,Yt,standardize=FALSE)
  model<- glmnet(Xt,Yt,family=family,standardize=FALSE)
  cvres<-cv.glmnet(Xt,Yt,family=family,standardize=FALSE)
  fits<-coef(model, s=cvres$lambda.min)
  fits<-as.numeric(fits)
  
  beta.fin[2,] =  fits
  S.fin[2,] = (fits[-1] != 0)*1
  
  
  
  
  
  
  ###get parameter estimate and SE using OLS after naive lasso
  # results=summary(lm(Yt~Xt[,which(fits[-1] != 0)]));
  
  results=summary(glm(Yt~Xt[,which(fits[-1] != 0)],family=family0));
  
  # naive.betahat[c(1,1+which(tmp != 0)),2]=paste(round(results$coefficients[,1],digits=3)," (",
  #                                               round(results$coefficients[,2],digits=3),")",sep="");
  # 
  
  ###naive estimate using only complete observations and adaptive lasso ###
  all = complete.cases(X)
  Yt=as.matrix(Y[all==1])
  Xt=as.matrix(X[all==1,])
  
  cv=dim(Xt)[1]
  Y.fit = adalasso(Xt,Yt)
  tmp = Y.fit$coefficients.adalasso
  if (is.null(Y.fit$intercept.adalasso)) intercept=0 else intercept=Y.fit$intercept.adalasso;
  beta.fin[3,] = c(intercept,tmp)
  S.fin[3,] = (tmp != 0)*1
  
  
  
  
  
  #####################################################
  ###proposed methods
  ### Step 1. bootstrap data and then impute missing X's ###
  Bimp=array(0,dim=c(M,N,K+1));
  for (i in 1: M) {
    
    ###bootstrap data
    data.bs=data.frame(X,Y)[sort(sample(1:N,size=N,replace=TRUE)),];
    
    ###two imputation methods
    if(MI.method=="mice")
    {#tmp = complete(mice(data.bs,m=1,method=imp.methods,printFlag=FALSE),1);
      tmp = complete(mice(data.bs,m=1,printFlag=FALSE),1);
      Bimp[i,,] = as.matrix(tmp);
    }
    
    if(MI.method=="blasso"){
      if(MI.method=="blasso"){
        Bimp[i,,] = as.matrix(data.bs);
        ###using blasso.vs for bayesian lasso imputation to impute X1
        ### where both Y and X are mean standardized, though X is not scaled by its SD;
        
        Bimp[i,,missing_col]<-MIdurr(data=Bimp[i,,],method="blasso",family=family,m=1)$impute
        
        
      }#end of if blasso;
    } # end of m loop
  }
  
  
  
  ###stability selection with randomized lasso, equation (7) in in meinshausen and Buhlmann (2010);
  S.fin0 = rep(0,K);
  S2.fin0 = array(0, dim=c(nsteps,K));
  S = array(0,dim=c(M,K)) # variable selection indicator;
  S2 = array(0,dim=c(M,nsteps,K)) # variable selection indicator for stability selection with randomized lasso;
  beta1 = beta2= array(0,dim=c(M,1+K)) # estimated parameters (include intercept)
  
  ###### Step 2-4 LASSO or adaptive lasso variable selection for each imputed dataset ###
  for (m in 1:M) {
    
    ### Step 2 LASSO or stability selection  variable selection for each imputed dataset ###
    Yt=Bimp[m,,K+1];
    Xt=Bimp[m,,1:K];
    
    ###lasso
    # model<- glmnet(Xt,Yt,standardize=FALSE)
    # cvres<-cv.glmnet(Xt,Yt,standardize=FALSE)
    # print(dim(Xt))
    # print(Xt)
    model<- glmnet(Xt,Yt,family=family,standardize=FALSE)
    cvres<-cv.glmnet(Xt,Yt,family=family,standardize=FALSE)
    fits<-coef(model, s=cvres$lambda.min)
    fits<-as.numeric(fits)
    
    # print(length(fits))
    # print(length(S[m,]))
    beta1[m,] = fits
    S[m,] = (fits[-1] != 0)*1
    
    
    
    ###stability selection with randomized lasso in meinshausen and Buhlmann (2010)
    random.factor=runif(K,RLasso.alpha,1);  #random number in (0.5,1)
    Xs <- t(t(Xt)*random.factor);
    
    
    # Y.fit = lars(Xs,Yt,type="lasso", max.steps=nsteps,normalize=FALSE)
    # tmp = coef.lars(Y.fit)[-1,]
    # S2[m,,] = (tmp != 0)*1
    
    model<- glmnet(Xs,Yt,family=family,standardize=FALSE)
    cvres<-cv.glmnet(Xs,Yt,family=family,standardize=FALSE)
    fits<-coef(model, s=cvres$lambda.min)
    fits<-as.numeric(fits)
    S2[m,,] = (fits[-1] != 0)*1
    
    
    
    
    ### Step 3. for computing the intersection of all M imputations ###
    S.fin0 = S.fin0 + S[m,]
    S2.fin0 = S2.fin0 + S2[m,,]
    
  } # end of m loop
  
  ###stability selection with randomized lasso, equation (7) in in meinshausen and Buhlmann (2010);
  S2.fin0=apply(S2.fin0,2,max);
  
  
  ### Step 4a. calculate the final estimate of beta using mean of beta's from the previous step ###
  for (i in 1:length(pct)) {
    ###lasso
    S.fin[q+i,] = (S.fin0 >= pct[i] * M)
    ### NOTE: did not exclude zero estimates ###
    
    beta.fin[q+i,] = apply(beta1,2,mean,na.rm=TRUE) * c(1,S.fin[q+i,])
    
    ###stability selection with randomized lasso
    S.fin[q+length(pct)+i,] = (S2.fin0 >= pct[i] * M)
    beta.fin[q+length(pct)+i,] = apply(beta1,2,mean,na.rm=TRUE) * c(1,S.fin[q+i,])
    
  }
  ###obtain the list of selected genes;
  ### Step 4b. alternative calculation of beta by rerun the regression analysis using only selected X's###
  S.fin[q+2*length(pct)+(1:(2*length(pct))),] = S.fin[q+(1:(2*length(pct))),]
  
  for (i in 1:length(pct)) {
    betat1=betat2 = vart1 = vart2 = array(0,dim=c(M,1+K)) 
    
    for (m in 1:M) {
      Yt=Bimp[m,,K+1];
      Xt=Bimp[m,,1:K];
      
      ###Bolasso
      
      if (sum(S.fin[q+2*length(pct)+i,]==1)==0) Y.fit = lm(Yt ~ 1) else Y.fit = lm(Yt ~ Xt[,S.fin[q+2*length(pct)+i,]==1])
      betat1[m,c(1,1+which(S.fin[q+2*length(pct)+i,]==1))] = Y.fit$coef
      #vart[m,c(1,1+which(S.fin[q+2*length(pct)+i,]==1))] = diag(vcov(Y.fit))
      
      ###stability selection with randomized lasso                        
      if (sum(S.fin[q+3*length(pct)+i,]==1)==0) Y.fit = lm(Yt ~ 1) else Y.fit = lm(Yt ~ Xt[,S.fin[q+3*length(pct)+i,]==1]);
      betat2[m,c(1,1+which(S.fin[q+3*length(pct)+i,]==1))] = Y.fit$coef
      
    }
    beta.fin[q+2*length(pct)+i,] = apply(betat1,2,mean,na.rm = TRUE)
    var.fin[q+2*length(pct)+i,] =  apply(betat1,2,var,na.rm = TRUE) #computing var based on bootstrap Imputation.	
    beta.fin[q+3*length(pct)+i,] = apply(betat2,2,mean,na.rm = TRUE)
    var.fin[q+3*length(pct)+i,] =  apply(betat2,2,var,na.rm = TRUE) #computing var based on bootstrap Imputation.	
  } # end of i loop
  
  

  S1.prob=S.fin0/M
  S2.prob=S2.fin0/M
  beta.se = sqrt(var.fin); ##only available for method b
  
  output=list(BIMP = Bimp, S.FIN = S.fin, BETA.FIN = beta.fin, S1.PROB = S1.prob,S2.PROB = S2.prob)
  return(output)
}
