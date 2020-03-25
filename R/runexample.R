#' @title runexample
#' @description Fit elastic net models (without/with batch effects) for a generated example dataset 
#' @details The function do following steps: 1) simulate dataset 2) Split data into train/test 3) For two comparing models we run cross-valid for calibration and obtain performance using test data
#' 
#' @param TestTrainProp Proportion of data to be used for training (remaining used as test)
#' @param ntot total number of observations (including train and test data)
#' @param nsites number of sites
#' @param nbatches number of batches
#' @param propAgeAs #proportion of sites with assosiation to age
#' @param sd_batch standard dev of batch effects (normal distr effects)
#' @param sd_signal standard dev of signal (norm distr noise)
#' @param sd_bparam sd of age coefficient 
#' @param ageRange Age range of considered ages
#' @param seed The user can set seed for reproducibility if wanted
#' @param verbose boolean, show progress Default: TRUE
#' @export
#' @examples
#' \dontrun{ 
#'  runexample(ntot=500,nsites=1000,nbatches=20,verbose=TRUE,seed=1, sd_batch = 1, sd_signal = 1, sd_bparam = 0.1)
#' }

runexample = function(TestTrainProp=0.6, ntot = 500, nsites = 1000, nbatches = 20, propAgeAs=0.05, sd_batch = 1, sd_signal = 1, sd_bparam = 0.1 , ageRange = c(10,35), seed=NULL,verbose=TRUE) {
  if(!is.null(seed)) set.seed(seed) #set seed if provided
  
# library(lme1en); TestTrainProp = 0.6;seed=1;ntot = 100;nsites = 500;nbatches = 5;propAgeAs=0.05;sd_batch = 1;sd_signal = 0.1;sd_bparam = 0.1;ageRange = c(10,35);verbose=TRUE
  dat = genData(ntot , nsites , nbatches , propAgeAs, sd_batch , sd_signal , sd_bparam  , ageRange , seed ,verbose) #create dataset (default settings)

  ########################
  #Divide train/test data#
  ########################
  # Proportion of data as test
  nTrain = ceiling(TestTrainProp*ntot) #data to train model
  nTest = ntot - nTrain #test data
  trainInd = sample(1:ntot,nTrain,replace = FALSE) #get train data
  testInd = setdiff(1:ntot,trainInd) #get test data

  #################################
  #COMPARING TWO MODELS (I and II)#
  #################################
  #Error function to use; MAD=median absolute devaiation
  errorfun0 = function(x,y) median(abs(x-y))  #error function to use in CV (MAD)
  #errorfun0 is equivalent to type.measure="mae" in cv.glmnet
      
  #########################################
  #MODEL I (elasitic net, no batch effect)#
  #########################################
  if(verbose) print("Running crossvalidation on glmnet model")
  
  #Tuningparameter vectors:
  alphav = 0.5 #seq(0.1,1,by=0.1) #steps
  lamv = c(10^(-c(5:2)),0.05,0.1,1:2)  # lambda penalties
  
  glmnetCVList = list()
  for(alpha in alphav) { #for each alpha
    glmnetCVList[[paste0(alpha)]] = glmnet::cv.glmnet(dat$X[trainInd,],  dat$y[trainInd], family="gaussian", alpha=alpha, lambda=lamv, type.measure="mae",nfolds=10,standardize=FALSE,intercept = FALSE, parallel=FALSE)
    #plot(glmnetCVList[[paste0(alpha)]])
  }
  
  #OBTAIN BEST MODEL (CALIBRATED)
  cvm = sapply(glmnetCVList,function(x) x$cvm) #get mean errors
  minInd = which(cvm==min(cvm), arr.ind = TRUE) #get best models
  alphahat = alphav[minInd[2]] #get best alpha
  lamhat = rev(lamv)[minInd[1]] #get best lambda
  coefhat = coef(glmnetCVList[[minInd[2]]]) #get coefficients
  #a0 = coefhat[1]
  bhat1 = coefhat[-1]
    
  showStat = function(selSites0,met="") {
    print(paste0("Number of sites found with ",met,": ",length(selSites0)))
    print(paste0("Number of true sites observed/number of true=(",sum(selSites0%in%names(dat$coefEffect)),"/",length(dat$coefEffect),")"))
    print(paste0("Number of false positive sites observed=",sum(!selSites0%in%names(dat$coefEffect))))
  }  
  showStat(rownames(coefhat)[-1][bhat1!=0],"glmnet") #get selected sites and show stats
  
  #PREDICT AND SHOW RESULTS
  age_glmnet = dat$X%*%bhat1 #predict age for train data
  showResult( dat$y ,age_glmnet,trainInd,testInd,errorfun0) 
  #plot(bi,bhat)
  
  ##################################################
  #MODEL II (elastic net, random effect for batches#
  ##################################################

  if(verbose) print("Running crossvalidation on lme1en model")
  
  #CV settings:
  CViter = 20 #number of CV samples
  Kfold=10 #Selection of K-fold
    
  #Hyperparameters to traverse:
  alphv = 0.5# alpha is fixed
  rhov = c(0.1,0.25,0.5) #.3 #high correlation within batches
  lamv = 10^(-c(4,2,1)) 
  CVpar = expand.grid(list(rhov,lamv,alphv))
  nCVpar = nrow(CVpar) #number of hyper params to traverse
  errorMat = matrix(nrow=nCVpar,ncol=3)
  colnames(errorMat) = c("Mean","2.5%","97.5%")
  rownames(errorMat) = paste0("rho=",CVpar[,1],",lambda=",CVpar[,2],",alpha=",CVpar[,3])
  for(i in 1:nrow(CVpar)) {
    rho0 = CVpar[i,1]
    lam0 = CVpar[i,2]
    alph0 = CVpar[i,3]
    cvlist = cv.lme1en(y=dat$y[trainInd],X=dat$X[trainInd,],batch=dat$batch[trainInd], rho=rho0,lambda=lam0,alpha=alph0,parallel=TRUE,seed=seed,errorFun = errorfun0, CViter = CViter, Kfold=Kfold)

    #calculating mean and 95% CI of population mean
    errorMat[i,1] <- exp <- mean(cvlist[,1]) #get average
    sd = sd(cvlist[,1]) #get empirical sd
    dev = qt(0.975,CViter-1)*sd/sqrt(CViter) #get deviation
    errorMat[i,-1] = c(exp-dev,exp+dev)  #get confidence interval based on t-distr

    if(verbose) print(paste0("Iterating set of hyperparameters: ",round(i/nrow(CVpar)*100),"% compelete"))
  }
  bestCVpar = CVpar[which.min(errorMat[,1]),] #get optimal set
  rho0 = as.numeric(bestCVpar[1])
  lam0 = as.numeric(bestCVpar[2])
  alph0 = as.numeric(bestCVpar[3])
#  y=dat$y[trainInd];X=dat$X[trainInd,];batch=dat$batch[trainInd];rho=rho0;lambda=lam0;alpha=alph0
  bhat2 = lme1en(y=dat$y[trainInd],X=dat$X[trainInd,],batch=dat$batch[trainInd], rho=rho0,lambda=lam0,alpha=alph0)
  
  showStat(rownames(coefhat)[-1][bhat2!=0],"lme1en") #get selected sites and show stats
  
  #PREDICT AND SHOW RESULTS
  age_lme1en = dat$X%*%bhat2 #predict age for train data
  showResult( dat$y ,age_lme1en,trainInd,testInd,errorfun0) 
} #end runexample function

