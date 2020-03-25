#' @title cv.lme1en
#' @description Performing crossvalidation of lme1en (extension of inputs)
#' @details The function performs (possibly parallelized) K-fold cross validation for given alpha,lambda,rho
#' The function assume that the response y and design matrix X are centralized (no intercept returned from function) 
#' @param y response vector. Must be normalized before input.
#' @param X Design matrix belonging to fixed effects coefficients (beta). Must be normalized before input.
#' @param batch factor with batch effect names (vector for each observations)
#' @param rho proportion of variation explained by batch effect
#' @param lambda numeric, penalty levels for fixed effects betas (component 1) 
#' @param alpha numeric, penalty levels for fixed effects betas (balancing LASSO/RIDGE erros) 
#' @param beta numeric, initial values of the beta coefficients (using glmnet or marginal estimates if not provided)  
#' @param glmnetPenalty boolean, whether to use the original peanalty (FALSE) or the glmnet penalty (TRUE)
#' @param glmnetWarmup boolean, whether to use glmnet beta-estimates as warmup (if not marginals are used)
#' @param maxit maximum number of iterations (i.e. forloops) in the coordinate decent algorithm 
#' @param toler tolerance level of beta changes for each iterations (similar to 'thresh' in glmnet)
#' @param verbose boolean, show progress Default: FALSE
#' @param errorFun function, The user can provide own error function to perform cross validation with (Default is MSE if function not provided). The function must be able to take equally long two vectors and return a number
#' @param Kfold The train-test fold split: The training proportion becomes '(Kfold-1)/Kfold'
#' @param CViter Number of iterations of the cross-validation
#' @param parallel Boolean of whether parallelization in CV process should be applied (utilizing all CPUs)
#' @param seed The user can set seed for reproducibility if wanted
#' @return matrix with error values (test and train) for each iterations
#' @export
#' @examples
#' \dontrun{ 
#' dat = genData(seed=1)
#' CVerrors = cv.lme1en(y=dat$y,X=dat$X,batch=dat$batch, rho=0.3,lambda=0.1,alpha=0.5,verbose=TRUE,seed=1)
#' }
#errorFun = NULL;Kfold=10;CViter=10;parallel=TRUE
cv.lme1en = function(y,X,batch, rho=0, lambda=1, alpha=0.5, beta=NULL, glmnetPenalty=TRUE,glmnetWarmup=TRUE, maxit = 10000, toler = 1e-5 ,verbose=FALSE, errorFun = NULL, Kfold=10, CViter=10, parallel=TRUE,seed=NULL) {
  if(!is.null(seed)) set.seed(seed) #set seed if provided
  set.seed(seed) #MUST BE SAME FOR EACH M
  
  if(is.null(errorFun)) { #error function returning value between true and pred value 
    errorFun = function(x,y) mean( (x-y)^2 ) #use mean squared error by default
  } 

  n = length(y) #number of data
  nTrainF = (Kfold-1)/Kfold #fraction of train inds 1/(1-0.5)=2 fold
  nTrain = ceiling(n*nTrainF) #number of train data
  
  #TrainTest simulation:
  isTrain = matrix(FALSE,nrow=CViter,ncol=n) #matrix which indicates which is train etc
  for(m in 1:CViter) isTrain[m,sample(1:n,nTrain,replace=FALSE)] = TRUE #randomize which is train

  if(parallel) { #prepare for parallel run
    ncores <- parallel::detectCores() #number of physical cores (parallel)=number of chains
    if(verbose) print(paste0("Number of paralell chains will be ",ncores ))
    doParallel::registerDoParallel(cores=ncores) #do paralle
  }
  
  #traverse CV-iterations
  forList = foreach::foreach(m = 1:CViter, .packages = c("lme1en","glmnet")) %dopar% {
    trainind = which(isTrain[m,])
    testind = which(!isTrain[m,])
    
    bhat = lme1en(y=y[trainind],X=X[trainind,,drop=FALSE],batch=batch[trainind], rho=rho,lambda=lambda,alpha=alpha,glmnetPenalty=glmnetPenalty,glmnetWarmup=glmnetWarmup,toler=toler,maxit=maxit,verbose=FALSE)
    
    #PREDICT AND return error results
    yPred = X%*%bhat  #predict response for all data (both train/test). Assume centralized y and X data (intercept not used)
    
    errvals = c( test=errorFun( yPred[testind],y[testind]) , train=errorFun( yPred[trainind],y[trainind]) ) #get test and train errors
    return(errvals)
  }
  errorList = t(matrix(unlist(forList),nrow=2))
  colnames(errorList) = names(forList[[1]]) 
  return(errorList) #return matrix with error estimates
  
} #end function


