#' @title lme1en
#' @description Fitting Elastic Net with random intercept 
#' @details The model is fitted using the coordinate decent algorithm (exact solutions)
#' The extended model includes the rho parameter (proportion of total variation) as argument. 
#' Calls a algorithm (iterate) in C which uses exact solutions of marginal beta's based on formula x = sgn(c/a)(|c/a| - b/a)_+, where ax + b*sgn(x) = c
#' The function assume that the response y and design matrix X are centralized (no intercept returned from function) 
#' @param y response vector. Should be standardized before input (use scale). 
#' @param X Design matrix belonging to fixed effects coefficients (beta). Should be standardized before input (use scale).
#' @param batch factor with batch effect names (vector for each observations)
#' @param rho proportion of variation explained by batch effect
#' @param lambda numeric, penalty levels for fixed effects betas
#' @param alpha numeric, penalty levels for fixed effects betas (balancing LASSO/RIDGE erros) 
#' @param beta numeric, initial values of the beta coefficients (using glmnet or marginal estimates if not provided)  
#' @param glmnetPenalty boolean, whether to use the original peanalty (FALSE) or the glmnet penalty (TRUE)
#' @param glmnetWarmup boolean, whether to use glmnet beta-estimates as warmup (if not marginals are used)
#' @param maxit maximum number of iterations (i.e. forloops) in the coordinate decent algorithm 
#' @param toler tolerance level of beta changes for each iterations (similar to 'thresh' in glmnet)
#' @param verbose boolean, show progress Default: FALSE
#' @return fitted beta values
#' @export
#' @examples
#' \dontrun{ 
#' dat = genData(seed=1)
#' bhat = lme1en(y=dat$y,X=dat$X,batch=dat$batch, rho=0.3,lambda=0.1,alpha=0.5,verbose=TRUE)
#' }
lme1en = function(y,X,batch, lambda=1, alpha=0.5,  rho=0, beta=NULL, glmnetPenalty=TRUE,glmnetWarmup=TRUE, maxit = 10000, toler = 1e-3 ,verbose=FALSE) {
  #NB; BE SURE TO CENTRALIZE y and X data (no intercept returned)
  
  #CHECK DATA INPUT:
  if(any(rho<0) || any(rho>1)) stop("rho was not within [0,1]")
  if(any(alpha<0) || any(alpha>1)) stop("alpha was not within [0,1]")
  if(any(lambda<0) ) stop("lambda cannot be less than zero")
  if(maxit<0 ) stop("maxit cannot be less than zero")

  if(length(rho)>1 || length(alpha)>1 || length(lambda)>1) {
    if(verbose) print("rho/alpha/lambda cannot be a vectors. Only first element will be used..")
  }
  lambda = as.numeric(lambda)[1]
  alpha = as.numeric(alpha)[1]
  rho = as.numeric(rho)[1]
  
  #Obtain dimensions:  
  batches = levels(batch) #batch must be factor
  if(is.null(batches)) stop("The batch input must be a factor vector")
  
  ni = tabulate(batch) #number of samples per batch
  ntot = sum(ni) #total number of observations
  nbatches = length(ni) #number of batches
  p = ncol(X) #dimension of fixed design matrix
  if(is.null(p)) stop("The X input must be a numerical matrix")
  #CHECK dimensions:
  if(length(y)!=ntot || nrow(X)!=ntot) stop("The dimensions of y,X and batch did not match")
  
  #Fit conventional EN REGRESSION to obtain init values of beta params (avoids zero betas) Better to use EN?
  if(is.null(beta)) { #if beta not provided
    if(glmnetWarmup) { #if glmnet installed
      if(verbose) print("Fitting glm-model for init start values of beta")
      fit = glmnet::glmnet(x = X,y = y, alpha=alpha, lambda=lambda, family="gaussian",standardize=FALSE,intercept = FALSE)
      beta = coef(fit)[-1] #obtain init betas (remove intercept which is zero)
    } else { #if glmnet not installed (get marginal beta estimates as warmups)
      beta = c(y%*%X)/colSums(X^2) #use marginal beta estimates (without rho correction)
    }
  } else {
    if(p!=length(beta)) stop("The provided initital coefficients did not have correct length!")
  }
  
  #converting X to Matrix type if "matrix"  (this speeds up crossproduct )
  if(class(X)=="matrix") X = Matrix::Matrix(X) #quite fast to do

  #lambda rescaled with n to avoid using mean instead of sum
  lambda = lambda*ntot  
  
  if(verbose) print("Preparing variables for C-call...")
  #Prepare variable to be used in calculations:
  YXsum <- rep(0,p) #suff.stat of cross prod of response and covars (vector)
  XsqSumRidge <- rep(0,p) # suff.stat of covar (vector), sum([X,-j]*Xj) and Sum(Xj^2) (cross prod and squared variant)
  bXXjSum = 0 #calculate #cross prod of beta-j,X-j,Xj (betaX_datavec = beta-j%*%X-j already)

  startInd_Batch = as.integer(c(0,cumsum(ni))) #get start index for vectorized batch vector
  betaX_datavec = rep(NA,ntot) #used to store data vector for each batch (since n<<p)
  invCXvec = rep(NA,ntot*p) #very long vector: n x p (vectorized)
  #following could be implemented in C code:
  #  betaXcrossjSum <- rep(0,p) #dynamic constant changing with j iterations: a vector before multiplying with beta
  for (i in 1:nbatches) { #run over all batches
    bat = batches[i] #get batch name
    if(ni[i]==0) next #skip if no data in batch (this is fine since all is vectorized)
    COVAR =matrix(rho,ncol=ni[i],nrow=ni[i]) #Obtain correlation matrix: insert proportion rho in (0,1)
    diag(COVAR) = 1 #total variance is factorized out
    invC = Matrix::Matrix(solve(t(chol(COVAR)))) #Convert to matrix type after inverting cholesky matrix (recognize as lower triangular)

    #Multiply cholesky matrix with ind. bsed
    ind = which(bat==batch) #get index of inds in batch
    
    #linear transform covars (ni x p) in dimension (VERY SLOW FOR LARGE p)
    X2 = invC%*%X[ind,] #Note: %*% is faster than crossprod and fastest using Matrix
    y2 = invC%*%y[ind] #linear transform response FAST

    ################################
    #prepare suff vars (summing up)#  
    ################################
    
    #TIME CONSUMING FOR LARGE p
    YXsum <- YXsum + Matrix::t(y2)%*%X2 #crossprod X,y (1 x p)

    #Update covar: 
    XsqSumRidge <- XsqSumRidge + Matrix::colSums(X2^2) #sum up vector
    
    #store data vector for each batch, to be multiplied by X,j later (note init beta, and index=1 not used  )
    betaX_data = Matrix::t(X2[,-1]%*%beta[-1])  #j=1 (ni x (p-1))
    bXXjSum = bXXjSum + sum(betaX_data*X2[,1]) #calculate constant sum 'on-the-fly' (j=1)
    
    #store long vectorized matrices:
    indvec1 = startInd_Batch[i] + 1:ni[i]
    betaX_datavec[indvec1] = betaX_data #insert in vectorized over batches
    
    indvec2 = startInd_Batch[i]*p + 1:(ni[i]*p) #get data indices to use at batch i
    invCXvec[indvec2] = as.numeric(X2) #insert in vectorized over batches
  	if(verbose) print(paste0("Done preprosessing with ",i,"/",nbatches," batches"))
  }    
 # rm(X,y,X2,y2,betaX_data);gc()
# invCXvec2=invCXvec;betaX_datavec2=betaX_datavec;bXXjSum2=bXXjSum;XsqSumRidge2=XsqSumRidge;YXsum2=YXsum
#  max(abs(invCXvec2-invCXvec));max(abs(betaX_datavec2-betaX_datavec2));max(abs(bXXjSum2-bXXjSum));max(abs(XsqSumRidge2-XsqSumRidge));max(abs(t(YXsum2)-YXsum))

  if(glmnetPenalty) {
    XsqSumRidge <- XsqSumRidge + lambda*(1-alpha) #whether to use glmnet penalty (divide Ridge penalty by 2)
  } else {
    XsqSumRidge <- XsqSumRidge + 2*lambda*(1-alpha) #whether to use ordinary penalty
  }
  
  #DATA INPUT FOR LOOP ITERATIONS: 3*(p-long vectors) + 1*(n-long vector)  beta, XsqSumRidge + YXmean  , betaX_datavec
  #beta,YXsum,XsqSumRidge: p-long vectors
  #betaX_datavec A vectorized list with nbatches batch list-elements each containing ni long vectors (n in total length)
  #invCXvec A vectorized list with nbatches batch list-elements each containing 'ni x p' large matriced ('n x p' in total size) 
  
  if(verbose) print("Running C code. This may take a while for very large p....")
  time = system.time({
   beta = .C("iterate",beta,p,nbatches, ni, startInd_Batch, as.numeric(YXsum), as.numeric(XsqSumRidge), betaX_datavec, invCXvec, bXXjSum, as.numeric(lambda), as.numeric(alpha), as.numeric(toler), as.integer(maxit),PACKAGE="lme1en")[[1]]
  })[3]
  if(verbose) print(paste0("Time=",round(time,2),"s"))
	return(beta)
} #end function


