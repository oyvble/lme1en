#' @title genData
#' @description Generating a example dataset (age estimation)
#' @details The function generates a dataset based on normal distribution where age (uniformly distributed) is a underlying explanatory variable
#' 
#' @param ntot total number of observations
#' @param nsites number of sites
#' @param nbatches number of batches
#' @param propAgeAs #proportion of sites with assosiation to age
#' @param sd_batch standard dev of batch effects (normal distr effects)
#' @param sd_signal standard dev of signal (norm distr noise)
#' @param sd_bparam sd of age coefficient 
#' @param ageRange Age range of considered ages
#' @param seed The user can set seed if wanted
#' @param verbose boolean, show progress Default: FALSE
#' @return A list with variables y (response), X (design matrix belonging to fixed effects coefficients), batch (factor with batch effect names), age (true ages)
#' @export
#' @examples
#' \dontrun{ 
#' dat = genData(seed=1)
#' }
#ntot = 1000;nsites = 1000;nbatches = 5;propAgeAs=0.05;sd_batch = 0.5;sd_signal = 0.1;sd_bparam = 0.03;ageRange = c(10,35);seed=NULL;verbose=TRUE
genData = function(ntot = 1000, nsites = 1000, nbatches = 5, propAgeAs=0.05, sd_batch = 0.5, sd_signal = 0.1, sd_bparam = 0.03 , ageRange = c(10,35), seed=NULL,verbose=FALSE) {
  if(!is.null(seed)) set.seed(seed) #set seed if provided

   #GENERATE AGES
  ageMid=sum(ageRange)/2 #get mid age (to centralize age when generating data)
  agei = runif(ntot,ageRange[1],ageRange[2]) #get age vector for each inds

  #GENERATE Batch levels and random effects
  batchInd = sort(sample( 1:nbatches,ntot,replace=TRUE)) #get batch index for each inds
  batchRE = rnorm(nbatches,0,sd_batch) #simulate random batch effect (per batch)
  
  ##Simulate age effect (site specific)
  nAgeAs = ceiling(propAgeAs*nsites)    #number of sites with assosiation to age
  AgeAsind = sort(sample(1:nsites,size=nAgeAs,replace=FALSE)) #get index of selected sites (to be assosiated with age)
  bj = rnorm(nAgeAs,0,sd_bparam) #simulate effect on age #site specific
  
  #Generate noise (all sites and observations)
  if(verbose) print("Generating noise...")
  epsij = Matrix::Matrix(rnorm(ntot*nsites,0,sd_signal ),nrow=ntot,ncol=nsites) #random noise per individal per marker

  #generate data on each site (X covariate)
  if(verbose) print("Including batch effect...")
  X = batchRE[batchInd] + epsij #include batch effect and noise (include batch effect on all sites)
  
  #include age effects on selected sites 
  if(verbose) print("Including age effect...")
  tmp = t(bj%*%t(agei-ageMid))
  X[,AgeAsind] = X[,AgeAsind] + tmp
#  for(j in AgeAsind)  X[,j] = X[,j] + bj[j==AgeAsind]*(agei-ageMid) #add age effect (SLOW FOR LARGE nAgeAs)
  colnames(X) = paste0("Site ",1:nsites)
  rownames(X) = paste0("ID",1:ntot) 
  
# max(abs(X[,AgeAsind[1]] - (bj[1]*(agei-ageMid) + batchRE[batchInd] +epsij[,AgeAsind[1]])))

  #SCALE AND STORE DATA:
  if(verbose) print("Standardizing data...")
  X = scale(X,center=TRUE,scale=TRUE) #normalize coeffs
  agei = scale(agei,center=TRUE,scale=TRUE) #scale age

  coefAgeAs = bj #extract age related coeffs
  names(coefAgeAs) = paste0("Site ",AgeAsind)
  
  dat = list(y=agei,X=X,batch=as.factor(batchInd), coefEffect=coefAgeAs) #store data in dataframe
  #dat = data.frame(age=agej,batch=as.factor(batchInd),beta=beta) #store data in dataframe

  return(dat)
}
  