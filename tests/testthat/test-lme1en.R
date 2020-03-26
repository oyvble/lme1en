#devtools::test("C:\\Users\\oyvbl\\Dropbox\\git\\lme1en")

test_that("generated data are reproducible", {
  dat = genData(ntot = 20, nsites = 10, nbatches = 3, propAgeAs=0.8, sd_batch = 0.5, sd_signal = 0.1, sd_bparam = 0.03 , ageRange = c(10,35), seed=1)
  
  expect_equal(sum(dat$coefEffect) , -0.005359324) #sum of coefficients 
  
  expect_equal( as.integer(table(dat$batch)) , c(9,8,3)) #batch levels 

  expect_equal( attr(dat$y,"scaled:center") ,  23.879177832 ) #mean of responsens
                
  expect_equal( attr(dat$y,"scaled:scale") ,  7.1529477905 ) #sd of responses
   
  expect_equal( as.numeric(attr(dat$X,"scaled:center")[1]) ,  0.0849312230835395 ) #mean of 1st covariate
  
  expect_equal( as.numeric(attr(dat$X,"scaled:scale")[1]) ,  0.448171693458162 ) #sd of 1st covariate
  
})

test_that("algorithm gives correct estimated coefficient", {
  dat = genData(ntot = 20, nsites = 10, nbatches = 3, propAgeAs=0.8, sd_batch = 0.5, sd_signal = 0.1, sd_bparam = 0.03 , ageRange = c(10,35), seed=1)
  
  bhat1 = lme1en(y=dat$y,X=dat$X,batch=dat$batch, lambda=1, alpha=0.5,  rho=0.3, beta=NULL, glmnetPenalty=TRUE,glmnetWarmup=TRUE, maxit = 10000, toler = 1e-3 ,verbose=FALSE)
  expect_equal( bhat1[2],0.0815136919) #use glmnetPenalty

  bhat2 = lme1en(y=dat$y,X=dat$X,batch=dat$batch, lambda=1, alpha=0.5,  rho=0.3, beta=NULL, glmnetPenalty=FALSE,glmnetWarmup=TRUE, maxit = 10000, toler = 1e-3 ,verbose=FALSE)
  expect_equal( bhat2[2],0.0566187871) #dont use glmnetPenalty
  
})