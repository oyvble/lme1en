#include<cmath>
using namespace std;

extern "C" {

//p,ni, beta,YXsum,XsqSumRidge, betaX_datavec, invCXvec,bXXjSum, lambda,alpha ,toler,maxiters
void iterate(double *beta, int *p, int *nbatches, int *ni,int *startInd_Batch, double *YXsum, double *XsqSumRidge, double *betaX_datavec, double *invCXvec, double *bXXjSum, double *lambda, double *alpha, double *toler, int *maxiters) {
	//beta = initiated coefficients (fixed)
	//p = number of coefficients
	//nbatches = number of batches
	//ni = number of observations inside batches	
	double betaj_prev,RESj,ascale; //used to store a copy of betaj
	double bshift = (*lambda)*(*alpha); //values to shift: this is constant b in formula
	int i,j,jnext,batchind; //used for calculating next index
	bool anyIsLarge; //used to check for non-convergence (difference is too large)
	int shift1,shift_j,shift_jnext; //used to store index of large vector (when traversing batches)
	double betajXjAdd,betajXjSub; //temporary vars
	
	int iter = 0; //counting number of loop-traversing
	bool done = false; //whether convergence has been reached
	while(!done) {
		anyIsLarge = false; //any still false after running through all this mean convergence
		for (j = 0; j < *p; j++) {	//for each variables	
		  RESj = YXsum[j] - *bXXjSum; //calculate (conditional) residual where var j is removed: this is constant c in formula
		  ascale = XsqSumRidge[j]; //get scaling constant a in formula (the Ridge penalty 2*lambda*(1-alpha) is already included)
		  RESj = RESj/ascale; //scaling residuals
      
		  //Using following formula to get betahat: x = sgn(c/a)*(|c/a| - b/a)_+
		  betaj_prev = beta[j]; //store val of beta
		  beta[j] = abs(RESj) - bshift/ascale; //obtain estimate
		  if( beta[j]< 0) { //if negative
			beta[j] = 0.0; //set to zero
		  } else if( RESj <0 ) { //else not zero and if the RESj variable was negative
			  beta[j] = -beta[j]; //change sign of betaj 
		  }
		  //Check if bhat is different from prev, set true if yes (must run loop more)
		  if( abs(beta[j]-betaj_prev) > *toler) {
			  anyIsLarge = true; 
		  }
		  
		  //POST CALCS; UPDATE bXXjSum and betaX_datavec: Addidative for index=j and subtract for index=j+1 (to be handled next)
		  jnext = j+1;
		  if(jnext==*p) {
		    jnext = 0; //start at index zero again
		  }
		  *bXXjSum = 0.0; //Used for calculating cross prod of beta-j,X-j,Xj (betaX_datavec = beta-j%*%X-j already calc)
		  for (batchind=0; batchind < *nbatches; batchind++)  { //run through all data (each batch)
			  
			  //Obtain shifts for accessing data in the long vectors
			  shift1 = startInd_Batch[batchind]; //get data indices to use at batch i
			  shift_j = startInd_Batch[batchind]*(*p) + ni[batchind]*j; //get data indicies for variable j at batch i
			  shift_jnext = startInd_Batch[batchind]*(*p) + ni[batchind]*jnext; //get data indicies for variable j at batch i
			
			  for (i=0; i < ni[batchind]; i++)  { //run through all observations within given batch
				betajXjAdd = beta[j]*invCXvec[shift_j+i]; //#vector to be added (current) [vector scaling]
				betajXjSub = beta[jnext]*invCXvec[shift_jnext+i]; //vector to be subtracted (next) [vector scaling]
				betaX_datavec[shift1+i] += betajXjAdd - betajXjSub; //calculate difference [vector additive]
				*bXXjSum += betaX_datavec[shift1+i]*invCXvec[shift_jnext+i]; //calc dot product (over data vector)
			  }
		  }
		  
		} //end for each j variables
		iter++; //adding iteration
		if(!anyIsLarge || iter> *maxiters ) { //check if done
			done = true; 
		}
	}
} //end main function

} //end external