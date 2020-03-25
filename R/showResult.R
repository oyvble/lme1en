
#' @title showResult
#' @description For plotting true and predicted values for train/test data
#' @details Allow for scaled variables (will scale back to original scale)
#' @param true True values (x)
#' @param pred Predicted values (y)
#' @param train index of train data
#' @param test index of test data
#' @param errorFun A error function (must return av constant when sending two equally long vectors)
#' @export

showResult = function(true,pred,train,test,errorFun) {
  y_mean = attr(true,"scaled:center") #get mean-
  y_sd = attr(true,"scaled:scale") #get sd
  y_mean = ifelse(is.null(y_mean),0,y_mean) 
  y_sd = ifelse(is.null(y_sd),1,y_sd) 
  
  true = true*y_sd + y_mean #scale back
  pred = pred*y_sd + y_mean #scale back
  txt1 = paste0("error(train)=", round( errorFun(true[train],pred[train]),1) ) 
  txt2 = paste0("error(test)=", round( errorFun(true[test],pred[test]),1) ) 
  txt = paste0(c(txt1,txt2),collapse="\n")
  plotdata(true[train],pred[train],yr=range(c(true,pred)),txt = txt)
  points( true[test],pred[test],col=2,pch=19)
}

