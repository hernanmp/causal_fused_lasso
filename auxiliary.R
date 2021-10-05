
fused_lasso_bic = function(y,lambda_grid =  10^seq(0,6,length=500))
{
  n =length(y)
  est_var =  mean(diff(y)^2)*n/(2*(n-1))
  #mean(diff(y)^2)
  
  
  temp =  glmgen::trendfilter(y,k=0,lambda= lambda_grid )
  BIC =  rep(0, length(lambda_grid))
  for(j in 1:length(lambda_grid))
  {
    est =  temp$beta[,j]
    
    
    df= length(which(abs(diff(est))>10^-5))
    
    BIC[j] = 0.5*sum((y- est )^2)   +   est_var*df*(log(n))
    #(log(n))*df
  }
  best_j =  which.min(BIC)
  plot(temp$beta[,best_j])
  #plot(BIC[50:100])
  #print(best_j)
  est =  temp$beta[,best_j]
  return(est)
}
############################################################################

