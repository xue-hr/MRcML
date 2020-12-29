MLE_theta_bg <- function(b_exp,b_out,
                         se_exp,se_out,
                         K, initial_theta = 0, 
                         max_ite = 100)
{
  # Using iterative method to get the MLE of theta and bg's
  
  # b_exp: vector of estimated effect sizes on X
  # b_out: vector of estimated effect sizes on Y
  # se_exp: vector of standard errors of b_exp
  # se_out: vector of standard errors of b_out
  # K: TLP constraint paramater, number of non-zero bg's
  # initial_theta: initial beta value at step 0
  # max_ite: maximum number of iteration
  p = length(b_exp)
  
  theta_old = initial_theta-1
  theta_new = initial_theta
  ite_ind = 0
  while ( (abs(theta_old - theta_new) > 1e-5) & (ite_ind<max_ite) ) {
    ite_ind = ite_ind + 1
    theta_old = theta_new
    if(K == 0)
    {
      new_bg = rep(0,p)
    } else{
      new_bg = estimate_bg(b_exp,b_out,
                           se_exp,se_out,
                           K,theta_new)
    }

    theta_new = estimate_thetam(b_exp,b_out,
                                se_exp,se_out,
                                new_bg)
    #cat(theta_new,"\n")
  }
  
  ### one more step for positive K
  if(K!=0)
  {
    new_bg = estimate_bg(b_exp,b_out,
                         se_exp,se_out,
                         K,theta_new)
  }

  
  return(list(MLE_theta = theta_new,
              MLE_bg = new_bg))
}