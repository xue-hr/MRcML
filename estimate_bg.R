estimate_bg <- function(b_exp,b_out,
                        se_exp,se_out,
                        K,thetam)
{
  # Given thetam from last iteration, get new bg's
  
  # b_exp: vector of estimated effect sizes on X
  # b_out: vector of estimated effect sizes on Y
  # se_exp: vector of standard errors of b_exp
  # se_out: vector of standard errors of b_out
  # K: TLP constraint paramater, number of non-zero bg's
  # thetam: estimated theta value from last iteration
  
  v_importance = (b_out - thetam*b_exp)^2/(se_out^2 + thetam^2*se_exp^2)
  nonzero_bg_ind = (order(v_importance,decreasing = T))[1:K]
  v_bg = rep(0,length(b_exp))
  v_bg[nonzero_bg_ind] = (b_out - thetam*b_exp)[nonzero_bg_ind]
  
  return(v_bg)
}