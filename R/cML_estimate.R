estmate_theta_bg_mu <- function(b_exp,b_out,
                                se_exp,se_out,
                                K,initial_theta = 0,
                                initial_mu = rep(0,length(b_exp)),
                                maxit = 100)
{
  p = length(b_exp)
  ### initialize
  theta = initial_theta
  theta_old = theta-1
  mu_vec = initial_mu

  ###
  ite_ind = 0
  while( (abs(theta_old - theta) > 1e-7) & (ite_ind<maxit))
  {
    theta_old = theta
    ite_ind = ite_ind + 1

    ### first, update v_bg
    if(K>0)
    {
      v_importance = (b_out - theta*mu_vec)^2 / se_out^2
      nonzero_bg_ind = sort((order(v_importance,decreasing = T))[1:K])
      v_bg = rep(0,p)
      v_bg[nonzero_bg_ind] = (b_out - theta*mu_vec)[nonzero_bg_ind]
    } else{
      v_bg = rep(0,p)
    }


    ### second, update mu_vec
    mu_vec =
      (b_exp / se_exp^2 + theta*(b_out - v_bg) / se_out^2) /
      (1 / se_exp^2 + theta^2 / se_out^2)

    ### third, update theta
    theta =
      sum((b_out - v_bg)*mu_vec / se_out^2) /
      sum(mu_vec^2 / se_out^2)

  }
  ### one more step for v_bg and mu
  if(K>0)
  {
    nonzero_ind = which(v_bg!=0)
    mu_vec[nonzero_ind] = b_exp[nonzero_ind]
    v_bg[nonzero_ind] = (b_out - theta*mu_vec)[nonzero_ind]
  }

  return(list(theta = theta,
              mu_vec = mu_vec,
              v_bg = v_bg))

}
