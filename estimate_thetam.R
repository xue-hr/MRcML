estimate_thetam <- function(b_exp,b_out,
                            se_exp,se_out,
                            v_bg)
{
  # Given v_bg from last iteration, get new thetam
  
  # b_exp: vector of estimated effect sizes on X
  # b_out: vector of estimated effect sizes on Y
  # se_exp: vector of standard errors of b_exp
  # se_out: vector of standard errors of b_out
  # v_bg: estimated bg's value from last iteration
  b_out = b_out - v_bg
  
  profile.loglike <- function(beta) {
    -(1/2) * sum((b_out - b_exp * beta)^2/(se_out^2 + se_exp^2 * 
                                             beta^2))
  }
  bound <- 5
  beta.hat <- optimize(profile.loglike, bound * c(-1, 1), maximum = TRUE, 
                       tol = .Machine$double.eps^0.5)$maximum
  while (abs(beta.hat) > 0.95 * bound) {
    bound <- bound * 2
    beta.hat <- optimize(profile.loglike, bound * c(-1, 1), 
                         maximum = TRUE, tol = .Machine$double.eps^0.5)$maximum
  }
  return(beta.hat)
}
