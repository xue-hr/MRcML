#' Estimate with Regular Likelihood Using Multiple Random Start Points
#'
#' Get estimated theta, se of estimated theta
#' and negative log-likelihood,
#' using multiple random starting points.
#'
#' @param b_exp Vector of estimated effects for exposure.
#' @param b_out Vector or estimated effects for outcome.
#' @param se_exp Vector of standard errors for exposure.
#' @param se_out Vector of standard errors for outcome.
#' @param K Constraint parameter, number of invalid IVs.
#' @param random_start Number of random starting points, default is 0.
#' @param maxit Maximum number of iteration.
#'
#' @return A list contains: theta is the estimate causal effect,
#' se is standard error of estimated theta,
#' l is negative log-likelihood,
#' r_est is estimated r vector.
#' @export
#'
#' @examples
cML_estimate_random <- function(b_exp, b_out,
                                se_exp, se_out,
                                K,random_start = 0,
                                maxit = 100)
{
  p = length(b_exp)
  min_theta_range = min(b_out/b_exp)
  max_theta_range = max(b_out/b_exp)

  theta_v_RandomCandidate = NULL
  sd_v_RandomCandidate = NULL
  l_v_RandomCandidate = NULL
  invalid_RandomCandidate = NULL

  for(random_ind in 1:(1+random_start))
  {
    if(random_ind == 1)
    {
      initial_theta = 0
      initial_mu = rep(0,p)
    } else {
      initial_theta = runif(1,min = min_theta_range,max = max_theta_range)
      initial_mu = rnorm(p,mean = b_exp,sd = se_exp)
    }


    MLE_result =
      cML_estimate(b_exp,b_out,
                   se_exp,se_out,
                   K = K,initial_theta = initial_theta,
                   initial_mu = initial_mu,
                   maxit = maxit)

    Neg_l =
      sum( (b_exp - MLE_result$b_vec)^2 / (2*se_exp^2) ) +
      sum( (b_out - MLE_result$theta * MLE_result$b_vec - MLE_result$r_vec)^2 /
             (2*se_out^2))

    sd_theta = cML_SdTheta(b_exp,b_out,
                           se_exp,se_out,
                           MLE_result$theta,
                           MLE_result$b_vec,
                           MLE_result$r_vec)

    theta_v_RandomCandidate = c(theta_v_RandomCandidate,MLE_result$theta)
    sd_v_RandomCandidate = c(sd_v_RandomCandidate,sd_theta)
    l_v_RandomCandidate = c(l_v_RandomCandidate,Neg_l)
    invalid_RandomCandidate = rbind(invalid_RandomCandidate,
                                    as.numeric(MLE_result$r_vec))

  }

  min_neg_l = which.min(l_v_RandomCandidate)

  theta_est = theta_v_RandomCandidate[min_neg_l]
  sd_est = sd_v_RandomCandidate[min_neg_l]
  l_est = l_v_RandomCandidate[min_neg_l]
  r_est = invalid_RandomCandidate[min_neg_l,]

  return(list(theta = theta_est,
              se = sd_est,
              l = l_est,
              r_est = r_est
  )
  )
}
