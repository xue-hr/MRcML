#' Estimate with With Data Perturbation
#'
#' With multiple perturbed data, get estimated theta,
#' se of estimated theta
#' and negative log-likelihood,
#' using multiple random starting points.
#'
#' @param b_exp Vector of estimated effects for exposure.
#' @param b_out Vector or estimated effects for outcome.
#' @param se_exp Vector of standard errors for exposure.
#' @param se_out Vector of standard errors for outcome.
#' @param K Constraint parameter, number of invalid IVs.
#' @param num_pert Number of perturbation, default is 200.
#' @param random_start Number of random starting points, default is 0.
#' @param maxit Maximum number of iteration.
#'
#' @return A list contains:
#' theta_v is the vector estimated thetas,
#' se_v is vector of standard errors of estimated thetas,
#' l_v is vector of negative log-likelihood. Vectors all have length
#' num_pert.
#' @export
#'
#' @examples
cML_estimate_DP <- function(b_exp,b_out,
                            se_exp,se_out,
                            K,num_pert = 200,
                            random_start = 0,
                            maxit = 100)
{
  p = length(b_exp)
  theta_v = NULL
  sd_v = NULL
  l_v = NULL
  for(pt_ind in 1:num_pert)
  {
    b_exp_new = b_exp + rnorm(p,0,1)*se_exp
    b_out_new = b_out + rnorm(p,0,1)*se_out
    MLE_result = cML_estimate_random(b_exp = b_exp_new,
                                     b_out = b_out_new,
                                     se_exp = se_exp,
                                     se_out = se_out,
                                     K = K,
                                     random_start = random_start,
                                     maxit = maxit)

    theta_v = c(theta_v,MLE_result$theta)
    sd_v = c(sd_v,MLE_result$se)
    l_v = c(l_v,MLE_result$l)
  }

  return(list(theta_v = theta_v,
              se_v = sd_v,
              l_v = l_v))
}
