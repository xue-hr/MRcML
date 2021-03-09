#' Estimate with Regular Likelihood Based on Multiple Permutations
#'
#' Estimate theta, standard error of theta, negative likelihood
#' with constrained maximum likelihood, averaged from multiple permutations.
#'
#' @param b_exp Vector of estimated effects for exposure.
#' @param b_out Vector or estimated effects for outcome.
#' @param se_exp Vector of standard errors for exposure.
#' @param se_out Vector of standard errors for outcome.
#' @param K Constraint parameter, number of invalid IVs.
#' @param num_permutation Number of permutations.
#' @param initial_theta Starting point for theta.
#' @param initial_mu Starting point for mu.
#' @param maxit Maximum number of iteration.
#'
#' @return A list contains: theta, estimated causal effect;
#' se_theta, standard error of theta;
#' Neg_l, minus log likelihood.
#' @export
#'
#' @examples
cML_estimate_pt <- function(b_exp,b_out,
                            se_exp,se_out,
                            K,num_permutation = 100,
                            initial_theta = 0,
                            initial_mu = rep(0, length(b_exp)),
                            maxit = 100)
{
  p = length(b_exp)
  theta_v = NULL
  sd_v = NULL
  l_v = NULL
  for(pt_ind in 1:num_permutation)
  {
    b_exp_new = b_exp + rnorm(p,0,1)*se_exp
    b_out_new = b_out + rnorm(p,0,1)*se_out
    MLE_result = cML_estimate(b_exp = b_exp_new,
                              b_out = b_out_new,
                              se_exp = se_exp,
                              se_out = se_out,
                              K = K,
                              initial_theta = initial_theta,
                              initial_mu = initial_mu,
                              maxit = maxit)
    sd_theta = cML_SdTheta(b_exp = b_exp_new,
                           b_out = b_out_new,
                           se_exp = se_exp,
                           se_out = se_out,
                           MLE_result$theta, MLE_result$b_vec, MLE_result$r_vec)
    theta_v = c(theta_v,MLE_result$theta)
    sd_v = c(sd_v,sd_theta)
    l_v = c(l_v,
            sum((b_exp_new - MLE_result$b_vec)^2/(2 * se_exp^2)) +
              sum((b_out_new - MLE_result$theta * MLE_result$b_vec -
                     MLE_result$r_vec)^2/(2 * se_out^2)))
  }
  weight_vec = 1/num_permutation
  Combined_theta = sum(theta_v * weight_vec)
  Combined_se = sum(weight_vec * sqrt(sd_v^2 + (theta_v - Combined_theta)^2),
                    na.rm = TRUE)
  return(list(theta = Combined_theta,
              se_theta = Combined_se,
              Neg_l = mean(l_v)))
}
