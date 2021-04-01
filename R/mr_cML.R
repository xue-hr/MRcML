#' MRcML method
#'
#' This is the main function of MRcML method, without data perturbation.
#'
#' @param b_exp Vector of estimated effects for exposure.
#' @param b_out Vector or estimated effects for outcome.
#' @param se_exp Vector of standard errors for exposure.
#' @param se_out Vector of standard errors for outcome.
#' @param K_vec Sets of candidate K's, the constraint parameter representing number of invalid IVs.
#' @param random_start Number of random start points for cML, default is 0.
#' @param maxit Maximum number of iterations for each optimization.
#' @param random_seed Random seed, an integer. Default is
#' 0, which does not set random seed; user could specify a positive integer
#' as random seed to get replicable results.
#' @param n Sample size.
#'
#' @return  A list contains full results of cML methods.
#' MA_BIC_theta, MA_BIC_se, MA_BIC_p:
#' Estimate of theta,
#' its standard error and p-value from cML-MA-BIC.
#' Similarly for BIC_theta, BIC_se, BIC_p from cML-BIC;
#' for MA_AIC_theta, MA_AIC_se, MA_AIC_p from cML-MA-AIC;
#' for AIC_DP_theta, AIC_DP_se, AIC_DP_p from cML-AIC.
#' BIC_invalid is the set of invalid IVs selected by cML-BIC,
#' AIC_invalid is the set of invalid IVs selected by cML-AIC.
#' @export
#'
#' @examples
mr_cML <- function(b_exp,b_out,
                   se_exp,se_out,
                   K_vec = 0:(length(b_exp) - 2),
                   random_start = 0,
                   maxit = 100,
                   random_seed = 0,
                   n)
{
  if(random_seed)
  {
    set.seed(random_seed)
  }

  rand_theta = NULL
  rand_sd = NULL
  rand_l = NULL
  invalid_mat = NULL
  for(K_value in K_vec)
  {
    rand_res = cML_estimate_random(b_exp = b_exp,
                                   b_out = b_out,
                                   se_exp = se_exp,
                                   se_out = se_out,
                                   K = K_value,
                                   random_start = random_start,
                                   maxit = maxit)
    rand_theta = c(rand_theta,rand_res$theta)
    rand_sd = c(rand_sd,rand_res$se)
    rand_l = c(rand_l,rand_res$l)
    invalid_mat = rbind(invalid_mat,rand_res$r_est)
  }

  ### get result
  theta_v = rand_theta
  sd_v = rand_sd
  l_v = rand_l

  # cML-MA-BIC
  BIC_vec = log(n) * K_vec + 2 * l_v
  BIC_vec = BIC_vec - min(BIC_vec)
  weight_vec = exp(-1/2 * BIC_vec)
  weight_vec = weight_vec/sum(weight_vec)
  MA_BIC_theta = sum(theta_v * weight_vec)
  MA_BIC_se = sum(weight_vec * sqrt(sd_v^2 + (theta_v - MA_BIC_theta)^2),
                  na.rm = TRUE)
  MA_BIC_p = pnorm(-abs(MA_BIC_theta/MA_BIC_se))*2

  # cML-BIC
  BIC_vec = log(n) * K_vec + 2 * l_v
  BIC_vec = BIC_vec - min(BIC_vec)
  min_ind = which.min(BIC_vec)
  BIC_theta = theta_v[min_ind]
  BIC_se = sd_v[min_ind]
  BIC_p = pnorm(-abs(BIC_theta/BIC_se))*2
  BIC_invalid = which(invalid_mat[min_ind,]!=0)

  # cML-MA-AIC
  AIC_vec = 2*K_vec + 2*l_v
  AIC_vec = AIC_vec - min(AIC_vec)
  weight_vec = exp(-1/2 * AIC_vec)
  weight_vec = weight_vec/sum(weight_vec)
  MA_AIC_theta = sum(theta_v * weight_vec)
  MA_AIC_se = sum(weight_vec * sqrt(sd_v^2 + (theta_v - MA_AIC_theta)^2),
                  na.rm = TRUE)
  MA_AIC_p = pnorm(-abs(MA_AIC_theta/MA_AIC_se))*2

  # cML-AIC
  AIC_vec = 2*K_vec + 2*l_v
  AIC_vec = AIC_vec - min(AIC_vec)
  min_ind = which.min(AIC_vec)
  AIC_theta = theta_v[min_ind]
  AIC_se = sd_v[min_ind]
  AIC_p = pnorm(-abs(AIC_theta/AIC_se))*2
  AIC_invalid = which(invalid_mat[min_ind,]!=0)


  return(list(MA_BIC_theta = MA_BIC_theta,
              MA_BIC_se = MA_BIC_se,
              MA_BIC_p = MA_BIC_p,
              BIC_theta = BIC_theta,
              BIC_se = BIC_se,
              BIC_p = BIC_p,
              BIC_invalid = BIC_invalid,
              MA_AIC_theta = MA_AIC_theta,
              MA_AIC_se = MA_AIC_se,
              MA_AIC_p = MA_AIC_p,
              AIC_theta = AIC_theta,
              AIC_se = AIC_se,
              AIC_p = AIC_p,
              AIC_invalid = AIC_invalid)
         )
}
