#' Mendelian randomization with cML method Using Permutation
#'
#' Main function for MR cML methods using random start points. Always use 0 as a
#' start point, plus random_start number of random start points. Estimates
#' are averaged from multiple permutations.
#'
#' @param b_exp Vector of estimated effects for exposure.
#' @param b_out Vector or estimated effects for outcome.
#' @param se_exp Vector of standard errors for exposure.
#' @param se_out Vector of standard errors for outcome.
#' @param n Sample size.
#' @param K_vec Sets of candidate K's, the constraint parameter representing number of invalid IVs.
#' @param num_permutation Number of permutations.
#' @param random_start Number of random start points, default is 0.
#' @param maxit Maximum number of iterations for each optimization.
#' @param random_seed Random seed, an integer. Default is
#' 0, which does not set random seed; user could specify a positive integer
#' as random seed to get replicable results.
#'
#' @return A list contains full results of cML methods.
#' AIC_K: K with minimum AIC;
#' AIC_theta, AIC_se, AIC_p: Estimate of theta, its standard error and p-value from AIC selected model;
#' Similarly for BIC_K, BIC_theta, BIC_se, BIC_p.
#' MA_AIC_theta, MA_AIC_se, MA_AIC_p: Estimate of theta, its standard error and p-value from model averaging with AIC-based weights.
#' Similarly for MA_BIC_theta, MA_BIC_se, MA_BIC_p.
#' @export
#'
#' @examples
mr_cML_pt =function (b_exp, b_out,
                     se_exp, se_out,
                     n, K_vec = 0:(length(b_exp) - 2),
                     num_permutation = 100,
                     random_start = 0, maxit = 100, random_seed = 0)
{
  if (random_seed) {
    set.seed(random_seed)
  }
  p = length(b_exp)
  min_theta_range = min(b_out/b_exp)
  max_theta_range = max(b_out/b_exp)
  theta_v = NULL
  sd_v = NULL
  l_v = NULL
  for (K_value in K_vec) {
    theta_v_RandomCandidate = NULL
    sd_v_RandomCandidate = NULL
    l_v_RandomCandidate = NULL

    for (random_ind in 1:(1 + random_start)) {
      if (random_ind == 1) {
        initial_theta = 0
        initial_mu = rep(0, p)
      }
      else {
        initial_theta = runif(1, min = min_theta_range,
                              max = max_theta_range)
        initial_mu = rnorm(p, mean = b_exp, sd = se_exp)
      }
      MLE_result = cML_estimate_pt(b_exp, b_out,
                                   se_exp, se_out,
                                   K = K_value,
                                   initial_theta = initial_theta,
                                   initial_mu = initial_mu,
                                   maxit = maxit,
                                   num_permutation = num_permutation)
      Neg_l = MLE_result$Neg_l
      sd_theta = MLE_result$se_theta
      theta_v_RandomCandidate = c(theta_v_RandomCandidate,
                                  MLE_result$theta)
      sd_v_RandomCandidate = c(sd_v_RandomCandidate, sd_theta)
      l_v_RandomCandidate = c(l_v_RandomCandidate, Neg_l)

    }
    min_neg_l = which.min(l_v_RandomCandidate)
    theta_v = c(theta_v, theta_v_RandomCandidate[min_neg_l])
    sd_v = c(sd_v, sd_v_RandomCandidate[min_neg_l])
    l_v = c(l_v, l_v_RandomCandidate[min_neg_l])

  }

  AIC_vec = 2 * K_vec + 2 * l_v
  AIC_vec = AIC_vec - min(AIC_vec)
  weight_vec = exp(-1/2 * AIC_vec)
  weight_vec = weight_vec/sum(weight_vec)
  Combined_theta = sum(theta_v * weight_vec)
  Combined_se = sum(weight_vec * sqrt(sd_v^2 + (theta_v - Combined_theta)^2),
                    na.rm = TRUE)
  AIC_theta = Combined_theta
  AIC_se = Combined_se
  AIC_vec = 2 * K_vec + 2 * l_v
  AIC_min_ind = which.min(AIC_vec)
  AIC_select_theta = theta_v[AIC_min_ind]
  AIC_select_se = sd_v[AIC_min_ind]
  AIC_select_K = K_vec[AIC_min_ind]

  BIC_vec = log(n) * K_vec + 2 * l_v
  BIC_vec = BIC_vec - min(BIC_vec)
  weight_vec = exp(-1/2 * BIC_vec)
  weight_vec = weight_vec/sum(weight_vec)
  Combined_theta = sum(theta_v * weight_vec)
  Combined_se = sum(weight_vec * sqrt(sd_v^2 + (theta_v - Combined_theta)^2),
                    na.rm = TRUE)
  BIC_theta = Combined_theta
  BIC_se = Combined_se
  BIC_vec = log(n) * K_vec + 2 * l_v
  BIC_min_ind = which.min(BIC_vec)
  BIC_select_theta = theta_v[BIC_min_ind]
  BIC_select_se = sd_v[BIC_min_ind]
  BIC_select_K = K_vec[BIC_min_ind]

  pval_MA_AIC = pnorm(-abs(AIC_theta/AIC_se)) * 2
  pval_AIC_select = pnorm(-abs(AIC_select_theta/AIC_select_se)) *
    2
  pval_MA_BIC = pnorm(-abs(BIC_theta/BIC_se)) * 2
  pval_BIC_select = pnorm(-abs(BIC_select_theta/BIC_select_se)) *
    2
  return(list(AIC_theta = AIC_select_theta,
              AIC_se = AIC_select_se,
              AIC_p = pval_AIC_select,
              AIC_K = AIC_select_K,
              BIC_theta = BIC_select_theta,
              BIC_se = BIC_select_se,
              BIC_p = pval_BIC_select,
              BIC_K = BIC_select_K,
              MA_AIC_theta = AIC_theta,
              MA_AIC_se = AIC_se,
              MA_AIC_p = pval_MA_AIC,
              MA_BIC_theta = BIC_theta,
              MA_BIC_se = BIC_se,
              MA_BIC_p = pval_MA_BIC
              ))
}

