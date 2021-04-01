#' MRcML method with Data Perturbation
#'
#' This is the main function of MRcML method with data perturbation.
#'
#' @param b_exp Vector of estimated effects for exposure.
#' @param b_out Vector or estimated effects for outcome.
#' @param se_exp Vector of standard errors for exposure.
#' @param se_out Vector of standard errors for outcome.
#' @param K_vec Sets of candidate K's, the constraint parameter representing number of invalid IVs.
#' @param random_start Number of random start points for cML, default is 0.
#' @param random_start_pert Number of random start points for cML with data perturbation, default is 0.
#' @param maxit Maximum number of iterations for each optimization.
#' @param num_pert Number of perturbation, default is 200.
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
#' for MA_BIC_DP_theta, MA_BIC_DP_se, MA_BIC_DP_p from cML-MA-BIC-DP;
#' for BIC_DP_theta, BIC_DP_se, BIC_DP_p from cML-BIC-DP.
#' BIC_invalid is the set of invalid IVs selected by cML-BIC.
#' @export
#'
#' @examples
mr_cML_DP <- function(b_exp,b_out,
                      se_exp,se_out,
                      K_vec = 0:(length(b_exp) - 2),
                      random_start = 0,
                      random_start_pert = 0,
                      maxit = 100,
                      num_pert = 200,
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
  rand_pert_theta = NULL
  rand_pert_sd = NULL
  rand_pert_l = NULL
  for(K_value in K_vec)
  {
    rand_res = cML_estimate_random(b_exp = b_exp,
                                   b_out = b_out,
                                   se_exp = se_exp,
                                   se_out = se_out,
                                   K = K_value,
                                   random_start = random_start,
                                   maxit = maxit)
    rand_pert_res = cML_estimate_DP(b_exp = b_exp,
                                    b_out = b_out,
                                    se_exp = se_exp,
                                    se_out = se_out,
                                    K = K_value,
                                    num_pert = num_pert,
                                    random_start = random_start_pert,
                                    maxit = maxit)
    rand_theta = c(rand_theta,rand_res$theta)
    rand_sd = c(rand_sd,rand_res$se)
    rand_l = c(rand_l,rand_res$l)
    invalid_mat = rbind(invalid_mat,rand_res$r_est)
    rand_pert_theta = rbind(rand_pert_theta,rand_pert_res$theta_v)
    rand_pert_sd = rbind(rand_pert_sd,rand_pert_res$se_v)
    rand_pert_l = rbind(rand_pert_l,rand_pert_res$l_v)
  }

  ### get result
  theta_v = rand_theta
  sd_v = rand_sd
  l_v = rand_l
  theta_pt_v = rand_pert_theta
  sd_pt_v = rand_pert_sd
  l_pt_v = rand_pert_l

  var_mat = sd_pt_v^2
  K_vec = 0:(nrow(theta_pt_v) - 1)
  numer_perturb = ncol(theta_pt_v)
  l_pt_v = rowMeans(l_pt_v)
  sd_pt_v = sqrt(diag(var(t(theta_pt_v))))
  theta_pt_mat = theta_pt_v
  theta_pt_v = rowMeans(theta_pt_v)

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

  # cML-MA-BIC-DP
  BIC_vec = log(n) * K_vec + 2 * l_pt_v
  BIC_vec = BIC_vec - min(BIC_vec)
  weight_vec = exp(-1/2 * BIC_vec)
  weight_vec = weight_vec/sum(weight_vec)
  MA_BIC_DP_theta = sum(theta_pt_v * weight_vec)
  MA_BIC_DP_se = sum(weight_vec * sqrt(sd_pt_v^2 + (theta_pt_v - MA_BIC_DP_theta)^2),
                     na.rm = TRUE)
  MA_BIC_DP_p = pnorm(-abs(MA_BIC_DP_theta/MA_BIC_DP_se))*2

  # cML-BIC-DP
  BIC_vec = log(n) * K_vec + 2 * l_pt_v
  BIC_vec = BIC_vec - min(BIC_vec)
  min_ind = which.min(BIC_vec)
  BIC_DP_theta = theta_pt_v[min_ind]
  BIC_DP_se = sd_pt_v[min_ind]
  BIC_DP_p = pnorm(-abs(BIC_DP_theta/BIC_DP_se))*2

  # GOF 1
  BIC_vec = log(n) * K_vec + 2 * l_v
  BIC_vec = BIC_vec - min(BIC_vec)
  min_ind = which.min(BIC_vec)
  pt_sd_v = sd_pt_v[min_ind]
  origin_sd_v = sd_v[min_ind]
  more_var = var(var_mat[min_ind,])
  #n = numer_perturb
  x = theta_pt_mat[min_ind,]
  sd_x = sqrt((mean((x - mean(x))^4) -
                 (numer_perturb-3)/(numer_perturb-1)*var(x)^2)/numer_perturb +
                more_var)
  Tstat =
    (origin_sd_v^2 - pt_sd_v^2)/(sd_x)
  GOF1_p = pnorm(-abs(Tstat))*2

  # GOF 2
  BIC_vec = log(n) * K_vec + 2 * l_v
  BIC_vec = BIC_vec - min(BIC_vec)
  min_ind = which.min(BIC_vec)
  pt_sd_v = sd_pt_v[min_ind]
  origin_sd_v = sd_v[min_ind]
  more_var = var(var_mat[min_ind,])

  Tstat =
    (origin_sd_v^2 - pt_sd_v^2)/(sqrt(2/(numer_perturb-1)*pt_sd_v^4 +
                                        more_var))
  GOF2_p = pnorm(-abs(Tstat))*2

  return(list(MA_BIC_theta = MA_BIC_theta,
              MA_BIC_se = MA_BIC_se,
              MA_BIC_p = MA_BIC_p,
              BIC_theta = BIC_theta,
              BIC_se = BIC_se,
              BIC_p = BIC_p,
              BIC_invalid = BIC_invalid,
              MA_BIC_DP_theta = MA_BIC_DP_theta,
              MA_BIC_DP_se = MA_BIC_DP_se,
              MA_BIC_DP_p = MA_BIC_DP_p,
              BIC_DP_theta = BIC_DP_theta,
              BIC_DP_se = BIC_DP_se,
              BIC_DP_p = BIC_DP_p,
              GOF1_p = GOF1_p,
              GOF2_p = GOF2_p
              ))
}
