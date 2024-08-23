#' MRcML method for overlapping samples with Data Perturbation
#'
#' This is the main function of MRcML method for overlapping samples with data perturbation.
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
#' @param rho Correlation between GWAS summary statistics due to overlapping samples.
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

mr_cML_DP_Overlap <- function(b_exp,b_out,
                        se_exp,se_out,
                        K_vec = 0:(length(b_exp) - 2),
                        random_start = 0,
                        random_start_pert = 0,
                        maxit = 100,
                        num_pert = 100,
                        random_seed = 0,
                        n,
                        rho=0)
{
  if(random_seed)
  {
    set.seed(random_seed)
  }
  t=0
  c1=1
  c2=1
  theta_v = theta_MA_v = NULL
  if(c1<1){c1=1}
  if(c2<1){c2=1}
  se_exp = sqrt(c1)*se_exp
  se_out = sqrt(c2)*se_out 
  ind_filter = which(abs(b_exp/se_exp)>t)
  b_exp_used = b_exp[ind_filter]; se_exp_used=se_exp[ind_filter]; b_out_used=b_out[ind_filter]; se_out_used=se_out[ind_filter]
  cML_res = mr_cML_Overlap(b_exp = b_exp_used, b_out = b_out_used, se_exp = se_exp_used, se_out = se_out_used,
                     random_start = random_start, n = n, rho = rho, maxit = maxit)
  p = length(se_exp)
  sigma <- lapply(1:p,function(i){ matrix(c(se_exp[i]^2, se_exp[i]*se_out[i]*rho,
                                            se_exp[i]*se_out[i]*rho, se_out[i]^2),
                                          2)})
  
  for(pt_ind in 1:num_pert){
    epis = lapply(1:p,function(i){MASS::mvrnorm(1, mu = c(0,0), Sigma = sigma[[i]])})
    epis = matrix(unlist(epis),ncol=2,byrow=T)
    b_exp_new = b_exp + epis[,1]
    b_out_new = b_out + epis[,2]
    ind_filter = which(abs(b_exp_new/se_exp)>t)
    b_exp_new = b_exp_new[ind_filter]
    b_out_new = b_out_new[ind_filter]
    se_exp_new = se_exp[ind_filter]
    se_out_new = se_out[ind_filter]
    K_vec_new = 0:(length(b_exp_new) - 2)
    cML_res_b = mr_cML_Overlap(b_exp = b_exp_new, b_out = b_out_new, se_exp = se_exp_new, se_out = se_out_new,
                         K_vec = K_vec_new, random_start = random_start_pert, n = n, rho = rho, maxit = maxit)
    theta_MA_v = c(theta_MA_v,cML_res_b$MA_BIC_theta)
    theta_v = c(theta_v,cML_res_b$BIC_theta)
  }
  
  
  
  # cML-MA-BIC-DP
  MA_BIC_DP_theta = mean(theta_MA_v)
  MA_BIC_DP_se = sd(theta_MA_v)
  MA_BIC_DP_p = pnorm(-abs(MA_BIC_DP_theta/MA_BIC_DP_se))*2
  
  # cML-BIC-DP
  BIC_DP_theta = mean(theta_v)
  BIC_DP_se = sd(theta_v)
  BIC_DP_p = pnorm(-abs(BIC_DP_theta/BIC_DP_se))*2
  
  
  return(list(MA_BIC_theta = cML_res$MA_BIC_theta,
              MA_BIC_se = cML_res$MA_BIC_se,
              MA_BIC_p = cML_res$MA_BIC_p,
              BIC_theta = cML_res$BIC_theta,
              BIC_se = cML_res$BIC_se,
              BIC_p = cML_res$BIC_p,
              BIC_invalid = cML_res$BIC_invalid,
              MA_BIC_DP_theta = MA_BIC_DP_theta,
              MA_BIC_DP_se = MA_BIC_DP_se,
              MA_BIC_DP_p = MA_BIC_DP_p,
              BIC_DP_theta = BIC_DP_theta,
              BIC_DP_se = BIC_DP_se,
              BIC_DP_p = BIC_DP_p
  ))
}
