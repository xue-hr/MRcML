#' Title
#'
#' @param b_exp
#' @param b_out
#' @param se_exp
#' @param se_out
#' @param n
#' @param K_vec
#' @param random_start
#'
#' @return
#' @export
#'
#' @examples
mr_cML <- function(b_exp,b_out,
                   se_exp,se_out,
                   n,
                   K_vec = 0:(length(b_exp)-2),
                   random_start = 0
                   )
{
  p = length(b_exp)

  min_theta_range = min(b_out/b_exp)
  max_theta_range = max(b_out/b_exp)
  # perform TLP AIC ---------------------------------------------------------
  theta_v = NULL
  sd_v = NULL
  l_v = NULL
  nonzero_mat = NULL

  for(K_value in K_vec)
  {
    theta_v_RandomCandidate = NULL
    sd_v_RandomCandidate = NULL
    l_v_RandomCandidate = NULL
    nonzero_mat_RandomCandidate = NULL

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
                     K = K_value,initial_theta = initial_theta,
                     initial_mu = initial_mu,
                     maxit = 100)

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
      nonzero_mat_RandomCandidate = rbind(nonzero_mat_RandomCandidate,
                                          t(MLE_result$r_vec))

    }

    min_neg_l = which.min(l_v_RandomCandidate)

    theta_v = c(theta_v,theta_v_RandomCandidate[min_neg_l])
    sd_v = c(sd_v,sd_v_RandomCandidate[min_neg_l])
    l_v = c(l_v,l_v_RandomCandidate[min_neg_l])
    nonzero_mat = rbind(nonzero_mat,
                        nonzero_mat_RandomCandidate[min_neg_l,])
  }

  ### AIC
  AIC_vec = 2*K_vec + 2*l_v
  AIC_vec = AIC_vec - min(AIC_vec)
  weight_vec = exp(-1/2*AIC_vec)
  weight_vec = weight_vec/sum(weight_vec)
  Combined_theta = sum(theta_v * weight_vec)
  Combined_se = sum(weight_vec*sqrt(sd_v^2 +
                                      (theta_v -
                                         Combined_theta)^2
  ),na.rm = TRUE
  )
  AIC_theta = Combined_theta
  AIC_se = Combined_se

  ### AIC_select
  AIC_vec = 2*K_vec + 2*l_v
  AIC_min_ind = which.min(AIC_vec)
  AIC_select_theta = theta_v[AIC_min_ind]
  AIC_select_se = sd_v[AIC_min_ind]
  AIC_select_K = K_vec[AIC_min_ind]
  AIC_select_invalid = which(nonzero_mat[AIC_min_ind,]!=0)


  ### BIC
  BIC_vec = log(n)*K_vec + 2*l_v
  BIC_vec = BIC_vec - min(BIC_vec)
  weight_vec = exp(-1/2*BIC_vec)
  weight_vec = weight_vec/sum(weight_vec)
  Combined_theta = sum(theta_v * weight_vec)
  Combined_se = sum(weight_vec*sqrt(sd_v^2 +
                                      (theta_v -
                                         Combined_theta)^2
  ),na.rm = TRUE
  )
  BIC_theta = Combined_theta
  BIC_se = Combined_se

  ### BIC_select
  BIC_vec = log(n)*K_vec + 2*l_v
  BIC_min_ind = which.min(BIC_vec)
  BIC_select_theta = theta_v[BIC_min_ind]
  BIC_select_se = sd_v[BIC_min_ind]
  BIC_select_K = K_vec[BIC_min_ind]
  BIC_select_invalid = which(nonzero_mat[BIC_min_ind,]!=0)

  ### p-values
  pval_MA_AIC = pnorm(-abs(AIC_theta/AIC_se))*2
  pval_AIC_select = pnorm(-abs(AIC_select_theta/AIC_select_se))*2
  pval_MA_BIC = pnorm(-abs(BIC_theta/BIC_se))*2
  pval_BIC_select = pnorm(-abs(BIC_select_theta/BIC_select_se))*2

  ### return
  return(list(AIC_theta = AIC_select_theta,
              AIC_se = AIC_select_se,
              AIC_p = pval_AIC_select,
              AIC_K = AIC_select_K,
              AIC_invalid = AIC_select_invalid,
              BIC_theta = BIC_select_theta,
              BIC_se = BIC_select_se,
              BIC_p = pval_BIC_select,
              BIC_K = BIC_select_K,
              BIC_invalid = BIC_select_invalid,
              MA_AIC_theta = AIC_theta,
              MA_AIC_se = AIC_se,
              MA_AIC_p = pval_MA_AIC,
              MA_BIC_theta = BIC_theta,
              MA_BIC_se = BIC_se,
              MA_BIC_p = pval_MA_BIC
              ))

}
