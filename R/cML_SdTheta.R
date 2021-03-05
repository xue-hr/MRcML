#' Standard Error of Estimated Theta
#'
#' Get the standard error of estimated theta from constrained maximum likelihood.
#'
#' @param b_exp Vector of estimated effects for exposure.
#' @param b_out Vector or estimated effects for outcome.
#' @param se_exp Vector of standard errors for exposure.
#' @param se_out Vector of standard errors for outcome.
#' @param theta Estimated theta from cML.
#' @param b_vec Estimated vector of b from cML.
#' @param r_vec Estimated vector of r from cML.
#'
#' @return Standard error of theta.
#' @export
#'
#' @examples
cML_SdTheta <- function(b_exp,b_out,
                        se_exp,se_out,
                        theta,b_vec,r_vec)
{
  nonzero_ind = which(r_vec!=0)
  zero_ind = which(r_vec==0)

  VarTheta =
    1/(sum((b_vec^2/se_out^2)[zero_ind])
       -sum(
         (
           (2*b_vec*theta - b_out)^2/se_out^4*
             1/(1/se_exp^2 + theta^2/se_out^2)
         )[zero_ind]
       ))

  if(VarTheta<=0)
  {
    warning("Variance of theta is not positive,
            due to not converging to a minimum so Fisher Information
            Matrix is not positive definite. Try increasing
            number of iteration (maxit) in cML_estimate,
            or try a different start point.")
    return(NaN)
  } else {
    return(sqrt(VarTheta))
  }

}
