FisherInfoMat <- function(b_exp,b_out,
                          se_exp,se_out,
                          theta,mu_vec,v_bg)
{
  p = length(b_exp)
  ind_nonzero = which(v_bg!=0)
  num_nonzero = length(ind_nonzero)
  
  InfoMat = matrix(0,nrow = 1+p+num_nonzero,ncol = 1+p+num_nonzero)
  
  diag(InfoMat) = c(sum(mu_vec^2/se_out^2),
                    c(1/se_exp^2 + theta^2/se_out^2),
                    (1/se_out^2)[ind_nonzero]
                    )
  
  InfoMat[1,2:(1+p)] = (2*theta*mu_vec + v_bg - b_out) / se_out^2
  InfoMat[2:(1+p),1] = (2*theta*mu_vec + v_bg - b_out) / se_out^2
  
  if(num_nonzero>1)
  {
    InfoMat[1,(p+2):(1+p+num_nonzero)] = (mu_vec/se_out^2)[ind_nonzero]
    InfoMat[(p+2):(1+p+num_nonzero),1] = (mu_vec/se_out^2)[ind_nonzero]
    
    diag(InfoMat[(2:(1+p))[ind_nonzero],(p+2):(1+p+num_nonzero)]) = 
      (theta/se_out^2)[ind_nonzero]
    diag(InfoMat[(p+2):(1+p+num_nonzero),(2:(1+p))[ind_nonzero]]) = 
      (theta/se_out^2)[ind_nonzero]
  }
  
  if(num_nonzero == 1)
  {
    InfoMat[1,(p+2):(1+p+num_nonzero)] = (mu_vec/se_out^2)[ind_nonzero]
    InfoMat[(p+2):(1+p+num_nonzero),1] = (mu_vec/se_out^2)[ind_nonzero]
    
    InfoMat[(2:(1+p))[ind_nonzero],(p+2):(1+p+num_nonzero)] = 
      (theta/se_out^2)[ind_nonzero]
    InfoMat[(p+2):(1+p+num_nonzero),(2:(1+p))[ind_nonzero]] = 
      (theta/se_out^2)[ind_nonzero]
  }
  
  
  return(InfoMat)
}
