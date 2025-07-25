#include "cocons_types.h"

using namespace Rcpp;

//' smoothed-L1 penalization over the covariate-driven parameters
//'
//' @param x a vector of parameters to penalize
//' @param lambda penalization hyperparameter
//' @param alpha hyperparameter to approximate L1 penalization. By default 1e6 provides a decent approximation to L1, while also being a differentiable function
//' @return penalization
// [[Rcpp::export]]
double sumsmoothlone(NumericVector& x,
                      double lambda,
                      double alpha = 1e6){
   
   double sum = 0;
   
   for(int ww = 0; ww < x.length(); ww++){
     
     if(std::abs(x[ww]) > 1e-4){
       sum = sum + std::abs(x[ww]);
     } else{
       sum = sum + std::pow(alpha,-1) * (std::log(1 + std::exp(-alpha * x[ww])) + std::log(1 + exp(alpha * x[ww])));
     }
     
   }
   
   return lambda * sum;
   
}

//' Dense covariance function (difference parameterization)
//'
//' @param theta vector of parameters
//' @param locs a matrix with locations
//' @param x_covariates design data.frame
//' @param smooth_limits smooth limits
//' @return dense covariance matrix
// [[Rcpp::export]]
NumericMatrix cov_rns(List& theta, 
                      NumericMatrix& locs, 
                      NumericMatrix& x_covariates, 
                      NumericVector& smooth_limits){
   
   const double epsilon = std::numeric_limits<double>::epsilon();
   const int m_dim = locs.nrow();
   const NumericVector std_dev_vector = theta["std.dev"];
   NumericVector scale = theta["scale"];
   NumericVector scale_je = clone(scale); 
   
   const NumericVector aniso = theta["aniso"];
   const NumericVector tilt = theta["tilt"];
   const NumericVector smooth = theta["smooth"];
   const NumericVector nugget = theta["nugget"];
   
   NumericMatrix dist_matrix(m_dim, m_dim);
   
   NumericVector dif_s (2);
   
   double sigma11_ij, sigma22_ij, sigma12_ij, det_ij, smooth_s_Q_ij, smtns;

   double global_range = 1 / std::exp(- 2 * scale[0]);
   
   scale_je[0] = 0.0;
   
   const NumericVector sqrt_vector = 2 * scale_je + aniso;
   
   // initialize to prevent redundancy later on:
   
   std::vector<double> tilt_vector(m_dim);
   std::vector<double> range_det_vector(m_dim);
   std::vector<double> aniso_det_vector(m_dim);
   std::vector<double> dets_vector(m_dim);
   std::vector<double> sigma_vector(m_dim);
   std::vector<double> nugget_vector(m_dim);
   
   bool q_smooth_fix = allzeroelements(smooth);
   
   int smooth_switch = 0;
   
   double smooth_value;
   
   std::vector<double> smooth_vector(m_dim);
   
   if(q_smooth_fix && (smooth_limits[0] == smooth_limits[1])){
     
     smooth_value = smooth_limits[0];
     smooth_switch =  mapSmoothValue(smooth_value);
     
   } else{
     
     for(int ww = 0; ww < m_dim; ww++){
       smooth_vector[ww] = std::sqrt(Pexpfma_new_smoothness(smooth, x_covariates(ww,_), smooth_limits[0], smooth_limits[1]));
     }
     
   }

   for(int ww = 0; ww < m_dim; ww++){
     
     tilt_vector[ww] = newinvlogitfma(tilt, x_covariates(ww,_));
     range_det_vector[ww] = Pexpfma_new(2 * scale_je, x_covariates(ww,_));
     aniso_det_vector[ww] = Pexpfma_new(aniso, x_covariates(ww,_));
     dets_vector[ww] = Pexpfma_new(sqrt_vector, x_covariates(ww,_));
     sigma_vector[ww] = Pexpfma_new(0.5 * std_dev_vector, x_covariates(ww,_));
     nugget_vector[ww] = Pexpfma_new(nugget, x_covariates(ww,_));
     
   }
   
   // diag
   for(int ii = 0; ii < m_dim; ii++){
     dist_matrix(ii,ii) = Pexpfma_new(std_dev_vector, x_covariates(ii,_)) + nugget_vector[ii];
   }
   
   switch(smooth_switch){
   
   // smoothness 0.5
   case 1:
     
     for(int ii = 0; ii < m_dim; ii++){
       for(int jj = ii + 1; jj < m_dim; jj++){
         
         sigma11_ij = (range_det_vector[ii] + range_det_vector[jj]) * 0.5; 
         sigma22_ij = kahan(range_det_vector[ii], aniso_det_vector[ii] * aniso_det_vector[ii], 
                            - range_det_vector[jj], aniso_det_vector[jj] * aniso_det_vector[jj]) * 0.5; 
         
         sigma12_ij = kahan(range_det_vector[ii] * aniso_det_vector[ii],
                            std::cos(tilt_vector[ii]),
                            - 1 * range_det_vector[jj] * aniso_det_vector[jj], 
                                                                         std::cos(tilt_vector[jj])
         ) * 0.5;
         
         det_ij = kahan(sigma11_ij, sigma22_ij, sigma12_ij, sigma12_ij);
         
         dif_s = locs(ii,_) - locs(jj,_);
         
         smooth_s_Q_ij = std::sqrt(8 *  smooth_value / (global_range * det_ij )) * 
           std::sqrt(std::fma(kahan(sigma22_ij, 
                                    dif_s(0) * dif_s(0), 
                                    - sigma11_ij, 
                                    dif_s(1) * dif_s(1)), 1, 
                                    -2 * sigma12_ij * dif_s(0) * dif_s(1)));
           
         // smaller than eps?
         if(smooth_s_Q_ij <= epsilon){
           dist_matrix(ii,jj) = dist_matrix(jj,ii) = Pexpfma_new(std_dev_vector, x_covariates(ii,_)) + nugget_vector[ii];
           continue;
           
         } else{
           
           dist_matrix(ii,jj) = dist_matrix(jj,ii) = std::exp(-smooth_s_Q_ij) * 
             sigma_vector[ii] * sigma_vector[jj] * 
             std::sqrt(dets_vector[ii] * std::sin(tilt_vector[ii]) * 
             dets_vector[jj] * std::sin(tilt_vector[jj])) / std::sqrt(det_ij);
           
           continue;
         }
       }
     }
     
     break;
     
   // smoothness 1.5
   case 2:
     
     for(int ii = 0; ii < m_dim; ii++){
       for(int jj = ii + 1; jj < m_dim; jj++){
         
         sigma11_ij = (range_det_vector[ii] + range_det_vector[jj]) * 0.5; 
         sigma22_ij = kahan(range_det_vector[ii], aniso_det_vector[ii] * aniso_det_vector[ii], 
                            - range_det_vector[jj], aniso_det_vector[jj] * aniso_det_vector[jj]) * 0.5; 
         
         sigma12_ij = kahan(range_det_vector[ii] * aniso_det_vector[ii],
                            std::cos(tilt_vector[ii]),
                            - 1 * range_det_vector[jj] * aniso_det_vector[jj], 
                                                                         std::cos(tilt_vector[jj])
         ) * 0.5;
         
         det_ij = kahan(sigma11_ij, sigma22_ij, sigma12_ij, sigma12_ij);
         
         dif_s = locs(ii,_) - locs(jj,_);
         
         smooth_s_Q_ij = std::sqrt(8 *  smooth_value / (global_range * det_ij )) * 
           std::sqrt(std::fma(kahan(sigma22_ij, 
                                    dif_s(0) * dif_s(0), 
                                    - sigma11_ij, 
                                    dif_s(1) * dif_s(1)), 1, 
                                    -2 * sigma12_ij * dif_s(0) * dif_s(1)));
         
         // smaller than eps?
         if(smooth_s_Q_ij <= epsilon){
           dist_matrix(ii,jj) = dist_matrix(jj,ii) = Pexpfma_new(std_dev_vector, x_covariates(ii,_)) + nugget_vector[ii];
           continue;
           
         } else{
           
           dist_matrix(ii,jj) = dist_matrix(jj,ii) = (1 + smooth_s_Q_ij) * std::exp(- smooth_s_Q_ij) * 
             sigma_vector[ii] * sigma_vector[jj] * 
             std::sqrt(dets_vector[ii] * std::sin(tilt_vector[ii]) * 
             dets_vector[jj] * std::sin(tilt_vector[jj])) / std::sqrt(det_ij);
           
           continue;
         }
       }
     }
     
     break;
     
   // smoothness 2.5
   case 3:
     
     for(int ii = 0; ii < m_dim; ii++){
       for(int jj = ii + 1 ; jj < m_dim; jj++){
         
         sigma11_ij = (range_det_vector[ii] + range_det_vector[jj]) * 0.5; 
         sigma22_ij = kahan(range_det_vector[ii], aniso_det_vector[ii] * aniso_det_vector[ii], 
                            - range_det_vector[jj], aniso_det_vector[jj] * aniso_det_vector[jj]) * 0.5; 
         
         sigma12_ij = kahan(range_det_vector[ii] * aniso_det_vector[ii],
                            std::cos(tilt_vector[ii]),
                            - 1 * range_det_vector[jj] * aniso_det_vector[jj], 
                                                                         std::cos(tilt_vector[jj])
         ) * 0.5;
         
         det_ij = kahan(sigma11_ij, sigma22_ij, sigma12_ij, sigma12_ij);
         
         dif_s = locs(ii,_) - locs(jj,_);
         
         smooth_s_Q_ij = std::sqrt(8 *  smooth_value / (global_range * det_ij )) * 
           std::sqrt(std::fma(kahan(sigma22_ij, 
                                    dif_s(0) * dif_s(0), 
                                    - sigma11_ij, 
                                    dif_s(1) * dif_s(1)), 1, 
                                    -2 * sigma12_ij * dif_s(0) * dif_s(1)));
         
         // smaller than eps?
         if(smooth_s_Q_ij <= epsilon){
           dist_matrix(ii,jj) = dist_matrix(jj,ii) = Pexpfma_new(std_dev_vector, x_covariates(ii,_)) + nugget_vector[ii];
           continue;
           
         } else{
           
           dist_matrix(ii,jj) = dist_matrix(jj,ii) = (1 + smooth_s_Q_ij + smooth_s_Q_ij * smooth_s_Q_ij / 3) * std::exp(- smooth_s_Q_ij) * 
             sigma_vector[ii] * sigma_vector[jj] * 
             std::sqrt(dets_vector[ii] * std::sin(tilt_vector[ii]) * 
             dets_vector[jj] * std::sin(tilt_vector[jj])) / std::sqrt(det_ij);
           
           continue;
         }
       }
     }
     
     break;
   
   // o.w.
   default:
     
     for(int ii = 0; ii < m_dim; ii++){
       for(int jj = ii + 1; jj < m_dim; jj++){
           
           sigma11_ij = (range_det_vector[ii] + range_det_vector[jj]) * 0.5; 
           sigma22_ij = kahan(range_det_vector[ii], aniso_det_vector[ii] * aniso_det_vector[ii], 
                              - range_det_vector[jj], aniso_det_vector[jj] * aniso_det_vector[jj]) * 0.5; 
           
           sigma12_ij = kahan(range_det_vector[ii] * aniso_det_vector[ii],
                              std::cos(tilt_vector[ii]),
                              - 1 * range_det_vector[jj] * aniso_det_vector[jj], 
                                                                           std::cos(tilt_vector[jj])
           ) * 0.5;
           
           det_ij = kahan(sigma11_ij, sigma22_ij, sigma12_ij, sigma12_ij);
           
           dif_s = locs(ii,_) - locs(jj,_);
           
           smtns = smooth_vector[ii] * smooth_vector[jj];
           
           smooth_s_Q_ij = std::sqrt(8 *  smtns / (global_range * det_ij )) * 
             std::sqrt(std::fma(kahan(sigma22_ij, 
                                      dif_s(0) * dif_s(0), 
                                      - sigma11_ij, 
                                      dif_s(1) * dif_s(1)), 1, 
                                      -2 * sigma12_ij * dif_s(0) * dif_s(1)));
           
           // smaller than eps?
           if(smooth_s_Q_ij <= epsilon){
             dist_matrix(ii,jj) = dist_matrix(jj,ii) = Pexpfma_new(std_dev_vector, x_covariates(ii,_)) + nugget_vector[ii];
             continue;
             
           } else{
             
             
             if(smooth_s_Q_ij < 706.0){
               
               dist_matrix(ii,jj) = dist_matrix(jj,ii) = std::pow(2.0 , -(smtns - 1)) / std::tgamma(smtns) * 
                 std::pow(smooth_s_Q_ij, smtns) * boost::math::cyl_bessel_k(smtns, smooth_s_Q_ij) * 
                 sigma_vector[ii] * sigma_vector[jj] * 
                 std::sqrt(dets_vector[ii] * std::sin(tilt_vector[ii]) * 
                 dets_vector[jj] * std::sin(tilt_vector[jj])) / std::sqrt(det_ij);
               
             } else {
               
               dist_matrix(ii,jj) = dist_matrix(jj,ii) = dist_matrix(jj,ii) = std::pow(2.0 , -(smtns - 1)) / std::tgamma(smtns) * 
                 std::pow(smooth_s_Q_ij, smtns) * std::sqrt(M_PI / (2.0 * smooth_s_Q_ij)) * std::exp(-smooth_s_Q_ij) * 
                 sigma_vector[ii] * sigma_vector[jj] * 
                 std::sqrt(dets_vector[ii] * std::sin(tilt_vector[ii]) * 
                 dets_vector[jj] * std::sin(tilt_vector[jj])) / std::sqrt(det_ij);
               
             }
             
             continue;
             
           }
         }
       }
       
     break;
       
   }
   
   return dist_matrix;
   
 }
 

//' Dense covariance function
//'
//' @param theta vector of parameters
//' @param locs a matrix with locations
//' @param locs_pred a matrix with prediction locations
//' @param x_covariates design data.frame
//' @param x_covariates_pred design data.frame at prediction locations
//' @param smooth_limits smooth limits
//' @return dense covariance matrix
// [[Rcpp::export]]
NumericMatrix cov_rns_pred(List& theta, NumericMatrix& locs, 
                           NumericMatrix& locs_pred,
                           NumericMatrix& x_covariates,
                           NumericMatrix& x_covariates_pred,
                           NumericVector& smooth_limits){
  
  const double epsilon = std::numeric_limits<double>::epsilon();
  const int m_dim = locs.nrow();
  const int pred_dim = locs_pred.nrow();
  const NumericVector std_dev_vector = theta["std.dev"];
  const NumericVector scale = theta["scale"];
  const NumericVector aniso = theta["aniso"];
  const NumericVector tilt = theta["tilt"];
  const NumericVector smooth = theta["smooth"];
  const NumericVector nugget = theta["nugget"];
  
  NumericVector scale_je = clone(scale); 
  double global_range = 1 / std::exp(- 2 * scale[0]);
  scale_je[0] = 0.0;
  
  NumericMatrix dist_matrix(pred_dim, m_dim);
  
  NumericVector dif_s (2);
  
  double sigma11_ij, sigma22_ij, sigma12_ij, det_ij, smooth_s_Q_ij,smtns;
  // const double fix_gamma = std::tgamma(smtns);
  // const double some_cte = std::pow(2.0 , -(smtns - 1));
  
  const NumericVector sqrt_vector = 2 * scale_je + aniso;
  
  // initialize to prevent redundancy later on:
  
  std::vector<double> tilt_vector(m_dim); 
  std::vector<double> range_det_vector(m_dim); 
  std::vector<double> aniso_det_vector(m_dim); 
  std::vector<double> dets_vector(m_dim); 
  std::vector<double> sigma_vector(m_dim); 
  std::vector<double> smooth_vector(m_dim); 
  std::vector<double> nugget_vector(m_dim);
  
  for(int ww = 0; ww < m_dim; ww++){
    
    tilt_vector[ww] = newinvlogitfma(tilt, x_covariates(ww,_));
    range_det_vector[ww] = Pexpfma_new(2 * scale_je, x_covariates(ww,_));
    aniso_det_vector[ww] = Pexpfma_new(aniso, x_covariates(ww,_));
    dets_vector[ww] = Pexpfma_new(sqrt_vector, x_covariates(ww,_));
    sigma_vector[ww] = Pexpfma_new(0.5 * std_dev_vector, x_covariates(ww,_));
    smooth_vector[ww] = std::sqrt(Pexpfma_new_smoothness(smooth, x_covariates(ww,_), smooth_limits[0], smooth_limits[1]));
    nugget_vector[ww] = Pexpfma_new(nugget, x_covariates(ww,_));
    
  }
  
  std::vector<double> p_tilt_vector(pred_dim); 
  std::vector<double> p_range_det_vector(pred_dim); 
  std::vector<double> p_aniso_det_vector(pred_dim); 
  std::vector<double> p_dets_vector(pred_dim); 
  std::vector<double> p_sigma_vector(pred_dim); 
  std::vector<double> p_smooth_vector(pred_dim); 
  std::vector<double> p_nugget_vector(pred_dim);
  
  for(int ww = 0; ww < pred_dim; ww++){
    
    p_tilt_vector[ww] = newinvlogitfma(tilt, x_covariates_pred(ww,_));
    p_range_det_vector[ww] = Pexpfma_new(2 * scale_je, x_covariates_pred(ww,_));
    p_aniso_det_vector[ww] = Pexpfma_new(aniso, x_covariates_pred(ww,_));
    p_dets_vector[ww] = Pexpfma_new(sqrt_vector, x_covariates_pred(ww,_));
    p_sigma_vector[ww] = Pexpfma_new(0.5 * std_dev_vector, x_covariates_pred(ww,_));
    p_smooth_vector[ww] = std::sqrt(Pexpfma_new_smoothness(smooth, x_covariates_pred(ww,_), smooth_limits[0], smooth_limits[1]));

    p_nugget_vector[ww] = Pexpfma_new(nugget, x_covariates_pred(ww,_));
    
  }
  
  for(int ii = 0; ii < pred_dim; ii++){
    for(int jj = 0; jj < m_dim; jj++){
      
      if(locs_pred(ii,0) == locs(jj,0) && locs_pred(ii,1) == locs(jj,1)){
        
        dist_matrix(ii,jj) = Pexpfma_new(std_dev_vector, x_covariates_pred(ii,_)) + p_nugget_vector[ii];
        
        continue;
        
      } else {
        
        sigma11_ij = (p_range_det_vector[ii] + range_det_vector[jj]) * 0.5; 
        sigma22_ij = kahan(p_range_det_vector[ii], p_aniso_det_vector[ii] * p_aniso_det_vector[ii], - range_det_vector[jj], aniso_det_vector[jj] * aniso_det_vector[jj]) * 0.5; 
        
        sigma12_ij = kahan(p_range_det_vector[ii] * p_aniso_det_vector[ii],
                           std::cos(p_tilt_vector[ii]),
                           - 1 * range_det_vector[jj] * aniso_det_vector[jj], 
                                                                        std::cos(tilt_vector[jj])
        ) * 0.5;
        
        det_ij = kahan(sigma11_ij, sigma22_ij, sigma12_ij, sigma12_ij);
        
        dif_s = locs_pred(ii,_) - locs(jj,_);
        
        smtns = p_smooth_vector[ii] * smooth_vector[jj];
        
        smooth_s_Q_ij = std::sqrt(8 * smtns / (global_range * det_ij)) * 
          std::sqrt(std::fma(kahan(sigma22_ij, 
                                   dif_s(0) * dif_s(0), 
                                   - sigma11_ij, 
                                   dif_s(1) * dif_s(1)),1, -2 * sigma12_ij * dif_s(0) * dif_s(1)));
        
        // smaller than eps?
        if(smooth_s_Q_ij <= epsilon){
          dist_matrix(ii,jj) = Pexpfma_new(std_dev_vector, x_covariates_pred(ii,_)) + p_nugget_vector[ii];
          continue;
          
        } else{
          
          
          if(smooth_s_Q_ij < 706.0){
            
            dist_matrix(ii,jj) = std::pow(2.0 , -(smtns - 1)) / std::tgamma(smtns) * 
              std::pow(smooth_s_Q_ij, smtns) * boost::math::cyl_bessel_k(smtns, smooth_s_Q_ij) * 
              p_sigma_vector[ii] * sigma_vector[jj] * 
              std::sqrt(p_dets_vector[ii] * std::sin(p_tilt_vector[ii]) * 
              dets_vector[jj] * std::sin(tilt_vector[jj])) / std::sqrt(det_ij);
            
          } else {
            
            dist_matrix(ii,jj) = std::pow(2.0 , -(smtns - 1)) / std::tgamma(smtns) * 
              std::pow(smooth_s_Q_ij, smtns) * std::sqrt(M_PI / (2.0 * smooth_s_Q_ij)) * std::exp(-smooth_s_Q_ij) * 
              p_sigma_vector[ii] * sigma_vector[jj] * 
              std::sqrt(p_dets_vector[ii] * std::sin(p_tilt_vector[ii]) * 
              dets_vector[jj] * std::sin(tilt_vector[jj])) / std::sqrt(det_ij);
            
          }
          
        }
      }
    }
  }
  
  return dist_matrix;
}

//' Dense covariance function (classic parameterization)
//'
//' @param theta vector of parameters
//' @param locs a matrix with locations
//' @param x_covariates design data.frame
//' @return dense covariance matrix with classic parameterization
// [[Rcpp::export]]
NumericMatrix cov_rns_classic(List& theta, NumericMatrix& locs, NumericMatrix& x_covariates){
  
  const double epsilon = std::numeric_limits<double>::epsilon();
  const int m_dim = locs.nrow();
  const NumericVector std_dev_vector = theta["std.dev"];
  NumericVector scale = theta["scale"];
  NumericVector scale_je = clone(scale); 
  
  const NumericVector aniso = theta["aniso"];
  const NumericVector tilt = theta["tilt"];
  const NumericVector smooth = theta["smooth"];
  const NumericVector nugget = theta["nugget"];
  
  NumericMatrix dist_matrix(m_dim);
  
  NumericVector dif_s (2);
  
  double sigma11_ij, sigma22_ij, sigma12_ij, det_ij, smooth_s_Q_ij, smtns;
  // const double fix_gamma = std::tgamma(smtns);
  // const double some_cte = std::pow(2.0 , -(smtns - 1));
  
  double global_range = 1 / std::exp(- 2 * scale[0]);
  
  scale_je[0] = 0.0;
  
  const NumericVector sqrt_vector = 2 * scale_je + aniso;
  
  // initialize to prevent redundancy later on:
  
  std::vector<double> tilt_vector(m_dim);
  std::vector<double> range_det_vector(m_dim);
  std::vector<double> aniso_det_vector(m_dim);
  std::vector<double> dets_vector(m_dim);
  std::vector<double> sigma_vector(m_dim);
  std::vector<double> smooth_vector(m_dim);
  std::vector<double> nugget_vector(m_dim);
  
  for(int ww = 0; ww < m_dim; ww++){
    
    tilt_vector[ww] = newinvlogitfma(tilt, x_covariates(ww,_));
    range_det_vector[ww] = Pexpfma_new(2 * scale_je, x_covariates(ww,_));
    aniso_det_vector[ww] = Pexpfma_new(aniso, x_covariates(ww,_));
    dets_vector[ww] = Pexpfma_new(sqrt_vector, x_covariates(ww,_));
    sigma_vector[ww] = Pexpfma_new(0.5 * std_dev_vector, x_covariates(ww,_));
    smooth_vector[ww] = Pexpfma_new(smooth, x_covariates(ww,_));
    nugget_vector[ww] = Pexpfma_new(nugget, x_covariates(ww,_));
    
  }
  
  for(int ii = 0; ii < m_dim; ii++){
    for(int jj = ii; jj < m_dim; jj++){
      
      if(ii == jj){
        
        dist_matrix(ii,jj) = Pexpfma_new(std_dev_vector, x_covariates(ii,_)) + nugget_vector[ii];
        
        continue;
        
      } else {
        
        sigma11_ij = (range_det_vector[ii] + range_det_vector[jj]) * 0.5; 
        sigma22_ij = kahan(range_det_vector[ii], aniso_det_vector[ii] * aniso_det_vector[ii], 
                           - range_det_vector[jj], aniso_det_vector[jj] * aniso_det_vector[jj]) * 0.5; 
        
        sigma12_ij = kahan(range_det_vector[ii] * aniso_det_vector[ii],
                           std::cos(tilt_vector[ii]),
                           - 1 * range_det_vector[jj] * aniso_det_vector[jj], 
                                                                        std::cos(tilt_vector[jj])
        ) * 0.5;
        
        det_ij = kahan(sigma11_ij, sigma22_ij, sigma12_ij, sigma12_ij);
        
        dif_s = locs(ii,_) - locs(jj,_);
        
        smtns = (smooth_vector[ii] + smooth_vector[jj])/2;
        
        smooth_s_Q_ij = std::sqrt(8 *  smtns / (global_range * det_ij )) * 
          std::sqrt(std::fma(kahan(sigma22_ij, 
                                   dif_s(0) * dif_s(0), 
                                   - sigma11_ij, 
                                   dif_s(1) * dif_s(1)), 1, 
                                   -2 * sigma12_ij * dif_s(0) * dif_s(1)));
        
        // smaller than eps?
        if(smooth_s_Q_ij <= epsilon){
          dist_matrix(ii,jj) = dist_matrix(jj,ii) = Pexpfma_new(std_dev_vector, x_covariates(ii,_)) + nugget_vector[ii];
          continue;
          
        } else{
          
          if(smooth_s_Q_ij < 706.0){
            
            dist_matrix(ii,jj) = dist_matrix(jj,ii) = std::pow(2.0 , -(smtns - 1)) / std::tgamma(smtns) * 
              std::pow(smooth_s_Q_ij, smtns) * boost::math::cyl_bessel_k(smtns, smooth_s_Q_ij) * 
              sigma_vector[ii] * sigma_vector[jj] * 
              std::sqrt(dets_vector[ii] * std::sin(tilt_vector[ii]) * 
              dets_vector[jj] * std::sin(tilt_vector[jj])) / std::sqrt(det_ij);
            
          } else {
            
            dist_matrix(ii,jj) = dist_matrix(jj,ii) = dist_matrix(jj,ii) = std::pow(2.0 , -(smtns - 1)) / std::tgamma(smtns) * 
              std::pow(smooth_s_Q_ij, smtns) * std::sqrt(M_PI / (2.0 * smooth_s_Q_ij)) * std::exp(-smooth_s_Q_ij) * 
              sigma_vector[ii] * sigma_vector[jj] * 
              std::sqrt(dets_vector[ii] * std::sin(tilt_vector[ii]) * 
              dets_vector[jj] * std::sin(tilt_vector[jj])) / std::sqrt(det_ij);
            
          }
          
        }
      }
    }
  }
  
  return dist_matrix;
}
