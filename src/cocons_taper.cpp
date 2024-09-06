#include "cocons_types.h"

using namespace Rcpp;

//' Sparse covariance function
//'
//' @param theta vector of parameters
//' @param locs a matrix with locations
//' @param locs_pred a matrix with prediction locations
//' @param x_covariates design data.frame
//' @param x_covariates_pred design data.frame at prediction locations
//' @param colindices from spam object
//' @param rowpointers from spam object
//' @param smooth_limits smooth limits
//' @return sparse covariance matrix at locs 
// [[Rcpp::export]]
NumericVector cov_rns_taper_pred(List& theta, 
                                                    NumericMatrix& locs,
                                                    NumericMatrix& locs_pred,
                                                    NumericMatrix& x_covariates, 
                                                    NumericMatrix& x_covariates_pred, 
                                                    NumericVector& colindices, 
                                                    NumericVector& rowpointers, 
                                                    NumericVector& smooth_limits){
  
  const double epsilon = std::numeric_limits<double>::epsilon();
  
  const int m_dim = locs.nrow();
  const int pred_dim = locs_pred.nrow();
  
  const NumericVector std_dev_vector = theta["std.dev"];
  const NumericVector scale = theta["scale"];
  //const NumericVector aniso = theta["aniso"];
  //const NumericVector tilt = theta["tilt"];
  const NumericVector smooth = theta["smooth"];
  const NumericVector nugget = theta["nugget"];
  
  NumericVector dif_s (2);
  
  double smooth_s_Q_ij, smtns, prefactor, global_range;
  
  std::vector<double> scale_vector(m_dim); 
  std::vector<double> sigma_vector(m_dim); 
  std::vector<double> smooth_vector(m_dim); 
  
  std::vector<double> scale_vector_pred(pred_dim); 
  std::vector<double> sigma_vector_pred(pred_dim); 
  std::vector<double> smooth_vector_pred(pred_dim); 
  
  const int taper_dim = colindices.length();
  
  Rcpp::NumericVector dist_vector(taper_dim);
  
  // double global_range = 1 / std::exp(-2 * scale[0]);
  
  for(int ww = 0; ww < m_dim; ww++){
    
    scale_vector[ww] = Pexpfma_new(2 * scale, x_covariates(ww,_));
    sigma_vector[ww] = Pexpfma_new(0.5 * std_dev_vector, x_covariates(ww,_));
    smooth_vector[ww] = std::sqrt(Pexpfma_new_smoothness(smooth, x_covariates(ww,_), smooth_limits[0], smooth_limits[1]));
    //smooth_vector[ww] = std::sqrt((1 / std::exp(- smooth[0])) + Pexpfma_new_smoothness_two(smooth, x_covariates(ww,_), smooth_limits[0], smooth_limits[1]));
    
  }
  
  for(int ww = 0; ww < pred_dim; ww++){
    
    scale_vector_pred[ww] = Pexpfma_new(2 * scale, x_covariates_pred(ww,_));
    sigma_vector_pred[ww] = Pexpfma_new(0.5 * std_dev_vector, x_covariates_pred(ww,_));
    smooth_vector_pred[ww] = std::sqrt(Pexpfma_new_smoothness(smooth, x_covariates_pred(ww,_), smooth_limits[0], smooth_limits[1]));

  }
  
  colindices = colindices - 1;
  rowpointers = rowpointers - 1;
  
  int jj;
  // int to_to;
  int acumulating = 0;
  
  for(int ii = 0; ii < pred_dim; ii++){
    
    //to_to = rowpointers(ii);
    
    for(int ww = rowpointers(ii); ww < (rowpointers(ii + 1)); ww++){
      
      jj = colindices(ww);
      
      if(locs_pred(ii, 0) == locs(jj, 0) && locs_pred(ii, 1) == locs(jj, 1)){
        
        dist_vector(acumulating) = sigma_vector_pred[ii] * sigma_vector_pred[ii] + Pexpfma_new(nugget, x_covariates_pred(ii,_));
        
        acumulating = acumulating + 1;
        continue;
        
      } else {
        
        smtns = smooth_vector_pred[ii] * smooth_vector[jj];
        
        prefactor = (2 * std::pow(scale_vector_pred[ii], 0.5) * std::pow(scale_vector[jj], 0.5)) / ( scale_vector_pred[ii] + scale_vector[jj]);
        
        global_range = (scale_vector_pred[ii] + scale_vector[jj]) / 2;
        
        dif_s = locs_pred(ii,_) - locs(jj,_);
        
        smooth_s_Q_ij = std::sqrt(8 * smtns) * std::sqrt(std::pow(dif_s(0), 2) + std::pow(dif_s(1), 2)) / std::sqrt(global_range);
        
        // smaller than eps?
        if(smooth_s_Q_ij <= epsilon){
          dist_vector(acumulating) = sigma_vector_pred[ii] * sigma_vector_pred[ii] + Pexpfma_new(nugget, x_covariates_pred(ii,_));
          acumulating = acumulating + 1;
          continue;
          
        } else{
          
          if(smooth_s_Q_ij < 706.0){
            
            dist_vector(acumulating) = prefactor * std::pow(2.0 , -(smtns - 1)) / std::tgamma(smtns) * 
              std::pow(smooth_s_Q_ij, smtns) * boost::math::cyl_bessel_k(smtns, smooth_s_Q_ij) * 
              sigma_vector_pred[ii] * sigma_vector[jj];
            
          } else {
            dist_vector(acumulating) = 0.0;
          }
        }
        
      }
      
      acumulating = acumulating + 1;
      
    }
  }
  
  return dist_vector;
}

//' Sparse covariance function
//'
//' @param theta vector of parameters
//' @param locs a matrix with locations
//' @param x_covariates design data.frame
//' @param colindices from spam object
//' @param rowpointers from spam object
//' @param smooth_limits smooth limits
//' @return sparse covariance matrix between locs and pred_locs
// [[Rcpp::export]]
NumericVector cov_rns_taper(List& theta, 
                                            NumericMatrix& locs, 
                                            NumericMatrix& x_covariates, 
                                            NumericVector& colindices, 
                                            NumericVector& rowpointers, 
                                            NumericVector& smooth_limits){
  
  const double epsilon = std::numeric_limits<double>::epsilon();
  const int m_dim = locs.nrow();
  const NumericVector std_dev_vector = theta["std.dev"];
  const NumericVector scale = theta["scale"];
  const NumericVector aniso = theta["aniso"];
  const NumericVector tilt = theta["tilt"];
  const NumericVector smooth = theta["smooth"];
  const NumericVector nugget = theta["nugget"];
  
  NumericVector dif_s (2);
  
  double smooth_s_Q_ij, smtns, prefactor, global_range;
  
  std::vector<double> range_vector(m_dim); 
  std::vector<double> sigma_vector(m_dim); 
  std::vector<double> smooth_vector(m_dim); 
  std::vector<double> nugget_vector(m_dim);
  
  const int taper_dim = colindices.length();
  
  Rcpp::NumericVector dist_vector(taper_dim);
  
  // global parameters
  
  // double global_range = 1 / exp(-2 * scale[0]);
  
  for(int ww = 0; ww < m_dim; ww++){
    
    range_vector[ww] = Pexpfma_new(2 * scale, x_covariates(ww,_));
    sigma_vector[ww] = Pexpfma_new(0.5 * std_dev_vector, x_covariates(ww,_));
    smooth_vector[ww] = std::sqrt(Pexpfma_new_smoothness(smooth, x_covariates(ww,_), smooth_limits[0], smooth_limits[1]));
  }
  
  colindices = colindices - 1;
  rowpointers = rowpointers - 1;
  
  int jj;
  // int to_to;
  int acumulating = 0;
  
  for(int ii = 0; ii < m_dim; ii++){
    
    for(int ww = rowpointers(ii); ww < (rowpointers(ii + 1)); ww++){
      
      jj = colindices(ww);
      
      if(ii == jj){
        
        dist_vector(acumulating) = Pexpfma_new(std_dev_vector, x_covariates(ii,_)) + Pexpfma_new(nugget, x_covariates(ii,_));
        
        acumulating = acumulating + 1;
        continue;
        
      } else {
        
        smtns = smooth_vector[ii] * smooth_vector[jj];
        
        prefactor = (2 * std::pow(range_vector[ii],0.5) * std::pow(range_vector[jj],0.5))/ ( range_vector[ii] + range_vector[jj]);
        
        global_range = (range_vector[ii] + range_vector[jj]) / 2;
        
        dif_s = locs(ii,_) - locs(jj,_);
        
        smooth_s_Q_ij = std::sqrt(8 * smtns) * std::sqrt(std::pow(dif_s(0), 2) + 
          std::pow(dif_s(1), 2)) / std::sqrt(global_range);
        
        // smaller than eps?
        if(smooth_s_Q_ij <= epsilon){
          dist_vector(acumulating) = Pexpfma_new(std_dev_vector, x_covariates(ii,_)) + Pexpfma_new(nugget, x_covariates(ii,_));
          acumulating = acumulating + 1;
          continue;
          
        } else{
          
          if(smooth_s_Q_ij < 706.0){
            
            dist_vector(acumulating) = prefactor *  std::pow(2.0 , -(smtns - 1)) / std::tgamma(smtns) * 
              std::pow(smooth_s_Q_ij, smtns) * boost::math::cyl_bessel_k(smtns, smooth_s_Q_ij) * 
              sigma_vector[ii] * sigma_vector[jj];
            
          } else {
            dist_vector(acumulating) = 0.0;
          }
        }
        
      }
      
      acumulating = acumulating + 1;
      
    }
  }
  
  return dist_vector;
}
