#include <RcppArmadillo.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
//' Get marginal functional covariance matrix
//' 
//' Integrate over the longitudinal dimension to obtain an nt x nt marginal 
//' covariance matrix
//' 
//' @param cov An ns*nt by ns*nt dimensional covariance matrix
//' @param ns Number of longitudinal points
//' @param nt Number of functional points
//' @export
//' @return An nt x nt marginal functional covariance matrix
//' @examples 
//' See Root/Simulation
// [[Rcpp::export]]
arma::mat getMarginalFunc(arma::mat &cov, int ns, int nt){
  arma::mat marginalT(nt, nt);
  //omp_set_num_threads(12);
  //#pragma omp parallel
  //{
    double sum;
    //#pragma omp for schedule(static)
    for(int j = 0; j < nt; j++){
      for(int k = 0; k < nt; k++){
        sum = 0;
        for(int kk = 0; kk < ns; kk++){
          //for(int jj = 0; jj < ns; jj++){
            sum = sum + cov(j + kk * nt, k + kk * nt);
          //}
        }
        marginalT(j, k) = sum / (ns*ns);
      }
    }
  //}
  return(marginalT);
}

//' Get marginal longitudinal covariance matrix
//' 
//' Integrate over the functional dimension to obtain an ns x ns marginal 
//' covariance matrix
//' 
//' @param cov An ns*nt by ns*nt dimensional covariance matrix
//' @param ns Number of longitudinal points
//' @param nt Number of functional points
//' @export
//' @return An ns x ns marginal longitudinal covariance
//' @examples 
//' See Root/Simulation
// [[Rcpp::export]]
arma::mat getMarginalLong(arma::mat cov, int ns, int nt){
  arma::mat marginalS(ns, ns);
  //omp_set_num_threads(12);
  //#pragma omp parallel for schedule(static)
  for(int j = 0; j < ns; j++){
    for(int k = 0; k < ns; k++){
      marginalS(j,k) = arma::accu(cov.submat(j * nt, k * nt, (j + 1) * 
        nt - 1, (k + 1) * nt - 1).diag()) / nt;
    }
  }
  return marginalS;
}
