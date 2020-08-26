#ifndef utility_h
#define utility_h
#include <RcppArmadillo.h>

arma::vec integrated_latent(arma::mat latent, arma::vec times){
  arma::uword latent_dim = latent.n_cols;
  arma::vec latent_int(latent_dim);
  for(arma::uword i = 0; i < latent_dim; i++){
    latent_int(i) = arma::as_scalar(arma::trapz(times,
                                    latent.col(i) % latent.col(i)));
  }
  return(latent_int);
}

arma::mat integrated(arma::mat spline, arma::vec times){
  arma::uword spline_dim = spline.n_cols;
  arma::mat spline_int(spline_dim, spline_dim);
  for(arma::uword i = 0; i < spline_dim; i++){
    for(arma::uword j = 0; j < spline_dim; j++){
      spline_int(i, j) = arma::as_scalar(arma::trapz(times, spline.col(i) % spline.col(j)));
    }
  }
  return(spline_int);
}

// [[Rcpp::export]]
arma::vec integrated_latent55(arma::mat latent, arma::vec times){
  arma::uword latent_dim = latent.n_cols;
  arma::vec latent_int(latent_dim);
  for(arma::uword i = 0; i < latent_dim; i++){
    latent_int(i) = arma::as_scalar(trapz(times, latent.col(i) % latent.col(i)));
  }
  return(latent_int);
}

// [[Rcpp::export]]
arma::mat integrated55(arma::mat spline, arma::vec times){
  arma::uword spline_dim = spline.n_cols;
  arma::mat spline_int(spline_dim, spline_dim);
  for(arma::uword i = 0; i < spline_dim; i++){
    for(arma::uword j = 0; j < spline_dim; j++){
      spline_int(i, j) = arma::as_scalar(trapz(times, spline.col(i) % spline.col(j)));
    }
  }
  return(spline_int);
}

// [[Rcpp::export]]
Rcpp::List extract_eigenfn(arma::mat latent, arma::mat S, arma::mat H, arma::mat spline,
                             arma::mat spline_int_sqrt, arma::mat spline_int_sqrt_inv,
                             arma::mat spline_int, arma::vec latent_trapz, arma::uword numeig){
  arma::uword spline_dim = spline.n_cols;
  arma::uword time_size = spline.n_rows;
  arma::mat HD = arma::diagmat(H * latent_trapz);
  arma::mat cov_latent = spline_int_sqrt * (arma::diagmat(S * spline_int.diag()) + latent * HD * latent.t()) * spline_int_sqrt;
  arma::mat eigenfn;
  arma::vec eigval;
  arma::eig_sym(eigval, eigenfn, cov_latent);
  arma::mat eigen_spline(time_size, numeig);
  for(arma::uword i = 0; i < numeig; i++){
    eigen_spline.col(i) = spline * spline_int_sqrt_inv * eigenfn.col(spline_dim - 1 - i);
  }
  
  return(Rcpp::List::create(Rcpp::Named("eigenfn", eigen_spline), Rcpp::Named("eigenval", eigval)));
}

#endif
