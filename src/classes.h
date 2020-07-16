#ifndef NUM_H
#define NUM_H 
#include <RcppArmadillo.h>

arma::vec integrated_latent55(arma::mat latent, arma::vec times);
arma::mat integrated55(arma::mat spline, arma::vec times);

class EigenStruct {
private:
  int counter = 0;
  arma::mat basis_spline;
  arma::mat basis_spline_other_dim;
  arma::mat eigval_mat;
  arma::vec eigval_temp;
  arma::vec eigval_mean;
  arma::vec eigval_upper;
  arma::vec eigval_lower;
  arma::vec eigval_temp_rev;
  arma::cube eigvec_cube;
  arma::mat eigvec_temp;
  arma::mat eigvec_temp_rev;
  arma::mat eigvec_mean;
  arma::mat eigvec_upper;
  arma::mat eigvec_lower;
  arma::mat eigvec_sd;
  arma::mat eigvec_norm;
  arma::vec time_points;
  arma::mat psi;
  arma::mat basis_spline_trapz;
  arma::vec psi_trapz;
  arma::mat basis_spline_int;
  arma::mat basis_spline_int_sqrt;
  arma::mat basis_spline_int_sqrt_inv;
  arma::mat h;
  arma::mat sigma;
  std::string dimension;
  arma::mat cov_latent;
  arma::uword numeig;
  arma::uword num_draws;
  arma::uword nchains;
public:
  EigenStruct(
    arma::mat& b,
    arma::mat& b2,
    arma::vec& t,
    arma::vec& t2,
    arma::uword ne, 
    arma::uword nd,
    arma::uword nc,
    std::string z
  );
  void preprocess(arma::mat& loading, arma::mat& x, arma::mat& y);
  void get_eigenfunction(arma::mat& loading);
  void set_h_sigma(arma::mat& x, arma::mat& y);
  void normalize();
  Rcpp::List get_eigen_bands();
};

class MeanStruct {
private:
  int counter = 0;
  arma::uword dim;
  arma::uword num_draws;
  arma::uword nchains;
  arma::mat functional_spline;
  arma::mat longitudinal_spline;
  arma::mat mean_mat;
  arma::vec mean_temp;
  arma::vec mean_upper;
  arma::vec mean_mean;
  arma::vec mean_lower;
  arma::vec mean_sd;
  arma::vec mean_norm;
public:
  MeanStruct(arma::mat& bs,
             arma::mat& bt,
             arma::uword num_draws,
             arma::uword nchains);
  void update_mean(arma::mat&, arma::mat&, arma::mat&);
  void normalize();
  Rcpp::List get_mean();
};

#endif
