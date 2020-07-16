#include <RcppArmadillo.h>
#include "classes.h"
#include "utility.h"

EigenStruct::EigenStruct(
  arma::mat& b,
  arma::mat& b2,
  arma::vec& t,
  arma::vec& t2,
  arma::uword ne,
  arma::uword nd,
  arma::uword nc,
  std::string z
) {
  
  basis_spline = b;
  basis_spline_other_dim = b2;
  time_points = t;
  eigval_mat = arma::mat(b.n_cols, nd);
  eigval_temp = arma::vec(b.n_cols);
  eigval_mean = arma::vec(b.n_cols);
  eigval_lower = arma::vec(b.n_cols);
  eigval_upper = arma::vec(b.n_cols);
  
  eigvec_cube = arma::cube(b.n_rows, ne, nd * nc);
  eigvec_temp = arma::mat(b.n_rows, ne);
  eigvec_mean = arma::mat(b.n_rows, ne, arma::fill::zeros);
  eigvec_upper = arma::mat(b.n_rows, ne);
  eigvec_lower = arma::mat(b.n_rows, ne);
  eigvec_sd = arma::mat(b.n_rows, ne);
  eigvec_norm = arma::mat(nd * nc, ne);
  basis_spline_int = integrated(b2, t2);
  basis_spline_int_sqrt = arma::sqrtmat_sympd(integrated(basis_spline, t));
  basis_spline_int_sqrt_inv = arma::inv_sympd(basis_spline_int_sqrt);
  
  numeig = ne;
  nchains = nc;
  num_draws = nd;
  dimension = z;
}

void EigenStruct::preprocess(arma::mat& loading, arma::mat& x,
                             arma::mat& y) {
  psi = basis_spline_other_dim * loading;
  psi_trapz = integrated_latent(psi, time_points);
  if (dimension == "longitudinal") {
    h = x.t();
    sigma = y.t();
  } else {
    h = x;
    sigma = y;
  }
  
}

void EigenStruct::normalize() {
  arma::vec alpha_quantile(numeig);
  for(arma::uword j = 0; j < numeig; j++) {
    for(arma::uword i = 0; i < basis_spline.n_rows; i++) {
      eigvec_sd(i, j) = arma::stddev(arma::vec(eigvec_cube.tube(i, j)));
    }
  }
  arma::vec quant(1);
  quant(0) = .95;
  for(arma::uword j = 0; j < numeig; j++) {
    for(arma::uword i = 0; i < num_draws * nchains; i++) {
      eigvec_norm(i, j) = arma::max(arma::abs(eigvec_cube.slice(i).col(j) - 
        eigvec_mean.col(j)) / eigvec_sd.col(j));
    } 
    alpha_quantile(j) = arma::as_scalar(arma::quantile(eigvec_norm.col(j),
                                        quant));
    eigvec_upper.col(j) = eigvec_mean.col(j) + alpha_quantile(j) *
      eigvec_sd.col(j);
    eigvec_lower.col(j) = eigvec_mean.col(j) - alpha_quantile(j) *
      eigvec_sd.col(j);
  }
  
  arma::vec eigval_quant(3);
  eigval_quant = {.025, .5, .975};
  arma::mat eigval_quant_mat = arma::quantile(eigval_mat, eigval_quant, 1);
  eigval_lower = eigval_quant_mat.col(0);
  eigval_mean = eigval_quant_mat.col(1);
  eigval_upper = eigval_quant_mat.col(2);
}

Rcpp::List EigenStruct::get_eigen_bands() {
  return(Rcpp::List::create(Rcpp::Named("lower", eigvec_lower),
                            Rcpp::Named("mean", eigvec_mean),
                            Rcpp::Named("upper", eigvec_upper),
                            Rcpp::Named("eigenvalue_lower", eigval_lower),
                            Rcpp::Named("eigenvalue_mean", eigval_mean),
                            Rcpp::Named("eigenvalue_upper", eigval_upper)));
}




MeanStruct::MeanStruct(arma::mat& bs, arma::mat& bt, arma::uword nd,
                       arma::uword nc) {
  functional_spline = bt;
  longitudinal_spline = bs;
  dim = bt.n_rows * bs.n_rows;
  nchains = nc;
  num_draws = nd;
  mean_mat = arma::mat(dim, nc * nd);
  mean_mean = arma::vec(dim);
  mean_temp = arma::vec(dim);
  mean_upper = arma::vec(dim);
  mean_lower = arma::vec(dim);
  mean_norm = arma::vec(nc * nd);
  mean_sd = arma::vec(dim);
}

void MeanStruct::update_mean(arma::mat& loading_functional,
                             arma::mat& loading_longitudinal,
                             arma::mat& beta_coef) {
  mean_temp = arma::vectorise(functional_spline * loading_functional * 
    beta_coef * loading_longitudinal.t() * longitudinal_spline.t());
  mean_mat.col(counter) = mean_temp;
  mean_mean = (mean_temp - mean_mean) / (counter + 1) + mean_mean;
  counter++;
}

Rcpp::List MeanStruct::get_mean() {
  return(Rcpp::List::create(Rcpp::Named("mean", mean_mean),
                            Rcpp::Named("lower", mean_lower),
                            Rcpp::Named("upper", mean_upper)));
}
void MeanStruct::normalize() {
  arma::vec alpha_quantile(1);
  alpha_quantile(0) = 0.95;
  mean_sd = arma::stddev(mean_mat, 0, 1);
  for(arma::uword i = 0; i < num_draws * nchains; i++) {
    mean_norm(i) = arma::max(arma::abs(mean_mat.col(i) - mean_mean) / mean_sd);
  }
  mean_lower = mean_mean - alpha_quantile(0) * mean_sd;
  mean_upper = mean_mean + alpha_quantile(0) * mean_sd;
}

void EigenStruct::get_eigenfunction(arma::mat& loading) {
  arma::uword chain = floor(counter / num_draws);
  arma::uword spline_dim = basis_spline.n_cols;
  arma::uword time_size = time_points.n_elem;
  arma::mat hd = arma::diagmat(h * psi_trapz);
  cov_latent = basis_spline_int_sqrt * (arma::diagmat(sigma *
    basis_spline_int.diag()) + loading * hd * loading.t()) *
    basis_spline_int_sqrt;
  arma::eig_sym(eigval_temp_rev, eigvec_temp_rev, cov_latent);
  for(arma::uword i = 0; i < numeig; i++) {
    eigvec_temp.col(i) = basis_spline * basis_spline_int_sqrt_inv * 
      eigvec_temp_rev.col(spline_dim - 1 - i);
  }
  for(arma::uword i = 0; i < spline_dim; i++) {
    eigval_temp(i) = eigval_temp_rev(spline_dim - 1 - i);
  }
  if(counter > 1) {
    for(arma::uword j = 0; j < numeig; j++) {
      if(arma::norm(eigvec_temp.col(j)*(-1)-eigvec_mean.col(j)) <
        arma::norm(eigvec_temp.col(j)-eigvec_mean.col(j))){
        eigvec_temp.col(j) = eigvec_temp.col(j) * (-1);
      }
    }
  }
  eigvec_mean = (eigvec_temp - eigvec_mean) / (counter + 1) + eigvec_mean;
  eigval_mat.col(counter + chain * num_draws) = eigval_temp;
  eigvec_cube.slice(counter + chain * num_draws) = eigvec_temp;
  counter++;
}