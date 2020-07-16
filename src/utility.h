#include <RcppArmadillo.h>

arma::vec integrated_latent(arma::mat latent, arma::vec times);
arma::mat integrated(arma::mat spline, arma::vec times);
void extract_eigenfn(arma::mat latent, arma::mat S, arma::mat H,
                     arma::mat spline, arma::mat spline_int_sqrt,
                     arma::mat spline_int_sqrt_inv, arma::mat spline_int,
                     arma::vec latent_trapz, arma::uword numeig);

arma::mat integrated55(arma::mat spline, arma::vec times);
arma::vec integrated_latent55(arma::mat latent, arma::vec times);
Rcpp::List extract_eigenfn55(arma::mat latent, arma::mat S, arma::mat H,
                             arma::mat spline, arma::mat spline_int_sqrt,
                             arma::mat spline_int_sqrt_inv, 
                             arma::mat spline_int, arma::vec latent_trapz,
                             arma::uword numeig);