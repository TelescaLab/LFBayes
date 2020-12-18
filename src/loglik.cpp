#include <RcppArmadillo.h>
//' Calculate log-likelihood
//'
//' Calculates log-likelihood after obtaining draws from the joint posterior
//' distribution. Use after running run_mcmc.
//'
//' @param y A list of of length n containing responses
//' @param X An n x p design matrix
//' @param Bs Basis matrix for longitudinal direction
//' @param Bt Basis matrix for functional direction
//' @param missing A list of length n containing missing indices for each
//' response
//' @param Theta Posterior draws of Theta, arranges in a cube
//' @param Varphi Posterior draws of Varphi, vector format
//' @param iter Number of total samples
//' @param burnin Number of samples to use as burnin
//' @export loglik
//' @return A Matrix of size (iter - burnin) x number of observed time points
//' over all subjects containing log-likelihood values
// [[Rcpp::export]]
arma::mat loglik(arma::vec y, arma::mat X, arma::mat Bs, arma::mat Bt,
                 arma::field<arma::vec> missing,
                 arma::field<arma::cube> Theta,
                 arma::field<arma::vec> Varphi,
                 arma::uword iter,
                 arma::uword burnin){
  arma::uword ntot = X.n_rows * Bs.n_rows * Bt.n_rows;
  arma::vec missing_indicator = arma::zeros<arma::vec>(ntot);
  arma::uword nsub = X.n_rows;// [[Rcpp::export]]
  arma::uword ns = Bs.n_rows;
  arma::uword nt = Bt.n_rows;
  for(arma::uword i = 0; i < nsub; i++){
    arma::vec missing_subj = missing(i);
    for(arma::uword j = 0; j < missing_subj.n_elem; j++){
      missing_indicator.subvec(i * ns * nt + nt * (missing_subj(j) - 1),
                               i * ns * nt + nt * (missing_subj(j)) - 1) =
                                 arma::ones<arma::vec>(nt);
    }
  }
  arma::uvec obs = arma::find(missing_indicator == 0);
  arma::uword nsim = iter - burnin;
  arma::uword num_not_missing = ntot -
    arma::uword(arma::sum(missing_indicator));
  arma::vec fit(ntot);
  arma::mat logliks(nsim, num_not_missing);

  for(arma::uword i = 0; i < nsim; i++){
    for(arma::uword j = 0; j < nsub; j++){
      fit.subvec(ns * nt * j, ns * nt * (j + 1) - 1) =
        arma::vectorise(Bt * Theta(i).slice(j) * Bs.t());
    }
      logliks.row(i) = arma::trans(
        arma::log_normpdf(y, fit(obs),
                          arma::vec(obs.n_elem).fill(1.0 / Varphi(0)(i))));
  }
  return(logliks);
}
