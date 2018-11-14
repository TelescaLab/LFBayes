#include <RcppArmadillo.h>

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
Rcpp::List getStatistics(arma::mat spline, Rcpp::List mod, arma::rowvec x, arma::vec trueMean, int burnin, int cov) {
  arma::cube LambdaC = mod["Lambda"];
  arma::cube GammaC = mod["Gamma"];
  arma::cube betaC = mod["beta"];
  arma::cube HC = mod["HC"];
  arma::vec varphi = mod["varphi"];
  arma::mat sigma1 = mod["sigma1"];
  arma::mat sigma2 = mod["sigma2"];
  double IMSE = 0;
  double RASE;
  int iter = LambdaC.n_slices;
  arma::mat meanM(spline.n_rows, iter - burnin);
  arma::cube covC(spline.n_rows, spline.n_rows, iter - burnin);
  arma::mat coefM(GammaC.n_rows * LambdaC.n_rows, iter - burnin);
  arma::colvec delta(iter - burnin);
  arma::mat postcov = arma::zeros<arma::mat>(spline.n_rows, spline.n_rows);
  //omp_set_num_threads(12);
  //#pragma omp parallel for schedule(static)
  for(int i = 0; i < iter - burnin; i++){
    meanM.col(i) = spline * arma::kron(GammaC.slice(burnin + i),
              LambdaC.slice(burnin + i)) * arma::trans(betaC.slice(burnin + i)) * x;

    coefM.col(i) = arma::kron(GammaC.slice(burnin + i),
              LambdaC.slice(burnin + i)) * arma::trans(betaC.slice(burnin + i)) * x;
  }

  if(cov == 1){
    //omp_set_num_threads(12);
    //#pragma omp parallel for schedule(static)
    for(int i = 0; i < iter - burnin; i++){
      covC.slice(i) = spline * (arma::kron(GammaC.slice(burnin + i), LambdaC.slice(burnin + i)) * arma::inv(arma::diagmat(HC.slice(burnin + i))) *
        arma::trans(arma::kron(GammaC.slice(burnin + i), LambdaC.slice(burnin + i))) + arma::kron(arma::diagmat(1.0 / sigma2.col(burnin + i)), arma::diagmat(1.0 / sigma1.col(burnin + i)))) * arma::trans(spline) +
        1.0 / varphi(burnin + i) * arma::eye<arma::mat>(spline.n_rows, spline.n_rows);
      postcov = postcov + covC.slice(i);
    }
    postcov = postcov / (iter - burnin);
  }
  arma::colvec postcoef = arma::mean(coefM, 1);
  arma::colvec postmean = arma::mean(meanM, 1);
  arma::colvec postsd = arma::trans(arma::stddev(arma::trans(meanM), 1));
  /*for(int i = 0; i < iter - burnin; i++){
    delta[i] = (postmean - (spline * arma::kron(GammaC.slice(burnin + i),
              LambdaC.slice(burnin + i)) * arma::trans(betaC.slice(burnin + i)) * x) / postsd).max();
  }*/
  IMSE = arma::as_scalar(arma::sum(arma::square(postmean - trueMean))) / spline.n_rows;
  //RASE = ::pow(IMSE, 1.0 / 2.0);
  /*
  return Rcpp::List::create(Rcpp::Named("postmean", postmean),
                            Rcpp::Named("postsd", postsd),
                            Rcpp::Named("delta", delta),
                            Rcpp::Named("meanM", meanM),
                            Rcpp::Named("postcoef", postcoef),
                            Rcpp::Named("covC", covC),
                            Rcpp::Named("postcov", postcov));
   */
  return Rcpp::List::create(Rcpp::Named("postmean", postmean),
                            Rcpp::Named("covC", covC),
                            Rcpp::Named("postcov", postcov),
                            Rcpp::Named("IMSE", IMSE));
                            //Rcpp::Named("meanM", meanM),
                            //Rcpp::Named("postsd",postsd));
}
