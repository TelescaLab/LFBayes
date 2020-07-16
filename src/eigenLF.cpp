#include <RcppArmadillo.h>
#include "updateParam.h"
#include "classes.h"
#include <typeinfo>

// [[Rcpp::export]]
Rcpp::List LFB_post(arma::mat splineS, arma::mat splineT, Rcpp::List mod, 
              arma::uword numeig, arma::uword iter, arma::uword burnin,
              arma::uword nchains,
              arma::vec s, arma::vec t) {
 arma::field<arma::cube> LambdaF = mod["Lambda"];
 arma::field<arma::cube> GammaF = mod["Gamma"];
 arma::field<arma::cube> BetaF = mod["Beta"];
 arma::field<arma::cube> HF = mod["H"];
 arma::field<arma::cube> SigmaF = mod["Sigma"];
 arma::uword num_draws = iter - burnin;
 arma::uword q1 = LambdaF(0).n_cols;
 arma::uword q2 = GammaF(0).n_cols;
 arma::mat sigma, h, beta;
 
 EigenStruct func_struct(splineT, splineS, t,
                         s, numeig, num_draws, nchains, "functional");
 EigenStruct long_struct(splineS, splineT, s,
                         t, numeig, num_draws, nchains, "longitudinal");
 MeanStruct mean_struct(splineS, splineT, num_draws, nchains);
 for(arma::uword k = 0; k < nchains; k++) {
   for(arma::uword i = 0; i < num_draws; i++){
     if(int(i) % int(floor(double(num_draws) / 10.0)) == 0){
       Rcpp::Rcout << 100 * i / (iter - burnin) << '%' << std::endl;
     }
     sigma = 1 / SigmaF(k).slice(burnin + i);
     h = arma::reshape(1 / arma::diagvec(HF(k).slice(burnin + i)), q1, q2);
     beta = arma::reshape(BetaF(k).slice(burnin + i), q1, q2);
     func_struct.preprocess(GammaF(k).slice(burnin + i), h, sigma);
     func_struct.get_eigenfunction(LambdaF(k).slice(burnin + i));
     long_struct.preprocess(LambdaF(k).slice(burnin + i), h, sigma);
     long_struct.get_eigenfunction(GammaF(k).slice(burnin + i));
     mean_struct.update_mean(LambdaF(k).slice(burnin + i),
                             GammaF(k).slice(burnin + i), beta);
   }
 }
 func_struct.normalize();
 long_struct.normalize();
 mean_struct.normalize();
 Rcpp::List func_list = func_struct.get_eigen_bands();
 Rcpp::List long_list = long_struct.get_eigen_bands();
 Rcpp::List mean_list = mean_struct.get_mean();
 return(Rcpp::List::create(Rcpp::Named("functional", func_list),
                           Rcpp::Named("longitudinal", long_list),
                           Rcpp::Named("mean", mean_list)));
}


