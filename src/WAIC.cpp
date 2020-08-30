#include <RcppArmadillo.h>

double mydnorm(arma::vec X, arma::vec mu, arma::mat Sigma, int givelog);
double dmvnrm_arma(arma::vec x,  
                   arma::vec mean,  
                   arma::mat sigma,
                   bool logd);


Rcpp::List calculate_WAIC(arma::mat Y, arma::mat X,
                          Rcpp::List mod, arma::mat splineS,
                          arma::mat splineT, arma::uword burnin){
  double LPPD = 0;
  double P_1 = 0;
  double P_2 = 0;
  double a;
  double b;
  arma::field<arma::cube> Eta = mod["Eta"];
  arma::field<arma::vec> Varphi = mod["Varphi"];
  arma::field<arma::cube> Sigma = mod["Sigma"];
  arma::field<arma::cube> Lambda = mod["Lambda"];
  arma::field<arma::cube> Gamma = mod["Gamma"];
  arma::field<arma::cube> Beta = mod["Beta"];
  arma::field<arma::cube> H = mod["H"];
  arma::uword nsubject = arma::size(Eta(0))(2);
  arma::uword q1 = Lambda(0).slice(0).n_cols;
  arma::uword q2 = Gamma(0).slice(0).n_cols;
  arma::uword nchain = Varphi.n_elem;
  arma::uword iter = Varphi(0).n_elem;
  arma::mat varmat;
  arma::vec meanvec;
  arma::vec meantemp;
  arma::mat L((iter - burnin) * nchain, nsubject);
  for(arma::uword i = burnin; i < iter; i = i + 10){
    Rcpp::Rcout << i << std::endl;
    for(arma::uword k = 0; k < nchain; k++){
      meantemp = arma::trans(Beta(k).slice(i)) * X.row(0);
      meanvec = arma::vectorise(splineT * Lambda(k).slice(i) *
        arma::reshape(meantemp, q1, q2) *
        arma::trans(Gamma(k).slice(i)) * arma::trans(splineS));
      varmat = arma::diagmat(arma::ones<arma::vec>(Y.col(0).n_elem) * .1) + 
        arma::kron(splineS, splineT) * 
        arma::inv(arma::diagmat(arma::vectorise(Sigma(k).slice(i)))) *
        arma::kron(splineS.t(), splineT.t()) + 
        arma::kron(splineS, splineT) * 
        arma::kron(Gamma(k).slice(i), Lambda(k).slice(i)) *
        arma::inv(arma::diagmat(H(k).slice(i))) * 
        arma::kron(arma::trans(Gamma(k).slice(i)), 
                   arma::trans(Lambda(k).slice(i))) *
        arma::kron(splineS.t(), splineT.t());
      for(arma::uword s = 0; s < nsubject; s++){
        L(i - (k + 1) * burnin + k * iter, s) = dmvnrm_arma(Y.col(s),
          meanvec, varmat, 0);
      }
    }
  }
  for(arma::uword s = 0; s < nsubject; s++){
    LPPD  = LPPD + log(arma::mean(L.col(s)));
    a = log(arma::mean(L.col(s)));
    b = arma::mean(arma::log(L.col(s)));
    
    P_1 = P_1 + 2 * (a - b);
    P_2 = P_2 + arma::var(arma::log(L.col(s)));
  }

  double WAIC_1 = -2 * (LPPD - P_1);
  double WAIC_2 = -2 * (LPPD - P_2);
  return(Rcpp::List::create(Rcpp::Named("WAIC_1", WAIC_1),
                            Rcpp::Named("WAIC_2", WAIC_2),
                            Rcpp::Named("P_1", P_1),
                            Rcpp::Named("P_2", P_2),
                            Rcpp::Named("L",L)));
}


Rcpp::List calculate_BIC(arma::mat Y, arma::mat X,
                         Rcpp::List mod, arma::mat splineS,
                         arma::mat splineT, arma::uword burnin,
                         arma::uword thin){
  arma::field<arma::cube> Lambda = mod["Lambda"];
  arma::field<arma::cube> Gamma = mod["Gamma"];
  arma::field<arma::cube> Beta = mod["Beta"];
  arma::field<arma::cube> H = mod["H"];
  arma::field<arma::cube> Sigma = mod["Sigma"];
  arma::field<arma::vec> Varphi = mod["Varphi"];
  arma::mat spline = arma::kron(splineS, splineT);
  arma::mat cov(spline.n_rows, spline.n_rows);
  arma::vec meantemp1, meantemp2;
  arma::vec meanvec = arma::zeros<arma::vec>(spline.n_rows);
  arma::mat postcov = arma::zeros<arma::mat>(spline.n_rows, spline.n_rows);
  arma::uword iter = Varphi(0).n_elem;
  arma::uword nchains = Varphi.n_elem;
  arma::uword nsubject = Y.n_cols;
  arma::uword q1 = Lambda(0).n_cols;
  arma::uword q2 = Gamma(0).n_cols;
  arma::uword p1 = splineT.n_rows;
  arma::uword p2 = splineS.n_rows;
  arma::uword counter = 0;
  double llik = 0;
  double BIC1 = 0;
  double BIC2 = 0;
  for(arma::uword k = 0; k < nchains; k++){
    for(arma::uword i = burnin; i < iter; i = i + thin){
      meantemp1 = arma::trans(Beta(k).slice(i)) * X.row(0);
      meantemp2 = arma::vectorise(splineT * Lambda(k).slice(i) *
        arma::reshape(meantemp1, q1, q2) *
        arma::trans(Gamma(k).slice(i)) * arma::trans(splineS));
      cov = spline * (arma::kron(Gamma(k).slice(i), Lambda(k).slice(i)) *
        arma::inv(arma::diagmat(H(k).slice(i))) *
        arma::trans(arma::kron(Gamma(k).slice(i), Lambda(k).slice(i)))) *
        spline.t() +
        spline * arma::solve(
            arma::diagmat(arma::vectorise(Sigma(k).slice(i))), spline.t()) + 
        arma::diagmat(arma::ones<arma::vec>(Y.col(0).n_elem) *
        1.0 / Varphi(k)(i));
      meanvec = meanvec + meantemp2;
      postcov = postcov + cov;
      counter++;
    }
  }
  postcov = postcov / counter;
  meanvec = meanvec / counter;
  
  for(arma::uword s = 0; s < nsubject; s++){
    llik = llik + dmvnrm_arma(Y.col(s), meanvec, postcov, 1);
  }
  
  BIC1 = log(nsubject) * ((p1 * q1) + (p2 * q2) + (q1 * q2) + 
    (q1 * q2) + (p1 * p2) + 1) - 2 * llik;
  BIC2 = log(nsubject * splineS.n_rows * splineT.n_rows) * 
    ((p1 * q1) + (p2 * q2) + (q1 * q2) + 1) + 
    log(nsubject) * ((p1 * p2) + (q1 * q2)) - 2 * llik;
  return(Rcpp::List::create(Rcpp::Named("BIC1", BIC1),
                            Rcpp::Named("BIC2", BIC2), 
                            Rcpp::Named("LogLik", llik), 
                            Rcpp::Named("postcov", postcov),
                            Rcpp::Named("meanvec", meanvec)));
}

Rcpp::List calculate_DIC(arma::mat Y, arma::mat X,
                         Rcpp::List mod, arma::mat splineS,
                         arma::mat splineT, arma::uword burnin,
                         arma::uword thin){
  arma::field<arma::cube> Lambda = mod["Lambda"];
  arma::field<arma::cube> Gamma = mod["Gamma"];
  arma::field<arma::cube> Beta = mod["Beta"];
  arma::field<arma::cube> H = mod["H"];
  arma::field<arma::cube> Sigma = mod["Sigma"];
  arma::field<arma::vec> Varphi = mod["Varphi"];
  arma::mat spline = arma::kron(splineS, splineT);
  arma::mat cov(spline.n_rows, spline.n_rows);
  arma::vec meantemp1, meantemp2;
  arma::vec meanvec = arma::zeros<arma::vec>(spline.n_rows);
  arma::mat postcov = arma::zeros<arma::mat>(spline.n_rows, spline.n_rows);
  arma::uword iter = Varphi(0).n_elem;
  arma::uword nchains = Varphi.n_elem;
  arma::uword nsubject = Y.n_cols;
  arma::uword q1 = Lambda(0).n_cols;
  arma::uword q2 = Gamma(0).n_cols;
  //arma::uword p1 = splineT.n_rows;
  //arma::uword p2 = splineS.n_rows;
  arma::uword counter = 0;
  double loglik = 0;
  double loglikiter = 0;
  double DIC = 0;
  double L = 0;
  double P = 0;
  for(arma::uword k = 0; k < nchains; k++){
    for(arma::uword i = burnin; i < iter; i = i + thin){
      Rcpp::Rcout << i << std::endl;
      meantemp1 = arma::trans(Beta(k).slice(i)) * X.row(0);
      meantemp2 = arma::vectorise(splineT * Lambda(k).slice(i) *
        arma::reshape(meantemp1, q1, q2) *
        arma::trans(Gamma(k).slice(i)) * arma::trans(splineS));
      cov = spline * (arma::kron(Gamma(k).slice(i), Lambda(k).slice(i)) *
        arma::inv(arma::diagmat(H(k).slice(i))) *
        arma::trans(arma::kron(Gamma(k).slice(i),
                               Lambda(k).slice(i)))) * spline.t() +
        spline * arma::solve(
            arma::diagmat(arma::vectorise(Sigma(k).slice(i))), spline.t()) + 
        arma::diagmat(arma::ones<arma::vec>(
            Y.col(0).n_elem) * 1.0 / Varphi(k)(i));
      meanvec = meanvec + meantemp2;
      postcov = postcov + cov;
      counter++;
      for(arma::uword s = 0; s < nsubject; s++){
        loglikiter = loglikiter + dmvnrm_arma(Y.col(s), meantemp2, cov, 1);
      }
    }
  }
  loglikiter = loglikiter / counter;
  postcov = postcov / counter;
  meanvec = meanvec / counter;
  
  for(arma::uword s = 0; s < nsubject; s++){
    loglik = loglik + dmvnrm_arma(Y.col(s), meanvec, postcov, 1);
  }
  P = 2 * (loglik - loglikiter);
  L = loglik;
  DIC = -2 * (L - P);
  return(Rcpp::List::create(Rcpp::Named("DIC", DIC),
                            Rcpp::Named("LogLik", loglik),
                            Rcpp::Named("postcov", postcov),
                            Rcpp::Named("meanvec", meanvec)));
}

Rcpp::List calculate_BIC_Missing(arma::field<arma::vec> Y,
                                 arma::mat X,
                                 arma::field<arma::uvec> observed,
                                 Rcpp::List mod,
                                 arma::mat splineS,
                                 arma::mat splineT,
                                 arma::uword burnin,
                                 arma::uword thin){
  arma::field<arma::cube> Lambda = mod["Lambda"];
  arma::field<arma::cube> Gamma = mod["Gamma"];
  arma::field<arma::cube> Beta = mod["Beta"];
  arma::field<arma::cube> H = mod["H"];
  arma::field<arma::cube> Sigma = mod["Sigma"];
  arma::field<arma::vec> Varphi = mod["Varphi"];
  arma::mat spline = arma::kron(splineS, splineT);
  arma::mat cov(spline.n_rows, spline.n_rows);
  arma::vec meantemp1, meantemp2;
  arma::vec meanvec = arma::zeros<arma::vec>(spline.n_rows);
  arma::mat postcov = arma::zeros<arma::mat>(spline.n_rows, spline.n_rows);
  arma::uword iter = Varphi(0).n_elem;
  arma::uword nchains = Varphi.n_elem;
  arma::uword nsubject = Y.n_elem;
  arma::uword q1 = Lambda(0).n_cols;
  arma::uword q2 = Gamma(0).n_cols;
  arma::uword p1 = splineT.n_rows;
  arma::uword p2 = splineS.n_rows;
  arma::uword counter = 0;
  int numelem = 0;
  double llik = 0;
  double BIC1 = 0;
  double BIC2 = 0;
  for(arma::uword k = 0; k < nchains; k++){
    for(arma::uword i = burnin; i < iter; i = i + thin){
      Rcpp::Rcout << i << std::endl;
      meantemp1 = arma::trans(Beta(k).slice(i)) * X.row(0);
      meantemp2 = arma::vectorise(splineT * 
        Lambda(k).slice(i) * arma::reshape(meantemp1, q1, q2) *
        arma::trans(Gamma(k).slice(i)) * arma::trans(splineS));
      cov = spline * (arma::kron(Gamma(k).slice(i), Lambda(k).slice(i)) *
        arma::inv(arma::diagmat(H(k).slice(i))) *
        arma::trans(arma::kron(Gamma(k).slice(i), Lambda(k).slice(i)))) *
        spline.t() +
        spline * arma::solve(arma::diagmat(arma::vectorise(Sigma(k).slice(i))),
                             spline.t()) + 
        arma::diagmat(arma::ones<arma::vec>(splineT.n_rows * 
        splineS.n_rows) * 1.0 / Varphi(k)(i));
      meanvec = meanvec + meantemp2;
      postcov = postcov + cov;
      counter++;
    }
  }
  postcov = postcov / counter;
  meanvec = meanvec / counter;
  
  for(arma::uword s = 0; s < nsubject; s++){
    arma::vec mymeanpost = meanvec.elem(observed(s));
    arma::mat mycovpost = postcov.submat(observed(s), observed(s));
    llik = llik + dmvnrm_arma(Y(s), mymeanpost, mycovpost, 1);
    numelem = numelem + observed(s).n_elem;
  }
  
  BIC1 = log(nsubject) * ((p1 * q1) + (p2 * q2) + (q1 * q2) + 
    (q1 * q2) + (p1 * p2) + 1) - 2 * llik;
  BIC2 = log(numelem) * ((p1 * q1) + (p2 * q2) + (q1 * q2) + 1) + 
    log(nsubject) * ((p1 * p2) + (q1 * q2)) - 2 * llik;
  return(Rcpp::List::create(Rcpp::Named("BIC1", BIC1),
                            Rcpp::Named("BIC2", BIC2), 
                            Rcpp::Named("LogLik", llik), 
                            Rcpp::Named("postcov", postcov),
                            Rcpp::Named("meanvec", meanvec)));
}


Rcpp::List calculate_DIC_Missing(arma::field<arma::vec> Y,
                                 arma::mat X,
                                 arma::field<arma::uvec> observed,
                                 Rcpp::List mod,
                                 arma::mat splineS,
                                 arma::mat splineT,
                                 arma::uword burnin,
                                 arma::uword thin){
  arma::field<arma::cube> Lambda = mod["Lambda"];
  arma::field<arma::cube> Gamma = mod["Gamma"];
  arma::field<arma::cube> Beta = mod["Beta"];
  arma::field<arma::cube> H = mod["H"];
  arma::field<arma::cube> Sigma = mod["Sigma"];
  arma::field<arma::vec> Varphi = mod["Varphi"];
  arma::mat spline = arma::kron(splineS, splineT);
  arma::mat cov(spline.n_rows, spline.n_rows);
  arma::vec meantemp1, meantemp2;
  arma::vec meanvec = arma::zeros<arma::vec>(spline.n_rows);
  arma::mat postcov = arma::zeros<arma::mat>(spline.n_rows, spline.n_rows);
  arma::uword iter = Varphi(0).n_elem;
  arma::uword nchains = Varphi.n_elem;
  arma::uword nsubject = Y.n_elem;
  arma::uword q1 = Lambda(0).n_cols;
  arma::uword q2 = Gamma(0).n_cols;
  arma::uword counter = 0;
  double loglik = 0;
  double loglikiter = 0;
  double DIC = 0;
  double L = 0;
  double P = 0;
  for(arma::uword k = 0; k < nchains; k++){
    for(arma::uword i = burnin; i < iter; i = i + thin){
      Rcpp::Rcout << i << std::endl;
      meantemp1 = arma::trans(Beta(k).slice(i)) * X.row(0);
      meantemp2 = arma::vectorise(splineT * Lambda(k).slice(i) *
        arma::reshape(meantemp1, q1, q2) *
        arma::trans(Gamma(k).slice(i)) * arma::trans(splineS));
      cov = spline * (arma::kron(Gamma(k).slice(i), Lambda(k).slice(i)) *
        arma::inv(arma::diagmat(H(k).slice(i))) *
        arma::trans(arma::kron(Gamma(k).slice(i), Lambda(k).slice(i)))) *
        spline.t() +
        spline * arma::solve(arma::diagmat(arma::vectorise(Sigma(k).slice(i))),
                             spline.t()) + 
        arma::diagmat(arma::ones<arma::vec>(splineT.n_rows * splineS.n_rows) *
        1.0 / Varphi(k)(i));
      meanvec = meanvec + meantemp2;
      postcov = postcov + cov;
      counter++;
      for(arma::uword s = 0; s < nsubject; s++){
        arma::vec mymean = meantemp2.elem(observed(s));
        arma::mat mycov = cov.submat(observed(s), observed(s));
        loglikiter = loglikiter + dmvnrm_arma(Y(s), mymean, mycov, 1);
      }
    }
  }
  loglikiter = loglikiter / counter;
  postcov = postcov / counter;
  meanvec = meanvec / counter;
  
  for(arma::uword s = 0; s < nsubject; s++){
    arma::vec mymeanpost = meanvec.elem(observed(s));
    arma::mat mycovpost = postcov.submat(observed(s), observed(s));
    loglik = loglik + dmvnrm_arma(Y(s), mymeanpost, mycovpost, 1);
  }
  P = 2 * (loglik - loglikiter);
  L = loglik;
  DIC = -2 * (L - P);
  return(Rcpp::List::create(Rcpp::Named("DIC", DIC),
                            Rcpp::Named("LogLik", loglik),
                            Rcpp::Named("postcov", postcov),
                            Rcpp::Named("meanvec", meanvec)));
}


double dmvnrm_arma(arma::vec x,  
                   arma::vec mean,  
                   arma::mat sigma, 
                   bool logd = false) { 
  double xlen = x.n_elem;
  double out;
  double log2pi = std::log(2.0 * M_PI);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xlen)/2.0) * log2pi;
  
  arma::vec z = rooti * (x - mean);    
  out = constants - 0.5 * arma::sum(z%z) + rootisum;
  
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

