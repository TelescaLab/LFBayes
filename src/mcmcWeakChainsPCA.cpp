#include <RcppArmadillo.h>
#include <cmath>
#include <typeinfo>
#include "updateParam.h"
#include <omp.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List mcmcWeakChainsPCA(arma::field<arma::vec> y, arma::field<arma::vec> missing, arma::mat X, arma::mat splineS, arma::mat splineT, int q1, int q2, int iter, int thin, int burnin, int nchains){
  // Allocate memory for parameters
  int p1 = splineT.n_cols;
  int p2 = splineS.n_cols;
  arma::field<arma::cube> LambdaF(nchains);
  arma::field<arma::cube> GammaF(nchains);
  arma::field<arma::cube> Phi1F(nchains);
  arma::field<arma::cube> Phi2F(nchains);
  arma::field<arma::mat> DM1F(nchains);
  arma::field<arma::mat> DM2F(nchains);
  arma::field<arma::cube> SigmaF(nchains);
  arma::field<arma::vec> varphiF(nchains);
  arma::field<arma::vec> a1LF(nchains);
  arma::field<arma::vec> a2LF(nchains);
  arma::field<arma::vec> a1GF(nchains);
  arma::field<arma::vec> a2GF(nchains);
  arma::field<arma::cube> BetaF(nchains);
  arma::field<arma::cube> HF(nchains);
  arma::field<arma::cube> EF(nchains);
  arma::field<arma::mat> A1F(nchains);
  arma::field<arma::mat> A2F(nchains);
  for(arma::uword i = 0; i < nchains; i++){
    LambdaF(i) = arma::cube(p1, q1, iter);
    GammaF(i) = arma::cube(p2, q2, iter);
    Phi1F(i) = arma::cube(p1, q1, iter);
    Phi2F(i) = arma::cube(p2, q2, iter);
    DM1F(i) = arma::mat(q1, iter);
    DM2F(i) = arma::mat(q2, iter);
    A1F(i) = arma::mat(q1, iter);
    A2F(i) = arma::mat(q2, iter);
    SigmaF(i) = arma::cube(p1, p2, iter);
    varphiF(i) = arma::vec(iter);
    a1LF(i) = arma::vec(iter);
    a2LF(i) = arma::vec(iter);
    a1GF(i) = arma::vec(iter);
    a2GF(i) = arma::vec(iter);
    BetaF(i) = arma::cube(X.n_cols, q1 * q2, iter);
    HF(i) = arma::cube(q1 * q2, q1 * q2, iter);
    EF(i) = arma::cube(q1 * q2, X.n_cols, iter);
    
  }
  arma::mat Lambda(p1, q1);
  arma::mat Gamma(p2, q2);
  arma::mat Phi1(p1, q1);
  arma::mat Phi2(p2,q2);
  arma::vec DM1(q1);
  arma::vec DM2(q2);
  arma::vec A1(q1);
  arma::vec A2(q2);
  arma::mat Sigma(p1,p2);
  double a1Ld = 1;
  double a2Ld = 1;
  double a1Gd = 1;
  double a2Gd = 1;
  double V;
  arma::mat Beta = arma::zeros(X.n_cols, q1 * q2);
  arma::mat H(q1 * q2, q1 * q2);
  arma::mat E(q1 * q2, X.n_cols);
  arma::field<arma::cube> etaF(nchains, iter);
  arma::field<arma::cube> thetaF(nchains, iter);
  for(int i = 0; i < nchains; i++){
    for(int j = 0; j < iter; j++){
      etaF(i,j) = arma::cube(q1, q2, X.n_rows);
      thetaF(i,j) = arma::cube(p1, p2, X.n_rows);
    }
    etaF(i, 0).randn();
    thetaF(i, 0).randn();
  }
  arma::cube Eta(q1, q2, X.n_rows);
  arma::cube Theta(p1, p2, X.n_rows);
  
  arma::mat imputedY(splineS.n_rows * splineT.n_rows, X.n_rows);
  
  for(arma::uword i = 0; i < X.n_rows; i++){
    imputedY.col(i) = initializeY(y(i), missing(i), splineS.n_rows, splineT.n_rows);
  }
  arma::mat initialY = imputedY;
  
  // Set initial values
  Lambda.randn();
  Gamma.randn();
  H.eye();
  V = R::rgamma(1,1);
  Eta.randn();
  Theta.randn();
  Beta.randn();
  E.randu();
  Phi1.ones();//
  Phi2.ones();//
  Sigma.fill(.1);
  Beta.randn();
  E.ones();
  DM1.ones();
  DM2.ones();
  A1.ones();
  A2.ones();
  Rcpp::Rcout << "Starting MCMC..." << std::endl;
  for(int k = 0; k < nchains; k++){
    for(int i = 0; i < iter; i++){
      //if(i % 1000 == 0){
      //Rcpp::Rcout << i << std::endl;
      //}
      for(int j = 0; j < thin; j++){
        
        Gamma = updateGammaPCA(Eta, Lambda, Sigma, DM2, A2, Theta);
        A2 = updateAlphaPCA(Gamma, DM2, A2);
        DM2 = updateDeltaProdPCA(Gamma, DM2, A2, 1.0, 2.0);
        Lambda = updateLambdaPCA(Eta, Gamma, Sigma, DM1, A1, Theta);
        A1 = updateAlphaPCA(Lambda, DM1, A1);
        DM1 = updateDeltaProdPCA(Lambda, DM1, A1, 1.0, 2.0);
        a1Ld = updatea1(DM1, a1Ld);
        a2Ld = updatea2(DM1, a2Ld);
        a1Gd = updatea1(DM2, a1Gd);
        a2Gd = updatea2(DM2, a2Gd);
        Sigma = updateSigma(Lambda, Gamma, Theta, Eta);
        V = updateVarphi(Theta, splineS, splineT, imputedY);
        //updateEta3Sig(Gamma, Lambda, Sigma, Theta, H, X, Beta, Eta);
        updateEtaSig(Lambda, Gamma, Sigma, H, splineS, splineT, imputedY, V, Beta, X, Eta);
        if(i % 25 == 0){
          //Theta = updateThetaSig(Lambda, Gamma, Sigma, Eta, splineS, splineT, imputedY, V);
        }
        H = updateH(Eta, Beta, X);
        updateBeta(Eta, H, E, X, Beta);
        E = updateE(Beta);
        updateMissing(imputedY, missing, Theta, V, splineS, splineT, X.n_rows);
        
      }
      LambdaF(k).slice(i) = Lambda;
      GammaF(k).slice(i) = Gamma;
      SigmaF(k).slice(i) = Sigma;
      HF(k).slice(i) = H;
      BetaF(k).slice(i) = Beta;
      
      A1F(k).col(i) = A1;
      A2F(k).col(i) = A2;
      DM1F(k).col(i) = DM1;
      DM2F(k).col(i) = DM2;
      Phi1F(k).slice(i) = Phi1;
      Phi2F(k).slice(i) = Phi2;
      varphiF(k)(i) = V;
      etaF(k, i) = Eta;
      thetaF(k, i) = Theta;
      
    }
  }
  
  
  Rcpp::Rcout << "All done!";
  
  Rcpp::List mod = Rcpp::List::create(Rcpp::Named("Lambda", LambdaF),
                                      Rcpp::Named("Gamma", GammaF),
                                      Rcpp::Named("H", HF), Rcpp::Named("Sigma", SigmaF),
                                      Rcpp::Named("Beta", BetaF), Rcpp::Named("Theta", thetaF),
                                      Rcpp::Named("Delta1", DM1F), Rcpp::Named("Delta2", DM2F),
                                      Rcpp::Named("Alpha1", A1F), Rcpp::Named("Alpha2", A2F),
                                      Rcpp::Named("Eta", etaF), Rcpp::Named("Phi1", Phi1F),
                                      Rcpp::Named("Phi2", Phi2F), Rcpp::Named("initialY", initialY),
                                      Rcpp::Named("Varphi", varphiF));
  
  return(mod);
  //return(eigenLFChains(splineS, splineT, mod, 3, iter, burnin, nchains));
  //return(eigenLFC)
  /*
  //Rcpp::Named("Delta1", Delta1),
  //Rcpp::Named("Delta2", Delta2),
  //Rcpp::Named("Tau1", Tau1),
  //Rcpp::Named("Tau2", Tau2),
  //Rcpp::Named("imputedY", imputedY), Rcpp::Named("initialY", initialY),
  //Rcpp::Named("a1L", a1L), Rcpp::Named("a2L", a2L),
  //Rcpp::Named("a1G", a1G), Rcpp::Named("a2G", a2G));
  //Rcpp::Rcout << arma::size(imputedY);
  //return( eigenLF(splineS, splineT, mod, 3, burnin));
  */
  
}

