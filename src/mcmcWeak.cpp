#include <RcppArmadillo.h>
#include <cmath>
#include <typeinfo>
#include "updateParam.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List mcmcWeak(arma::field<arma::vec> y, arma::field<arma::vec> missing, arma::mat X, arma::mat splineS, arma::mat splineT, int q1, int q2, int iter, int thin, int burnin){
  // Allocate memory for parameters
  int p1 = splineT.n_cols;
  int p2 = splineS.n_cols;
  arma::cube LambdaC(p1, q1, iter);
  arma::mat Lambda(p1, q1);
  arma::cube GammaC(p2, q2, iter);
  arma::mat Gamma(p2, q2);
  arma::cube Phi1C(p1, q1, iter);//
  arma::mat Phi1(p1, q1);
  arma::cube Phi2C(p2, q2, iter);//
  arma::mat Phi2(p2,q2);

  arma::mat Delta1(q1, iter);
  arma::vec D1(q1);
  arma::mat Delta2(q2, iter);
  arma::vec D2(q2);
  arma::mat DeltaM(q1, iter);//
  arma::mat DeltastarM(q2, iter);//
  arma::vec DM1(q1);
  arma::vec DM2(q2);
  arma::mat Tau1(q1, iter);
  arma::mat Tau2(q2, iter);
  arma::mat sigma1M(p1, iter);
  arma::vec S1(p1);
  arma::mat sigma2M(p2, iter);
  arma::vec S2(p2);
  arma::vec varphiV(iter);
  arma::vec a1L(iter);
  arma::vec a2L(iter);
  arma::vec a1G(iter);
  arma::vec a2G(iter);
  arma::vec bv(iter);
  double a1Ld = 1;
  double a2Ld = 1;
  double a1Gd = 1;
  double a2Gd = 1;
  double V;
  //arma::mat Delta(q1 + q2 - 1, iter);
  arma::field<arma::cube> etaF(iter);
  for(int i = 0; i < iter; i++){
    etaF(i) = arma::cube(LambdaC.n_cols, GammaC.n_cols, X.n_rows);
  }
  arma::cube Eta(Lambda.n_cols, Gamma.n_cols, X.n_rows);
  arma::field<arma::cube> thetaF(iter);
  for(int i = 0; i < iter; i++){
    thetaF(i) = arma::cube(p1, p2, X.n_rows);
  }
  arma::cube Theta(p1, p2, X.n_rows);
  arma::cube betaC = arma::zeros(X.n_cols, q1 * q2, iter);
  arma::mat Beta = arma::zeros(X.n_cols, q1 * q2);
  arma::cube HC(q1 * q2, q1 * q2, iter);
  arma::mat H(q1 * q2, q1 * q2);
  H.zeros();
  arma::cube EC(q1 * q2, X.n_cols, iter);
  arma::mat E(q1 * q2, X.n_cols);
  arma::mat imputedY(splineS.n_rows * splineT.n_rows, X.n_rows);
  for(arma::uword i = 0; i < X.n_rows; i++){
    imputedY.col(i) = initializeY(y(i), missing(i), splineS.n_rows, splineT.n_rows);
  }

  arma::mat initialY = imputedY;

  // Set initial values
  Lambda.randn();
  Gamma.randn();
  D1.ones();
  D2.ones();
  S1.ones();
  S2.ones();
  //H.diag() = arma::randu<arma::vec>(q1 * q2);
  H.diag() = arma::square(arma::randn<arma::vec>(q1 * q2));
  V = R::rgamma(1,1);
  Eta.randn();
  Theta.randn();
  Beta.randn();
  E.randu();
  LambdaC.slice(0).randn();
  GammaC.slice(0).randn();
  Tau1.ones();//
  Tau2.ones();//
  //Delta.ones();
  DeltaM.ones();//
  DeltastarM.ones();//
  Delta1.ones();
  Delta2.ones();
  Phi1.ones();//
  Phi2.ones();//
  sigma1M.fill(.01);
  sigma2M.fill(.01);
  varphiV(0) = 1;
  etaF(0).randn();
  thetaF(0).ones();
  betaC.slice(0).ones();
  HC.slice(0).eye();
  EC.slice(0).ones();
  DM1.ones();
  DM2.ones();
  //Rcpp::Rcout << H;

  /*
  for(int i = 0; i < iter - 1; i++){
  EC.slice(i+1) = updateE(betaC.slice(0));

  }
  */

  Rcpp::Rcout << "Starting MCMC..." << std::endl;
  for(int i = 0; i < iter; i++){
    if(i % 10000 == 0){
      Rcpp::Rcout << i << std::endl;
    }
    for(int j = 0; j < thin; j++){
    /*
      updateLambda(Eta, Gamma, DM1,  Phi1,
                    S1, S2, Theta, Lambda);
      Phi1 = updatePhi(Lambda, DM1);
      DM1 = updateDelta(Lambda, Phi1, DM1, a1Ld, a2Ld);

      updateGamma(Eta, Lambda, DM2, Phi2, S1,
              S2, Theta, Gamma);
      Phi2 = updatePhistar(Gamma, DM2);
      DM2 = updateDelta(Gamma, Phi2, DM2, a1Ld, a2Ld);
    */

      Lambda = updateLambdaSmoothD(Eta, Gamma, S1, S2, D1, Theta);
      Gamma = updateGammaSmoothD(Eta, Lambda, S1, S2, D2, Theta);
      D1 = updateDeltaProdTemp(Lambda, D1, a1Ld, a2Ld);
      a1Ld = updatea1(D1, a1Ld);
      a2Ld = updatea2(D1, a2Ld);
      D2 = updateDeltaProdTemp(Gamma, D2, a1Gd, a2Gd);
      a1Gd = updatea1(D2, a1Gd);
      a2Gd = updatea2(D2, a2Gd);
      S1 = updateSigma1(Eta, Theta, Lambda, Gamma, S2);
      S2 = updateSigma2(Eta, Theta, Lambda, Gamma, S1);
      V = updateVarphi(Theta, splineS, splineT, imputedY);
      updateEta(Lambda, Gamma, S1, S2, H, splineS, splineT, imputedY, V, Beta, X, Eta);
      Theta = updateTheta(Lambda, Gamma, S1, S2, Eta, splineS, splineT, imputedY, V);
      H = updateH(Eta, Beta, X);
      //H.eye();
      updateBeta(Eta, H, E, X, Beta);
      E = updateE(Beta);
      updateMissing(imputedY, missing, Theta, splineS, splineT, X.n_rows);

    }
    LambdaC.slice(i) = Lambda;
    GammaC.slice(i) = Gamma;
    a1L(i) = a1Ld;
    a2L(i) = a2Ld;
    a1G(i) = a1Gd;
    a2G(i) = a2Gd;
    //Delta1.col(i) = D1;
    //Delta2.col(i) = D2;
    Delta1.col(i) = DM1;
    Delta2.col(i) = DM2;
    Phi1C.slice(i) = Phi1;
    Phi2C.slice(i) = Phi2;
    sigma1M.col(i) = S1;
    sigma2M.col(i) = S2;
    varphiV(i) = V;
    etaF(i) = Eta;
    thetaF(i) = Theta;
    HC.slice(i) = H;
    betaC.slice(i) = Beta;
    EC.slice(i) = E;

    /*

    //updateLambda(etaF(i), GammaC.slice(i), DeltaM.col(i),  PhiC.slice(i),
      //           sigma1M.col(i), sigma2M.col(i), thetaF(i), LambdaC.slice(i+1));
    //LambdaC.slice(i+1) = updateLambdaSmooth(etaF(i), GammaC.slice(i), sigma1M.col(i), sigma2M.col(i), Tau1.col(i), thetaF(i));
    //LambdaC.slice(i+1) = Lambda;
    //updateGamma(etaF(i), LambdaC.slice(i+1),
      //          DeltastarM.col(i), PhistarC.slice(i), sigma1M.col(i),
        //        sigma2M.col(i), thetaF(i), GammaC.slice(i+1));
    LambdaC.slice(i+1) = updateLambdaSmoothD(etaF(i), GammaC.slice(i), sigma1M.col(i), sigma2M.col(i), Delta1.col(i), thetaF(i));
    GammaC.slice(i+1) = updateGammaSmoothD(etaF(i), LambdaC.slice(i+1), sigma1M.col(i), sigma2M.col(i), Delta2.col(i), thetaF(i));
    //GammaC.slice(i+1) = updateGammaSmooth(etaF(i), LambdaC.slice(i+1), sigma1M.col(i), sigma2M.col(i), Tau2.col(i), thetaF(i));
    //Tau1.col(i+1) = updateTau(LambdaC.slice(i+1));
    //Tau2.col(i+1) = updateTau(GammaC.slice(i+1));
    Delta1.col(i+1) = updateDeltaProdTemp(LambdaC.slice(i+1), Delta1.col(i));
    Delta2.col(i+1) = updateDeltaProdTemp(GammaC.slice(i+1), Delta2.col(i));
    //GammaC.slice(i+1) = Gamma;
    //PhiC.slice(i + 1) = updatePhi(LambdaC.slice(i + 1), DeltaM.col(i));
    //PhistarC.slice(i + 1) = updatePhistar(GammaC.slice(i + 1), DeltastarM.col(i));
    //DeltaM.col(i + 1) = updateDelta(LambdaC.slice(i + 1), PhiC.slice(i + 1), DeltaM.col(i), 1, 2);
    //DeltastarM.col(i + 1) = updateDelta(GammaC.slice(i + 1), PhistarC.slice(i + 1), DeltastarM.col(i), 1, 2);
    sigma1M.col(i+1) = updateSigma1(etaF(i), thetaF(i), LambdaC.slice(i+1), GammaC.slice(i+1), sigma2M.col(i));
    //sigma1M.col(i+1) = sigma1;
    sigma2M.col(i+1) = updateSigma2(etaF(i), thetaF(i), LambdaC.slice(i+1), GammaC.slice(i+1), sigma1M.col(i+1));
    //sigma2M.col(i+1) = sigma2;

    varphiV(i+1) = updateVarphi(thetaF(i), splineS, splineT, imputedY);
    //varphiV(i+1) = varphi;


    //updateEtaProd(LambdaC.slice(i+1), GammaC.slice(i+1), sigma1M.col(i+1),
      //        sigma2M.col(i+1), Delta.col(i), splineS, splineT, imputedY, varphiV(i+1),
        //      betaC.slice(i), X, etaF(i+1));


    updateEta(LambdaC.slice(i+1), GammaC.slice(i+1), sigma1M.col(i+1), sigma2M.col(i+1), HC.slice(i), splineS, splineT, imputedY, varphiV(i+1), betaC.slice(i), X, etaF(i+1));
    //etaF(i+1) = etaF(0);
    thetaF(i+1) = updateTheta(LambdaC.slice(i+1), GammaC.slice(i+1), sigma1M.col(i+1),
           sigma2M.col(i+1), etaF(i+1), splineS, splineT, imputedY, varphiV(i+1));
    //thetaF(i+1) = thetaF(0);
    HC.slice(i+1) = updateH(etaF(i+1), betaC.slice(i), X);

    //if(i > 6000){
      //Delta.col(i + 1) = updateHProd(etaF(i + 1), betaC.slice(i), X, Delta.col(i));
    //}

    //HC.slice(i+1) = HC.slice(0);
    //updateBetaProd(etaF(i+1), Delta.col(i + 1), EC.slice(i), X, betaC.slice(i+1));
    updateBeta(etaF(i+1), HC.slice(i+1), EC.slice(i), X, betaC.slice(i+1));
    //betaC.slice(i+1) = beta;
    EC.slice(i+1) = updateE(betaC.slice(i+1));
    updateMissing(imputedY, missing, thetaF(i), splineS, splineT, X.n_rows);
    //meanM.col(i+1) = spline * arma::kron(GammaC.slice(i+1), LambdaC.slice(i+1)) * arma::trans(betaC.slice(i+1)) * arma::trans(X.row(0));
     */

  }
  Rcpp::Rcout << "All done!";

  Rcpp::List mod = Rcpp::List::create(Rcpp::Named("Lambda", LambdaC),
                                      Rcpp::Named("Gamma", GammaC),
                                      Rcpp::Named("sigma1", sigma1M), Rcpp::Named("sigma2", sigma2M),
                                      Rcpp::Named("varphi", varphiV),
                                      Rcpp::Named("eta", etaF),
                                      Rcpp::Named("theta", thetaF),
                                      Rcpp::Named("HC", HC),
                                      Rcpp::Named("beta", betaC), Rcpp::Named("Delta1", Delta1),
                                      Rcpp::Named("Delta2", Delta2),Rcpp::Named("a1L", a1L), Rcpp::Named("a2L", a2L),
                                      Rcpp::Named("a1G", a1G), Rcpp::Named("a2G", a2G),
                                      Rcpp::Named("bv", bv));
  //return(mod);
  //Rcpp::Named("Delta1", Delta1),
  //Rcpp::Named("Delta2", Delta2),
  //Rcpp::Named("Tau1", Tau1),
  //Rcpp::Named("Tau2", Tau2),
  //Rcpp::Named("imputedY", imputedY), Rcpp::Named("initialY", initialY),
  //Rcpp::Named("a1L", a1L), Rcpp::Named("a2L", a2L),
  //Rcpp::Named("a1G", a1G), Rcpp::Named("a2G", a2G));
  return( eigenLF(splineS, splineT, mod, 3, burnin));

  /*
  return Rcpp::List::create(Rcpp::Named("Lambda", LambdaC),
                            Rcpp::Named("Gamma", GammaC),
                            Rcpp::Named("sigma1", sigma1M), Rcpp::Named("sigma2", sigma2M),
                            Rcpp::Named("varphi", varphiV),
                            Rcpp::Named("eta", etaF),
                            Rcpp::Named("theta", thetaF),
                            Rcpp::Named("HC", HC),
                            Rcpp::Named("beta", betaC), //Rcpp::Named("Delta", Delta),
                            //Rcpp::Named("Delta1", Delta1),
                            //Rcpp::Named("Delta2", Delta2),
                            //Rcpp::Named("Tau1", Tau1),
                            //Rcpp::Named("Tau2", Tau2),
                            //Rcpp::Named("imputedY", imputedY), Rcpp::Named("initialY", initialY),
                            Rcpp::Named("a1L", a1L), Rcpp::Named("a2L", a2L),
                            Rcpp::Named("a1G", a1G), Rcpp::Named("a2G", a2G));
   */

}

