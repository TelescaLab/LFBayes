#ifndef MODSTRING_H
#define MODSTRING_H

#include <RcppArmadillo.h>

void updateLambda(arma::cube &eta, arma::mat &Gamma, arma::vec Delta, arma::mat &Phi, arma::vec sigma1, arma::vec sigma2, arma::cube &theta, arma::mat &Lambda);
arma::mat updateLambda2(arma::mat Theta, arma::mat eta, arma::vec Sigma, arma::vec Tau);
arma::mat updateLambda3(arma::mat Theta, arma::mat eta, arma::vec Sigma, arma::mat Phi, arma::vec Delta);
arma::mat updatePhi(arma::mat &Lambda, arma::vec delta);
arma::vec updateDelta(arma::mat &Lambda, arma::mat &Phi, arma::vec Delta, double a1, double a2);
arma::vec updateDelta3(arma::mat &Lambda, arma::mat &Phi, arma::vec Delta, double a1, double a2);
void updateGamma(arma::cube &eta, arma::mat &Lambda, arma::vec Deltastar,
                      arma::mat &Phistar, arma::vec sigma1, arma::vec sigma2,
                      arma::cube &theta, arma::mat &Gamma);
arma::mat updatePhistar(arma::mat &Gamma, arma::vec deltastar);
arma::vec updateDeltastar(arma::mat &Gamma, arma::mat &Phistar, arma::vec Deltastar, double a1star, double a2star);
arma::vec updateSigma1(arma::cube &eta, arma::cube &theta, arma::mat &Lambda, arma::mat &Gamma, arma::vec sigma2);
arma::vec updateSigma2(arma::cube &eta, arma::cube &theta, arma::mat &Lambda, arma::mat &Gamma,arma::vec sigma1);
arma::vec updateSigma2(arma::mat Theta, arma::mat Lambda, arma::mat Eta);
double updateVarphi(arma::cube &theta, arma::mat &splineS, arma::mat &splineT, arma::mat &y);
double updateVarphi2(arma::mat Data, arma::mat Theta, arma::mat Basis);
double updateVarphiSparse(arma::cube &theta, arma::field<arma::mat> &splineS, arma::mat &splineT, arma::field<arma::vec> &y);
void updateEta(arma::mat &Lambda, arma::mat &Gamma, arma::vec sigma1,
                     arma::vec sigma2, arma::mat H, arma::mat &splineS,
                     arma::mat &splineT, arma::mat &y, double varphi,
                     arma::mat &beta, arma::mat &X, arma::cube &eta);
void updateEtaSig(arma::mat &Lambda, arma::mat &Gamma, arma::mat Sigma,
                  arma::mat H, arma::mat &splineS,
                  arma::mat &splineT, arma::mat &y, double varphi,
                  arma::mat &beta, arma::mat &X, arma::cube &eta);
arma::mat updateEta2(arma::mat Lambda, arma::mat Basis, arma::vec Sigma, arma::mat Data, arma::mat Beta, arma::mat X, double Varphi);
void updateEta3(arma::mat Gamma, arma::mat Lambda, arma::vec sigma2, arma::vec sigma1, arma::cube Theta, arma::mat H, arma::mat X, arma::mat Beta, arma::cube &eta);
void updateEtaProd(arma::mat &Lambda, arma::mat &Gamma, arma::vec sigma1,
                   arma::vec sigma2, arma::vec Delta, arma::mat &splineS,
                   arma::mat &splineT, arma::mat &y, double varphi,
                   arma::mat &beta, arma::mat &X, arma::cube &eta);
void updateEtaSparse(arma::mat &Lambda, arma::mat &Gamma, arma::vec sigma1,
                     arma::vec sigma2, arma::mat &H, arma::field<arma::mat> &splineS,
                     arma::mat &splineT, arma::field<arma::vec> &y, double varphi,
                     arma::mat &beta, arma::mat &X, arma::cube &eta);
arma::cube updateTheta(arma::mat &Lambda, arma::mat &Gamma, arma::vec sigma1,
                       arma::vec sigma2, arma::cube &eta, arma::mat &splineS,
                       arma::mat &splineT, arma::mat &y, double varphi);
arma::mat updateTheta2(arma::mat Lambda, arma::vec Sigma, arma::mat eta, arma::mat Data, double Varphi, arma::mat Basis);
arma::cube updateThetaSparse(arma::mat &Lambda, arma::mat &Gamma, arma::vec sigma1,
                             arma::vec sigma2, arma::cube &eta, arma::field<arma::mat> &splineS,
                             arma::mat &splineT, arma::field<arma::vec> &y, double varphi);
void updateBeta(arma::cube &eta, arma::mat &H, arma::mat &E, arma::mat &X, arma::mat &beta);
arma::mat updateBeta2(arma::mat E, arma::mat Eta, arma::mat X);
void updateBetaProd(arma::cube &eta, arma::vec Delta, arma::mat &E, arma::mat &X, arma::mat &beta);

arma::mat updateH(arma::cube &eta, arma::mat &beta, arma::mat &X);
arma::vec updateHProdOne(arma::mat Eta, arma::vec Delta);
arma::mat updateE(arma::mat beta);
arma::mat get_marginal_func(arma::mat &cov, int ns, int nt);
arma::mat get_marginal_long(arma::mat cov, int ns, int nt);
arma::vec initializeY(arma::vec x, arma::vec missing, int s, int t);
void updateMissing(arma::mat &y, arma::field<arma::vec> &missing, arma::cube &theta, double Varphi, arma::mat &splineS, arma::mat &splineT, int n);
arma::mat getPenalty(arma::uword n);
arma::mat updateLambdaSmooth(arma::cube Eta, arma::mat Gamma, arma::vec Sigma1, arma::vec Sigma2, arma::vec Tau, arma::cube Theta);
arma::mat updateLambdaSmoothD(arma::cube Eta, arma::mat Gamma, arma::vec Sigma1, arma::vec Sigma2, arma::vec Delta, arma::cube Theta);
arma::mat updateLambdaSmoothDSig(arma::cube Eta, arma::mat Gamma, arma::mat Sigma, arma::vec Delta, arma::cube Theta);
arma::mat updateLambda4(arma::mat Theta, arma::mat eta, arma::vec Sigma, arma::vec Delta);
arma::mat updateGammaSmooth(arma::cube Eta, arma::mat Lambda, arma::vec Sigma1, arma::vec Sigma2, arma::vec Tau, arma::cube Theta);
arma::mat updateGammaSmoothD(arma::cube Eta, arma::mat Lambda, arma::vec Sigma1, arma::vec Sigma2, arma::vec Delta, arma::cube Theta);
arma::mat updateGammaSmoothDSig(arma::cube Eta, arma::mat Lambda, arma::mat Sigma, arma::vec Delta, arma::cube Theta);
arma::vec updateTau(arma::mat Lambda);
arma::vec updateDeltaProd(arma::mat Lambda, arma::vec Delta);
arma::vec updateHProd(arma::cube Eta, arma::mat beta, arma::mat X, arma::vec Delta);
arma::mat convertToPrecision(arma::vec Delta, arma::uword q1, arma::uword q2);
arma::vec updateDeltaProdTemp(arma::mat Lambda, arma::vec Delta, double a1, double a2);
double updatea1(arma::vec Delta, double a1);
double updatea2(arma::vec Delta, double a2);
double updateHb(arma::mat H);
Rcpp::List getStatistics(arma::mat spline, Rcpp::List mod, arma::rowvec x, arma::vec trueMean, int burnin, int cov);
Rcpp::List eigenLF(arma::mat splineS, arma::mat splineT, Rcpp::List mod, arma::uword numeig, int burnin);
Rcpp::List eigenLFChains(arma::mat splineS, arma::mat splineT, Rcpp::List mod, arma::uword numeig, int iter, int burnin, int nchains);
arma::mat FuncProcess(arma::mat Y, arma::mat splineS, arma::mat splineT);
arma::mat LongProcess(arma::mat Y, arma::mat splineS, arma::mat splineT);
void updateEta3Sig(arma::mat Gamma, arma::mat Lambda, arma::mat Sigma, arma::cube Theta, arma::mat H, arma::mat X, arma::mat Beta, arma::cube &eta);
void updateGammaSig(arma::cube &eta, arma::mat &Lambda, arma::vec Deltastar,
                    arma::mat &Phistar, arma::mat Sigma,
                    arma::cube &theta, arma::mat &Gamma);
void updateLambdaSig(arma::cube &eta, arma::mat &Gamma, arma::vec Delta,
                     arma::mat &Phi, arma::mat Sigma, arma::cube &theta, arma::mat &Lambda);
arma::cube updateThetaSig(arma::mat &Lambda, arma::mat &Gamma, arma::mat Sigma, arma::cube &eta, arma::mat &splineS,
                          arma::mat &splineT, arma::mat &y, double varphi);
arma::mat updateSigma(arma::mat Lambda, arma::mat Gamma, arma::cube Theta, arma::cube Eta);
arma::vec updateDeltaProdPCA(arma::mat Lambda, arma::vec Delta, arma::vec Alpha, double a1, double a2);
arma::vec updateAlphaPCA(arma::mat Lambda, arma::vec Delta, arma::vec Alpha);
arma::mat updateLambdaPCA(arma::cube Eta, arma::mat Gamma, arma::mat Sigma, arma::vec Delta, arma::vec Alpha, arma::cube Theta);
arma::mat updateGammaPCA(arma::cube Eta, arma::mat Lambda, arma::mat Sigma, arma::vec Delta, arma::vec Alpha, arma::cube Theta);

extern int p;
extern int k;
#endif
