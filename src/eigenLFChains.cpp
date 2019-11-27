#include <RcppArmadillo.h>
#include "updateParam.h"
#include <typeinfo>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]

Rcpp::List eigenLFChains(arma::mat splineS, arma::mat splineT, Rcpp::List mod, arma::uword numeig, int iter, int burnin, int nchains){
  arma::field<arma::cube> LambdaF = mod["Lambda"];
  arma::field<arma::cube> GammaF = mod["Gamma"];
  arma::field<arma::cube> BetaF = mod["Beta"];
  arma::field<arma::cube> HF = mod["H"];
  arma::field<arma::cube> SigmaF = mod["Sigma"];
  arma::mat spline = arma::kron(splineS, splineT);
  arma::mat cov(spline.n_rows, spline.n_rows);
  arma::mat postcov = arma::zeros<arma::mat>(spline.n_rows, spline.n_rows);
  //arma::mat theta_postcov = arma::zeros<arma::mat>(spline.n_cols, spline.n_cols);
  arma::mat eigvalFunc(splineT.n_rows, nchains * (iter - burnin));
  arma::cube eigvecFunc(splineT.n_rows, numeig, nchains * (iter - burnin));
  arma::vec eigvalFunc_temp;
  arma::mat eigvecFunc_temp;
  arma::mat eigvecFuncmean = arma::zeros<arma::mat>(splineT.n_rows, numeig);
  arma::mat eigvalLong(splineS.n_rows, nchains * (iter - burnin));
  arma::cube eigvecLong(splineS.n_rows, numeig, nchains * (iter - burnin));
  arma::vec eigvalLong_temp;
  arma::mat eigvecLong_temp;
  arma::mat eigvecLongmean = arma::zeros<arma::mat>(splineS.n_rows,numeig);
  arma::mat marginalFunc;
  arma::mat marginalLong;
  arma::mat meanM(spline.n_rows, nchains * (iter - burnin));
  arma::vec mmean(nchains * (iter- burnin));
  //arma::mat thetacov;
  /*
  cov = spline * (arma::kron(GammaC.slice(burnin), LambdaC.slice(burnin)) * arma::inv(arma::diagmat(HC.slice(burnin))) *
  arma::trans(arma::kron(GammaC.slice(burnin), LambdaC.slice(burnin))) + arma::kron(arma::diagmat(1.0 / sigma2.col(burnin)),
              arma::diagmat(1.0 / sigma1.col(burnin)))) * arma::trans(spline);
  */
  //cov = spline * (arma::kron(GammaF.slice(burnin), LambdaC.slice(burnin)) * arma::inv(arma::diagmat(HC.slice(burnin))) *
    //arma::trans(arma::kron(GammaC.slice(burnin), LambdaC.slice(burnin)))) * spline.t() + spline * arma::inv(arma::diagmat(arma::vectorise(Sigma.slice(burnin)))) * spline.t();
  /*
  cov = spline * (arma::kron(GammaC.slice(burnin), LambdaC.slice(burnin)) * arma::inv(convertToPrecision(Delta.col(burnin), LambdaC.n_cols, GammaC.n_cols)) *
  arma::trans(arma::kron(GammaC.slice(burnin), LambdaC.slice(burnin))) + arma::kron(arma::diagmat(1.0 / sigma2.col(burnin)),
              arma::diagmat(1.0 / sigma1.col(burnin)))) * arma::trans(spline);
  */
  marginalFunc = getMarginalFunc(cov, splineS.n_rows, splineT.n_rows);
  marginalLong = getMarginalLong(cov, splineS.n_rows, splineT.n_rows);
  arma::eig_sym(eigvalFunc_temp, eigvecFunc_temp, marginalFunc);
  arma::eig_sym(eigvalLong_temp, eigvecLong_temp, marginalLong);
  
  eigvalFunc.col(0) = eigvalFunc_temp / arma::accu(eigvalFunc_temp);
  eigvalLong.col(0) = eigvalLong_temp / arma::accu(eigvalLong_temp);
  eigvecFunc.slice(0) = eigvecFunc_temp.cols(splineT.n_rows-1-numeig+1, splineT.n_rows-1);
  
  eigvecLong.slice(0) = eigvecLong_temp.cols(splineS.n_rows - 1 - numeig+1, splineS.n_rows - 1);
  //eigvec = eigvec_temp.cols(splineT.n_rows-1-numeig+1, splineT.n_rows-1);
  //eigvecsq = arma::square(eigvec);
  meanM.col(0) = spline * arma::kron(GammaF(0).slice(burnin),
            LambdaF(0).slice(burnin)) * arma::trans(BetaF(0).slice(burnin));
  for(int k = 0; k < nchains; k++){
    for(int i = 0; i < iter - burnin; i++){
      meanM.col(i + k * (iter - burnin)) = spline * arma::kron(GammaF(k).slice(burnin + i), LambdaF(k).slice(burnin + i)) * arma::trans(BetaF(k).slice(burnin + i));
      cov = spline * (arma::kron(GammaF(k).slice(burnin + i), LambdaF(k).slice(burnin + i)) * arma::inv(arma::diagmat(HF(k).slice(burnin + i))) *
        arma::trans(arma::kron(GammaF(k).slice(burnin + i), LambdaF(k).slice(burnin + i)))) * spline.t() + spline * arma::inv(arma::diagmat(arma::vectorise(SigmaF(k).slice(burnin + i)))) * spline.t();
      postcov = postcov + cov;
      marginalFunc = getMarginalFunc(cov, splineS.n_rows, splineT.n_rows);
      marginalLong = getMarginalLong(cov, splineS.n_rows, splineT.n_rows);
      arma::eig_sym(eigvalFunc_temp, eigvecFunc_temp, marginalFunc);
      arma::eig_sym(eigvalLong_temp, eigvecLong_temp, marginalLong);
      if(i % 100 == 0){
        Rcpp::Rcout << i << std::endl;
      }
      if(i > 1){
        for(arma::uword j = 0; j < numeig; j++){
          if(arma::norm(eigvecFunc_temp.col(splineT.n_rows-1-numeig+j+1)*(-1)-eigvecFuncmean.col(j)) <
            arma::norm(eigvecFunc_temp.col(splineT.n_rows-1-numeig+j+1)-eigvecFuncmean.col(j))){
            eigvecFunc_temp.col(splineT.n_rows-1-numeig+j+1) = eigvecFunc_temp.col(splineT.n_rows-1-numeig+j+1) * (-1);
          }
          if(arma::norm(eigvecLong_temp.col(splineS.n_rows - 1-numeig+j+1)*(-1)-eigvecLongmean.col(j)) <
            arma::norm(eigvecLong_temp.col(splineS.n_rows - 1-numeig+j+1)-eigvecLongmean.col(j))){
            eigvecLong_temp.col(splineS.n_rows - 1-numeig+j+1) = eigvecLong_temp.col(splineS.n_rows - 1-numeig+j+1) * (-1);
          }
        }
        
        eigvecFuncmean = eigvecFunc_temp.cols(splineT.n_rows-1-numeig+1, splineT.n_rows-1) + eigvecFuncmean;
        eigvecLongmean = eigvecLong_temp.cols(splineS.n_rows - 1-numeig+1, splineS.n_rows - 1) + eigvecLongmean;
      }
      eigvalFunc_temp = eigvalFunc_temp / arma::accu(eigvalFunc_temp);
      eigvalLong_temp = eigvalLong_temp / arma::accu(eigvalLong_temp);
      eigvalFunc.col(i + k * (iter - burnin)) = eigvalFunc_temp;
      eigvalLong.col(i + k * (iter - burnin)) = eigvalLong_temp;
      eigvecFunc.slice(i + k * (iter - burnin)) = eigvecFunc_temp.cols(splineT.n_rows-1-numeig+1,splineT.n_rows-1);
      eigvecLong.slice(i + k * (iter - burnin)) = eigvecLong_temp.cols(splineS.n_rows - 1-numeig+1,splineS.n_rows - 1);
    }
  }
  postcov = postcov / (nchains * (iter - burnin));
  
  arma::colvec postmean = arma::mean(meanM, 1);
  eigvecFuncmean = arma::mean(eigvecFunc, 2);
  eigvecLongmean = arma::mean(eigvecLong, 2);
  arma::colvec postsd = arma::trans(arma::stddev(arma::trans(meanM), 1));
  
  arma::mat eigvecFuncsd(splineT.n_rows,numeig);
  arma::mat eigvecLongsd(splineS.n_rows,numeig);
  arma::vec vtempFunc(nchains * (iter - burnin));
  arma::vec vtempLong(nchains * (iter - burnin));
  
  for(arma::uword i = 0; i < splineT.n_rows; i++){
    for(arma::uword j = 0; j < numeig; j++){
      vtempFunc = eigvecFunc.tube(i, j);
      eigvecFuncsd(i, j) = stddev(vtempFunc);
    }
  }
  for(arma::uword i = 0; i < splineS.n_rows; i++){
    for(arma::uword j = 0; j < numeig; j++){
      vtempLong = eigvecLong.tube(i, j);
      eigvecLongsd(i, j) = stddev(vtempLong);
    }
  }
  arma::mat eigvecFuncm(eigvecFunc.n_slices,numeig);
  arma::mat eigvecLongm(eigvecLong.n_slices,numeig);
  for(arma::uword i = 0; i < eigvecFunc.n_slices; i++){
    mmean(i) = arma::max(arma::abs(meanM.col(i) - postmean) / postsd);
    for(arma::uword j = 0; j < numeig; j++){
      eigvecFuncm(i, j) = arma::max(arma::abs(eigvecFunc.slice(i).col(j)-eigvecFuncmean.col(j)) / eigvecFuncsd.col(j));
      eigvecLongm(i, j) = arma::max(arma::abs(eigvecLong.slice(i).col(j)-eigvecLongmean.col(j)) / eigvecLongsd.col(j));
    }
  }
  Rcpp::Environment base("package:stats");
  Rcpp::Function quantile_r = base["quantile"];
  Rcpp::NumericVector a =  quantile_r(mmean, .975);
  arma::mat eigvecFuncupper(splineT.n_rows, numeig);
  arma::mat eigvecFunclower(splineT.n_rows, numeig);
  arma::mat eigvecLongupper(splineS.n_rows, numeig);
  arma::mat eigvecLonglower(splineS.n_rows, numeig);
  
  for(arma::uword j = 0; j < numeig; j++){
    a = quantile_r(eigvecFuncm.col(j), .95);
    eigvecFuncupper.col(j) = eigvecFuncmean.col(j) + a[0] * eigvecFuncsd.col(j);
    eigvecFunclower.col(j) = eigvecFuncmean.col(j) - a[0] * eigvecFuncsd.col(j);
    a = quantile_r(eigvecLongm.col(j), .95);
    eigvecLongupper.col(j) = eigvecLongmean.col(j) + a[0] * eigvecLongsd.col(j);
    eigvecLonglower.col(j) = eigvecLongmean.col(j) - a[0] * eigvecLongsd.col(j);
  }
  arma::vec upper = postmean + a[0] * postsd;
  arma::vec lower = postmean - a[0] * postsd;
  
  return Rcpp::List::create(//Rcpp::Named("meanM", meanM),
    Rcpp::Named("postmean", postmean),
    //                        Rcpp::Named("postsd", postsd),
                          Rcpp::Named("upper", upper),
                        Rcpp::Named("lower", lower),
    //                    Rcpp::Named("eigvecFunc", eigvecFunc),
    Rcpp::Named("eigvecFuncmean", eigvecFuncmean),
    //                Rcpp::Named("eigvecFuncsd", eigvecFuncsd),
                  Rcpp::Named("eigvecFunclower", eigvecFunclower),
                Rcpp::Named("eigvecFuncupper", eigvecFuncupper),
    Rcpp::Named("eigvecLongmean", eigvecLongmean),
              Rcpp::Named("eigvecLonglower", eigvecLonglower),
            Rcpp::Named("eigvecLongupper", eigvecLongupper),
    //Rcpp::Named("eigvecLong", eigvecLong),
    Rcpp::Named("eigvalFunc", eigvalFunc),
    Rcpp::Named("eigvalLong", eigvalLong),
    Rcpp::Named("postcov", postcov));
  //Rcpp::Named("postcov_SE", postcov_SE));
}

