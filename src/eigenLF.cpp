#include <RcppArmadillo.h>
#include "updateParam.h"
#include <typeinfo>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]

Rcpp::List eigenLF(arma::mat splineS, arma::mat splineT, Rcpp::List mod, arma::uword numeig, int burnin){
  arma::cube LambdaC = mod["Lambda"];
  arma::cube GammaC = mod["Gamma"];
  arma::cube betaC = mod["Beta"];
  //arma::mat Delta = mod["Delta"];
  //arma::field<arma::cube> theta = mod["theta"];
  arma::cube HC = mod["HC"];
  //arma::vec varphi = mod["varphi"];
  //arma::mat sigma1 = mod["sigma1"];
  //arma::mat sigma2 = mod["sigma2"];
  arma::cube Sigma = mod["Sigma"];
  int iter = LambdaC.n_slices;
  arma::mat spline = arma::kron(splineS, splineT);
  arma::mat cov(spline.n_rows, spline.n_rows);
  arma::mat postcov = arma::zeros<arma::mat>(spline.n_rows, spline.n_rows);
  //arma::mat theta_postcov = arma::zeros<arma::mat>(spline.n_cols, spline.n_cols);
  arma::mat eigvalFunc(splineT.n_rows, iter - burnin);
  arma::cube eigvecFunc(splineT.n_rows, numeig, iter-burnin);
  arma::vec eigvalFunc_temp;
  arma::mat eigvecFunc_temp;
  arma::mat eigvecFuncmean = arma::zeros<arma::mat>(splineT.n_rows, numeig);
  arma::mat eigvalLong(splineS.n_rows, iter - burnin);
  arma::cube eigvecLong(splineS.n_rows, numeig, iter - burnin);
  arma::vec eigvalLong_temp;
  arma::mat eigvecLong_temp;
  arma::mat eigvecLongmean = arma::zeros<arma::mat>(splineS.n_rows,numeig);
  arma::mat marginalFunc;
  arma::mat marginalLong;
  arma::mat meanM(spline.n_rows, iter - burnin);
  arma::vec mmean(iter- burnin);
  //arma::mat thetacov;
  
  cov = spline * (arma::kron(GammaC.slice(burnin), LambdaC.slice(burnin)) * arma::inv(arma::diagmat(HC.slice(burnin))) *
    arma::trans(arma::kron(GammaC.slice(burnin), LambdaC.slice(burnin))) + arma::inv(arma::diagmat(arma::vectorise(Sigma.slice(burnin))))) * arma::trans(spline);
  /*
  cov = spline * (arma::kron(GammaC.slice(burnin), LambdaC.slice(burnin)) * arma::inv(arma::diagmat(HC.slice(burnin))) *
    arma::trans(arma::kron(GammaC.slice(burnin), LambdaC.slice(burnin)))) * spline.t() + spline * arma::inv(arma::diagmat(arma::vectorise(Sigma.slice(burnin)))) * spline.t();
  */
  
  // CHANGE BURNIN INDEXES FIX THIS CRAP
  /*
   cov = spline * (arma::kron(GammaC.slice(burnin + i), LambdaC.slice(burnin + i)) * arma::inv(convertToPrecision(Delta.col(burnin), LambdaC.n_cols, GammaC.n_cols)) *
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

  meanM.col(0) = spline * arma::kron(GammaC.slice(burnin),
            LambdaC.slice(burnin)) * arma::trans(betaC.slice(burnin));
  for(int i = 0; i < iter - burnin; i++){
    meanM.col(i) = spline * arma::kron(GammaC.slice(burnin + i),
              LambdaC.slice(burnin + i)) * arma::trans(betaC.slice(burnin + i));
    /*
    cov = spline * (arma::kron(GammaC.slice(burnin + i), LambdaC.slice(burnin + i)) * arma::inv(arma::diagmat(HC.slice(burnin + i))) *
      arma::trans(arma::kron(GammaC.slice(burnin + i), LambdaC.slice(burnin + i))) + arma::kron(arma::diagmat(1.0 / sigma2.col(burnin + i)),
                  arma::diagmat(1.0 / sigma1.col(burnin + i)))) * arma::trans(spline);
    */
    cov = spline * (arma::kron(GammaC.slice(burnin + i), LambdaC.slice(burnin + i)) * arma::inv(arma::diagmat(HC.slice(burnin + i))) *
      arma::trans(arma::kron(GammaC.slice(burnin + i), LambdaC.slice(burnin + i)))) * spline.t() + spline * arma::inv(arma::diagmat(arma::vectorise(Sigma.slice(burnin + i)))) * spline.t();
    
    //thetacov = arma::cov(arma::trans(arma::reshape( arma::mat(theta(burnin + i).memptr(), theta(burnin + i).n_elem, 1, false), theta(0).n_rows * theta(0).n_cols, 30)));
    
    postcov = postcov + cov;
    //theta_postcov = theta_postcov + thetacov;
   /*
    cov = spline * (arma::kron(GammaC.slice(burnin + i), LambdaC.slice(burnin + i)) * arma::inv(convertToPrecision(Delta.col(burnin + i), LambdaC.n_cols, GammaC.n_cols)) *
      arma::trans(arma::kron(GammaC.slice(burnin + i), LambdaC.slice(burnin + i))) + arma::kron(arma::diagmat(1.0 / sigma2.col(burnin + i)),
                  arma::diagmat(1.0 / sigma1.col(burnin + i)))) * arma::trans(spline);
     */
    marginalFunc = getMarginalFunc(cov, splineS.n_rows, splineT.n_rows);
    marginalLong = getMarginalLong(cov, splineS.n_rows, splineT.n_rows);
    arma::eig_sym(eigvalFunc_temp, eigvecFunc_temp, marginalFunc);
    arma::eig_sym(eigvalLong_temp, eigvecLong_temp, marginalLong);
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

    eigvalFunc.col(i) = eigvalFunc_temp;
    eigvalLong.col(i) = eigvalLong_temp;
    eigvecFunc.slice(i) = eigvecFunc_temp.cols(splineT.n_rows-1-numeig+1,splineT.n_rows-1);
    eigvecLong.slice(i) = eigvecLong_temp.cols(splineS.n_rows - 1-numeig+1,splineS.n_rows - 1);
  }
  postcov = postcov / (iter - burnin);
  //theta_postcov = theta_postcov / (iter - burnin);
  //arma::mat postcov_SE = spline * theta_postcov * arma::trans(spline);
  //eigvecFunc = eigvecFunc.tail_slices(eigvecFunc.n_slices - 100);
  //eigvecLong = eigvecLong.tail_slices(eigvecLong.n_slices - 100);

  arma::colvec postmean = arma::mean(meanM, 1);
  eigvecFuncmean = arma::mean(eigvecFunc, 2);
  eigvecLongmean = arma::mean(eigvecLong, 2);
  //eigvecFuncmean.zeros();
  //eigvecLongmean.zeros();
  //for(arma::uword i = 0; i < eigvecFunc.n_slices; i++){
   // eigvecFuncmean = eigvecFuncmean + eigvecFunc.slice(i) / eigvecFunc.n_slices;
   // eigvecLongmean = eigvecLongmean + eigvecLong.slice(i) / eigvecLong.n_slices;
 // }
  //eigvecFuncmean = eigvecFuncmean / (iter - burnin);
  //eigvecLongmean = eigvecLongmean / (iter - burnin);
  arma::colvec postsd = arma::trans(arma::stddev(arma::trans(meanM), 1));

  arma::mat eigvecFuncsd(splineT.n_rows,numeig);
  arma::mat eigvecLongsd(splineS.n_rows,numeig);
  arma::vec vtempFunc(iter - burnin);
  arma::vec vtempLong(iter - burnin);
  
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
      //                      Rcpp::Named("upper", upper),
        //                    Rcpp::Named("lower", lower),
        //                    Rcpp::Named("eigvecFunc", eigvecFunc),
                            Rcpp::Named("eigvecFuncmean", eigvecFuncmean),
            //                Rcpp::Named("eigvecFuncsd", eigvecFuncsd),
              //              Rcpp::Named("eigvecFunclower", eigvecFunclower),
                //            Rcpp::Named("eigvecFuncupper", eigvecFuncupper),
                            Rcpp::Named("eigvecLongmean", eigvecLongmean),
                  //          Rcpp::Named("eigvecLonglower", eigvecLonglower),
                    //        Rcpp::Named("eigvecLongupper", eigvecLongupper),
                            //Rcpp::Named("eigvecLong", eigvecLong),
                            //Rcpp::Named("eigvalFunc", eigvalFunc),
                            //Rcpp::Named("eigvalLong", eigvalLong),
                            Rcpp::Named("postcov", postcov));
                            //Rcpp::Named("postcov_SE", postcov_SE));
}

