#include <RcppArmadillo.h>
#include "updateParam.h"
#include <typeinfo>
#include "utility.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List eigenLFChains(arma::mat splineS, arma::mat splineT, Rcpp::List mod, arma::uword numeig, int iter, int burnin, int nchains, arma::vec s, arma::vec t, double alpha){

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
  arma::uword q1s = LambdaF(0).n_cols;
  arma::uword q2s = GammaF(0).n_cols;
  arma::mat H, S;
  arma::mat Func_cov;
  arma::mat Long_cov;
  arma::mat Func_weight;
  arma::mat Long_weight;
  //Rcpp::Environment pracma("package:pracma");
  //Rcpp::Function trapz = pracma["trapz"];
  Rcpp::Rcout << "2" << std::endl;

  arma::vec splineT_trapz(splineT.n_cols);

  arma::vec splineS_trapz(splineS.n_cols);

  arma::vec Psi_trapz(GammaF(0).n_cols);

  arma::vec Phi_trapz(LambdaF(0).n_cols);

  arma::mat Psi(splineS.n_rows, GammaF(0).n_cols);
  arma::mat Phi(splineT.n_rows, LambdaF(0).n_cols);

  arma::mat splineS_int(splineS.n_cols, splineS.n_cols);
  arma::mat splineT_int(splineT.n_cols, splineT.n_cols);

  splineS_int = integrated(splineS, s);
  splineT_int = integrated(splineT, t);

  arma::mat splineS_int_sqrt = arma::sqrtmat_sympd(splineS_int);
  arma::mat splineS_int_sqrt_inv = arma::inv_sympd(splineS_int_sqrt);
  arma::mat splineT_int_sqrt = arma::sqrtmat_sympd(splineT_int);
  arma::mat splineT_int_sqrt_inv = arma::inv_sympd(splineT_int_sqrt);
  Rcpp::List List_func, List_long;
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
  /*
  marginalFunc = getMarginalFunc(cov, splineS.n_rows, splineT.n_rows);
  marginalLong = getMarginalLong(cov, splineS.n_rows, splineT.n_rows);
  arma::eig_sym(eigvalFunc_temp, eigvecFunc_temp, marginalFunc);
  arma::eig_sym(eigvalLong_temp, eigvecLong_temp, marginalLong);
  
  eigvalFunc.col(0) = eigvalFunc_temp / arma::accu(eigvalFunc_temp);
  eigvalLong.col(0) = eigvalLong_temp / arma::accu(eigvalLong_temp);
  eigvecFunc.slice(0) = eigvecFunc_temp.cols(splineT.n_rows-1-numeig+1, splineT.n_rows-1);
  
  eigvecLong.slice(0) = eigvecLong_temp.cols(splineS.n_rows - 1 - numeig+1, splineS.n_rows - 1);
  */
  //eigvec = eigvec_temp.cols(splineT.n_rows-1-numeig+1, splineT.n_rows-1);
  //eigvecsq = arma::square(eigvec);
  meanM.col(0) = spline * arma::kron(GammaF(0).slice(burnin),
            LambdaF(0).slice(burnin)) * arma::trans(BetaF(0).slice(burnin));
  Rcpp::Rcout << "Starting post-processing..." << std::endl;
  for(int k = 0; k < nchains; k++){
    for(int i = 0; i < iter - burnin; i++){
      if(int(i) % int(floor(double(iter - burnin) / 10.0)) == 0){
        Rcpp::Rcout << 100 * i / (iter - burnin) << '%' << std::endl;
      }

      Psi = splineS * GammaF(k).slice(burnin + i);
      Phi = splineT * LambdaF(k).slice(burnin + i);

      Psi_trapz = integrated_latent(Psi, s);
      Phi_trapz = integrated_latent(Phi, t);

      meanM.col(i + k * (iter - burnin)) = arma::vectorise(Phi * arma::reshape(BetaF(k).slice(burnin+i), Phi.n_cols, Psi.n_cols) * Psi.t());
      //cov = spline * (arma::kron(GammaF(k).slice(burnin + i), LambdaF(k).slice(burnin + i)) * arma::inv(arma::diagmat(HF(k).slice(burnin + i))) *
      //arma::trans(arma::kron(GammaF(k).slice(burnin + i), LambdaF(k).slice(burnin + i)))) * spline.t() + spline * arma::inv(arma::diagmat(arma::vectorise(SigmaF(k).slice(burnin + i)))) * spline.t();
      postcov = postcov + cov;
      //marginalFunc = getMarginalFunc(cov, splineS.n_rows, splineT.n_rows);
      //marginalLong = getMarginalLong(cov, splineS.n_rows, splineT.n_rows);
      //arma::eig_sym(eigvalFunc_temp, eigvecFunc_temp, marginalFunc);
      //arma::eig_sym(eigvalLong_temp, eigvecLong_temp, marginalLong);
      H = arma::reshape(1 / arma::diagvec(HF(k).slice(burnin + i)), q1s, q2s);
      S = 1 / SigmaF(k).slice(burnin + i);
      /* 
      H_func = arma::diagmat(H * Psi_trapz);
      H_long = arma::diagmat(H.t() * Phi_trapz);
      Func_cov = splineT_int_sqrt * (arma::diagmat(1 / SigmaF(k).slice(burnin + i) * splineS_int.diag()) + LambdaF(k).slice(burnin + i) * arma::diagmat(H_func) * arma::trans(LambdaF(k).slice(burnin + i))) * splineT_int_sqrt;
      Long_cov = splineS_int_sqrt * (arma::diagmat(1 / SigmaF(k).slice(burnin + i).t() * splineT_int.diag()) + GammaF(k).slice(burnin + i) * arma::diagmat(H_long) * arma::trans(GammaF(k).slice(burnin + i))) * splineS_int_sqrt;
      arma::eig_sym(eigvalFunc_temp, Func_weight, Func_cov);
      arma::eig_sym(eigvalLong_temp, Long_weight, Long_cov);
      
      eigvecFunc_temp = splineT * splineT_int_sqrt_inv * Func_weight.cols(splineT.n_cols - numeig, splineT.n_cols - 1);
      eigvecLong_temp = splineS * splineS_int_sqrt_inv * Long_weight.cols(splineS.n_cols - numeig, splineS.n_cols - 1);
      */
      List_func = extract_eigenfn(LambdaF(k).slice(burnin + i), S, H, splineT, splineT_int_sqrt, splineT_int_sqrt_inv, splineS_int, Psi_trapz, numeig);
      List_long = extract_eigenfn(GammaF(k).slice(burnin + i), S.t(), H.t(), splineS, splineS_int_sqrt, splineS_int_sqrt_inv, splineT_int, Phi_trapz, numeig);
      eigvecFunc_temp = Rcpp::as<arma::mat>(List_func["eigenfn"]);
      eigvalFunc_temp = Rcpp::as<arma::vec>(List_func["eigenval"]);
      eigvecLong_temp = Rcpp::as<arma::mat>(List_long["eigenfn"]);
      eigvalLong_temp = Rcpp::as<arma::vec>(List_long["eigenval"]);
      eigvecFunc_temp = arma::normalise(eigvecFunc_temp, 2, 0);
      eigvecLong_temp = arma::normalise(eigvecLong_temp, 2, 0);
      
      if(i > 1){
        for(arma::uword j = 0; j < numeig; j++){
          /*
          if(arma::norm(eigvecFunc_temp.col(splineT.n_rows-1-numeig+j+1)*(-1)-eigvecFuncmean.col(j)) <
            arma::norm(eigvecFunc_temp.col(splineT.n_rows-1-numeig+j+1)-eigvecFuncmean.col(j))){
          eigvecFunc_temp.col(splineT.n_rows-1-numeig+j+1) = eigvecFunc_temp.col(splineT.n_rows-1-numeig+j+1) * (-1);
          }
          */
          
          if(arma::norm(eigvecFunc_temp.col(j)*(-1)-eigvecFuncmean.col(j)) <
            arma::norm(eigvecFunc_temp.col(j)-eigvecFuncmean.col(j))){
            eigvecFunc_temp.col(j) = eigvecFunc_temp.col(j) * (-1);
          }
          
          /*
          if(arma::norm(eigvecLong_temp.col(splineS.n_rows - 1-numeig+j+1)*(-1)-eigvecLongmean.col(j)) <
            arma::norm(eigvecLong_temp.col(splineS.n_rows - 1-numeig+j+1)-eigvecLongmean.col(j))){
          eigvecLong_temp.col(splineS.n_rows - 1-numeig+j+1) = eigvecLong_temp.col(splineS.n_rows - 1-numeig+j+1) * (-1);
            }
          */
          if(arma::norm(eigvecLong_temp.col(j)*(-1)-eigvecLongmean.col(j)) <
            arma::norm(eigvecLong_temp.col(j)-eigvecLongmean.col(j))){
            eigvecLong_temp.col(j) = eigvecLong_temp.col(j) * (-1);
          }
        }
        eigvecFuncmean = (eigvecFunc_temp - eigvecFuncmean) / i + eigvecFuncmean;
        eigvecLongmean = (eigvecLong_temp - eigvecLongmean) / i + eigvecLongmean;
        //eigvecFuncmean = eigvecFunc_temp.cols(splineT.n_rows-1-numeig+1, splineT.n_rows-1) + eigvecFuncmean;
        //eigvecLongmean = eigvecLong_temp.cols(splineS.n_rows - 1-numeig+1, splineS.n_rows - 1) + eigvecLongmean;
      }
      eigvalFunc_temp = eigvalFunc_temp / arma::accu(eigvalFunc_temp);
      
      eigvalLong_temp = eigvalLong_temp / arma::accu(eigvalLong_temp);
      
      eigvalFunc.col(i + k * (iter - burnin)).head(splineT.n_cols) = eigvalFunc_temp;
      
      eigvalLong.col(i + k * (iter - burnin)).head(splineS.n_cols) = eigvalLong_temp;
      
      eigvecFunc.slice(i + k * (iter - burnin)) = eigvecFunc_temp;
      
      eigvecLong.slice(i + k * (iter - burnin)) = eigvecLong_temp;
      //eigvecFunc.slice(i + k * (iter - burnin)) = eigvecFunc_temp.cols(splineT.n_rows-1-numeig+1,splineT.n_rows-1);
      //eigvecLong.slice(i + k * (iter - burnin)) = eigvecLong_temp.cols(splineS.n_rows - 1-numeig+1,splineS.n_rows - 1);
      
      //H_func.zeros();
      //H_long.zeros();
      
    }
  }
  postcov = postcov / (nchains * (iter - burnin));
  
  arma::colvec postmean = arma::mean(meanM, 1);
  //eigvecFuncmean = arma::mean(eigvecFunc, 2);
  //eigvecLongmean = arma::mean(eigvecLong, 2);
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
  //Rcpp::Environment base("package:stats");
  //Rcpp::Function quantile_r = base["quantile"];
  arma::vec quantiles(1);
  quantiles(0) = 1 - alpha;
  //Rcpp::NumericVector a = arma::quantile(mmean, quantiles);
  double a = arma::as_scalar(arma::quantile(mmean, quantiles));
  arma::mat eigvecFuncupper(splineT.n_rows, numeig);
  arma::mat eigvecFunclower(splineT.n_rows, numeig);
  arma::mat eigvecLongupper(splineS.n_rows, numeig);
  arma::mat eigvecLonglower(splineS.n_rows, numeig);
  arma::vec upper = postmean + a * postsd;
  arma::vec lower = postmean - a * postsd;
  
  for(arma::uword j = 0; j < numeig; j++){
    a = arma::as_scalar(arma::quantile(eigvecFuncm.col(j), quantiles));
    eigvecFuncupper.col(j) = eigvecFuncmean.col(j) + a * eigvecFuncsd.col(j);
    eigvecFunclower.col(j) = eigvecFuncmean.col(j) - a * eigvecFuncsd.col(j);
    a = arma::as_scalar(arma::quantile(eigvecLongm.col(j), quantiles));
    eigvecLongupper.col(j) = eigvecLongmean.col(j) + a * eigvecLongsd.col(j);
    eigvecLonglower.col(j) = eigvecLongmean.col(j) - a * eigvecLongsd.col(j);
  }
  
  Rcpp::Rcout << "All done!" << std::endl;
  
  return(Rcpp::List::create(//Rcpp::Named("meanM", meanM),
    Rcpp::Named("postmean", postmean),
    Rcpp::Named("postsd", postsd),
    Rcpp::Named("upper", upper),
    Rcpp::Named("lower", lower),
    Rcpp::Named("eigvecFunc", eigvecFunc),
    Rcpp::Named("eigvecFuncmean", eigvecFuncmean),
    Rcpp::Named("eigvecFuncsd", eigvecFuncsd),
    Rcpp::Named("eigvecFunclower", eigvecFunclower),
    Rcpp::Named("eigvecFuncupper", eigvecFuncupper),
    Rcpp::Named("eigvecLongmean", eigvecLongmean),
    Rcpp::Named("eigvecLonglower", eigvecLonglower),
    Rcpp::Named("eigvecLongupper", eigvecLongupper),
    Rcpp::Named("eigvecLong", eigvecLong),
    Rcpp::Named("eigvecLongsd", eigvecLongsd),
    Rcpp::Named("eigvalFunc", eigvalFunc),
    Rcpp::Named("eigvalLong", eigvalLong),
    Rcpp::Named("postcov", postcov)));
}

