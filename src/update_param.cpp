#include <RcppArmadillo.h>
#include <cmath>
#include <typeinfo>
#ifdef _OPENMP
#include <omp.h>
#endif

double gam_trunc_left(double a, double b,  double cut);
void updateBeta(arma::cube &eta, arma::mat &H, arma::mat &E,
                arma::mat &X, arma::mat &beta){
  int q = eta.n_rows;
  int qstar = eta.n_cols;
  // X is design matrix. Each row for each patient. 
  // Number of columns for number of covariates.
  int d = X.n_cols;
  
  // Prior mean
  
  // Thinking of M&m formula
  arma::mat Minv = arma::zeros(d, d);
  arma::mat mychol;
  arma::vec m = arma::zeros(d);
  int idx = 0;
  arma::vec b, w, mu, z, v;
  
  for(arma::uword(i) = 0; i < arma::uword(q); i++){
    for(arma::uword(j) = 0; j < arma::uword(qstar); j++){
      idx = arma::sub2ind(arma::size(eta.slice(0)), i, j);
      Minv = H(idx,idx)*(arma::trans(X) * X) + arma::diagmat(E.row(idx));
      mychol = arma::trans(arma::chol(Minv));
      b = H(idx,idx)*arma::trans(X) * arma::vec(eta.tube(i, j));
      w = arma::solve(arma::trimatl(mychol), b);
      mu = arma::solve(arma::trimatu(arma::trans(mychol)), w);
      z = arma::randn<arma::vec>(d);
      v = arma::solve(arma::trimatu(arma::trans(mychol)), z);
      beta.col(idx) = mu + v;
    }
  }
}


arma::vec updateDelta(arma::mat &Lambda, arma::mat &Phi,
                      arma::vec Delta, double a1, double a2){
  arma::vec output;
  output.copy_size(Delta);
  double p = Lambda.n_rows;
  double k = Lambda.n_cols;
  
  arma::uword num = Delta.n_elem;
  arma::umat A(num-1,num-2);
  for(arma::uword col = 0; col < num-2; ++col){
    for(arma::uword row = 0; row < num - 1; ++row){
      if(row > col){
        A(row, col) = col + 1;
      }
      else {
        A(row, col) = col + 2;
      }
    }
  }
  A.insert_cols(0, arma::zeros<arma::uvec>(num-1));
  
  arma::rowvec colsums = sum(Phi % arma::square(Lambda), 0);
  arma::vec tauH = arma::cumprod(Delta);
  
  // temporary for delta
  arma::vec tempD = tauH;
  tempD(0) = 1.0;
  Delta(0) = R::rgamma(a1 + p*k / 2.0, 1.0 /
    (1.0 + 1.0 / 2.0 * arma::dot(tempD, colsums)));
  // temporary for other delta, reusing tempD
  for(arma::uword i = 1; i < Delta.n_elem; ++i){
    tempD = cumprod(Delta(A.row(i - 1)));
    double tempSum = arma::dot(tempD.tail(Delta.n_elem - i),
                               colsums.tail(Delta.n_elem - i));
    Delta(i) = gam_trunc_left(a2 + p*(k-i)/2.0, 1.0 /
      (1.0 + 1.0/2.0 * tempSum), 1.0);
  }
  return Delta;
}
arma::mat updateE(arma::mat beta){
  double a = 1.0 / 2.0;
  double b = 1.0 / 2.0;
  int q = beta.n_cols;
  int d = beta.n_rows;
  arma::mat E(q, d);
  for(int i = 0; i < q; i++){
    for(int j = 0; j < d; j++){
      E(i, j) = R::rgamma(a + 1.0 / 2.0, 1.0 / (b + ::pow(beta(j,i), 2.0)));
    }
  }
  
  return E;
  
}






void updateEta3Sig(arma::mat Gamma, arma::mat Lambda,
                   arma::mat Sigma, arma::cube Theta,
                   arma::mat H, arma::mat X,
                   arma::mat Beta, arma::cube &eta){
  int q = arma::size(Lambda, 1);
  int qstar = arma::size(Gamma, 1);
  arma::mat sigmat = arma::diagmat(arma::vectorise(Sigma));
  arma::mat precision = arma::kron(Gamma.t(), Lambda.t()) *
    arma::diagmat(sigmat) * arma::kron(Gamma, Lambda) + H;
  arma::mat mychol = arma::trans(arma::chol(precision));
  arma::vec mu;
  arma::vec b;
  arma::vec a;
  arma::vec w;
  arma::vec z;
  arma::vec v;
  for(arma::uword i = 0; i < arma::size(Theta, 2); i++){
    a = sigmat * arma::vectorise(Theta.slice(i));
    b = H * arma::trans(Beta) * arma::trans(X.row(i)) +
      arma::vectorise(Lambda.t() *
      arma::reshape(a, Lambda.n_rows, Gamma.n_rows) * Gamma);
    z = arma::randn<arma::vec>(arma::uword(q * qstar));
    w = arma::solve(arma::trimatl(mychol), b);
    mu = arma::solve(arma::trimatu(arma::trans(mychol)), w);
    v = arma::solve(arma::trimatu(arma::trans(mychol)), z);
    eta.slice(i) = arma::reshape(mu + v, q, qstar);
    
  }
}
void updateGammaSig(arma::cube &eta, arma::mat &Lambda, arma::vec Deltastar,
                    arma::mat &Phistar, arma::mat Sigma,
                    arma::cube &theta, arma::mat &Gamma){
  
  arma::rowvec cumDelta = 
    arma::conv_to<arma::rowvec>::from(arma::cumprod(Deltastar));
  arma::mat muMat = arma::zeros<arma::mat>(Phistar.n_rows, Phistar.n_cols);
  arma::mat mychol;
  arma::mat precision(eta.n_cols, eta.n_cols);
  precision.zeros();
  arma::vec b, w, mu, z, v;
  arma::uword subject;
  for(arma::uword i = 0; i < Phistar.n_rows; i++){
    precision = arma::diagmat(cumDelta % Phistar.row(i));
    for( subject = 0; subject < eta.n_slices; subject++){
      precision = precision + eta.slice(subject).t() *
        Lambda.t() * arma::diagmat(Sigma.col(i)) * Lambda * eta.slice(subject);
      
      muMat.row(i) = muMat.row(i) + 
        theta.slice(subject).col(i).t() * 
        arma::diagmat(Sigma.col(i)) * Lambda * eta.slice(subject);
    }
    mychol = arma::trans(arma::chol(precision));
    b = arma::trans(muMat.row(i));
    w = arma::solve(arma::trimatl(mychol), b);
    mu = arma::solve(arma::trimatu(arma::trans(mychol)), w);
    z = arma::randn<arma::vec>(Phistar.n_cols);
    v = arma::solve(arma::trimatu(arma::trans(mychol)), z);
    Gamma.row(i) = arma::trans(mu + v);
    precision.zeros();
  }
}

arma::mat updateH(arma::cube &eta, arma::mat &beta, arma::mat &X){
  // Gives the inverse of H
  double n = eta.n_slices;
  double a = 1;
  double b = 1;
  int q = eta.n_rows;
  int qstar = eta.n_cols;
  arma::mat H(q, qstar);
  arma::vec z;
  for(arma::uword i = 0; i < arma::uword(q); i++){
    for(arma::uword j = 0; j < arma::uword(qstar); j++){
      z = arma::vec(eta.tube(i,j)) - X * 
        beta.col(arma::sub2ind(arma::size(eta.slice(0)), i, j));
      H(i, j) = R::rgamma(a + 2*i + 2*j + 
        n/2.0, 1.0 / (b + 1.0/2.0 * arma::sum(z % z)));
    }
  }
  return arma::diagmat(arma::vectorise(H));
  
}
void updateLambdaSig(arma::cube &eta, arma::mat &Gamma,
                     arma::vec Delta,arma::mat &Phi,
                     arma::mat Sigma, arma::cube &theta,
                     arma::mat &Lambda){
  
  arma::rowvec cumDelta = 
    arma::conv_to<arma::rowvec>::from(arma::cumprod(Delta));
  arma::mat precision(eta.n_rows, eta.n_rows);
  arma::mat muMat = arma::zeros<arma::mat>(Phi.n_rows, Phi.n_cols);
  arma::mat mychol;
  arma::vec b, w, mu, z, v;
  for(arma::uword i = 0; i < Phi.n_rows; ++i){
    precision = arma::diagmat(cumDelta % Phi.row(i));
    for(arma::uword subject = 0; subject < eta.n_slices; ++subject){
      muMat.row(i) = muMat.row(i) + theta.slice(subject).row(i) *
        arma::diagmat(Sigma.row(i)) * Gamma * eta.slice(subject).t();
      precision = precision + eta.slice(subject) * Gamma.t() * 
        arma::diagmat(Sigma.row(i)) * Gamma * eta.slice(subject).t();
    }
    mychol = arma::trans(arma::chol(precision));
    b = arma::trans(muMat.row(i));
    w = arma::solve(arma::trimatl(mychol), b);
    mu = arma::solve(arma::trimatu(arma::trans(mychol)), w);
    z = arma::randn<arma::vec>(Phi.n_cols);
    v = arma::solve(arma::trimatu(arma::trans(mychol)), z);
    Lambda.row(i) = arma::trans(mu + v);
    precision.zeros();
  }
}







arma::mat updatePhi(arma::mat &Lambda, arma::vec delta){
  double v = 5;
  arma::vec tauH = arma::cumprod(delta);
  arma::mat Phi(Lambda.n_rows, Lambda.n_cols);
  for(arma::uword i = 0; i < Phi.n_cols; ++i){
    for(arma::uword j = 0; j < Phi.n_rows; ++j){
      Phi(j,i) = R::rgamma((v + 1.0) / 2.0, 1.0 / (v /
        2.0 + 1.0 / 2.0 * tauH(i) * ::pow(Lambda(j,i), 2.0)));
    }
  }
  return Phi;
}


arma::mat updatePhistar(arma::mat &Gamma, arma::vec deltastar){
  double v = 5;
  arma::vec tauH = arma::cumprod(deltastar);
  arma::mat Phistar(Gamma.n_rows, Gamma.n_cols);
  for(arma::uword i = 0; i < Phistar.n_cols; ++i){
    for(arma::uword j = 0; j < Phistar.n_rows; ++j){
      Phistar(j,i) = R::rgamma((v + 1.0) / 2.0, 1.0 / (v /
        2.0 + 1.0/2.0 * tauH(i) * ::pow(Gamma(j,i), 2.0)));
    }
  }
  return Phistar;
}
arma::mat updateSigma(arma::mat Lambda, arma::mat Gamma,
                      arma::cube Theta, arma::cube Eta){
  double asig = .5;
  double bsig = .25;
  arma::uword p1 = arma::size(Lambda,0);
  arma::uword p2 = arma::size(Gamma, 0);
  arma::mat Sigma(p1,p2);
  double ztz = 0;
  for(arma::uword i=0; i < p2; i++){
    for(arma::uword j=0; j < p1;j++){
      for(arma::uword subject=0; subject<Theta.n_slices;subject++){
        ztz = ztz + ::pow(arma::as_scalar(Theta.slice(subject).col(i).row(j) - Lambda.row(j) * Eta.slice(subject) * arma::trans(Gamma.row(i))),2);
      }
      Sigma(j,i) = R::rgamma(asig + double(Theta.n_slices)/2.0, 1.0 / (bsig + 1.0/2.0 * ztz));
      ztz = 0;
    }
  }
  return(Sigma);
}

arma::cube updateThetaSig(arma::mat &Lambda, arma::mat &Gamma,
                          arma::mat Sigma, arma::cube &eta,
                          arma::mat &splineS,arma::mat &splineT,
                          arma::mat &y, double varphi){
  arma::uword p = Sigma.n_rows;
  arma::uword pstar = Sigma.n_cols;
  arma::uword n = arma::size(y, 1);
  arma::mat precision;
  //arma::mat cov;
  arma::cube theta(p, pstar, n);
  arma::mat sigmat = arma::diagmat(arma::vectorise(Sigma));
  precision = varphi * arma::kron(arma::trans(splineS) *
    splineS, arma::trans(splineT) * splineT) + arma::diagmat(sigmat);
  arma::mat mychol = arma::chol(precision, "lower");
  arma::vec b, w, mu, z, v;
  for(arma::uword i = 0; i < n; i++){
    b = (varphi * arma::vectorise(splineT.t() * 
      arma::reshape(y.col(i), splineT.n_rows, splineS.n_rows) * splineS) +
      sigmat * arma::vectorise(Lambda * eta.slice(i) * Gamma.t()));
    w = arma::inv(arma::trimatl(mychol)) * b;
    mu = arma::inv(arma::trimatu(arma::trans(mychol))) * w;
    z = arma::randn<arma::vec>(arma::uword(p * pstar));
    v = arma::inv(arma::trimatu(arma::trans(mychol))) * z;
    theta.slice(i) = arma::reshape(mu + v, p, pstar);
  }
  return theta;
}

double updateVarphi(arma::cube &theta, arma::mat &splineS, arma::mat &splineT, arma::mat &y){
  double aVarphi = .0001;
  double bVarphi = .0001;
  arma::uword n = theta.n_slices;
  double ztz = 0;
  double totalN = 0;

  arma::vec z;
  for(arma::uword i = 0; i < arma::uword(n); i++){
    arma::vec z;
    z = y.col(i) - arma::vectorise(splineT * theta.slice(i) * splineS.t());
    ztz = ztz + arma::sum(z % z);
    totalN = totalN + splineS.n_rows * splineT.n_rows;
  }
  return R::rgamma(aVarphi + totalN/2, 1.0 / (bVarphi + 1.0/2.0 * ztz));
}

arma::vec initializeY(arma::vec x, arma::vec missing, int s, int t){
  arma::vec y(s*t, arma::fill::zeros);
  int counter = 0;
  for(int i = 0; i < s; i++){
    if(arma::any(missing == i+1)){
      y.subvec(t * i, t * i + t - 1).zeros();
    } else {
      y.subvec(t * i, t * i + t - 1) = x.subvec(t *
        counter,  t * counter + t - 1);
      counter = counter + 1;
    }
  }
  return(y);
}

void updateMissing(arma::mat &y, arma::field<arma::vec> &missing, arma::cube &theta, double Varphi, arma::mat &splineS, arma::mat &splineT, int n){
  for(int i = 0; i < n; i++){
    
    arma::vec currentMissing = missing(i);
    for(arma::uword m = 0; m < currentMissing.n_elem; m++){
      y.col(i).subvec(splineT.n_rows * (currentMissing(m)-1), 
            splineT.n_rows * (currentMissing(m)-1) + splineT.n_rows - 1) =
              arma::vectorise(splineT * theta.slice(i) * 
              arma::trans(splineS.row(currentMissing(m)-1))) +
        1.0/Varphi * arma::randn<arma::vec>(splineT.n_rows);
    }
  }
}


double updatea1(arma::vec Delta, double a1){
  double a = 1;
  double b = 1;
  double proposal;
  while(true){
    proposal = a1 + R::rnorm(0,1);
    if(proposal > 0){
      break;
    }
  }
  double A = log(R::dgamma(Delta(0), proposal, 1/b, 0)) +
    log(R::dgamma(proposal, a, 1/b, 0)) +
    log(R::pnorm(a1, 0.0, 1.0, 1, 0)) - 
    log(R::dgamma(Delta(0), a1, 1/b, 0)) - log(R::dgamma(a1, a, 1/b, 0)) -
    log(R::pnorm(proposal, 0.0, 1.0, 1, 0));
  if(R::runif(0.0, 1.0) < exp(A)){
    a1 = proposal;
  }
  return(a1);
}

double updatea2(arma::vec Delta, double a2){
  double a = 2.0;
  double b = 1.0;
  double proposal;
  int q = arma::size(Delta, 0);
  while(true){
    proposal = a2 + R::rnorm(0,1);
    if(proposal > 0){
      break;
    }
  }
  Rcpp::NumericVector D = Rcpp::wrap(Delta.tail(q - 1));
  double A = Rcpp::sum(Rcpp::log(Rcpp::dgamma(D, proposal, b, false))) +
    log(R::dgamma(proposal, a, 1/b, 0)) +
    log(R::pnorm(a2, 0.0, 1.0, 1, 0)) - 
    Rcpp::sum(Rcpp::log(Rcpp::dgamma(D, a2, b, false))) - 
    log(R::dgamma(a2, a, 1/b, 0)) - log(R::pnorm(proposal, 0.0, 1.0, 1, 0));
  if(R::runif(0.0, 1.0) < exp(A)){
    a2 = proposal;
  }
  return(a2);
}



arma::mat FuncProcess(arma::mat Y, arma::mat splineS, arma::mat splineT){
  arma::mat result(Y.n_cols * splineS.n_rows, splineT.n_rows);
  for(int i = 0; i < Y.n_cols; i++){
    for(int s = 0; s < splineS.n_rows; s++){
      result.row(s + i * splineS.n_rows) = (Y.submat(s * 
        splineT.n_rows,i,(s+1) * splineT.n_rows - 1,i)).t();
    }
  }
  return(result);
}

arma::mat LongProcess(arma::mat Y, arma::mat splineS, arma::mat splineT){
  arma::mat result(Y.n_cols * splineT.n_rows, splineS.n_rows);
  return(arma::reshape(Y,Y.n_cols * splineT.n_rows, splineS.n_rows));
}



double gam_trunc_left(double a, double b,  double cut){ 
  double u, pcut, y; 
  
  pcut = R::pgamma(cut,a,b,1,0);
  if(pcut>0.99){
    return(cut);
  } 
  u = R::runif(0.0,1.0); 
  u = pcut + (1-pcut)*u; 
  y = R::qgamma(u, a, b, 1, 0);
  return y; 
} 