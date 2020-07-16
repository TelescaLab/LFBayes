# Generate loading matrix for brownian bridge
Loading.Brown.Bridge <- function(t, p, k){
  B <- bs(t, df = p, intercept = TRUE)
  Loading <- matrix(nrow = p, ncol = k)
  for(i in 1:k){
    eigval <- 1/(i^2*pi^2)
    psi <- sqrt(2)*sin(i*pi*t)
    Loading[,i] <- sqrt(eigval)*solve(t(B)%*%B)%*%t(B)%*%psi
  }
  return(Loading)
}

# Matern covariance function
Matern.Cov <- function(s){
  rho <- 0.5
  sigmasq <- 1
  Matern.Cov <- matrix(nrow = length(s), ncol = length(s))
  for(i in 1:length(s)){
    for(j in 1:length(s)){
      d <- abs(s[i] - s[j])
      Matern.Cov[i, j] <- sigmasq * (1 + sqrt(3) * d / rho) * exp(- sqrt(3) * d / rho)
      #Matern.Cov[i, j] <- (1 + sqrt(5) * d / rho + 5 * d^2 / (3*rho^2)) * exp(-sqrt(5) * d / rho)
    }
  }
  Matern.Cov
}

# Generate loading matrix for matern covariance
Loading.Matern <- function(s, p, k, B){
  rho <- 0.5
  sigmasq <- 1
  Matern.Cov <- matrix(nrow = length(s), ncol = length(s))
  for(i in 1:length(s)){
    for(j in 1:length(s)){
      d <- abs(s[i] - s[j])
      Matern.Cov[i, j] <- sigmasq * (1 + sqrt(3) * d / rho) * exp(- sqrt(3) * d / rho)
    }
  }
  Loading <- matrix(nrow = p, ncol = k)
  evec <- eigen(Matern.Cov)$vectors
  eval <- eigen(Matern.Cov)$values
  for(i in 1:k){
    Loading[,i] <- sqrt(eval[i]) * solve(t(B)%*%B)%*%t(B)%*%evec[,i]
  }
  Loading
}

# Generate H values. If H is the identity matrix, the model is weakly separable
GenerateH <- function(q1,q2){
  H <- matrix(0,q1, q2)
  for(i in 1:q1){
    for(j in 1:q2){
      H[i,j] <- exp(-(sqrt(.01*i) + sqrt(.1*j)))
    }
  }
  H <- diag(c(H))
  H
}

# Mean function
GenerateMu1 <- function(s,t){
  mu <- matrix(0, nrow= length(t),ncol=length(s))
  for(i in 1:length(t)){
    for(j in 1:length(s)){
      mu[i,j] <- sqrt(1/(5*sqrt(s[j]+1)))*sin(5*t[i])
    }
  }
  c(mu)
}