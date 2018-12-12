op <- par(mfrow = c(2,1), mgp = c(2,.8,0), mar = .1+c(3,3,3,1))
n <- 9
x <- 1:n
y <- rnorm(n)
plot(x, y, main = paste("spline[fun](.) through", n, "points"))
lines(spline(x, y))
lines(spline(x, y, n = 201), col = 2)


y <- (x-6)^2
plot(x, y, main = "spline(.) -- 3 methods")
lines(spline(x, y, n = 201), col = 2)
lines(spline(x, y, n = 201, method = "natural"), col = 3)
lines(spline(x, y, n = 201, method = "periodic"), col = 4)
legend(6,25, c("fmm","natural","periodic"), col=2:4, lty=1)

install.packages("multisensi")
library(multisensi)
knots <- c(.2,.5,.8)
bspline(x = seq(0, 1, len = 101), k = knots, i = 1, m = 2)
?bspline


t <- seq(from = 0, to = 3, length.out = 30)
t[11] <- 3
B <- bs(t, intercept = TRUE)
y <- 3* B[,1] -1 *B[,2] + 2 * B[,3] -5 * B[,4]
plot(t,y, type = "l")
lines(t,y, col = "red")
predict(B, .95)
B

nugget <- function(s1, s2, t1, t2, beta){
  h <- s2 - s1
  u <- t2 - t1
  alpha <- .722
  if(h == 0){
    return((.901*abs(u)^(2*alpha)+1)^(-1))
  }
  return(.968*(.901*abs(u)^(2*alpha)+1)^(-1)*exp(-.134*(abs(h))/(.901*abs(u)^(2*alpha)+1)^(beta/2)))
}
nugget2 <- function(s1, s2, t1, t2){
  h <- s2 - s1
  u <- t2 - t1
  alpha <- .722
  c <- 1.5
  sig <- 1
  delta <- .5
  d <- 1
  beta <- 1
  a <- .8
  first <- sig/(a*abs(u)^(2*alpha)+1)^(delta+beta*d/2)
  second <- 1 + c*abs(h)/(a*abs(u)^(2*alpha)+1)^(beta/2)
  third <- exp(-c *abs(h)/(a*abs(u)^(2*alpha)+1)^(beta/2))
  return(first*second*third)
}

s <- seq(from = 0, to = 3, length.out = 10)
t <- seq(from = 0, to = 3, length.out = 20)
Cov <- array(dim = c(200,200))
for(s1 in 1:10){
  for(s2 in 1:10){
    for(t1 in 1:20){
      for(t2 in 1:20){
        Cov[20*(s1-1) + t1, 20*(s2-1) + t2] <- nugget2(s[s1],s[s2],t[t1],t[t2])
      }
    }
  }
}
image(Cov[,200:1], col = heat.colors(200))
Cov[200:200]
Cov
View(Cov)
