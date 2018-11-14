m <- matrix(1:30, ncol=6)
colnames(m) <- paste("C", 1:6, sep="")
rownames(m) <- paste("R", 1:5, sep="")
m
image(1:ncol(m), 1:nrow(m), t(m), col = heat.colors(60), axes = TRUE)
image(1:ncol(m), 1:nrow(m), t(m), axes = FALSE)
image(1:ncol(t(m)), 1:nrow(t(m)), m, col = heat.colors(60), axes = FALSE)

axis(1, 1:ncol(m), colnames(m))
axis(2, 1:nrow(m), rownames(m))
for (x in 1:ncol(m))
  for (y in 1:nrow(m))
    text(x, y, m[y,x])

require(fields)
# Make a 10x10 matrix
m = matrix(rnorm(100), nrow=10)
image.plot(m)
for (x in 1:10)
  for (y in 1:10)
    text((x-1)/9, (y-1)/9, sprintf("%0.2f", m[x,y]))
image.plot(m)


require(grDevices) # for colours
x <- y <- seq(-4*pi, 4*pi, len = 27)
r <- sqrt(outer(x^2, y^2, "+"))
image(z = z <- cos(r^2)*exp(-r/6), col  = gray((0:32)/32))
image(z, axes = FALSE, main = "Math can be beautiful ...",
      xlab = expression(cos(r^2) * e^{-r/6}))
contour(z, add = TRUE, drawlabels = FALSE)


# Volcano data visualized as matrix. Need to transpose and flip
# matrix horizontally.
image(t(volcano)[ncol(volcano):1,])
dim(volcano)
image(t(m), col = heat.colors(60))
m
data = data.frame(
  x = rep( c(0.1, 0.2, 0.3, 0.4, 0.5), each=5),
  y = rep( c(1, 2, 3, 4, 5), 5)
)

data$z = runif(
  25,
  min = (data$x*data$y - 0.1 * (data$x*data$y)),
  max = (data$x*data$y + 0.1 * (data$x*data$y))
)
m = matrix(runif(100),10,10)
par(mar=c(0, 0, 0, 0))
image(m, useRaster=TRUE, axes=FALSE)
library(lattice)
levelplot(t(m[5:1,]))
n <- t(m)[6:1,]
levelplot(n)
n
m
dim(Brown.Motion.Cov)
levelplot(t(Brown.Bridge.Cov[20:1,]))
image(t(Brown.Bridge.Cov[,20:1]), col = heat.colors(100))
image(t(cov2cor(MaternCov)[,20:1]), col = heat.colors(60))

#image(t(cov2cor(Brown.Bridge.Cov)[20:1,]), col = heat.colors(60))
image(t(cov2cor(Brown.Motion.Cov))[,20:1], col = heat.colors(60))
#image(t(cov2cor(Brown.Motion.Cov)[20:1,]), col = heat.colors(60))
#levelplot(t(Brown.Motion.Cov[20:1,]))
#plot(diag(Brown.Motion.Cov))
Cov.Strong <- kronecker(Brown.Motion.Cov, Brown.Bridge.Cov)
mylist <- list()
mylist[[1]] <- Smooth_scaled_cov
mylist[[5]] <- mcmc$postcov
mylist[[6]] <- EmpCov
mylist[[7]] <- resMarginalCov
mylist[[8]] <- resProductCov
#mylist[[9]] <- EmpCovBayes
colnum <- 200
image(t((Smooth_scaled_cov))[1:colnum,1:colnum][,colnum:1], zlim = c(min(unlist(mylist)), max(unlist(mylist))), col = heat.colors(100))
image(t(mcmc$postcov)[1:colnum,1:colnum][,colnum:1],zlim = c(min(unlist(mylist)), max(unlist(mylist))), col = heat.colors(100))
image(t(resPACE$fittedCov)[1:colnum,1:colnum][,colnum:1], zlim = c(min(unlist(mylist)), max(unlist(mylist))),col = heat.colors(100))
image(t((EmpCov))[1:colnum,1:colnum][,colnum:1], zlim = c(min(unlist(mylist)), max(unlist(mylist))),col = heat.colors(100))
image(t(resMarginalCov)[1:colnum,1:colnum][,colnum:1], zlim = c(min(unlist(mylist)), max(unlist(mylist))),col = heat.colors(100))
image(t(resProductCov)[1:colnum,1:colnum][,colnum:1], zlim = c(min(unlist(mylist)), max(unlist(mylist))),col = heat.colors(100))
#image(t((EmpCovBayes))[1:colnum,1:colnum][,colnum:1], zlim = c(min(unlist(mylist)), max(unlist(mylist))), col = heat.colors(100))


image(t(cov2cor(Cov.Weak))[1:400,1:400][,1:400], zlim = c(-1,1), col = heat.colors(100))
image(t(cov2cor(D$postcov))[1:400,1:400][,1:400],zlim = c(-1,1), col = heat.colors(100))
image(t(cov2cor(resPACE$fittedCov))[1:400,1:400][,1:400], zlim = c(-1,1),col = heat.colors(100))
image(t(cov2cor(EmpCov))[1:400,1:400][,1:400], zlim = c(-1,1),col = heat.colors(100))
image(t(cov2cor(resMarginalCov))[1:400,1:400][,1:400], zlim = c(-1,1),col = heat.colors(100))
image(t(resProductCov)[1:400,1:400][400:1,], zlim = c(min(unlist(mylist)), max(unlist(mylist))),col = heat.colors(100))

image(t(S[1:370,1:370])[370:1,], col = heat.colors(60))
sum((Cov.Weak - aa)^2)/400
sum((Cov.Weak - bb)^2)/400
sum((Cov.Weak - EmpCov)^2)/400
sum((Cov.Weak - D$postcov)^2)/400
sum((Cov.Weak - FPCADense$fittedCov)^2)/400
