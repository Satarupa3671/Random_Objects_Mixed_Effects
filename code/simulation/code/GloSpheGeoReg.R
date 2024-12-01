### Global Fréchet Regression for Spherical Data with respect to the geodesic distance.
### Input: xin An \eqn{n}-by-\eqn{p} matrix with input measurement points.
###        yin An \eqn{n}-by-\eqn{m} matrix holding the spherical data, of which the sum of squares of elements within each row is 1.
###        xout A vector of length \eqn{p} with output measurement points;
### Output: A vector of length \eqn{m} holding the fitted responsee, which is a spherical vector, corresponding to each element in \code{xout}.
### Example:
###   n <- 101
###   xin <- seq(-1,1,length.out = n)
###   xin = as.matrix(xin, ncol = 1)
###   theta_true <- rep(pi/2,n)
###   phi_true <- (xin + 1) * pi / 4
###   ytrue <- apply( cbind( 1, phi_true, theta_true ), 1, pol2car )
###   yin <- t( ytrue )
###   xout <- mean(xin)
###   res <- GloSpheGeoReg(xin=xin, yin=yin, xout=xout)
###   references: \cite{Petersen, A., & Müller, H.-G. (2019). "Fréchet regression for random objects with Euclidean predictors." The Annals of Statistics, 47(2), 691--719.}

### Important: using trust package and perturbation for initial value


GloSpheGeoReg <- function(xin, yin, xout) {
  k = length(xout)
  n = nrow(xin)
  m = ncol(yin)
  
  xbar <- colMeans(xin)
  Sigma <- cov(xin) * (n-1) / n
  invSigma <- solve(Sigma)
  
  s <- 1 + t(t(xin) - xbar) %*% invSigma %*% (xout - xbar)
  s <- as.vector(s)
  
  # initial guess
  y0 = colMeans(yin*s)
  y0 = y0 / l2norm(y0)
  if (sum(sapply(1:n, function(i) sum(yin[i,]*y0)) > 1-1e-8)){
    #if (sum( is.infinite (sapply(1:n, function(i) (1 - sum(yin[i,]*y0)^2)^(-0.5) )[ker((xout[j] - xin) / bw)>0] ) ) +
    #   sum(sapply(1:n, function(i) 1 - sum(yin[i,] * y0)^2 < 0)) > 0){
    # return(y0)
    y0 = y0 + rnorm(3) * 1e-3
    y0 = y0 / l2norm(y0)
  }
  objFctn = function(y){
    y <- y / l2norm(y)
    f = mean(s * sapply(1:n, function(i) SpheGeoDist(yin[i,], y)^2))
    if (abs(l2norm(y)-1) > 1e-15) {
      return(list(value = Inf))
    }
    g = 2 * colMeans(t(sapply(1:n, function(i) SpheGeoDist(yin[i,], y) * SpheGeoGrad(yin[i,], y))) * s)
    res = sapply(1:n, function(i){
      grad_i = SpheGeoGrad(yin[i,], y)
      return((grad_i %*% t(grad_i) + SpheGeoDist(yin[i,], y) * SpheGeoHess(yin[i,], y)) * s[i])
    }, simplify = "array")
    h = 2 * apply(res, 1:2, mean)
    return(list(value=f, gradient=g, hessian=h))
  }
  res = trust::trust(objFctn, y0, 0.1, 1e5)
  # res = trust::trust(objFctn, y0, 0.1, 1)
  return(res$argument / l2norm(res$argument))
}
