## Geodesic distance on spheres.
# y1,y2 Two unit vectors, i.e., with \eqn{L^2} norm equal to 1, of the same length.
# returns A scalar holding the geodesic distance between \code{y1} and \code{y2}.
## Example
#d <- 3
#y1 <- rnorm(d)
#y1 <- y1 / sqrt(sum(y1^2))
#y2 <- rnorm(d)
#y2 <- y2 / sqrt(sum(y2^2))
#dist <- SpheGeoDist(y1,y2)


SpheGeoDist <- function(y1,y2) {
  tol <- 1e-6 # Distance that the L2 norms of y1 and y2 are allowed to be away from 1.
  if (abs(length(y1) - length(y2)) > 0) {
    stop("y1 and y2 should be of the same length.")
  }
  if (abs(l2norm(y1)-1) > tol) {
    stop("y1 is not a unit vector.")
  }
  if (abs(l2norm(y2)-1) > tol) {
    stop("y2 is not a unit vector.")
  }
  y1 = y1 / l2norm(y1)
  y2 = y2 / l2norm(y2)
  if (sum(y1 * y2) > 1){
    return(0)
  } else if (sum(y1*y2) < -1){
    return(pi)
  } else return(acos(sum(y1 * y2)))
}
