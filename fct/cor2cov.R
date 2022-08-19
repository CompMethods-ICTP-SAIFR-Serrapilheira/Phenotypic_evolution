#Function to make covariance matrix
cor2cov <- function(cor=cor.mat, sd=c(2,2)){
  require(mvtnorm)
  sigma <- sd ### variance to sample the effects in alleles
  cov.1 <- sigma %*% t(sigma) ### transform the matrix of correlation in a covariance matrix
  cov <- cor * cov.1 ## estimating the matrices with a given standard deviation
  return(cov)
}
