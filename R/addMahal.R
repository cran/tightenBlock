addMahal <-
function (costmatrix, z, X)
{
  stopifnot(is.vector(z))
  stopifnot(all((z == 0) | (z == 1)))
  stopifnot(is.matrix(costmatrix))
  stopifnot((dim(costmatrix)[1])==sum(z))
  stopifnot((dim(costmatrix)[2])==sum(1-z))
  if (is.vector(X))
    X <- matrix(X, length(X), 1)
  if (is.data.frame(X))
    X <- as.matrix(X)
  stopifnot(is.matrix(X))
  stopifnot(length(z) == (dim(X)[1]))

  # Check that z and costmatrix have compatible names, if any.
  if (!is.null(names(z))) nm<-names(z)
  else nm<-1:length(z)
  stopifnot(all(rownames(costmatrix)==nm[z==1]))
  stopifnot(all(colnames(costmatrix)==nm[z==0]))


  n <- dim(X)[1]
  rownames(X) <- 1:n
  k <- dim(X)[2]
  m <- sum(z)
  for (j in 1:k) X[, j] <- rank(X[, j])
  cv <- stats::cov(X)
  vuntied <- stats::var(1:n)
  rat <- sqrt(vuntied/diag(cv))
  cv <- diag(rat) %*% cv %*% diag(rat)
  if (dim(X)[2]==1){
    out<-outer(X[z==1,1],X[z==0,1],"-")
    out<-(out^2)/cv[1,1]
    costmatrix+out
  }
  else{
    out <- matrix(NA, m, n - m)
    Xc <- X[z == 0, ]
    Xt <- X[z == 1, ]
    icov <- MASS::ginv(cv)
    for (i in 1:m) out[i, ] <- stats::mahalanobis(Xc, Xt[i, ],
                                                icov, inverted = TRUE)

    rownames(out)<-nm[z==1]
    colnames(out)<-nm[z==0]

    costmatrix+out
  }
}
