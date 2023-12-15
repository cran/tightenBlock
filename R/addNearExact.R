addNearExact <-
function(costmatrix,z,exact,penalty=1000){
  stopifnot(is.vector(z))
  stopifnot(all((z == 0) | (z == 1)))
  stopifnot((dim(costmatrix)[1])==sum(z))
  stopifnot((dim(costmatrix)[2])==sum(1-z))
  stopifnot(is.vector(exact)&(length(z)==length(exact)))
  if (!is.null(names(z))) nm<-names(z)
  else nm<-1:length(z)
  stopifnot(all(rownames(costmatrix)==nm[z==1]))
  stopifnot(all(colnames(costmatrix)==nm[z==0]))
  o<-outer(exact[z==1],exact[z==0],"!=")*penalty
  rownames(o)<-nm[z==1]
  colnames(o)<-nm[z==0]
  costmatrix+o
}
