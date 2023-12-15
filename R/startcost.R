startcost <-
function(z){
  # returns a treated x control cost matrix with zero costs
  stopifnot(is.vector(z))
  stopifnot(all((z==0)|(z==1)))
  if (!is.null(names(z))) nm<-names(z)
  else nm<-1:length(z)
  o<-matrix(0,sum(z),sum(1-z))
  rownames(o)<-nm[z==1]
  colnames(o)<-nm[z==0]
  o
}
