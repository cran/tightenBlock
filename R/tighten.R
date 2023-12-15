tighten<-function(dat, z, block, x=NULL, f=NULL, ncontrols=1,
                  subset=NULL, pspace=10, solver="rrelaxiv"){

  # Check input
  stopifnot((solver=="rrelaxiv")|(solver=="rlemon"))
  if (is.null(x)&is.null(f)) stop("Either x or f or both must not be NULL.")
  stopifnot(is.vector(pspace)&(length(pspace)==1)&(pspace>0))
  stopifnot(is.vector(ncontrols)&(length(ncontrols)==1)&(ncontrols>=1))
  stopifnot(is.matrix(dat)|is.data.frame(dat))
  stopifnot(is.vector(z)&(length(z)==(dim(dat)[1])))
  stopifnot(all((z==0)|(z==1)))
  stopifnot((min(z)==0)&(max(z)==1))
  stopifnot(is.vector(block)&(length(block)==length(z)))
  if (!is.null(x)) stopifnot(is.vector(x)|is.matrix(x)|is.data.frame(x))
  if ((!is.null(x))&(is.vector(x))) stopifnot(length(x)==length(z))
  if ((!is.null(x))&(!is.vector(x))) stopifnot(length(z)==(dim(x)[1]))

  if (is.null(names(z))) names(z)<-1:length(z)
  rownames(dat)<-names(z)


  # Check that this is a balanced block design
  blocksize<-table(block)
  if (min(blocksize)<2) stop("Each block must contain at least two people.")
  btable<-table(z,block)
  if (!all(btable[2,]==1)) stop("Every block must contain exactly one
                                treated individual")

  left<-startcost(z)
  right<-startcost(z)
  if (!is.null(x)) left<-addMahal(left,z,x)
  penalty<-(10+ceiling(max(as.vector(left))))*pspace

  if (!is.null(f)){
    stopifnot(is.vector(f)|is.factor(f)|is.matrix(f)|is.data.frame(f))
    if (is.vector(f)|is.factor(f)) stopifnot(length(f)==length(z))
    else stopifnot((dim(f)[1])==length(z))
    if (is.factor(f)) f<-as.integer(f)
    if (is.vector(f)) f<-matrix(f,length(f),1)
    for (j in 1:(dim(f)[2])) right<-addNearExact(right,z,f[,j],penalty=penalty)
    penalty<-(penalty+ceiling(max(right)))*pspace
  }


  left<-addNearExact(left,z,block,penalty=penalty)

  if (is.null(subset)) treatedcosts<-NULL
  else {
    stopifnot(is.vector(subset)&(length(subset)==1)&(subset>0))
    treatedcosts<-rep(subset,sum(z))
  }

  makematch(dat,left,right,ncontrols=ncontrols,treatedcosts=treatedcosts,
            solver=solver)
}
