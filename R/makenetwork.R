makenetwork <-
function(costL,costR,ncontrols=1,controlcosts=NULL,treatedcosts=NULL,large=100){
  stopifnot(is.matrix(costR))
  stopifnot(is.matrix(costL))
  stopifnot((dim(costR)[1])==(dim(costL)[1]))
  stopifnot((dim(costR)[2])==(dim(costL)[2]))
  stopifnot(is.vector(ncontrols)&(length(ncontrols==1))&(ncontrols>=1))
  stopifnot(ncontrols==round(ncontrols))
  if ((!is.null(treatedcosts))&(ncontrols>1)) {
    warning("Subset matching is not available with ncontrols>1.")
    stopifnot((is.null(treatedcosts))|(ncontrols==1))
  }
  if (!is.null(controlcosts)) {
    stopifnot((dim(costR)[2])==length(controlcosts))
    stopifnot(all(controlcosts>=0))
  }
  else controlcosts<-rep(0,dim(costR)[2])
  stopifnot(all(rownames(costR)==rownames(costL)))
  stopifnot(all(colnames(costR)==colnames(costL)))

  if (!is.null(treatedcosts)) {
    stopifnot((dim(costR)[1])==length(treatedcosts))
    stopifnot(all(treatedcosts>=0))
  }
  else {
    stopifnot(is.vector(large)&(length(large)==1)&(large>1))
    tmx<-max(c(as.vector(costL),as.vector(costR),controlcosts))
    treatedcosts<-rep(tmx*large,dim(costR)[1])
  }

  stopifnot(all(rownames(costR)==rownames(costL)))
  stopifnot(all(colnames(costR)==colnames(costL)))

  ntreated<-dim(costR)[1]
  ncontrol<-dim(costR)[2]

  treated<-1:ntreated
  control<-(ntreated+1):(ntreated+ncontrol)

  idtreated<-rownames(costR)
  idcontrol<-colnames(costR)

  startn<-rep(treated,ncontrol)
  endn<-as.vector(t(matrix(rep(control,ntreated),ncontrol,ntreated)))
  cost<-as.vector(costL)

  cc<-(ntreated+ncontrol+1):(ntreated+2*ncontrol) # duplicate controls

  startn<-c(startn,control) # control-to-duplicate-control edges
  endn<-c(endn,cc)
  cost<-c(cost,controlcosts)

  startn<-c(startn,as.vector(t(matrix(rep(cc,ntreated),ncontrol,ntreated))))
  tt<-(1+ntreated+2*ncontrol):(2*(ntreated+ncontrol)) # duplicate treated
  endn<-c(endn,rep(tt,ncontrol))
  cost<-c(cost,as.vector(costR))

  startn<-c(startn,treated) # treated-to-duplicate-treated subset edges
  endn<-c(endn,tt)
  cost<-c(cost,treatedcosts)


  ucap<-rep(1,length(startn))
  b<-c(rep(ncontrols,ntreated),rep(0,2*ncontrol),rep(-ncontrols,ntreated))

  net=list(startn=startn,endn=endn,ucap=ucap,cost=cost,b=b)
  list(idtreated=idtreated,idcontrol=idcontrol,net=net)
}
