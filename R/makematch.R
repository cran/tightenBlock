makematch <-
function(dat,costL,costR,ncontrols=1,controlcosts=NULL,
         treatedcosts=NULL,large=100,solver="rrelaxiv"){

  stopifnot((solver=="rrelaxiv")|(solver=="rlemon"))
  stopifnot(is.matrix(costR))
  stopifnot(is.matrix(costL))
  stopifnot((dim(costR)[1])==(dim(costL)[1]))
  stopifnot((dim(costR)[2])==(dim(costL)[2]))
  stopifnot(is.vector(ncontrols)&(length(ncontrols==1))&(ncontrols>=1))
  stopifnot(ncontrols==round(ncontrols))
  if (!is.null(controlcosts)) {
    stopifnot((dim(costR)[2])==length(controlcosts))
    stopifnot(all(controlcosts>=0))
  }
  else controlcosts<-rep(0,dim(costR)[2])
  stopifnot(all(rownames(costR)==rownames(costL)))
  stopifnot(all(colnames(costR)==colnames(costL)))
  stopifnot(is.data.frame(dat)|is.matrix(dat))
  stopifnot((dim(dat)[1])==sum(dim(costR)))
  if (!is.null(treatedcosts)) {
    stopifnot((dim(costR)[1])==length(treatedcosts))
    stopifnot(all(treatedcosts>=0))
    stopifnot(is.vector(large)&(length(large)==1)&(large>1))
  }

  net<-makenetwork(costL,costR,ncontrols=ncontrols,controlcosts=controlcosts,
                   treatedcosts=treatedcosts,large=large)$net

  result<-rcbalance::callrelax(net,solver=solver)
  if (result$crash==1) {
    warning("callrelax crashed")
    stop()
  }
  if (result$feasible==0) {
    warning("problem is infeasible")
    stop()
  }
  pairs<-result$x[1:(dim(costR)[1]*dim(costR)[2])]
  pairs<-matrix(pairs,dim(costR)[1],dim(costR)[2])
  rownames(pairs)<-rownames(costR)
  colnames(pairs)<-colnames(costR)

  o<-NULL
  notMatch<-NULL
  cmatch<-0 # count of matched sets
  for (i in 1:(dim(costR)[1])){
    if (sum(pairs[i,])==0) notMatch<-c(notMatch,rownames(pairs)[i])
    else {
      o<-rbind(o,dat[which(rownames(dat)==(rownames(pairs)[i])),])
      w<-colnames(pairs)[pairs[i,]==1]
      if (length(w)!=ncontrols) warning(paste(rownames(pairs)[i]," has fewer
                                              than ",ncontrols," controls."))
      stopifnot(length(w)==ncontrols)
      o<-rbind(o,dat[is.element(rownames(dat),w),])
      cmatch<-cmatch+1
    }
  }

  if (!is.null(notMatch))
    warning(paste("Number of treated individuals not matched: ",length(notMatch)))
  mset<-gl(cmatch,ncontrols+1)
  o<-cbind(o,mset)
  o
}
