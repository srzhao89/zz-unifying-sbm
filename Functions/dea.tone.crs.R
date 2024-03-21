require(Rglpk)

dea.tone.crs <- function(XOBS,YOBS,XREF=NULL,YREF=NULL) {
  
  if (is.null(XREF) | is.null(YREF)) {
    XREF=XOBS
    YREF=YOBS
  }
  
  n=nrow(XOBS)
  p=ncol(XOBS)
  q=ncol(YOBS)
  n.ref=nrow(XREF)
  res=rep(NA,n)
  
  for (i in 1:n) {
    xi=as.vector(XOBS[i,])
    yi=as.vector(YOBS[i,])
    #
    f.obj=c(rep(0,n.ref),(-1)/(p*xi),rep(0,q),1)
    #
    cons1=c(rep(0,n.ref),rep(0,p),1/(q*yi),1)
    xblock=rbind(XREF,diag(1,p,p),matrix(0,nrow=q,ncol=p),-xi)
    yblock=rbind(YREF,matrix(0,nrow=p,ncol=q),-1*diag(1,q,q),-yi)
    f.con=rbind(cons1,t(xblock),t(yblock))
    # Set unequality signs
    f.dir=c("==",rep("==",p),rep("==",q))
    #
    f.rhs=c(1,rep(0,p),rep(0,q))
    #
    results=Rglpk_solve_LP(f.obj, f.con, f.dir, f.rhs, max=FALSE)
    #
    res[i]=results$optimum
  }
  return(res)
}