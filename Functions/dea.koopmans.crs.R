require(Rglpk)

dea.koopmans.crs <- function(XOBS,YOBS,XREF=NULL,YREF=NULL,weights=weights) {
  
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
    f.obj=c(weights,rep(0,n.ref))
    #
    xblock=rbind(diag(-1,p,p),XREF)
    yblock=rbind(matrix(0,nrow=p,ncol=q),YREF)
    f.con=t(cbind(xblock,yblock))
    # Set unequality signs
    f.dir=c(rep("<=",p),rep(">=",q))
    #
    f.rhs=c(rep(0,p),yi)
    #
    results=Rglpk_solve_LP(f.obj, f.con, f.dir, f.rhs, max=FALSE)
    #
    res[i]=sum(weights*results$optimum)/sum(weights*xi)
  }
  return(res)
}