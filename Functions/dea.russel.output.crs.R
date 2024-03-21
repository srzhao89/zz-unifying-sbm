require(Rglpk)

dea.russel.output.crs <- function(XOBS,YOBS,XREF=NULL,YREF=NULL) {
  
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
    f.obj=c(rep(1/q,q),rep(0,n.ref))
    #
    xblock=rbind(matrix(0,nrow=q,ncol=p),XREF)
    yblock=rbind(-1*diag(yi,q,q),YREF)
    f.con=t(cbind(xblock,yblock))
    # Set unequality signs
    f.dir=c(rep("<=",p),rep(">=",q))
    #
    f.rhs=c(xi,rep(0,q))
    #
    bounds <- list(lower=list(ind = c(1:q), val = rep(1,q)))
    #
    results=Rglpk_solve_LP(f.obj, f.con, f.dir, f.rhs, bounds, max=TRUE)
    #
    res[i]=results$optimum
  }
  return(res)
}