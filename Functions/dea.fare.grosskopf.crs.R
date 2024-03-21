require(Rglpk)

dea.fare.grosskopf.crs <- function(XOBS,YOBS,XREF=NULL,YREF=NULL) {
  
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
    f.obj=c(rep(1,p+q),rep(0,n.ref))
    #
    xblock=rbind(diag(1,p,p),matrix(0,nrow=q,ncol=p),XREF)
    yblock=rbind(matrix(0,nrow=p,ncol=q),diag(-1,q,q),YREF)
    f.con=t(cbind(xblock,yblock))
    # Set unequality signs
    f.dir=c(rep("<=",p),rep(">=",q))
    #
    f.rhs=c(xi,yi)
    #
    results=Rglpk_solve_LP(f.obj, f.con, f.dir, f.rhs, max=TRUE)
    #
    res[i]=results$optimum
  }
  return(res)
}