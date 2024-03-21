###########################################################################################
## Title: Russell and Slack-Based Measures of Efficiency: A Unifying Framework
## Authors: Valentin Zelenyuk and Shirong Zhao
## Date: March 21, 2024
## It is free for academic use only with adequate citation and acknowledgments.
## For any other use, contact the authors. All rights reserved. 
##############
## The following codes will report the results for the numerical illustration
## Here we use CRS-DEA method
###########################################################################################

require(Rglpk)
require(FEAR)

if (exists(".Random.seed")) {
  save.seed=.Random.seed
  flag.seed=TRUE
} else {
  flag.seed=FALSE
}
set.seed(900001)
#####################################################
######################################################



x1=c(1,1,2,1,1,2)
x2=c(1,2,1,1,1,3)
y1=c(2,2,2,1,2,1)
y2=c(2,2,2,2,1,2)


x=cbind(x1,x2)
y=cbind(y1,y2)

#########################################
###### Various Russell Measure    #######
#########################################

### input-oriented, crs ###
source("./Functions/dea.russel.input.crs.R") 
russel.input.crs=dea.russel.input.crs(XOBS=x,YOBS=y,XREF=NULL,YREF=NULL)
round(russel.input.crs,3)


### output-oriented, crs ###
source("./Functions/dea.russel.output.crs.R") 
russel.output.crs=dea.russel.output.crs(XOBS=x,YOBS=y,XREF=NULL,YREF=NULL)
round(russel.output.crs,3)


### hyperbolic-oriented, crs ###

# the following results are obtained from
# russell_hyperbolic_crs.m in the folder of Matlab-codes
russel.hyperbolic.crs=c(1,0.875,0.875,0.875,0.875,0.559)
round(russel.hyperbolic.crs,3)



### additive version, crs ###
source("./Functions/dea.russel.additive.crs.R") 
russel.additive.crs=dea.russel.additive.crs(XOBS=x,YOBS=y,XREF=NULL,YREF=NULL)
round(russel.additive.crs,3)


### weighted input-oriented, crs ###
source("./Functions/dea.russel.input.weighted.crs.R") 

library("quadprog")
y12=sqrt(y1^2+y2^2)
Rinv <- solve(chol(t(log(x)) %*% log(x)));
C <- cbind(rep(1,2), diag(2))
b <- c(1,rep(0,2))
d <- t(log(y12)) %*% log(x)  
solve=solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq = 1)
weights=solve$solution
weights


russel.input.weighted.crs=dea.russel.input.weighted.crs(XOBS=x,YOBS=y,XREF=NULL,YREF=NULL,weights=weights)
round(russel.input.weighted.crs,3)


#########################################
#######   Koopmans Measures       #######
#########################################


# input-oriented, crs
# weights are the same as weighted input-oriented Russell Measure
source("./Functions/dea.koopmans.crs.R") 
koopmans.input.crs=dea.koopmans.crs(XOBS=x,YOBS=y,XREF=NULL,YREF=NULL,weights=weights)
round(koopmans.input.crs,3)


#########################################
############    Tone's SBM     ##########
#########################################

# input-oriented, crs
source("./Functions/dea.tone.crs.R") 
tone.crs=dea.tone.crs(XOBS=x,YOBS=y,XREF=NULL,YREF=NULL)
round(tone.crs,3)


##################################################################################
#######   Directional Slack-Based Measures (Fukuyama and Weber, 2009)      #######
##################################################################################

source("./Functions/dea.fukuyama.weber.crs.R") 
fukuyama.weber.crs=dea.fukuyama.weber.crs(XOBS=x,YOBS=y,XREF=NULL,YREF=NULL,XDIREC=x,YDIREC=y)
round(fukuyama.weber.crs,3)


#################################################################################################################
#######  Slack-Based Directional Distance Function (Fare and Grosskopf,2010; Fare et al., 2015)         #######
#################################################################################################################

source("./Functions/dea.fare.grosskopf.crs.R") 
fare.grosskopf.crs=dea.fare.grosskopf.crs(XOBS=x,YOBS=y,XREF=NULL,YREF=NULL)
round(fare.grosskopf.crs,3)


######################################################
#######  Directional Distance Function        #######
######################################################

# ddf, crs
ddf.crs=FEAR::dea.direc(XOBS=t(x), YOBS=t(y), XDIREC=t(x), YDIREC=t(y),
          RTS = 3, XREF = NULL, YREF = NULL,
          DISP = "strong", IS.EFF = NULL, errchk = TRUE)
round(ddf.crs,3)


###################################################################
##############    Debreu-Farrell-Shephard Measures    #############
###################################################################

# input-oriented, crs
input.crs=FEAR::dea(XOBS=t(x),YOBS=t(y),RTS=3,ORIENTATION=1,METRIC=2)
round(input.crs,3)

# output-oriented, crs
output.crs=FEAR::dea(XOBS=t(x),YOBS=t(y),RTS=3,ORIENTATION=2,METRIC=2)
round(output.crs,3)

# hyperbolic-oriented, crs
hyperbolic.crs=FEAR::dea(XOBS=t(x),YOBS=t(y),RTS=3,ORIENTATION=3,METRIC=2)
round(hyperbolic.crs,3)

#######################################################################
#######################################################################
############################  Summary     #############################
#######################################################################
eff.all=cbind(russel.input.crs, 
              russel.output.crs, 
              russel.hyperbolic.crs,
              russel.additive.crs,
              russel.input.weighted.crs,
              koopmans.input.crs,
              tone.crs,
              fukuyama.weber.crs,
              fare.grosskopf.crs, 
              ddf.crs,
              input.crs,
              output.crs,
              hyperbolic.crs)

colnames(eff.all)<-NULL

round(eff.all,3)


### construct the Table ###
id=1:nrow(eff.all)
tex=formatC(id,width=7,digits=0,format="f")

for (k in 1:ncol(eff.all)) {
  tex = paste(tex,"&",formatC(eff.all[,k],width=5,digits = 3,format = "f"))
}
tex = paste(tex,"\\\\")
write(tex,file="./Output/estimates.tex")





