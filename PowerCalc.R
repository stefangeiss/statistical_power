
### EXECUTE THIS FUNCTION CODE TO ESTABLISH THE FUNCTION IN YOUR R

PowerCalc <- function(E,Q,N,alpha){

	# Loads the equations estimated from the simulation data
	load("eq_ci_025_ED.RData")
	load("eq_ci_025_MD.RData")
	load("eq_median_ED.RData")
	load("eq_median_MD.RData")

	# Uses the equations to calculate the lower limit of the confidence interval
	ci.m <- predict(fitmodel.m.025x,newdata=data.frame(ceffsize=E,cquality=Q,csampsize=N))
	ci.e <- predict(fitmodel.s.025x,newdata=data.frame(ceffsize=E,cquality=Q,csampsize=N))
	mdn.m <- predict(fitmodel.m.500,newdata=data.frame(ceffsize=E,cquality=Q,csampsize=N))
	mdn.e <- predict(fitmodel.s.500,newdata=data.frame(ceffsize=E,cquality=Q,csampsize=N))
	LL.m <- mdn.m+ci.m
	LL.e <- mdn.e+ci.e

	# Computes the critical values the lower limit of the confidence interval needs to surpass
	critT <- qt(1-alpha,N-2)
	critR <- critT/((N-2+(E^2))^0.5)

	# Compares obtained lower limits with critical values
	powerOK.m <- LL.m>critR
	powerSurplus.m <- LL.m-critR
	powerOK.e <- LL.e>critR
	powerSurplus.e <- LL.e-critR

	# Returns the results
	return(list(MDPrior=list(PowerOK=powerOK.m,PowerSurplus=powerSurplus.m,MedianCorr=mdn.m,CorrConf=ci.m,LowerLimit=LL.m),
			EDPrior=list(PowerOK=powerOK.e,PowerSurplus=powerSurplus.e,MedianCorr=mdn.e,CorrConf=ci.e,LowerLimit=LL.e),
			CriticalValues=list(criticalR=critR,criticalT=critT)))
}

### USE THE FUNCTION TO CHECK WHETHER YOU HAVE SUFFICIENT POWER

PowerCalc(E=0.25,Q=0.667,N=500,alpha=.05) # Desired alpha level: p<.05; Coding accuracy=0.667; sample size=500; effect size=R=0.25

### ANNOTED OUTPUT

# Your output looks like this:
#
# $MDPrior					### THESE ARE THE RESULTS FOR MD-PRIORS
# $MDPrior$PowerOK				### IS THE POWER GREAT ENOUGH TO RELIABLY DISCOVER THE EFFECT YOU SPECIFIED?
# [1] TRUE 							# TRUE: Power is OK; FALSE: Power is too low. 
# 
# $MDPrior$PowerSurplus			### HOW MUCH GREATER IS THE EXPECTED LOWER BOUNDARY OF THE 95%CONFIDENCE INTERVAL COMPARED TO THE CRITICAL VALUE FOR R?
# [1] 0.05015891					# POSITIVE. Power is OK. NEGATIVE. Power is too low. 
# 
# $MDPrior$MedianCorr			### EXPECTED MEDIAN CORRELATION OBTAINED
# [1] 0.2055127
#  
# $MDPrior$CorrConf				### LOWER 95% CONFIDENCE RANGE FOR THE MEDIAN
# [1] -0.08011315
# 
# $MDPrior$LowerLimit			### LOWER BOUNDARY OF THE 95% CONFIDENCE REGION AROUND THE MEDIAN; SHOULD BE GREATER THAN THE CRITICAL R --> PowerOK, PowerSurplus
# [1] 0.1253995
#
#
# $EDPrior					### THESE ARE THE RESULTS FOR MD-PRIORS
# $EDPrior$PowerOK				### IS THE POWER GREAT ENOUGH TO RELIABLY DISCOVER THE EFFECT YOU SPECIFIED?
# [1] TRUE 							# TRUE: Power is OK; FALSE: Power is too low. 
# 
# $EDPrior$PowerSurplus			### HOW MUCH GREATER IS THE EXPECTED LOWER BOUNDARY OF THE 95%CONFIDENCE INTERVAL COMPARED TO THE CRITICAL VALUE FOR R?
# [1] 0.01767901					# POSITIVE. Power is OK. NEGATIVE. Power is too low. 
#
# $EDPrior$MedianCorr			### EXPECTED MEDIAN CORRELATION OBTAINED
# [1] 0.1730328
#
# $EDPrior$CorrConf				### LOWER 95% CONFIDENCE RANGE FOR THE MEDIAN
# [1] -0.08151337
# 
# $EDPrior$LowerLimit			### LOWER BOUNDARY OF THE 95% CONFIDENCE REGION AROUND THE MEDIAN; SHOULD BE GREATER THAN THE CRITICAL R --> PowerOK, PowerSurplus
# [1] 0.0915194
#
#
# $CriticalValues
# $CriticalValues$criticalR
# [1] 0.07384039 				### THIS IS THE R VALUE THAT MUST BE SURPASSED IN ANY SINGLE STUDY TO OBTAIN A STATISTICALLY SIGNIFICANT RESULT AT "ALPHA" --> PowerOK, PowerSurplus
# 
# $CriticalValues$criticalT
# [1] 1.647919					### THE CRITICAL T VALUE. NOT THAT IMPORTANT.

### SOME MORE APPLICATIONS OF THE FUNCTION

PowerCalc(E=0.20,Q=0.667,N=500,alpha=.05)
PowerCalc(E=0.15,Q=0.667,N=500,alpha=.05)

PowerCalc(E=0.15,Q=0.667,N=500,alpha=.05)
PowerCalc(E=0.15,Q=0.667,N=1000,alpha=.05)
PowerCalc(E=0.15,Q=0.667,N=2000,alpha=.05)

PowerCalc(E=0.15,Q=0.667,N=500,alpha=.05)
PowerCalc(E=0.15,Q=0.750,N=500,alpha=.05)
PowerCalc(E=0.15,Q=0.850,N=500,alpha=.05)
PowerCalc(E=0.15,Q=0.950,N=500,alpha=.05)
PowerCalc(E=0.15,Q=1.000,N=500,alpha=.05)
