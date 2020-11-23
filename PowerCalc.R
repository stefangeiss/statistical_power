######################################################
### EXECUTE THIS FUNCTION CODE TO ESTABLISH THE Powercalc FUNCTION IN YOUR R
######################################################

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

### Check out the PowerCalc function

PowerCalc(E=0.25,Q=0.667,N=500,alpha=.05) # Desired alpha level: p<.05; Coding accuracy=0.667; sample size=500; effect size=R=0.25

######################################################
### EXECUTE THESE 3 FUNCTION CODES TO ESTABLISH THE FindMinE, FindMinQ and FindMinN FUNCTIONS IN YOUR R
######################################################

FindMinN <- function(Q,E,alpha,STEPS){
	N_MD <- cbind(STEPS,rep(NA,times=length(STEPS)))
	N_ED <- cbind(STEPS,rep(NA,times=length(STEPS)))
	for (i in 1:length(STEPS))
	{
	N_MD[i,2] <- PowerCalc(E=E,Q=Q,N=STEPS[i],alpha=alpha)$MDPrior$PowerOK
	N_ED[i,2] <- PowerCalc(E=E,Q=Q,N=STEPS[i],alpha=alpha)$EDPrior$PowerOK
	}
	low_ED <- which(diff(N_ED[,2])==1)
	low_MD <- which(diff(N_MD[,2])==1)
	min_ED <- N_ED[low_ED+1,]
	min_MD <- N_MD[low_MD+1,]
	return(list(N_MD,N_ED,MinimumSampleSizeED=min_ED[1],MinimumSampleSizeMD=min_MD[1]))
	}

# This is a vector that lists all sample sizes you want to consider. Change manually according to your needs. 
# The default here is: 1-unit steps between 1 and 20; 5 unit steps between 20 and 50; 10 unit steps between 50 and 1000; 100 unit steps between 1000 and 4000. 1000 unit steps between 4000 and 10000. 10000 unit steps between 10000 and 100000. Sample sizes larger than 100 000 not considered.
nSTEPS <- c(seq(1,20,1),seq(25,50,5),seq(60,1000,10),seq(1100,4000,100),seq(5000,10000,1000),seq(20000,100000,10000))

FindMinN(STEPS=nSTEPS,Q=.667,E=.250,alpha=.05)

FindMinE <- function(Q,N,alpha){
	E_MD <- cbind(seq(0,1,.01),rep(NA,times=100))
	E_ED <- cbind(seq(0,1,.01),rep(NA,times=100))
	for (i in 0:100)
	{
	E_MD[i+1,2] <- PowerCalc(E=(i/100),Q=Q,N=N,alpha=alpha)$MDPrior$PowerOK
	E_ED[i+1,2] <- PowerCalc(E=(i/100),Q=Q,N=N,alpha=alpha)$EDPrior$PowerOK
	}
	low_ED <- which(diff(E_ED[,2])==1)
	low_MD <- which(diff(E_MD[,2])==1)
	min_ED <- E_ED[low_ED+1,]
	min_MD <- E_MD[low_MD+1,]
	return(list(E_MD,E_ED,MinimumEffectSizeED=min_ED[1],MinimumEffectSizeMD=min_MD[1]))
	}

### Check out the FindMinE function: It gives a minumum necessary effect size if Q, n, and alpha (type I error probability) are specified. 
FindMinE(Q=.667,N=500,alpha=.05)

FindMinQ <- function(E,N,alpha){
	Q_MD <- cbind(seq(0,1,.01),rep(NA,times=100))
	Q_ED <- cbind(seq(0,1,.01),rep(NA,times=100))
	for (i in 0:100)
	{
	Q_MD[i+1,2] <- PowerCalc(Q=(i/100),E=E,N=N,alpha=alpha)$MDPrior$PowerOK
	Q_ED[i+1,2] <- PowerCalc(Q=(i/100),E=E,N=N,alpha=alpha)$EDPrior$PowerOK
	}
	low_ED <- which(diff(Q_ED[,2])==1)
	low_MD <- which(diff(Q_MD[,2])==1)
	min_ED <- Q_ED[low_ED+1,]
	min_MD <- Q_MD[low_MD+1,]
	return(list(Q_MD,Q_ED,MinimumCodingAccuracyED=min_ED[1],MinimumCodingAccuracyMD=min_MD[1]))
	}
	
### Check out the FindMinQ function: It gives a minumum necessary coding accuracy if E(effect size), n, and alpha (type I error probability) are specified. 
FindMinQ(E=.250,N=500,alpha=.05)

FindMinAlpha <- function(Q,E,N,STEPS){
	Alpha_MD <- cbind(STEPS,rep(NA,times=length(STEPS)))
	Alpha_ED <- cbind(STEPS,rep(NA,times=length(STEPS)))
	for (i in 1:length(STEPS))
	{
	Alpha_MD[i,2] <- PowerCalc(E=E,Q=Q,N=N,alpha=STEPS[i])$MDPrior$PowerOK
	Alpha_ED[i,2] <- PowerCalc(E=E,Q=Q,N=N,alpha=STEPS[i])$EDPrior$PowerOK
	}
	low_ED <- which(diff(Alpha_ED[,2])==1)
	low_MD <- which(diff(Alpha_MD[,2])==1)
	min_ED <- Alpha_ED[low_ED+1,]
	min_MD <- Alpha_MD[low_MD+1,]
	return(list(Alpha_MD,Alpha_ED,MinimumAlphaED=min_ED[1],MinimumAlphaMD=min_MD[1]))
	}

# Manually define the alpha (type I error probability) levels to be considered. 
# Default here: 0.0001 steps between .0001 and .0010, .001 steps between .001 and .010; .01 steps between .01 and .10. Alphas greater than .10 are not considered.
alphaSTEPS <- c(seq(0.0001,0.0009,0.0001),seq(0.001,0.009,0.001),seq(0.01,0.09,0.01),seq(0.1,1.0,0.1))

### Check out the FindMinAlpha function: It gives a minumum possible alpha (type I error probability) if E(effect size), n, and Q are specified. 
FindMinAlpha(STEPS=alphaSTEPS,Q=.667,E=.250,N=500)



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

#########################################
### See the online appendix to the following article for an application and demo code: 
### Geiß, Stefan (forthcoming): Statistical Power in Content Analysis Designs: How Effect Size, Sample Size and Coding Accuracy Jointly Affect Hypothesis Testing – A Monte Carlo Simulation Approach. Computational Communication Research.

