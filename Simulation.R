# setwd()

##################################
##################################
### GETTING THE FUNCTIONS IN

install.packages("combinat")

library(combinat)
library(mvtnorm)
library(car)
library(reshape2)
library(stringr)

tval <- function(r,n){
	df <- n-2
	ser <- sqrt((1-r^2)/df)
	tva <- r/ser
	return(tva)
}

rval <- function(t,n){
	df <- n-2
	ser <- sqrt((1-r^2)/df)
	tva <- r/ser
	return(tva)
}




#######################################
#######################################
##### FUNCTIONS FOR CALCULATION INTER-CODER RELIABILITY COEFFICIENTS (BY KILEM GWET)
##### --> http://inter-rater-reliability.blogspot.com/2014/03/some-r-functions-for-calculating-chance.html
##### Now available under library(irrCAC)

identity.weights<-function(categ){
	weights<-diag(length(categ))
	return (weights)
}

quadratic.weights<-function(categ){
	q<-length(categ)
	weights <- diag(q)
	if (is.numeric(categ)) { 
	   categ.vec <- sort(categ)
	}
	else {
	   categ.vec<-1:length(categ)
	}
	xmin<-min(categ.vec)
	xmax<-max(categ.vec)
	for(k in 1:q){
	    for(l in 1:q){
		  weights[k,l] <- 1-(categ.vec[k]-categ.vec[l])^2/(xmax-xmin)^2 
	    }
      }
	return (weights)
}

linear.weights<-function(categ){
	q<-length(categ)
	weights <- diag(q)
	if (is.numeric(categ)) { 
	   categ.vec <- sort(categ)
	}
	else {
	   categ.vec<-1:length(categ)
	}
	xmin<-min(categ.vec)
	xmax<-max(categ.vec)
	for(k in 1:q){
	    for(l in 1:q){
		  weights[k,l] <- 1-abs(categ.vec[k]-categ.vec[l])/abs(xmax-xmin)
	    }
      }
	return (weights)
}
#--------------------------------
radical.weights<-function(categ){
	q<-length(categ)
	weights <- diag(q)
	if (is.numeric(categ)) { 
	   categ.vec <- sort(categ)
	}
	else {
	   categ.vec<-1:length(categ)
	}
	xmin<-min(categ.vec)
	xmax<-max(categ.vec)
	for(k in 1:q){
	    for(l in 1:q){
		  weights[k,l] <- 1-sqrt(abs(categ.vec[k]-categ.vec[l]))/sqrt(abs(xmax-xmin))
	    }
      }
	return (weights)
}

#--------------------------------
ratio.weights<-function(categ){
	q<-length(categ)
	weights <- diag(q)
	if (is.numeric(categ)) { 
	   categ.vec <- sort(categ)
	}
	else {
	   categ.vec<-1:length(categ)
	}
	xmin<-min(categ.vec)
	xmax<-max(categ.vec)
	for(k in 1:q){
	    for(l in 1:q){
		  weights[k,l] <- 1-((categ.vec[k]-categ.vec[l])/(categ.vec[k]+categ.vec[l]))^2 / ((xmax-xmin)/(xmax+xmin))^2
	    }
      }
	return (weights)
}

#--------------------------------
circular.weights<-function(categ){
	q<-length(categ)
	weights <- diag(q)
	if (is.numeric(categ)) { 
	   categ.vec <- sort(categ)
	}
	else {
	   categ.vec<-1:length(categ)
	}
	xmin<-min(categ.vec)
	xmax<-max(categ.vec)
	U = xmax-xmin+1
	for(k in 1:q){
	    for(l in 1:q){
		  weights[k,l] <- (sin(pi*(categ.vec[k]-categ.vec[l])/U))^2
	    }
      }
	weights <- 1-weights/max(weights)
	return (weights)
}

#--------------------------------
bipolar.weights<-function(categ){
	q<-length(categ)
	weights <- diag(q)
	if (is.numeric(categ)) { 
	   categ.vec <- sort(categ)
	}
	else {
	   categ.vec<-1:length(categ)
	}
	xmin<-min(categ.vec)
	xmax<-max(categ.vec)
	for(k in 1:q){
	    for(l in 1:q){
		  if (k!=l)
		  	weights[k,l] <- (categ.vec[k]-categ.vec[l])^2 / (((categ.vec[k]+categ.vec[l])-2*xmin)*(2*xmax-(categ.vec[k]+categ.vec[l])))
		  else weights[k,l] <- 0
	    }
      }
	weights <- 1-weights/max(weights)
	return (weights)
}


#--------------------------------
ordinal.weights<-function(categ){
	q<-length(categ)
	weights <- diag(q)
      categ.vec<-1:length(categ)
	for(k in 1:q){
	    for(l in 1:q){
		  nkl <- max(k,l)-min(k,l)+1
		  weights[k,l] <- nkl * (nkl-1)/2
	    }
      }
	weights <- 1-weights/max(weights)
	return (weights)
}


#								AGREE.COEFF3.RAW.R
#						 		 (March 26, 2014)
#Description: This script file contains a series of R functions for computing various agreement coefficients
#		  for multiple raters (2 or more) when the input data file is in the form of nxr matrix or data frame showing 
#             the actual ratings each rater (column) assigned to each subject (in row). That is n = number of subjects, and r = number of raters.
#             A typical table entry (i,g) represents the rating associated with subject i and rater g. 
#Author: Kilem L. Gwet, Ph.D.
#
#


# ==============================================================
# This is an r function for trimming leading and trealing blanks
# ==============================================================
trim <- function( x ) {
  gsub("(^[[:space:]]+|[[:space:]]+$)", "", x) }

#===========================================================================================
#gwet.ac1.raw: Gwet's AC1/Ac2 coefficient (Gwet(2008)) and its standard error for multiple raters when input 
#		   dataset is a nxr matrix of alphanumeric ratings from n subjects and r raters 
#-------------
#The input data "ratings" is a nxr data frame of raw alphanumeric ratings
#from n subjects and r raters. Exclude all subjects that are not rated by any rater.
#Bibliography:
#Gwet, K. L. (2008). ``Computing inter-rater reliability and its variance in the presence of high
#		agreement." British Journal of Mathematical and Statistical Psychology, 61, 29-48.
#======================================================================================
gwet.ac1.raw <- function(ratings,weights="unweighted",conflev=0.95,N=Inf,print=TRUE){ 
  ratings.mat <- as.matrix(ratings) 
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N # final population correction 

  # creating a vector containing all categories used by the raters
 
  categ.init <- unique(as.vector(ratings.mat))
  if (is.numeric(categ.init))
     categ <- sort(as.vector(na.omit(categ.init)))
  else{
     categ.init <- trim(categ.init) #trim vector elements to remove leading and trailing blanks
     categ <- categ.init[nchar(categ.init)>0]
  }
  q <- length(categ)

  # creating the weights matrix

  if (is.character(weights)){
     if (weights=="quadratic")
	  weights.mat<-quadratic.weights(categ)
     else if (weights=="ordinal")
	  weights.mat<-ordinal.weights(categ)
     else if (weights=="linear")
	  weights.mat<-linear.weights(categ)
     else if (weights=="radical")
	  weights.mat<-radical.weights(categ)
     else if (weights=="ratio")
	  weights.mat<-ratio.weights(categ)
     else if (weights=="circular")
	  weights.mat<-circular.weights(categ)
     else if (weights=="bipolar")
	  weights.mat<-bipolar.weights(categ)
     else weights.mat<-identity.weights(categ)
  }else weights.mat= as.matrix(weights)
  
  # creating the nxq agreement matrix representing the distribution of raters by subjects and category

  agree.mat <- matrix(0,nrow=n,ncol=q)
  for(k in 1:q){
	if (is.numeric(ratings.mat)){
          k.mis <-(ratings.mat==categ[k])
          in.categ.k <- replace(k.mis,is.na(k.mis),FALSE)
	    agree.mat[,k] <- in.categ.k%*%rep(1,r) 
      }else
          agree.mat[,k] <- (ratings.mat==categ[k])%*%rep(1,r)
  }
  agree.mat.w <- t(weights.mat%*%t(agree.mat))

  # calculating gwet's ac1 coefficient

  ri.vec <- agree.mat%*%rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more

  pi.vec <- t(t(rep(1/n,n))%*%(agree.mat/(ri.vec%*%t(rep(1,q)))))
  pe <- sum(weights.mat) * sum(pi.vec*(1-pi.vec)) / (q*(q-1))
  gwet.ac1 <- (pa-pe)/(1-pe)

  # calculating variance, stderr & p-value of gwet's ac1 coefficient
  
  den.ivec <- ri.vec*(ri.vec-1)
  den.ivec <- den.ivec - (den.ivec==0) # this operation replaces each 0 value with -1 to make the next ratio calculation always possible.
  pa.ivec <- sum.q/den.ivec

  pe.r2 <- pe*(ri.vec>=2)
  ac1.ivec <- (n/n2more)*(pa.ivec-pe.r2)/(1-pe)
  pe.ivec <- (sum(weights.mat)/(q*(q-1))) * (agree.mat%*%(1-pi.vec))/ri.vec
  ac1.ivec.x <- ac1.ivec - 2*(1-gwet.ac1) * (pe.ivec-pe)/(1-pe)
  
  var.ac1 <- ((1-f)/(n*(n-1))) * sum((ac1.ivec.x - gwet.ac1)^2)
  stderr <- sqrt(var.ac1)# ac1's standard error
  p.value <- 2*(1-pt(gwet.ac1/stderr,n-1))
  
  lcb <- gwet.ac1 - stderr*qt(1-(1-conflev)/2,n-1) # lower confidence bound
  ucb <- min(1,gwet.ac1 + stderr*qt(1-(1-conflev)/2,n-1)) # upper confidence bound
#  if(print==TRUE) {
#    cat("Gwet's AC1/AC2 Coefficient\n")
#    cat('==========================\n')	
#    cat('Percent agreement:',pa,'Percent chance agreement:',pe,'\n')
#    if (weights=="unweighted") {
#       cat('AC1 coefficient:',gwet.ac1,'Standard error:',stderr,'\n')
#    }
#    else {
#	 cat('AC2 coefficient:',gwet.ac1,'Standard error:',stderr,'\n')
#       if (!is.numeric(weights)) {
#	    cat('Weights: ', weights,'\n')
#	 }
#      else
#	    cat('Weights: Custom Weights\n')
#   }
#   cat(conflev*100,'% Confidence Interval: (',lcb,',',ucb,')\n')
#   cat('P-value: ',p.value,'\n')
# }
  return(list(pa=pa,pe=pe,value=gwet.ac1,se=stderr,p=p.value))
}

#=====================================================================================
#fleiss.kappa.raw: This function computes Fleiss' generalized kappa coefficient (see Fleiss(1971)) and 
#		   its standard error for 3 raters or more when input dataset is a nxr matrix of alphanumeric 
#		   ratings from n subjects and r raters.
#-------------
#The input data "ratings" is a nxr data frame of raw alphanumeric ratings
#from n subjects and r raters. Exclude all subjects that are not rated by any rater.
#Bibliography:
#Fleiss, J. L. (1981). Statistical Methods for Rates and Proportions. John Wiley & Sons.
#======================================================================================
fleiss.kappa.raw <- function(ratings,weights="unweighted",conflev=0.95,N=Inf,print=TRUE){ 
  ratings.mat <- as.matrix(ratings) 
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N # final population correction 

  # creating a vector containing all categories used by the raters
 
  categ.init <- unique(as.vector(ratings.mat))
  if (is.numeric(categ.init)){
     categ <- sort(as.vector(na.omit(categ.init)))
  }else{
     categ.init <- trim(categ.init) #trim vector elements to remove leading and trailing blanks
     categ <- categ.init[nchar(categ.init)>0]
  }
  q <- length(categ)

  # creating the weights matrix

  if (is.character(weights)){
     if (weights=="quadratic")
	  weights.mat<-quadratic.weights(categ)
     else if (weights=="ordinal")
	  weights.mat<-ordinal.weights(categ)
     else if (weights=="linear")
	  weights.mat<-linear.weights(categ)
     else if (weights=="radical")
	  weights.mat<-radical.weights(categ)
     else if (weights=="ratio")
	  weights.mat<-ratio.weights(categ)
     else if (weights=="circular")
	  weights.mat<-circular.weights(categ)
     else if (weights=="bipolar")
	  weights.mat<-bipolar.weights(categ)
     else weights.mat<-identity.weights(categ)
  }else weights.mat= as.matrix(weights)
  
  # creating the nxq agreement matrix representing the distribution of raters by subjects and category

  agree.mat <- matrix(0,nrow=n,ncol=q)
  for(k in 1:q){
	if (is.numeric(ratings.mat)){
          k.mis <-(ratings.mat==categ[k])
          in.categ.k <- replace(k.mis,is.na(k.mis),FALSE)
	    agree.mat[,k] <- in.categ.k%*%rep(1,r) 
      }else
          agree.mat[,k] <- (ratings.mat==categ[k])%*%rep(1,r)
  }
  agree.mat.w <- t(weights.mat%*%t(agree.mat))

  # calculating fleiss's generalized kappa coefficient

  ri.vec <- agree.mat%*%rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more

  pi.vec <- t(t(rep(1/n,n))%*%(agree.mat/(ri.vec%*%t(rep(1,q)))))
  pe <- sum(weights.mat * (pi.vec%*%t(pi.vec)))
  fleiss.kappa <- (pa-pe)/(1-pe)

  # calculating variance, stderr & p-value of gwet's ac1 coefficient
  
  den.ivec <- ri.vec*(ri.vec-1)
  den.ivec <- den.ivec - (den.ivec==0) # this operation replaces each 0 value with -1 to make the next ratio calculation always possible.
  pa.ivec <- sum.q/den.ivec

  pe.r2 <- pe*(ri.vec>=2)
  kappa.ivec <- (n/n2more)*(pa.ivec-pe.r2)/(1-pe)
  pi.vec.wk. <- weights.mat%*%pi.vec
  pi.vec.w.k <- t(weights.mat)%*%pi.vec
  pi.vec.w <- (pi.vec.wk. + pi.vec.w.k)/2

  pe.ivec <- (agree.mat%*%pi.vec.w)/ri.vec
  kappa.ivec.x <- kappa.ivec - 2*(1-fleiss.kappa) * (pe.ivec-pe)/(1-pe)
  
  var.fleiss <- ((1-f)/(n*(n-1))) * sum((kappa.ivec.x - fleiss.kappa)^2)
  stderr <- sqrt(var.fleiss)# kappa's standard error
  p.value <- 2*(1-pt(fleiss.kappa/stderr,n-1))
  
  lcb <- fleiss.kappa - stderr*qt(1-(1-conflev)/2,n-1) # lower confidence bound
  ucb <- min(1,fleiss.kappa + stderr*qt(1-(1-conflev)/2,n-1)) # upper confidence bound
#  if(print==TRUE){
#    cat("Fleiss' Kappa Coefficient\n")
#    cat('==========================\n')	
#    cat('Percent agreement:',pa,'Percent chance agreement:',pe,'\n')
#    cat('Fleiss kappa coefficient:',fleiss.kappa,'Standard error:',stderr,'\n')
#    if (weights!="unweighted") {
#       if (!is.numeric(weights)) {
#	    cat('Weights: ', weights,'\n')
#	 }
#       else
#	    cat('Weights: Custom Weights\n')
#    }
#    cat(conflev*100,'% Confidence Interval: (',lcb,',',ucb,')\n')
#    cat('P-value: ',p.value,'\n')
#  }
#  invisible(c(pa,pe,fleiss.kappa,stderr,p.value))
#}
  return(list(pa=pa,pe=pe,value=fleiss.kappa,se=stderr,p=p.value))
}

#=====================================================================================
#krippen.alpha.raw: This function computes Krippendorff's alpha coefficient (see Krippendorff(1970, 1980)) and 
#		   its standard error for 3 raters or more when input dataset is a nxr matrix of alphanumeric 
#		   ratings from n subjects and r raters.
#-------------
#The algorithm used to compute krippendorff's alpha is very different from anything that was published on this topic. Instead,
#it follows the equations presented by K. Gwet (2010)
#The input data "ratings" is a nxr data frame of raw alphanumeric ratings
#from n subjects and r raters. Exclude all subjects that are not rated by any rater.
#Bibliography:
#Gwet, K. (2012). Handbook of Inter-Rater Reliability: The Definitive Guide to Measuring the Extent of Agreement Among 
#	Multiple Raters, 3rd Edition.  Advanced Analytics, LLC; 3rd edition (March 2, 2012)
#Krippendorff (1970). "Bivariate agreement coefficients for reliability of data." Sociological Methodology,2,139-150
#Krippendorff (1980). Content analysis: An introduction to its methodology (2nd ed.), New-bury Park, CA: Sage.
#======================================================================================
krippen.alpha.raw <- function(ratings,weights="unweighted",conflev=0.95,N=Inf,print=TRUE){ 

  ratings.mat <- as.matrix(ratings) 
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N # final population correction 

  # creating a vector containing all categories used by the raters
 
  categ.init <- unique(as.vector(ratings.mat))
  if (is.numeric(categ.init)){
     categ <- sort(as.vector(na.omit(categ.init)))
  }else{
     categ.init <- trim(categ.init) #trim vector elements to remove leading and trailing blanks
     categ <- categ.init[nchar(categ.init)>0]
  }
  q <- length(categ)

  # creating the weights matrix

  if (is.character(weights)){
     if (weights=="quadratic")
	  weights.mat<-quadratic.weights(categ)
     else if (weights=="ordinal")
	  weights.mat<-ordinal.weights(categ)
     else if (weights=="linear")
	  weights.mat<-linear.weights(categ)
     else if (weights=="radical")
	  weights.mat<-radical.weights(categ)
     else if (weights=="ratio")
	  weights.mat<-ratio.weights(categ)
     else if (weights=="circular")
	  weights.mat<-circular.weights(categ)
     else if (weights=="bipolar")
	  weights.mat<-bipolar.weights(categ)
     else weights.mat<-identity.weights(categ)
  }else weights.mat= as.matrix(weights)
  
  # creating the nxq agreement matrix representing the distribution of raters by subjects and category

  agree.mat <- matrix(0,nrow=n,ncol=q)
  for(k in 1:q){
	if (is.numeric(ratings.mat)){
          k.mis <-(ratings.mat==categ[k])
          in.categ.k <- replace(k.mis,is.na(k.mis),FALSE)
	    agree.mat[,k] <- in.categ.k%*%rep(1,r) 
      }else
          agree.mat[,k] <- (ratings.mat==categ[k])%*%rep(1,r)
  }
  agree.mat.w <- t(weights.mat%*%t(agree.mat))

  # calculating krippendorff's alpha coefficient

  ri.vec <- agree.mat%*%rep(1,q)
  agree.mat<-agree.mat[(ri.vec>=2),]
  agree.mat.w <- agree.mat.w[(ri.vec>=2),]
  ri.vec <- ri.vec[(ri.vec>=2)]
  ri.mean <- mean(ri.vec)
  n <- nrow(agree.mat)
  epsi <- 1/sum(ri.vec)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  pa <- (1-epsi)* sum(sum.q/(ri.mean*(ri.vec-1)))/n + epsi

  pi.vec <- t(t(rep(1/n,n))%*%(agree.mat/ri.mean))
  pe <- sum(weights.mat * (pi.vec%*%t(pi.vec)))
  krippen.alpha <- (pa-pe)/(1-pe)

  # calculating variance, stderr & p-value of gwet's ac1 coefficient
  
  den.ivec <- ri.mean*(ri.vec-1)
  pa.ivec <- sum.q/den.ivec
  pa.v <- mean(pa.ivec)
  pa.ivec <- (1-epsi)*(pa.ivec-pa.v*(ri.vec-ri.mean)/ri.mean) + epsi

  krippen.ivec <- (pa.ivec-pe)/(1-pe)
  pi.vec.wk. <- weights.mat%*%pi.vec
  pi.vec.w.k <- t(weights.mat)%*%pi.vec

  pi.vec.w <- (pi.vec.wk. + pi.vec.w.k)/2

  pe.ivec <- (agree.mat%*%pi.vec.w)/ri.mean - sum(pi.vec) * (ri.vec-ri.mean)/ri.mean
  krippen.ivec.x <- krippen.ivec - (1-krippen.alpha) * (pe.ivec-pe)/(1-pe)
  
  var.krippen <- ((1-f)/(n*(n-1))) * sum((krippen.ivec.x - krippen.alpha)^2)
  stderr <- sqrt(var.krippen)# alpha's standard error
  p.value <- 2*(1-pt(krippen.alpha/stderr,n-1))
  
  lcb <- krippen.alpha - stderr*qt(1-(1-conflev)/2,n-1) # lower confidence bound
  ucb <- min(1,krippen.alpha + stderr*qt(1-(1-conflev)/2,n-1)) # upper confidence bound
  if(print==TRUE){
    cat("Krippendorff's Alpha Coefficient\n")
    cat('==========================\n')	
    cat('Percent agreement:',pa,'Percent chance agreement:',pe,'\n')
    cat('Krippendorff alpha coefficient:',krippen.alpha,'Standard error:',stderr,'\n')
    if (weights!="unweighted") {
       if (!is.numeric(weights)) {
	    cat('Weights: ', weights,'\n')
	 }
       else
	    cat('Weights: Custom Weights\n')
    }
    cat(conflev*100,'% Confidence Interval: (',lcb,',',ucb,')\n')
    cat('P-value: ',p.value,'\n')
  }
  invisible(c(pa,pe,krippen.alpha,stderr,p.value))
}


#===========================================================================================
#conger.kappa.raw: Conger's kappa coefficient (see Conger(1980)) and its standard error for multiple raters when input 
#		   dataset is a nxr matrix of alphanumeric ratings from n subjects and r raters 
#-------------
#The input data "ratings" is a nxr data frame of raw alphanumeric ratings
#from n subjects and r raters. Exclude all subjects that are not rated by any rater.
#Bibliography:
#Conger, A. J. (1980), ``Integration and Generalization of Kappas for Multiple Raters,"
#		Psychological Bulletin, 88, 322-328.
#======================================================================================
conger.kappa.raw <- function(ratings,weights="unweighted",conflev=0.95,N=Inf,print=TRUE){ 
  ratings.mat <- as.matrix(ratings) 
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N # final population correction 

  # creating a vector containing all categories used by the raters
 
  categ.init <- unique(as.vector(ratings.mat))
  if (is.numeric(categ.init))
     categ <- sort(as.vector(na.omit(categ.init)))
  else {
     categ.init <- trim(categ.init) #trim vector elements to remove leading and trailing blanks
     categ <- categ.init[nchar(categ.init)>0]
  }
  q <- length(categ)

  # creating the weights matrix

  if (is.character(weights)){
     if (weights=="quadratic")
	  weights.mat<-quadratic.weights(categ)
     else if (weights=="ordinal")
	  weights.mat<-ordinal.weights(categ)
     else if (weights=="linear")
	  weights.mat<-linear.weights(categ)
     else if (weights=="radical")
	  weights.mat<-radical.weights(categ)
     else if (weights=="ratio")
	  weights.mat<-ratio.weights(categ)
     else if (weights=="circular")
	  weights.mat<-circular.weights(categ)
     else if (weights=="bipolar")
	  weights.mat<-bipolar.weights(categ)
     else weights.mat<-identity.weights(categ)
  }else weights.mat= as.matrix(weights)
  
  # creating the nxq agreement matrix representing the distribution of raters by subjects and category

  agree.mat <- matrix(0,nrow=n,ncol=q)
  for(k in 1:q){
	if (is.numeric(ratings.mat)){
          k.mis <-(ratings.mat==categ[k])
          in.categ.k <- replace(k.mis,is.na(k.mis),FALSE)
	    agree.mat[,k] <- in.categ.k%*%rep(1,r) 
      }else
          agree.mat[,k] <- (ratings.mat==categ[k])%*%rep(1,r)
  }
  agree.mat.w <- t(weights.mat%*%t(agree.mat))

  # creating the rxq rater-category matrix representing the distribution of subjects by rater and category

  classif.mat <- matrix(0,nrow=r,ncol=q)
  for(k in 1:q){
	if (is.numeric(ratings.mat)){
	    with.mis <-(t(ratings.mat)==categ[k])
          without.mis <- replace(with.mis,is.na(with.mis),FALSE)
          classif.mat[,k] <- without.mis%*%rep(1,n)
	}else
	    classif.mat[,k] <- (t(ratings.mat)==categ[k])%*%rep(1,n)
  }

  # calculating conger's kappa coefficient

  ri.vec <- agree.mat%*%rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more

  ng.vec <- classif.mat%*%rep(1,q)
  pgk.mat <- classif.mat/(ng.vec%*%rep(1,q))
  p.mean.k <- (t(pgk.mat)%*%rep(1,r))/r 
  s2kl.mat <- (t(pgk.mat)%*%pgk.mat - r * p.mean.k%*%t(p.mean.k))/(r-1)
  pe <- sum(weights.mat * (p.mean.k%*%t(p.mean.k) -  s2kl.mat/r))
  conger.kappa <- (pa-pe)/(1-pe)

  # calculating variance, stderr & p-value of conger's kappa coefficient
  
  bkl.mat <- (weights.mat+t(weights.mat))/2
  pe.ivec1 <- r*(agree.mat%*%t(t(p.mean.k)%*%bkl.mat))
  pe.ivec2 = rep(0,n)
  if (is.numeric(ratings.mat)){
      for(l in 1:q){
          l.mis <-(ratings.mat==categ[l])
          delta.ig.mat <- replace(l.mis,is.na(l.mis),FALSE)
          pe.ivec2 <- pe.ivec2 + delta.ig.mat%*%(pgk.mat%*%bkl.mat[,l]) 
      }
      }else{
	    for(k in 1:q){
             delta.ig.mat <- (ratings.mat==categ[k])
             pe.ivec2 <- pe.ivec2 + delta.ig.mat%*%(pgk.mat%*%bkl.mat[,l])
          }
      }
  pe.ivec <- (pe.ivec1-pe.ivec2)/(r*(r-1)) 
  den.ivec <- ri.vec*(ri.vec-1)
  den.ivec <- den.ivec - (den.ivec==0) # this operation replaces each 0 value with -1 to make the next ratio calculation always possible.
  pa.ivec <- sum.q/den.ivec
  pe.r2 <- pe*(ri.vec>=2)
  conger.ivec <- (n/n2more)*(pa.ivec-pe.r2)/(1-pe) 
  conger.ivec.x <- conger.ivec - 2*(1-conger.kappa) * (pe.ivec-pe)/(1-pe)
  
  var.conger <- ((1-f)/(n*(n-1))) * sum((conger.ivec.x - conger.kappa)^2)
  stderr <- sqrt(var.conger)# conger's kappa standard error
  p.value <- 2*(1-pt(conger.kappa/stderr,n-1))
  
  lcb <- conger.kappa - stderr*qt(1-(1-conflev)/2,n-1) # lower confidence bound
  ucb <- min(1,conger.kappa + stderr*qt(1-(1-conflev)/2,n-1)) # upper confidence bound
  if(print==TRUE) {
    cat("Conger's Kappa Coefficient\n")
    cat('==========================\n')	
    cat('Percent agreement: ',pa,'Percent chance agreement: ',pe,'\n')
    cat("Conger's kappa coefficient: ",conger.kappa,'Standard error:',stderr,'\n')
    if (weights!="unweighted") {
       if (!is.numeric(weights)) {
          cat('Weights: ', weights,'\n')
       }
    else
	cat('Weights: Custom Weights\n')
    }
    cat(conflev*100,'% Confidence Interval: (',lcb,',',ucb,')\n')
    cat('P-value: ',p.value,'\n')
  }
  invisible(c(pa,pe,conger.kappa,stderr,p.value))
}

#===========================================================================================
#bp.coeff.raw: Brennan-Prediger coefficient (see Brennan & Prediger(1981)) and its standard error for multiple raters when input 
#		   dataset is a nxr matrix of alphanumeric ratings from n subjects and r raters 
#-------------
#The input data "ratings" is a nxr data frame of raw alphanumeric ratings
#from n subjects and r raters. Exclude all subjects that are not rated by any rater.
#Bibliography:
#Brennan, R.L., and Prediger, D. J. (1981). ``Coefficient Kappa: some uses, misuses, and alternatives."
#           Educational and Psychological Measurement, 41, 687-699.
#======================================================================================
bp.coeff.raw <- function(ratings,weights="unweighted",conflev=0.95,N=Inf,print=TRUE){ 
  ratings.mat <- as.matrix(ratings) 
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N # final population correction 

  # creating a vector containing all categories used by the raters
 
  categ.init <- unique(as.vector(ratings.mat))
  if (is.numeric(categ.init))
     categ <- sort(as.vector(na.omit(categ.init)))
  else{
     categ.init <- trim(categ.init) #trim vector elements to remove leading and trailing blanks
     categ <- categ.init[nchar(categ.init)>0]
  }
  q <- length(categ)

  # creating the weights matrix

  if (is.character(weights)){
     if (weights=="quadratic")
	  weights.mat<-quadratic.weights(categ)
     else if (weights=="ordinal")
	  weights.mat<-ordinal.weights(categ)
     else if (weights=="linear")
	  weights.mat<-linear.weights(categ)
     else if (weights=="radical")
	  weights.mat<-radical.weights(categ)
     else if (weights=="ratio")
	  weights.mat<-ratio.weights(categ)
     else if (weights=="circular")
	  weights.mat<-circular.weights(categ)
     else if (weights=="bipolar")
	  weights.mat<-bipolar.weights(categ)
     else weights.mat<-identity.weights(categ)
  }else weights.mat= as.matrix(weights)
  
  # creating the nxq agreement matrix representing the distribution of raters by subjects and category

  agree.mat <- matrix(0,nrow=n,ncol=q)
  for(k in 1:q){
	if (is.numeric(ratings.mat)){
          k.mis <-(ratings.mat==categ[k])
          in.categ.k <- replace(k.mis,is.na(k.mis),FALSE)
	    agree.mat[,k] <- in.categ.k%*%rep(1,r) 
      }else
          agree.mat[,k] <- (ratings.mat==categ[k])%*%rep(1,r)
  }
  agree.mat.w <- t(weights.mat%*%t(agree.mat))

  # calculating gwet's ac1 coefficient

  ri.vec <- agree.mat%*%rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more

  pi.vec <- t(t(rep(1/n,n))%*%(agree.mat/(ri.vec%*%t(rep(1,q)))))
  pe <- sum(weights.mat) / (q^2)
  bp.coeff <- (pa-pe)/(1-pe)

  # calculating variance, stderr & p-value of gwet's ac1 coefficient
  
  den.ivec <- ri.vec*(ri.vec-1)
  den.ivec <- den.ivec - (den.ivec==0) # this operation replaces each 0 value with -1 to make the next ratio calculation always possible.
  pa.ivec <- sum.q/den.ivec

  pe.r2 <- pe*(ri.vec>=2)
  bp.ivec <- (n/n2more)*(pa.ivec-pe.r2)/(1-pe)
  var.bp <- ((1-f)/(n*(n-1))) * sum((bp.ivec - bp.coeff)^2)
  stderr <- sqrt(var.bp)# BP's standard error
  p.value <- 2*(1-pt(bp.coeff/stderr,n-1))
  
  lcb <- bp.coeff - stderr*qt(1-(1-conflev)/2,n-1) # lower confidence bound
  ucb <- min(1,bp.coeff + stderr*qt(1-(1-conflev)/2,n-1)) # upper confidence bound
#  if(print==TRUE) {
#    cat("Brennan-Prediger Coefficient\n")
#    cat('============================\n')	
#    cat('Percent agreement:',pa,'Percent chance agreement:',pe,'\n')
#    cat('B-P coefficient:',bp.coeff,'Standard error:',stderr,'\n')
#    if (weights!="unweighted") {
#       if (!is.numeric(weights)) {
#	    cat('Weights: ', weights,'\n')
#	 }
#       else
#	    cat('Weights: Custom Weights\n')
#    }
#    cat(conflev*100,'% Confidence Interval: (',lcb,',',ucb,')\n')
#    cat('P-value: ',p.value,'\n')
#  }
#  invisible(c(pa,pe,bp.coeff,stderr,p.value))
#}
  return(list(pa=pa,pe=pe,value=bp.coeff,se=stderr,p=p.value))
}

gwets.ac2 <- function(codings,weights) {
	
	### CREATE RATING-MATRIX
	r <- sum(!is.na(codings))/dim(codings)[1]
	r_ <- matrix(0,nrow=dim(codings)[1],ncol=length(table(codings)))
	colnames(r_) <- as.numeric(names(table(codings)))
	codes <- table(codings)
	for (o in 1:dim(codings)[1])
	{
		for (p in 1:length(table(codings)))
		{
			r_[o,p] <- sum(codings[o,]==as.numeric(colnames(r_)[p]))
		}	
	}

	### CREATE WEIGHTING-MATRIX
if (weights=="quad")
{
	w <- matrix(0,nrow=length(table(codings)),ncol=length(table(codings)))
	colnames(w) <- as.numeric(names(table(codings)))
	rownames(w) <- as.numeric(names(table(codings)))
	for (i in 1:dim(w)[1])
	{
		for (j in 1:dim(w)[2])
		{	
			w[i,j] <- c(1	-	(	(as.numeric(names(table(codings)))[i]-as.numeric(names(table(codings)))[j]	)	/	max(as.numeric(names(table(codings)))) )^2	)
		}
	}
}

if (weights=="lin")
{
	w <- matrix(0,nrow=length(table(codings)),ncol=length(table(codings)))
	colnames(w) <- as.numeric(names(table(codings)))
	rownames(w) <- as.numeric(names(table(codings)))
	for (i in 1:dim(w)[1])
	{
		for (j in 1:dim(w)[2])
		{	
			w[i,j] <- c(1	-	(	abs((as.numeric(names(table(codings)))[i]-as.numeric(names(table(codings)))[j]	))	/	max(as.numeric(names(table(codings)))) )	)
		}
	}
}

if (weights=="ord")
{
	w <- matrix(0,nrow=length(table(codings)),ncol=length(table(codings)))
	colnames(w) <- as.numeric(names(table(codings)))
	rownames(w) <- as.numeric(names(table(codings)))
	for (i in 1:dim(w)[1])
	{
		for (j in 1:dim(w)[2])
		{	
			w[i,j] <- 1-abs((	i-j	)	/	length(names(table(codings)))  )
		}
	}
}















###################################
###################################
### DEFINING THE DATA

rellevels <- list(0,.05,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.60,.65,.70,.75,.80,.85,.90,.95,1.00)
sampsize <- list(10,20,30,40,50,100,200,300,400,500,700,1000)
means <- c(5,3,4.5,3.5,6)
sds <-c (1,1.5)
corlist <- c(.70,.46,.32,.29,.36,.25,.19,.18,.09,.07)
corlist2 <- c(.75,.50,.40,.30,.40,.30,.20,.20,.10,.05)

#real <- rbinom(n=1000,size=5,prob=.3)
#rand.bi <- rbinom(n=1000,size=5,prob=.3)
#rand.un <- runif(n=1000,min=0,max=5)
#rand.hi <- rbinom(n=1000,size=5,prob=.8)
#rand.lo <- rbinom(n=1000,size=5,prob=.1)

q <- seq(0,1,0.05)

c1 <- matrix(NA,ncol=21,nrow=1000)
c2 <- matrix(NA,ncol=21,nrow=1000)

data <- list()

for (i in 1:length(sampsize))
{
covarmat <- matrix(c(	1,	.75,	.50,	.4,	.3,
			  	.75,	1,	.4,	.3,	.2,
				.50,  .40,	1, 	.2,	.1,
				.40,	.30,	.20,	1,	.05,
				.30,	.20,	.10,	.05,	1),ncol=5,nrow=5)


x <- rmvnorm(n=sampsize[[i]],mean=means,sigma=covarmat,method="eigen")

data[[i]] <- x

}

real <- data


###################################
###################################
### FUNCTIONS FOR SIMULATED "CODING"

code <- function(level,truval,trudis,expag){
		if (expag=="scale") {guess <- sample(as.numeric(names(table(trudis))),1)}		
		if (expag=="marginal") {guess <- sample(as.numeric(names(table(trudis))),1,prob=prop.table(table(trudis)))}		
		score <- level*truval+(1-level)*guess
		}

code <- function(level,truval,trudis,expag){
		if (expag=="scale") {guess <- round(runif(n=1,min=min(trudis),max=max(trudis)))}		
		if (expag=="marginal") {guess <- sample(trudis,1)}		
		score <- level*truval+(1-level)*guess
		}



###################################
###################################
### RUNNING THE SIMULATION OF CODING

	arr.cor.s <- array(NA,dim=c(1000,length(rellevels),length(sampsize),dim(combn(5,2))[2]))
	arr.cor.m <- array(NA,dim=c(1000,length(rellevels),length(sampsize),dim(combn(5,2))[2]))

	rel.bp.s <- array(NA,dim=c(1000,length(rellevels),length(sampsize),dim(combn(5,2))[2]))
	rel.bp.m <- array(NA,dim=c(1000,length(rellevels),length(sampsize),dim(combn(5,2))[2]))
	rel.ka.s <- array(NA,dim=c(1000,length(rellevels),length(sampsize),dim(combn(5,2))[2]))
	rel.ka.m <- array(NA,dim=c(1000,length(rellevels),length(sampsize),dim(combn(5,2))[2]))
	rel.fk.s <- array(NA,dim=c(1000,length(rellevels),length(sampsize),dim(combn(5,2))[2]))
	rel.fk.m <- array(NA,dim=c(1000,length(rellevels),length(sampsize),dim(combn(5,2))[2]))
	rel.ck.s <- array(NA,dim=c(1000,length(rellevels),length(sampsize),dim(combn(5,2))[2]))
	rel.ck.m <- array(NA,dim=c(1000,length(rellevels),length(sampsize),dim(combn(5,2))[2]))
	
	
##### ATTENTION: the following loop is the core of the simulation and takes a LOT of time. For a regular PC, 1-2 complete days are not unusual. Process is finished when the count displayed reaches 1000.
	
for (z in 1:1000)
{
	coderes1s <- array(data=NA,dim=c(length(rellevels),length(sampsize),ncol(real[[1]]),nrow(real[[length(sampsize)]])))
	coderes2s <- array(data=NA,dim=c(length(rellevels),length(sampsize),ncol(real[[1]]),nrow(real[[length(sampsize)]])))
	coderes1m <- array(data=NA,dim=c(length(rellevels),length(sampsize),ncol(real[[1]]),nrow(real[[length(sampsize)]])))
	coderes2m <- array(data=NA,dim=c(length(rellevels),length(sampsize),ncol(real[[1]]),nrow(real[[length(sampsize)]])))


	for (m in 1:length(sampsize))
	{
	real[[m]] <- data[[12]][sample(1:1000,size=sampsize[[m]]),]
	real[[m]] <- round(real[[m]])
	}

	for (r in 1:length(rellevels))
	{
		for (i in 1:length(sampsize))
		{
			for (j in 1:ncol(real[[1]]))
			{
				for (k in 1:nrow(real[[i]]))	
				{
				
				# CODING
				
				coderes1s[r,i,j,k] <- code(level=rellevels[[r]],truval=real[[i]][k,j],trudis=real[[i]][,j],"scale")
				coderes2s[r,i,j,k] <- code(level=rellevels[[r]],truval=real[[i]][k,j],trudis=real[[i]][,j],"scale")
				coderes1m[r,i,j,k] <- code(level=rellevels[[r]],truval=real[[i]][k,j],trudis=real[[i]][,j],"marginal")
				coderes2m[r,i,j,k] <- code(level=rellevels[[r]],truval=real[[i]][k,j],trudis=real[[i]][,j],"marginal")

				}	
				
				# RELIABILITY COMPUTATION
				
				code.s <- cbind(coderes1s[r,i,j,],coderes2s[r,i,j,])
				code.s <- matrix(code.s[rowSums(is.na(code.s))==0],ncol=2)
				code.m <- cbind(coderes1m[r,i,j,],coderes2m[r,i,j,])
				code.m <- matrix(code.m[rowSums(is.na(code.m))==0],ncol=2)
				rel.bp.s[z,r,i,j] <- bp.coeff.raw(code.s,weights="ordinal")$value
				rel.bp.m[z,r,i,j] <- bp.coeff.raw(code.m,weights="ordinal")$value
				rel.ka.s[z,r,i,j] <- krippen.alpha.raw(code.s,weights="ordinal",print=FALSE)[3]
				rel.ka.m[z,r,i,j] <- krippen.alpha.raw(code.m,weights="ordinal",print=FALSE)[3]
				rel.fk.s[z,r,i,j] <- fleiss.kappa.raw(code.s,weights="ordinal")$value
				rel.fk.m[z,r,i,j] <- fleiss.kappa.raw(code.m,weights="ordinal")$value
				rel.ck.s[z,r,i,j] <- conger.kappa.raw(code.s,weights="ordinal",print=FALSE)[3]
				rel.ck.m[z,r,i,j] <- conger.kappa.raw(code.m,weights="ordinal",print=FALSE)[3]
			}	
		}
	}

###################################
###################################
### COMPUTE CORRELATION ESTIMATES

	corrs <- array(NA,dim=c(length(rellevels),length(sampsize),dim(combn(5,2))[2]))

	corrm <- array(NA,dim=c(length(rellevels),length(sampsize),dim(combn(5,2))[2]))

	varpacks <- combn(5,2)

	for (r in 1:dim(corrs)[1])
	{
		for (i in 1:dim(corrs)[2])
		{
			for (j in 1:dim(corrs)[3])
			{
			corrs[r,i,j] <- cor(coderes1s[r,i,varpacks[1,j],],coderes1s[r,i,varpacks[2,j],],use="complete.obs")
			corrm[r,i,j] <- cor(coderes1m[r,i,varpacks[1,j],],coderes1m[r,i,varpacks[2,j],],use="complete.obs")
			}
		}
	}

arr.cor.s[z,,,] <- corrs
arr.cor.m[z,,,] <- corrm

flush.console()
print(z)
}

##### This is the end of the long loop. Congratulations. The simulation has been run successfully.





###############################
###############################
### REARRANGE THE DATA FOR EASIER VISUALIZATION / TABULATION / ANALYSIS

reliability.test <- data.frame(matrix(NA,ncol=11,nrow=21))
names(reliability.test) <- c("bp.m","bp.s","ka.m","ka.s","ck.m","ck.s","fk.m","fk.s","gg.m","gg.s","accuracy")

for (i in 1:21)
{
reliability.test[i,"bp.m"] <- mean(rel.bp.m[,i,,],na.rm=T)
reliability.test[i,"bp.s"] <- mean(rel.bp.s[,i,,],na.rm=T)
reliability.test[i,"ka.m"] <- mean(rel.ka.m[,i,,],na.rm=T)
reliability.test[i,"ka.s"] <- mean(rel.ka.s[,i,,],na.rm=T)
reliability.test[i,"fk.m"] <- mean(rel.fk.m[,i,,],na.rm=T)
reliability.test[i,"fk.s"] <- mean(rel.fk.s[,i,,],na.rm=T)
reliability.test[i,"ck.m"] <- mean(rel.ck.m[,i,,],na.rm=T)
reliability.test[i,"ck.s"] <- mean(rel.ck.s[,i,,],na.rm=T)
reliability.test[i,"gg.m"] <- mean(rel.gg.m[,i,,],na.rm=T)
reliability.test[i,"gg.s"] <- mean(rel.gg.s[,i,,],na.rm=T)
}
reliability.test$accuracy <- seq(0,1,0.05)

for (i in 1:21)
{
reliability.test[i,"bp.m"] <- mean(rel.bp.m[,i,12,],na.rm=T)
reliability.test[i,"bp.s"] <- mean(rel.bp.s[,i,12,],na.rm=T)
reliability.test[i,"ka.m"] <- mean(rel.ka.m[,i,12,],na.rm=T)
reliability.test[i,"ka.s"] <- mean(rel.ka.s[,i,12,],na.rm=T)
reliability.test[i,"fk.m"] <- mean(rel.fk.m[,i,12,],na.rm=T)
reliability.test[i,"fk.s"] <- mean(rel.fk.s[,i,12,],na.rm=T)
reliability.test[i,"ck.m"] <- mean(rel.ck.m[,i,12,],na.rm=T)
reliability.test[i,"ck.s"] <- mean(rel.ck.s[,i,12,],na.rm=T)
reliability.test[i,"gg.m"] <- mean(rel.gg.m[,i,12,],na.rm=T)
reliability.test[i,"gg.s"] <- mean(rel.gg.s[,i,12,],na.rm=T)
}
reliability.test$accuracy <- seq(0,1,0.05)


df.reli <- melt(reliability.test,id.vars="accuracy")

df.reli[,c("coef","guess")] <- t(matrix(unlist(str_split(df.reli$variable,pattern="\\.")),nrow=2))

df.reli$Coefficient <- Recode(df.reli$coef,"'gg'='Geiss´ Gamma';'bp'='Kappa Brennan-Prediger';'ka'='Alpha Krippendorff';'fk'='Kappa Fleiss';'ck'='Kappa Cohen'")
df.reli$Guess <- Recode(df.reli$guess,"'s'='Scale based coder guessing';'m'='Marginal based coder guessing'")


siglim <- c(0.0005,0.005,0.025,0.500,0.975,0.995,0.9995)

arr.corr.s <- array(NA,c(7,21,12,10))

arr.corr.s <- matrix(NA,nrow=21*12*10,ncol=10)

colnames(arr.corr.s) <- c("p.0005","p.005","p.025","p.500","p.975","p.995","p.9995","quality","sample","effsize")

for (q in 1:21)
{
	for (s in 1:12)
	{
		for (c in 1:10)
		{
		arr.corr.s[(q-1)*120+(s-1)*10+c,1:7] <- quantile(arr.cor.s[,q,s,c],siglim,na.rm=T)
		arr.corr.s[(q-1)*120+(s-1)*10+c,"quality"] <- rellevels[[q]]
		arr.corr.s[(q-1)*120+(s-1)*10+c,"sample"] <- sampsize[[s]]
		arr.corr.s[(q-1)*120+(s-1)*10+c,"effsize"] <- corlist[[c]]
		}
	}
}

dat.corr.s <- data.frame(arr.corr.s)

dat.corr.s$quality <- factor(dat.corr.s$quality)
dat.corr.s$sample <- factor(dat.corr.s$sample)
dat.corr.s$effsize <- factor(dat.corr.s$effsize)

arr.corr.m <- matrix(NA,nrow=21*12*10,ncol=10)

colnames(arr.corr.m) <- c("p.0005","p.005","p.025","p.500","p.975","p.995","p.9995","quality","sample","effsize")

for (q in 1:21)
{
	for (s in 1:12)
	{
		for (c in 1:10)
		{
		arr.corr.m[(q-1)*120+(s-1)*10+c,1:7] <- quantile(arr.cor.m[,q,s,c],siglim,na.rm=T)
		arr.corr.m[(q-1)*120+(s-1)*10+c,"quality"] <- rellevels[[q]]
		arr.corr.m[(q-1)*120+(s-1)*10+c,"sample"] <- sampsize[[s]]
		arr.corr.m[(q-1)*120+(s-1)*10+c,"effsize"] <- corlist[[c]]
		}
	}
}

dat.corr.m <- data.frame(arr.corr.m)

dat.corr.m$quality <- factor(dat.corr.m$quality)
dat.corr.m$sample <- factor(dat.corr.m$sample)
dat.corr.m$effsize <- factor(dat.corr.m$effsize)

dat.corr.m$type <- "marginal"
dat.corr.s$type <- "scale"

dat.corr.m$cquality <- (as.numeric(dat.corr.m$quality)-1)/20
dat.corr.s$cquality <- (as.numeric(dat.corr.s$quality)-1)/20

dat.corr.m$ceffsize <- (as.numeric(as.character(dat.corr.m$effsize)))
dat.corr.s$ceffsize <- (as.numeric(as.character(dat.corr.s$effsize)))

dat.corr.m$csampsize <- (as.numeric(as.character(dat.corr.m$sample)))
dat.corr.s$csampsize <- (as.numeric(as.character(dat.corr.s$sample)))

dat.corr$ctype <- Recode(dat.corr$type,"'marginal'='MD-prior';'scale'='ED-prior'")

dat.corr <- rbind(dat.corr.s,dat.corr.m)

dat.corr$sig.05 <- as.numeric(as.character(Recode(dat.corr$sample,"10=.632;20=.444;30=.361;40=.312;50=.279;100=.197;200=.139;300=.113;400=.098;500=.088;700=.074;1000=.062")))
dat.corr$sig.01 <- as.numeric(as.character(Recode(dat.corr$sample,"10=.765;20=.561;30=.463;40=.403;50=.361;100=.256;200=.182;300=.149;400=.129;500=.115;700=.097;1000=.081")))
# dat.corr$sig.001








#######################################
#######################################
##### Equation-Fitting procedure

dat.corr.m <- subset(dat.corr,type=="marginal")

dat.corr.m$ci.975 <- dat.corr.m$p.975-dat.corr.m$p.500
dat.corr.m$ci.025 <- dat.corr.m$p.025-dat.corr.m$p.500

fitmodel.m.500 <- nls(p.500~I(1)*ceffsize / (1 + exp(-(b+d*ceffsize) * (cquality-c))),start=list(b=1,c=1,d=1),data=dat.corr.m,control=list(maxiter = 500))
dat.corr.m$pre.500 <- predict(fitmodel.m.500)

fitmodel.m.025 <- nls(ci.025~(b+1-cquality)*(1/(sqrt(csampsize)))+(d*(cquality^2)/csampsize)+c*cquality^2*csampsize,start=list(b=1,c=1,d=1),data=dat.corr.m,control=list(maxiter = 500))
fitmodel.m.975 <- nls(ci.975~(b+1-cquality)*(1/(sqrt(csampsize)))+(d*(cquality^2)/csampsize)+c*cquality^2*csampsize,start=list(b=1,c=1,d=1),data=dat.corr.m,control=list(maxiter = 500))
fitmodel.m.025 <- nls(ci.025~(b+1-cquality)*(1/(sqrt(csampsize)))+(d*(cquality^2)/csampsize)+c*cquality^2*csampsize,start=list(b=1,c=1,d=1),data=dat.corr.m,control=list(maxiter = 500))
fitmodel.m.975 <- nls(ci.975~(b+1-cquality)*(1/(sqrt(csampsize)))+(d*(cquality^2)/csampsize)+c*cquality^2*csampsize,start=list(b=1,c=1,d=1),data=dat.corr.m,control=list(maxiter = 500))
fitmodel.m.025 <- nls(ci.025~(b+1-cquality)*(a/(sqrt(csampsize)))+(d*(cquality^2)/csampsize)+c*cquality^2*csampsize,start=list(a=-1,b=1,c=1,d=1),data=dat.corr.m,control=list(maxiter = 500))
fitmodel.m.975 <- nls(ci.975~(b+1-cquality)*(a/(sqrt(csampsize)))+(d*(cquality^2)/csampsize)+c*cquality^2*csampsize,start=list(a=-1,b=1,c=1,d=1),data=dat.corr.m,control=list(maxiter = 500))

fitmodel.m.025x <- nls(ci.025~(b+1-cquality)*(a/(sqrt(csampsize)))+(d*(cquality^2)/csampsize),start=list(a=-1,b=1,d=1),data=dat.corr.m,control=list(maxiter = 500))
fitmodel.m.975x <- nls(ci.975~(b+1-cquality)*(a/(sqrt(csampsize)))+(d*(cquality^2)/csampsize),start=list(a=-1,b=1,d=1),data=dat.corr.m,control=list(maxiter = 500))

dat.corr.s$ci.975 <- dat.corr.s$p.975-dat.corr.s$p.500
dat.corr.s$ci.025 <- dat.corr.s$p.025-dat.corr.s$p.500

fitmodel.s.500 <- nls(p.500~I(1)*ceffsize / (1 + exp(-(b+d*ceffsize) * (cquality-c))),start=list(b=1,c=1,d=1),data=dat.corr.s,control=list(maxiter = 500))
dat.corr.s$pre.500 <- predict(fitmodel.s.500)

fitmodel.s.025 <- nls(ci.025~(b+1-cquality)*(1/(sqrt(csampsize)))+(d*(cquality^2)/csampsize)+c*cquality^2*csampsize,start=list(b=1,c=1,d=1),data=dat.corr.s,control=list(maxiter = 500))
fitmodel.s.975 <- nls(ci.975~(b+1-cquality)*(1/(sqrt(csampsize)))+(d*(cquality^2)/csampsize)+c*cquality^2*csampsize,start=list(b=1,c=1,d=1),data=dat.corr.s,control=list(maxiter = 500))
fitmodel.s.025 <- nls(ci.025~(b+1-cquality)*(1/(sqrt(csampsize)))+(d*(cquality^2)/csampsize)+c*cquality^2*csampsize,start=list(b=1,c=1,d=1),data=dat.corr.s,control=list(maxiter = 500))
fitmodel.s.975 <- nls(ci.975~(b+1-cquality)*(1/(sqrt(csampsize)))+(d*(cquality^2)/csampsize)+c*cquality^2*csampsize,start=list(b=1,c=1,d=1),data=dat.corr.s,control=list(maxiter = 500))
fitmodel.s.025 <- nls(ci.025~(b+1-cquality)*(a/(sqrt(csampsize)))+(d*(cquality^2)/csampsize)+c*cquality^2*csampsize,start=list(a=-1,b=1,c=1,d=1),data=dat.corr.s,control=list(maxiter = 500))
fitmodel.s.975 <- nls(ci.975~(b+1-cquality)*(a/(sqrt(csampsize)))+(d*(cquality^2)/csampsize)+c*cquality^2*csampsize,start=list(a=-1,b=1,c=1,d=1),data=dat.corr.s,control=list(maxiter = 500))

dat.corr.s$ci.pre.025 <- predict(fitmodel.s.025)
dat.corr.s$ci.pre.975 <- predict(fitmodel.s.975)
dat.corr.s$pre.025 <- dat.corr.s$pre.500+dat.corr.s$ci.pre.025
dat.corr.s$pre.975 <- dat.corr.s$pre.500+dat.corr.s$ci.pre.975
dat.corr.s$pre.025 <- replace(dat.corr.s$pre.025,dat.corr.s$pre.025<(-1),(-1))
dat.corr.s$pre.975 <- replace(dat.corr.s$pre.975,dat.corr.s$pre.975>(1),(1))


fitmodel.s.025x <- nls(ci.025~(b+1-cquality)*(a/(sqrt(csampsize)))+(d*(cquality^2)/csampsize),start=list(a=-1,b=1,d=1),data=dat.corr.s,control=list(maxiter = 500))
fitmodel.s.975x <- nls(ci.975~(b+1-cquality)*(a/(sqrt(csampsize)))+(d*(cquality^2)/csampsize),start=list(a=-1,b=1,d=1),data=dat.corr.s,control=list(maxiter = 500))

dat.corr.m$ci.pre.025 <- predict(fitmodel.m.025)
dat.corr.m$ci.pre.975 <- predict(fitmodel.m.975)
dat.corr.m$pre.025 <- dat.corr.m$pre.500+dat.corr.m$ci.pre.025
dat.corr.m$pre.975 <- dat.corr.m$pre.500+dat.corr.m$ci.pre.975
dat.corr.m$pre.025 <- replace(dat.corr.m$pre.025,dat.corr.m$pre.025<(-1),(-1))
dat.corr.m$pre.975 <- replace(dat.corr.m$pre.975,dat.corr.m$pre.975>(1),(1))

dat.table.m <- expand.grid(seq(0,1,0.05),seq(0,1,0.05),c(seq(10,90,10),seq(100,500,50),seq(600,1900,100),seq(2000,5000,500),seq(6000,100000,1000)))

dat.table.e <- expand.grid(seq(0,1,0.05),seq(0,1,0.05),c(seq(10,90,10),seq(100,500,50),seq(600,1900,100),seq(2000,5000,500),seq(6000,100000,1000)))

names(dat.table.m) <- c("ceffsize","cquality","csampsize")

names(dat.table.e) <- c("ceffsize","cquality","csampsize")

dat.table.m$median <- predict(fitmodel.m.500,newdata=dat.table.m)
dat.table.m$lower <- predict(fitmodel.m.025x,newdata=dat.table.m)
dat.table.m$upper <- predict(fitmodel.m.975x,newdata=dat.table.m)
dat.table.m$LL <- dat.table.m$median+dat.table.m$lower
dat.table.m$UL <- dat.table.m$median+dat.table.m$upper
dat.table.m$crit.T <- qt(0.975,dat.table.m$csampsize-2)
dat.table.m$crit.R <- with(dat.table.m,crit.T/(csampsize-2+(ceffsize^2))^(0.5))
dat.table.m$power.ok <- dat.table.m$LL>dat.table.m$crit.R
dat.table.m$power.surplus <- dat.table.m$LL-dat.table.m$crit.R

dat.table.e$median <- predict(fitmodel.s.500,newdata=dat.table.e)
dat.table.e$lower <- predict(fitmodel.s.025x,newdata=dat.table.e)
dat.table.e$upper <- predict(fitmodel.s.975x,newdata=dat.table.e)
dat.table.e$LL <- dat.table.e$median+dat.table.e$lower
dat.table.e$UL <- dat.table.e$median+dat.table.e$upper
dat.table.e$crit.T <- qt(0.975,dat.table.e$csampsize-2)
dat.table.e$crit.R <- with(dat.table.e,crit.T/(csampsize-2+(ceffsize^2))^(0.5))
dat.table.e$power.ok <- dat.table.e$LL>dat.table.e$crit.R
dat.table.e$power.surplus <- dat.table.e$LL-dat.table.e$crit.R

subset(dat.table.m,ceffsize==0.5&cquality==0.5)[,c("csampsize","power.surplus")]
subset(dat.table.m,ceffsize==0.0&cquality==0.0)[,c("csampsize","power.surplus")]

sampsizelist <- c(seq(10,90,10),seq(100,500,50),seq(600,1900,100),seq(2000,5000,500),seq(6000,100000,1000))
effsizelist <- seq(0,1,0.05)
qualitylist <- seq(0,1,0.05)

eq2 <- matrix(NA,ncol=length(effsizelist),nrow=length(qualitylist))

for (i in 1:length(effsizelist))
{
for (j in 1:length(qualitylist))
{
eq2[j,i] <- sampsizelist[min(which(subset(dat.table.e,ceffsize==effsizelist[i]&cquality==qualitylist[j])[,c("power.surplus")]>0))]
}
}

eq2.df <- data.frame(min.sample.size=c(eq2),quality=rep(seq(0,1,.05),times=21),effect.size=rep(seq(0,1,.05),each=21))

th2.df <- data.frame(	author=c(rep("Landis & Koch / Altman",times=6),rep("Fleiss, Levin & Paik",times=4),c(rep("Krippendorff",times=4))),
				y=c(rep(-.05,times=6),rep(-.025,times=4),rep(.00,times=4)),
				threshold=c(0.0,0.2,0.4,0.6,0.8,1.0,0.0,0.4,0.75,1.0,0.0,0.667,0.800,1.000),
				clr=(c("#e41a1c","#ff7f00","#ffff33","#377eb8","#4daf4a","#4daf4a","#e41a1c","#ffff33","#377eb8","#4daf4a","#e41a1c","#ffff33","#377eb8","#4daf4a")))



eq <- matrix(NA,ncol=length(effsizelist),nrow=length(qualitylist))

for (i in 1:length(effsizelist))
{
for (j in 1:length(qualitylist))
{
eq[j,i] <- sampsizelist[min(which(subset(dat.table.m,ceffsize==effsizelist[i]&cquality==qualitylist[j])[,c("power.surplus")]>0))]
}
}

eq.df <- data.frame(min.sample.size=c(eq),quality=rep(seq(0,1,.05),times=21),effect.size=rep(seq(0,1,.05),each=21))

th.df <- data.frame(	author=c(rep("Landis & Koch / Altman",times=6),rep("Fleiss, Levin & Paik",times=4),c(rep("Krippendorff",times=4))),
				y=c(rep(-.05,times=6),rep(-.025,times=4),rep(.00,times=4)),
				threshold=c(0.0,0.2,0.4,0.6,0.8,1.0,0.0,0.4,0.75,1.0,0.0,0.667,0.800,1.000),
				clr=(c("#e41a1c","#ff7f00","#ffff33","#377eb8","#4daf4a","#4daf4a","#e41a1c","#ffff33","#377eb8","#4daf4a","#e41a1c","#ffff33","#377eb8","#4daf4a")))

eq.df$diff.sample.size <- eq2.df$min.sample.size-eq.df$min.sample.size
eq.df$diff.sample.size <- ifelse(is.na(eq2.df$min.sample.size)&!is.na(eq.df$min.sample.size),4*eq.df$min.sample.size,eq.df$diff.sample.size)

colnames(eq) <- effsizelist
rownames(eq) <- qualitylist

dat.corr.s$sig.05 <- as.numeric(as.character(Recode(dat.corr.s$sample,"10=.632;20=.444;30=.361;40=.312;50=.279;100=.197;200=.139;300=.113;400=.098;500=.088;700=.074;1000=.062")))
dat.corr.s$sig.01 <- as.numeric(as.character(Recode(dat.corr.s$sample,"10=.765;20=.561;30=.463;40=.403;50=.361;100=.256;200=.182;300=.149;400=.129;500=.115;700=.097;1000=.081")))

dat.corr <- rbind(dat.corr.s,dat.corr.m)

dat.corr$issig05 <- ifelse(dat.corr$p.025-dat.corr$sig.05>0,dat.corr$p.500,NA)


dat.corr$delta.upper <- dat.corr$p.975-dat.corr$sig.05



save(fitmodel.s.025x,file="eq_ci_025_ED.RData")
save(fitmodel.s.500,file="eq_median_ED.RData")
save(fitmodel.s.975x,file="eq_ci_975_ED.RData")
save(fitmodel.m.025x,file="eq_ci_025_MD.RData")
save(fitmodel.m.500,file="eq_median_MD.RData")
save(fitmodel.m.975x,file="eq_ci_025_MD.RData")



###########################
###########################
### THE POWERCALC FUNCTION

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

PowerCalc(E=0.25,Q=0.667,N=500,alpha=.05) # Desired alpha level: p<.05; Coding accuracy=0.667; sample size=500; effect size=R=0.25

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


