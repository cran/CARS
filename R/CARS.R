#' The CARS procedure for controlling the false discovery rate
#'
#' @docType package
#' @name CARS
library(np);

#' CARS
#'
#' This function runs the CARS procedure, constructing the auxiliary variables,
#' computing the test statistics, choosing the cutoff and selecting the locations.
#'
#' @param X the first matrix or data frame of observation
#' @param Y the second matrix or date frame of observation
#' @param alpha targeted FDR (false discovery rate) level
#' @param tau the threshold for choosing interesting locations for density estimation, default is 0.5
#' @param variance for X and Y, default is NULL. If provided, in the form of a m*2 matrix, the columns are representing x and y's variance for each location
#'
#' @return A list containing the following components: 
#' \item{de}{decision for each location (0 or 1)}
#' \item{cars}{estimated CARS statistics}
#' \item{th}{threshold for CARS procedure}
#' 
#' @examples
#' X <- matrix(rnorm(1000),ncol=5,nrow=200);
#' Y <- matrix(rep(c(0,3),c(800,200))+rnorm(1000),ncol=5,nrow=200);
#' CARS(X,Y,0.05,tau=0.9);
#'
#' @importFrom np npudensbw
#' @importFrom stats density dt pnorm pt qnorm var rnorm
#' @export 
CARS <- function(X,Y,alpha,tau=0.9,variance){
	#Validate input types.
	if (is.data.frame(X)){
		X.names=names(X)
		X = as.matrix(X,rownames.force=F)
	} else if (is.matrix(X))
		X.names=colnames(X)
	else
		stop('Input X must be a matrix or data frame')

	if (is.data.frame(Y)){
		Y.names=names(Y)
		Y = as.matrix(Y,rownames.force=F)
	} else if (is.matrix(Y))
		Y.names=colnames(Y)
	else
		stop('Input Y must be a matrix or data frame')

		#validate input dimensions.
		n_x <- ncol(X); m <- nrow(X);
		stopifnot(nrow(Y)==m)
		n_y <- ncol(Y);
		n <- n_x+n_y;

		#Calculate the mean of X, Y for each location
		X.mean <- rowMeans(X);
		Y.mean <- rowMeans(Y);

	if (missing(variance)){
		# Calculate the estimated variance
		X.var <- apply(X,1,var);
		Y.var <- apply(Y,1,var);
		
		#Calculate the pooled variance and kappa
   		pool.sd <- sqrt(n_x/n*Y.var+n_y/n*X.var);
   		kappa <- n_y*X.var/(n_x*Y.var);

   		#Calculate primary and auxiliary statistics
   		t_1 <- (X.mean-Y.mean)/pool.sd*sqrt(n_x*n_y/n);
   		t_2 <- (X.mean+kappa*Y.mean)/(sqrt(kappa)*pool.sd)*sqrt(n_x*n_y/n);

   		#Estimate the non-null proportions of primary and auxiliary statistics
   		deg <- (X.var/n_x+Y.var/n_y)^2/((X.var/n_x)^2/(n_x-1)+(Y.var/n_y)^2/(n_y-1));
   		t_1.p.Est <- epsest.func(qnorm(pt(t_1,deg)),0,1);
   		t_2.p.Est <- epsest.func(t_2,0,1);
}


	else {
		#Assign corresponding variances for each location
		X.var <- variance[,1];
		Y.var <- variance[,2];

		deg <- (X.var/n_x+Y.var/n_y)^2/((X.var/n_x)^2/(n_x-1)+(Y.var/n_y)^2/(n_y-1));

		#Calculate kappa and pooled variance
		kappa <- n_y*X.var/(n_x*Y.var);
		pool.sd <- sqrt(n_y/n*X.var+n_x/n*Y.var);

		#Calculate primary and auxiliary statistics
   		t_1 <- (X.mean-Y.mean)/pool.sd*sqrt(n_x*n_y/n);
   		t_2 <- (X.mean+kappa*Y.mean)/(sqrt(kappa)*pool.sd)*sqrt(n_x*n_y/n);

   		#Estimate the bandwidths
   		bandwidth <- npudensbw(~t_1+t_2,bwmethod="normal-reference")$bw;
   		hx <- bandwidth[1];
   		ht <- bandwidth[2];

   		#Estimate the non-null proportions of primary and auxiliary statistics
   		t_1.p.Est <- epsest.func(t_1,0,1);
   		t_2.p.Est <- epsest.func(t_2,0,1);
}
	   		#Estimate the lfdrs
   		t_1.density.Est <- density(t_1,from=min(t_1)-10,to=max(t_1)+10,n=1000);
	   	t_1.density.Est <- lin.itp(t_1,t_1.density.Est$x,t_1.density.Est$y);
	   	t_1.Lfdr.Est <- (1-t_1.p.Est)*dt(t_1,n_x-1)/t_1.density.Est;
	   	t_1.Lfdr.Est[which(t_1.Lfdr.Est>1)] <- 1;

	   	t_2.density.Est <- density(t_2,from=min(t_2)-10,to=max(t_2)+10,n=1000);
	   	t_2.density.Est <- lin.itp(t_2,t_2.density.Est$x,t_2.density.Est$y);
	   	t_2.Lfdr.Est <- (1-t_2.p.Est)*dt(t_2,n_y-1)/t_2.density.Est;
	   	t_2.Lfdr.Est[which(t_2.Lfdr.Est>1)] <- 1;

	   	S <- which(t_1.Lfdr.Est<=0.8);
	   	bandwidth <- np::npudensbw(~t_1[S]+t_2[S],bwmethod="normal-reference")$bw;
	   	hx <- bandwidth[1];
	   	ht <- bandwidth[2];

	   	#Calculate the estimated CARS statistics
	   	cars.denominator <- np::npudens(~t_1+t_2,bws=bandwidth)$dens;

      #Estimate Correction 
      sample_null <- rnorm(10000);
      t_1.den.Est <- density(t_1,bw=hx,from=min(t_1)-10,to=max(t_1)+10,n=1000);
      sample_lfdr <- (1-t_1.p.Est)*dt(sample_null,deg)/lin.itp(sample_null,t_1.den.Est$x,t_1.den.Est$y);
      correction <- length(which(sample_lfdr>=tau))/10000;  

      t_1.pval <- 2*pnorm(-abs(t_1));
	   	T.tau <- which(t_1.Lfdr.Est>=tau);

	   	t_2.star.den.Est <- density(t_2[T.tau],bw=ht,from=min(t_2[T.tau])-10,to=max(t_2[T.tau])+10,n=1000);
	   	t_2.star.den.Est <- lin.itp(t_2,t_2.star.den.Est$x,t_2.star.den.Est$y);

	   	cars.numerator <- dt(t_1,deg)*length(T.tau)/m*t_2.star.den.Est/correction;

	   	cars.Est <- cars.numerator/cars.denominator;
      cars.Est[which(cars.Est>=1)] <- 1;



   	cars.sorted <- sort(cars.Est,decreasing=FALSE,index.return=TRUE);
   	cars.sorted.cumsum <- cumsum(cars.sorted$x);

   	decision <- rep(0,m);
   	threshold <- 0;
   	for (i in 1:m){
   		if(cars.sorted.cumsum[i]/i <= alpha & cars.sorted.cumsum[i+1]/(i+1) > alpha){
   			decision[cars.sorted$ix[1:i]] <- 1;
   			threshold <- cars.Est[cars.sorted$ix[i]];
   			break;
   		}
   	}

   	y <- list(de=decision,cars=cars.Est,th=threshold);
   	return(y);

}



#' Estimation of the non-null proportion
#'
#' Estimates the proportion of non-nulls.
#' 
#' @param x the corresponding vector to be estimated
#' @param u the mean of the null distribution
#' @param sigma the standard deviation of the null distribution
#' @return a value indicating the estimated non-null proportion
#'
#' @examples
#' X <- rep(c(0,2),c(800,200))+rnorm(1000);
#' epsest.func(X,0,1);
#'
#' @export
epsest.func <- function(x,u,sigma)
{
  z  = (x - u)/sigma
  xi = c(0:100)/100
  tmax=sqrt(log(length(x)))
  tt=seq(0,tmax,0.5)

  epsest=NULL

  for (j in 1:length(tt)) { 

    t=tt[j]
    f  = t*xi
    f  = exp(f^2/2)
    w  = (1 - abs(xi))
    co  = 0*xi

    for (i in 1:101) {
      co[i] = mean(cos(t*xi[i]*z));
    } 
    epshat = 1 - sum(w*f*co)/sum(w)
    epsest=c(epsest,epshat)
  }
  return(epsest=max(epsest))
}

#' Linear interpolation
#'
#' Interpolates desired vector given density estimation.
#' 
#' @param x the coordinates of points where the density needs to be interpolated
#' @param X the coordinates of the estimated densities
#' @param Y the values of the estimated densities
#' @return the interpolated densities
#' 
#' @examples
#' X <- seq(-10,10,length.out=20);
#' Y <- dnorm(X);
#' x <- seq(-10,10,length.out=100);
#' lin.itp(x,X,Y)
#'
#' @export
lin.itp<-function(x, X, Y){

  x.N<-length(x)
  X.N<-length(X)
  y<-rep(0, x.N)
  for (k in 1:x.N){
    i<-max(which((x[k]-X)>=0))
    if (i<X.N)
      y[k]<-Y[i]+(Y[i+1]-Y[i])/(X[i+1]-X[i])*(x[k]-X[i])
    else 
      y[k]<-Y[i]
  }
  return(y)
}









