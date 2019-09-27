#' Old version f24 function using fft
#'
#' This function tests if the data is 24h rhythmic
#'
#' @param x the data of time series
#' @param offset offset of time points
#' @param nf24 nb of replicates
#' @author Felix Naef
#'
f24=function(x, offset=0, nf24=2){
	#this works only when the length of x is even
	# x=2^x
	f=fft(as.double(x))
	phase=12/pi*atan2(-Im(f[nf24]), Re(f[nf24]))
	# cat(phase, "\n")
	if(phase<0) phase=phase+24
	if(phase>24) phase=phase-24
	phase=(phase+offset)%%24

	f2=abs(f)^2
	z=f2[nf24]/sum(f2[2:(length(x)/2)])
	if(is.nan(z)) z=0
	lu=length(unique(x))
	f3 = abs(f)
	#amp = f3[1]/length(x)
	rel.amp = f3[nf24]/length(x)*2/(f3[1]/length(x))

	#amp=max(x)-min(x)
	# if (amp<0.5 | lu != 6) z=NA
	# if (lu != 6) z=NA

	# pv=pbeta(z, 1, length(x)/2-2, lower.tail = FALSE, log.p = FALSE)

	c(lu=lu, mean=mean(x), amp=max(x)-min(x), rel.amp=rel.amp, phase=phase, pvalue=(1-z)^(length(x)/2-2))
}

f24_R2_lm=function(x, t=2*(0:(length(x)-1)), period=24, offset=0)
{
	n=length(x)
	# mu=mean(x)
	sig2=var(x)

	#
	c=cos(2*pi*t/period)
	s=sin(2*pi*t/period)

	# x=x-mean(x)
	# a=mean(x*(c-mean(c))
	# b=mean(x*(s-mean(s))


	# a=mean(x*c)
	# b=mean(x*s)
	# x.hat=mu+2*a*c+2*b*s
	# sig2.1=var(x-x.hat)

	fit = lm(x~c+s)

	a=fit$coef[2]
	b=fit$coef[3]

	R2=0
	if(sig2>0) R2 =1.0-sum(fit$residuals^2)/(n-1)/sig2

	# R2=0
	# if(sig2>0) R2=1-sig2.1/sig2
	# http://www.combustion-modeling.com/downloads/beta-distribution-for-testing-r-squared.pdf
	# I checked that it works
	p=3

#	if(n<3)
#	{
#
#	}
#	else
#	{
	pv=pbeta(R2, (p-1)/2, (n-p)/2, lower.tail = FALSE, log.p = FALSE)

	lu=length(unique(x))
	amp=max(x)-min(x)
	phase=period/(2*pi)*atan2(b, a)
	if(phase<0) phase=phase+period
	if(phase>period) phase=phase-period
	phase=(phase+offset)%%period
	#x=2^x
	c(lu=lu, mean=mean(x), amp=2*sqrt(a^2+b^2),relamp=sqrt(a^2+b^2)/(fit$coef[1]),phase=phase, pval=pv)
#	}
}

#An alternative function of f24_R2
f24_R2_alt=function(x, t=2*(0:(length(x)-1)), period=24, offset=0)
{
	n=length(x)
	#mu=mean(x)
	lu=length(unique(x))
	if(n<4)
	{
		if(n==0) c(lu=lu, mean=NA, amp=NA, relamp=NA,phase=NA,pval=NA)
		else
		{
			c(lu=lu, mean=mean(x), amp=NA, relamp=NA,phase=NA,pval=NA)
		}
	}
	else
	{
		sig2=var(x)
		c=cos(2*pi*t/period)
		s=sin(2*pi*t/period)
		A = mean(x*c)-mean(x)*mean(c)
		B = mean(x*s)-mean(x)*mean(s)
		c1 = mean(c^2)-mean(c)^2
		c2 = mean(c*s)-mean(c)*mean(s)
		c3 = mean(s^2)-mean(s)^2
		b = (A*c2-B*c1)/(c2^2-c1*c3)
		a = (A-b*c2)/c1
		mu = mean(x)-a*mean(c)-b*mean(s)
		#	b=2*mean(x*s)
		x.hat=mu+a*c+b*s
		sig2.1=var(x-x.hat)
		if(is.na(a)||is.na(b)) {c(lu=lu, mean=mean(x), amp=NA, relamp=NA,phase=NA,pval=NA)}
		else
		{
			p=3
			R2=0
			if(sig2>0) R2=1-sig2.1/sig2
			# http://www.combustion-modeling.com/downloads/beta-distribution-for-testing-r-squared.pdf
			# I checked that it works

			k = 1
			amp=max(x)-min(x)
			phase=period/(2*pi)*atan2(b, a)
			if(phase<0) phase=phase+period
			if(phase>period) phase=phase-period
			phase=(phase+offset)%%period
			if(n>p)
			{
				pval = pbeta(R2, (p-1)/2, (n-p)/2, lower.tail = FALSE, log.p = FALSE)

			}
			else
			{
				pval = NA

			}
			#print("here")
			c(lu=lu, mean =mean(x), amp=2*sqrt(a^2+b^2),relamp=sqrt(a^2+b^2)/(mu),phase=phase, pval=pval)
		}
	}
}

#' latest version f24 function using lm
#'
#' This function tests if the data is 24h rhythmic
#'
#' @param x the data of time series
#' @param t the corresponding time points
#' @param period the period to test
#' @param offset the offset of time points (not used)
#' @return A vector of nb of time points, mean of data, amplitude and relative amplitude of fitted result,  phase and pvalue
#' @author Jingkui Wang
#' @details
#' the function can deal with missing time points
#' if there are less 4 time points, no amplitude, relative amplitude, phase and pvalue were calculated
#' @seealso \code{lm}\
#'
#' @export
#' @import stats
f24_R2_ls=function(x, t=2*(0:(length(x)-1)), period=24, offset=0)
{
	kk = which(!is.na(x)==TRUE)
	x = x[kk]
	t = t[kk]
	n=length(x)
	#mu=mean(x)
	nb.timepoints=length(x)
	if(n<4)
	{
		if(n==0) c(nb.timepoints=nb.timepoints, mean=NA, amp=NA, relamp=NA,phase=NA,pval=NA)
		else
		{
			c(nb.timepoints=nb.timepoints, mean=mean(x), amp=NA, relamp=NA,phase=NA,pval=NA)
		}
	}
	else
	{
		sig2=var(x)
		c=cos(2*pi*t/period)
		s=sin(2*pi*t/period)
		A = mean(x*c)-mean(x)*mean(c)
		B = mean(x*s)-mean(x)*mean(s)
		c1 = mean(c^2)-mean(c)^2
		c2 = mean(c*s)-mean(c)*mean(s)
		c3 = mean(s^2)-mean(s)^2
		b = (A*c2-B*c1)/(c2^2-c1*c3)
		a = (A-b*c2)/c1
		mu = mean(x)-a*mean(c)-b*mean(s)
		#	b=2*mean(x*s)
		x.hat=mu+a*c+b*s
		sig2.1=var(x-x.hat)
		if(is.na(a)||is.na(b)) {c(nb.timepoints=nb.timepoints, mean=mean(x), amp=NA, relamp=NA,phase=NA,pval=NA)}
		else
		{
			p=3
			R2=0
			if(sig2>0) R2=1-sig2.1/sig2
			# http://www.combustion-modeling.com/downloads/beta-distribution-for-testing-r-squared.pdf
			# I (Felix) checked that it works
			amp=max(x)-min(x)
			phase=period/(2*pi)*atan2(b, a)
			if(phase<0) phase=phase+period
			if(phase>period) phase=phase-period
			phase=(phase+offset)%%period
			pval = pbeta(R2, (p-1)/2, (n-p)/2, lower.tail = FALSE, log.p = FALSE)

			c(nb.timepoints=nb.timepoints, mean =mean(x), amp=2*sqrt(a^2+b^2),relamp=sqrt(a^2+b^2)/(mu),phase=phase, pval=pval)
		}
	}
}


f24_R2_log=function(x, t=2*(0:(length(x)-1)), period=24, offset=0)
{
	n=length(x)
	mu=mean(x)
	sig2=var(x)

	c=cos(2*pi*t/period)
	s=sin(2*pi*t/period)
	a=2*mean(x*c)
	b=2*mean(x*s)
	x.hat=mu+a*c+b*s
	sig2.1=var(x-x.hat)

	p=3
	R2=0
	if(sig2>0) R2=1-sig2.1/sig2
	# http://www.combustion-modeling.com/downloads/beta-distribution-for-testing-r-squared.pdf
	# I checked that it works
	p=pbeta(R2, (p-1)/2, (n-p)/2, lower.tail = FALSE, log.p = FALSE)

	lu=length(unique(x))
	amp=max(x)-min(x)
	phase=period/(2*pi)*atan2(b, a)
	if(phase<0) phase=phase+period
	if(phase>period) phase=phase-period
	phase=(phase+offset)%%period
	#x=2^x
	c(lu=lu, mean =mean(x), amp=2*sqrt(a^2+b^2), relamp=(max(x)-min(x))/(2.0*mean(x)), phase=phase, pval=p)
}

##########
########## likelihood ratio test
##########
f24_ratio=function(x, t=2*(0:(length(x)-1)), period=24, offset=0)
{
  n=length(x)
  # mu=mean(x)
  sig2=var(x)
  c=cos(2*pi*t/period)
  s=sin(2*pi*t/period)

  fit = lm(x~c+s)

  a=fit$coef[2]
  b=fit$coef[3]

  R2=0
  if(sig2>0) R2 =1.0-sum(fit$residuals^2)/(n-1)/sig2
  # http://www.combustion-modeling.com/downloads/beta-distribution-for-testing-r-squared.pdf
  # I checked that it works
  p=3
  pv=pbeta(R2, (p-1)/2, (n-p)/2, lower.tail = FALSE, log.p = FALSE)
  lu=length(unique(x))
  amp=max(x)-min(x)
  phase=period/(2*pi)*atan2(b, a)
  if(phase<0) phase=phase+period
  if(phase>period) phase=phase-period
  phase=(phase+offset)%%period

  ############### compare two models
  logf1 = -sum((x-mean(x))^2)
  logf2 = -sum(fit$residuals^2)
  D = -2*logf1+2*logf2
  p_ratio = dchisq(D, df = 2)


  c(lu=lu, mean=mean(x), amp=2*sqrt(a^2+b^2),relamp=sqrt(a^2+b^2)/(fit$coef[1]),phase=phase, pval=pv,p.ratio=p_ratio)
  #	}
}

#######################################################################
####################################################################### some homemade statistic tests
#######################################################################
#' function to compute q value
#'
#' This function implement the Benjamin-Hochberg (BH) method for q value calculation
#'
#' @param pv a vector of p values
#' @return A vector of q values
#' @author Laura Symul
#' @details
#' the function can deal with missing time points
#'
#' @export
qvals=function(pv)
{
  PV=pv
  if(any(is.na(PV))) pv=pv[!is.na(PV)]
  # print(any(is.na(pv)))

  r=rank(-pv)
  o=order(-pv)
  lpv=length(pv)
  q=double(lpv)

  pv.s=pv[o]

  for (n in 1:lpv)
  {
    q[n] = pv.s[n] * lpv / (lpv-n+1)
    if(n==1) qmin=q[1]
    if(n>1 & q[n]>qmin) {q[n]=qmin}
    qmin=q[n]
  }
  Q=rep(NA,length(PV))
  if(any(is.na(PV))) Q[!is.na(PV)]=q[r]
  else Q=q[r]
  Q
}

#y1 and y2 are two vectors
#Kolmogorov-Smirnov D statistic
D.stat <- function(y1,y2) {
    x1 <- ecdf(y1)
    x2 <- ecdf(y2)
    z <- sort(c(y1,y2))
    D.all <- abs( x1(z) - x2(z) )
    return(max(D.all))
}

#Kuiper V statistic
V.stat <- function( y1,y2 ) {
    x1 <- ecdf(y1)
    x2 <- ecdf(y2)
    z <- sort(c(y1,y2))
    D.all <- x1(z) - x2(z)
    D.pos <- D.all[D.all > 0]
    D.neg <- D.all[D.all <= 0]
    a <- 0
    if (length(D.pos) > 0) {a <- max(D.pos) }
    b <- 0
    if (length(D.neg) > 0) {b <- max(abs(D.neg)) }
    V <- a + b
    return(V)
}

#Kuiper counterpart to ks.boot routine
kuiper.test <- function(y1,y2,nboots=10) {
#Calculate the Kuiper V-eff statistic
    x1 <- ecdf(y1)
    x2 <- ecdf(y2)
    z <- sort(c(y1,y2))
    n1 <- length(y1)
    n2 <- length(y2)
    n.eff <- n1*n2 / (n1+n2)
    V.ku <- V.stat(y1,y2)
    z <- V.ku * sqrt(n.eff)
    summ <- 0
    for (j in 1:100) { summ <- summ + (4*j^2*z^2 - 1) * exp(-2*j^2*z^2) }
    p.value <- 2*summ
    if( z < 0.4 ) {p.value <- 1}
#Bootstrap KS Test using Abadie algorithm
    y <- c(y1,y2)
    m <- nboots
    V.b <- numeric(m)
    for( i in 1:m ) {
        z <- sample(y,(n1+n2),replace=TRUE)
        z1 <- z[1:n1]
        z2 <- z[(n1+1):(n1+n2)]
        V.b[i] <- V.stat(z1,z2)
    }
    p.value.b <- length(V.b[V.b > V.ku])/m
    return(list(V=V.ku,V.eff=z,p.value=p.value,p.value.b=p.value.b,nboots=nboots))
}
############################################
############################################
Perm.test <- function(x, y, R=999, testfun=ks.test) {
	   z <- c(x, y)  # pooled sample
	   myfun <- function(a, b) suppressWarnings(unname(testfun(a, b)$statistic))
	   DoIt <- function() {
		     i <- sample(length(z), length(x))
		     myfun(z[i], z[-i])  # z[-i] is everything *except* the "i" elements of z
		   }
	   pstats <- replicate(R, DoIt())
	   stat <- myfun(x, y)
	   c("p-value" = mean(c(stat, pstats) >= stat))
}

####file://localhost/Users/jiwang/Dropbox/gachonprot/figure_plot_final/ex_MicroArray/JTK_Cycle/run_JTK_CYCLE%20(Example2).R
JTK_cycle=function(data.matrix, nb_timepoints, interval_timepoints,nb_replicats=1, periods=24)
{
	#source('lib/stats/big_stats.R', chdir=T)
	source("~/R_folder/JTK_CYCLE.R")
	options(stringsAsFactors=FALSE)
	data = data.matrix
	jtkdist(nb_timepoints, nb_replicats)  # total time points,  numbers of replicates per time point
	periods = periods/interval_timepoints  # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
	#periods = 6
	jtk.init(periods,interval_timepoints)  # 4 is the number of hours between time points

	#cat("JTK analysis started on",date(),"\n")
	#flush.console()

	#st <- system.time({
	res = apply(data,1,function(z){
				jtkx(z)
				c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
				})
	res <- as.data.frame(t(res))
	bhq <- p.adjust(unlist(res[,1]),"BH")
	res <- cbind(bhq,res)
	colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
	res
	#results <- cbind(annot,res,data)
	#results <- results[order(res$ADJ.P,-res$AMP),]

	#)
	#print(st)

	#save(results,file=paste("JTK",project,"rda",sep="."))
	#write.table(results,file=paste("JTK",project,"txt",sep="."),row.names=F,col.names=T,quote=F,sep="\t")
}
