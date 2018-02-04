rm(list=ls())

## This lab builds on the week 2 lab, but now includes the Aronow variance

### Code from Aronow, Green, Lee, Annals of Statistics, 2014

sharp.var <- function(yt,yc,N=length(c(yt,yc)),upper=TRUE) {
  m <- length(yt)
  n <- m + length(yc)
  FPvar <- function(x,N) (N-1)/(N*(length(x)-1))*sum((x - mean(x))^2)
  yt <- sort(yt)
  if(upper == TRUE) yc <- sort(yc) else
    yc <- sort(yc,decreasing=TRUE)
  p_i <- unique(sort(c(seq(0,n-m,1)/(n-m),seq(0,m,1)/m))) -.Machine$double.eps^.5
  p_i[1] <- .Machine$double.eps^.5
  yti <- yt[ceiling(p_i*m)]
  yci <- yc[ceiling(p_i*(n-m))]
  p_i_minus <- c(NA,p_i[1: (length(p_i)-1)])
  return(((N-m)/m * FPvar(yt,N) + (N-(n-m))/(n-m) * FPvar(yc,N) + 2*sum(((p_i-p_i_minus)*yti*yci)[2:length(p_i)]) - 2*mean(yt)*mean(yc))/(N-1))
}


## (1) Neyman-based inference; binary treatment and response

## Data
x <- c(rep(0,20), rep(1,20))
y <- c(0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,
       1,1,1,1,0,0,0,1,1,0,0,1,1,0,0,1,1,0,1,1)
n <- length(x)
k <- sum(x)  # size of treatment group

mean.ctrl <- mean(y[x==0]) 
mean.trt <- mean(y[x==1])
hat.ace <-  mean.trt - mean.ctrl


## Confidence interval:
var.ctrl <- var(y[x==0])
var.trt <- var(y[x==1]) 

hat.var.hat.ace <- (var.trt/k + var.ctrl/(n-k))

hat.var.hat.ace

##95% asymptotic confidence interval:

hat.ace - 1.96*sqrt(hat.var.hat.ace)
hat.ace + 1.96*sqrt(hat.var.hat.ace)

## Neyman's proposal:

## Note the formula below has n in the denominator, not (n-1) as in lecture note
## This is because the sample variance estimator var(.) in R does not include a finite pop correction factor
## of (n-1)/n.

hat.var.hat.ace.neyman <- (var.trt*(n-k)/k + var.ctrl*k/(n-k) + 2*(var.trt*var.ctrl)^0.5)/n
hat.var.hat.ace.neyman

hat.ace - 1.96*sqrt(hat.var.hat.ace.neyman)
hat.ace + 1.96*sqrt(hat.var.hat.ace.neyman)


## Aronow calculation

hat.var.hat.ace.aronow <- sharp.var(y[x==1],y[x==0])
hat.var.hat.ace.aronow

hat.ace - 1.96*sqrt(hat.var.hat.ace.aronow)
hat.ace + 1.96*sqrt(hat.var.hat.ace.aronow)

## Ratio of estimated variances
hat.var.hat.ace.aronow / hat.var.hat.ace   




#######################################################

# (2) Neyman-based inference; binary treatment and 
#     continuous response

x <- c(rep(0,20), rep(1,20))
y <- c(9.87, 12.14,9.62,8.63,11.40,7.40,6.88,10.00,7.39,11.31,
       11.56, 8.78, 3.70, 10.36, 10.68, 6.60, 11.11, 7.94, 10.97, 10.43,
       14.90, 17.06, 20.62, 18.18, 15.29, 12.38, 18.27, 16.00, 13.45, 16.29,
       14.93, 11.35, 17.02, 16.60, 13.99, 16.59, 10.76, 16.79, 13.42, 15.30)
n <- length(x)
k <- sum(x)

mean.ctrl <- mean(y[x==0])
mean.trt <- mean(y[x==1])
hat.ace <-  mean.trt - mean.ctrl


## Confidence interval:
var.ctrl <- var(y[x==0])
var.trt <- var(y[x==1])

hat.var.hat.ace <- var.trt/k + var.ctrl/(n-k)
hat.var.hat.ace

##95% asymptotic confidence interval:

hat.ace - 1.96*sqrt(hat.var.hat.ace)
hat.ace + 1.96*sqrt(hat.var.hat.ace)


## Neyman's proposal:

## Note the formula below has n in the denominator, not (n-1) as in lecture note
## This is because the sample variance estimator var(.) in R does not include a finite pop correction factor
## of (n-1)/n.

hat.var.hat.ace.neyman <- (var.trt*(n-k)/k + var.ctrl*k/(n-k) + 2*(var.trt*var.ctrl)^0.5)/n
hat.var.hat.ace.neyman

hat.ace - 1.96*sqrt(hat.var.hat.ace.neyman)
hat.ace + 1.96*sqrt(hat.var.hat.ace.neyman)

## Aronow calculation

hat.var.hat.ace.aronow <- sharp.var(y[x==1],y[x==0])
hat.var.hat.ace.aronow

hat.ace - 1.96*sqrt(hat.var.hat.ace.aronow)
hat.ace + 1.96*sqrt(hat.var.hat.ace.aronow)

## Ratio of estimated variances
hat.var.hat.ace.aronow / hat.var.hat.ace   


#######################################################

# (3) Neyman-based inference; binary treatment and 
#     continuous response
#     Investigating coverage properties and additivity
#     of naive, Neyman and Aronow


## (3a) Simulating potential outcomes under additivity
n <- 40  #finite population size
k <- 20  #treatment group size

ace.true <- 3
yctrl <- rnorm(n,10,2) # simulating Y(0) [from N(10,2)]
ytrt <- yctrl + ace.true  # additivity
mean(ytrt) - mean(yctrl) # true ACE
nsims <- 10000
covered <- rep(NA,nsims) # vector for coverage indicators
covered.aronow <- rep(NA,nsims)
hat.ace.vec <- rep(NA,nsims)
hat.var.hat.ace.vec <- rep(NA,nsims)
hat.var.hat.ace.vec.aronow <- rep(NA,nsims)

for(i in 1:nsims){
  trt.grp <- sample((1:n),k)  #choose our treatment group
  x <- rep(0,n)
  x[trt.grp] <- 1  #construct treatment vector
  y <- yctrl*(x==0) + ytrt*(x==1)  #build observed data vector
  ### Neyman procedure and Aronow
  mean.ctrl <- mean(y[x==0])
  mean.trt <- mean(y[x==1])
  hat.ace <-  mean.trt - mean.ctrl
  hat.ace.vec[i] <- hat.ace	# store the ACE estimate
  var.ctrl <- var(y[x==0])  # variance of control grp
  var.trt <- var(y[x==1])   # variance of trt grp
  hat.var.hat.ace <- var.trt/k + var.ctrl/(n-k)
  hat.var.hat.ace.aronow <- sharp.var(y[x==1],y[x==0])
  covered[i] <- (((hat.ace - 1.96*sqrt(hat.var.hat.ace)) < ace.true) &
                   ((hat.ace + 1.96*sqrt(hat.var.hat.ace)) > ace.true))
  covered.aronow[i] <- (((hat.ace - 1.96*sqrt(hat.var.hat.ace.aronow)) < ace.true) &
                          ((hat.ace + 1.96*sqrt(hat.var.hat.ace.aronow)) > ace.true))
  hat.var.hat.ace.vec[i] <- hat.var.hat.ace # store var(hat ACE) estimate
  hat.var.hat.ace.vec.aronow[i] <- hat.var.hat.ace.aronow # store Aronow var(hat ACE) estimate
}
sum(covered)/nsims  ## about 95%
sum(covered.aronow)/nsims 
mean(hat.ace.vec)  # are we unbiased ?
var(hat.ace.vec)   # true variability
mean(hat.var.hat.ace.vec) # are we unbiased?
mean(hat.var.hat.ace.vec.aronow) # what about Aronow's estimate?

par(mfrow=c(1,2))
plot(hat.ace.vec,hat.var.hat.ace.vec, main="hat(ACE) vs. Naive est of ACE",
     ylab="hat(var(hat(ACE)))", xlab="hat(ACE)" )

plot(hat.ace.vec,hat.var.hat.ace.vec.aronow, main="hat(ACE) vs. Aronow Var Est of ACE",
     ylab="hat(var(hat(ACE)))^Aro", xlab="hat(ACE)" )

#######
par(mfrow=c(1,2))

hist(hat.var.hat.ace.vec,main="Naive Estimate")
abline(v=(mean(hat.var.hat.ace.vec)),col="red") #mean of distribution of estimates hat.var(ace.hat)
abline(v=(var(hat.ace.vec)),col="blue")  #true variance var(ace.hat)

hist(hat.var.hat.ace.vec.aronow,main="Aronow Estimate")
abline(v=(mean(hat.var.hat.ace.vec.aronow)),col="red") #mean of distribution
abline(v=(var(hat.ace.vec)),col="blue")  #true variability


par(mfrow=c(1,1))

plot(hat.var.hat.ace.vec,hat.var.hat.ace.vec.aronow, main="Naive vs. Aronow",
     ylab="Aronow", xlab="Naive" )
abline(a=0,b=1,col="red")
abline(h=(var(hat.ace.vec)),col="blue")
abline(v=(var(hat.ace.vec)),col="blue")



## (3a.2) Simulating potential outcomes under additivity with larger sample size
n <- 200  #finite population size
k <- 100  #treatment group size

ace.true <- 3
yctrl <- rnorm(n,10,2) # simulating Y(0) [from N(10,2)]
ytrt <- yctrl + ace.true  # additivity
mean(ytrt) - mean(yctrl) # true ACE
nsims <- 10000
covered <- rep(NA,nsims) # vector for coverage indicators
covered.aronow <- rep(NA,nsims)
hat.ace.vec <- rep(NA,nsims)
hat.var.hat.ace.vec <- rep(NA,nsims)
hat.var.hat.ace.vec.aronow <- rep(NA,nsims)

for(i in 1:nsims){
  trt.grp <- sample((1:n),k)  #choose our treatment group
  x <- rep(0,n)
  x[trt.grp] <- 1  #construct treatment vector
  y <- yctrl*(x==0) + ytrt*(x==1)  #build observed data vector
  ### Neyman procedure and Aronow
  mean.ctrl <- mean(y[x==0])
  mean.trt <- mean(y[x==1])
  hat.ace <-  mean.trt - mean.ctrl
  hat.ace.vec[i] <- hat.ace	# store the ACE estimate
  var.ctrl <- var(y[x==0])  # variance of control grp
  var.trt <- var(y[x==1])   # variance of trt grp
  hat.var.hat.ace <- var.trt/k + var.ctrl/(n-k)
  hat.var.hat.ace.aronow <- sharp.var(y[x==1],y[x==0])
  covered[i] <- (((hat.ace - 1.96*sqrt(hat.var.hat.ace)) < ace.true) &
                   ((hat.ace + 1.96*sqrt(hat.var.hat.ace)) > ace.true))
  covered.aronow[i] <- (((hat.ace - 1.96*sqrt(hat.var.hat.ace.aronow)) < ace.true) &
                          ((hat.ace + 1.96*sqrt(hat.var.hat.ace.aronow)) > ace.true))
  hat.var.hat.ace.vec[i] <- hat.var.hat.ace # store var(hat ACE) estimate
  hat.var.hat.ace.vec.aronow[i] <- hat.var.hat.ace.aronow # store Aronow var(hat ACE) estimate
}
sum(covered)/nsims  ## about 95%
sum(covered.aronow)/nsims 
mean(hat.ace.vec)  # are we unbiased for the ACE?
var(hat.ace.vec)   # true variability
mean(hat.var.hat.ace.vec) # are we unbiased?
mean(hat.var.hat.ace.vec.aronow) # what about Aronow's estimate?

par(mfrow=c(1,2))
plot(hat.ace.vec,hat.var.hat.ace.vec, main="hat(ACE) vs. Naive est of ACE",
     ylab="hat(var(hat(ACE)))", xlab="hat(ACE)" )

plot(hat.ace.vec,hat.var.hat.ace.vec.aronow, main="hat(ACE) vs. Aronow Var Est of ACE",
     ylab="hat(var(hat(ACE)))^Aro", xlab="hat(ACE)" )

#######
par(mfrow=c(1,2))

hist(hat.var.hat.ace.vec,main="Naive Estimate")
abline(v=(mean(hat.var.hat.ace.vec)),col="red") #mean of distribution of estimates hat.var(ace.hat)
abline(v=(var(hat.ace.vec)),col="blue")  #true variance var(ace.hat)

hist(hat.var.hat.ace.vec.aronow,main="Aronow Estimate")
abline(v=(mean(hat.var.hat.ace.vec.aronow)),col="red") #mean of distribution
abline(v=(var(hat.ace.vec)),col="blue")  #true variability


par(mfrow=c(1,1))

plot(hat.var.hat.ace.vec,hat.var.hat.ace.vec.aronow, main="Naive vs. Aronow",
     ylab="Aronow", xlab="Naive" )
abline(a=0,b=1,col="red")
abline(h=(var(hat.ace.vec)),col="blue")
abline(v=(var(hat.ace.vec)),col="blue")



## (3b) Simulating potential outcomes under non-additivity

n <- 200  #finite population size
k <- 100  #treatment group size

ace.true <- 3
yctrl <- rnorm(n,10,2) # simulating from N(10,2)
heterog <- rnorm(n,0,3) # building disturbance vector
heterog <- heterog - mean(heterog) # ensuring sums to 0
ytrt <- yctrl + ace.true + heterog # building ytrt
mean(ytrt) - mean(yctrl) # True ACE
nsims <- 10000
covered <- rep(NA,nsims) # vector for coverage indicators
covered.aronow <- rep(NA,nsims)
hat.ace.vec <- rep(NA,nsims)
hat.var.hat.ace.vec <- rep(NA,nsims)
hat.var.hat.ace.vec.aronow <- rep(NA,nsims)

for(i in 1:nsims){
  trt.grp <- sample((1:n),k)  #choose our treatment group
  x <- rep(0,n)
  x[trt.grp] <- 1  #construct treatment vector
  y <- yctrl*(x==0) + ytrt*(x==1)  #build observed data vector
  ### Neyman procedure and Aronow
  mean.ctrl <- mean(y[x==0])
  mean.trt <- mean(y[x==1])
  hat.ace <-  mean.trt - mean.ctrl
  hat.ace.vec[i] <- hat.ace	# store the ACE estimate
  var.ctrl <- var(y[x==0])  # variance of control grp
  var.trt <- var(y[x==1])   # variance of trt grp
  hat.var.hat.ace <- var.trt/k + var.ctrl/(n-k)
  hat.var.hat.ace.aronow <- sharp.var(y[x==1],y[x==0])
  covered[i] <- (((hat.ace - 1.96*sqrt(hat.var.hat.ace)) < ace.true) &
                   ((hat.ace + 1.96*sqrt(hat.var.hat.ace)) > ace.true))
  covered.aronow[i] <- (((hat.ace - 1.96*sqrt(hat.var.hat.ace.aronow)) < ace.true) &
                          ((hat.ace + 1.96*sqrt(hat.var.hat.ace.aronow)) > ace.true))
  hat.var.hat.ace.vec[i] <- hat.var.hat.ace # store var(hat ACE) estimate
  hat.var.hat.ace.vec.aronow[i] <- hat.var.hat.ace.aronow # store Aronow var(hat ACE) estimate
}
sum(covered)/nsims  ## about 95%
sum(covered.aronow)/nsims 
mean(hat.ace.vec)  # are we unbiased for the ACE?
var(hat.ace.vec)   # true variability
mean(hat.var.hat.ace.vec) # are we unbiased?
mean(hat.var.hat.ace.vec.aronow) # what about Aronow's estimate?

par(mfrow=c(1,2))
plot(hat.ace.vec,hat.var.hat.ace.vec, main="hat(ACE) vs. Naive est of ACE",
     ylab="hat(var(hat(ACE)))", xlab="hat(ACE)" )

plot(hat.ace.vec,hat.var.hat.ace.vec.aronow, main="hat(ACE) vs. Aronow Var Est of ACE",
     ylab="hat(var(hat(ACE)))^Aro", xlab="hat(ACE)" )

#######
par(mfrow=c(1,2))

hist(hat.var.hat.ace.vec,main="Naive Estimate")
abline(v=(mean(hat.var.hat.ace.vec)),col="red") #mean of distribution of estimates hat.var(ace.hat)
abline(v=(var(hat.ace.vec)),col="blue")  #true variance var(ace.hat)

hist(hat.var.hat.ace.vec.aronow,main="Aronow Estimate")
abline(v=(mean(hat.var.hat.ace.vec.aronow)),col="red") #mean of distribution
abline(v=(var(hat.ace.vec)),col="blue")  #true variability


par(mfrow=c(1,1))

plot(hat.var.hat.ace.vec,hat.var.hat.ace.vec.aronow, main="Naive vs. Aronow",
     ylab="Aronow", xlab="Naive" )
abline(a=0,b=1,col="red")
abline(h=(var(hat.ace.vec)),col="blue")
abline(v=(var(hat.ace.vec)),col="blue")




#######################################################
#######################################################

## (4) Fisher's test of the null hypothesis of no effect

## Vesikari data (from lecture)

x <- c(rep(1,100),rep(0,100))
y   <- c(rep(1,63),rep(0,37),rep(1,48),rep(0,52))
table(y,x)

obs.test.stat <- sum(x*(1-y))

n <- length(x)
k <- sum(x)

nsims <- 5000
test.stat.vec <- rep(NA,nsims)

for(i in 1:nsims){
  trt.grp <- sample((1:n),k)  #choose our treatment group
  xsim <- rep(0,n)
  xsim[trt.grp] <- 1  #construct treatment vector
  test.stat <- sum(xsim*(1-y))
  test.stat.vec[i] <- test.stat
}


sum(test.stat.vec < (obs.test.stat+0.5))/nsims

## Compare to the true value
phyper(37,100,100,89)





set.seed(17)
#######################################################
##  (5) Fisher's randomization test continuous response
##  using difference in means as statistic
##
##  C is control; T is treatment
##
#data 
trt<-c(rep("C",10),rep("T",10))
y<-c(rnorm(10,0.75,3),rnorm(10,2.25,1)) # simulated data; different vars; ACE = 1.5

yt<-y[trt=="T"]
yc<-y[trt=="C"]


#descriptive analysis
par(mfrow=c(1,2))
hist(yt) ; hist(yc)

par(mfrow=c(1,1))
boxplot(y~as.factor(trt))

mean(yc)   
mean(yt)   

sd(yc)
sd(yt)    

## Using coin package
library(coin)
oneway_test(y~as.factor(trt), distribution="exact")                              


# Via simulation 

diff.obs<- (mean(yt)-mean(yc)) 

nsim<-10000
diff.null<-rep(0,nsim)
p.value<-rep(0,nsim)      


for(i in 1:nsim) {
  
  trt.sim<-sample(trt) 
  yc.sim<-y[trt.sim=="C"]
  yt.sim<-y[trt.sim=="T"]
  
  
  diff.sim<- (mean(yt.sim)-mean(yc.sim))
  
  diff.null[i]<-diff.sim
  
  p.value[i]<-mean(abs(diff.null[1:i])>=abs(diff.obs))
  #can you figure out what the above line does?
} 

#has the *estimated* p-value converged?
plot(p.value)
#note that p.value[nsim] will converge to the 
#correct p-value as nsim -> infinity. To decide
#if your nsim is big enough, use the above plot 
#and decide if it has converged

#ok, I think it has converged
hist(diff.null,prob=T, main="Difference in means")    # note prob=T makes the height equal to a proportion
abline(v= c( -abs(diff.obs), abs(diff.obs)), lty=2,col="blue")
p.value.rdist<-mean( abs(diff.null)>=abs(diff.obs) )
p.value.rdist

# compare to a normal distribution

sd.diff <- sqrt(var(yc)/(length(yc)) + var(yt)/(length(yt)))

x<-seq(-3*sd.diff,3*sd.diff,length=100)
lines(x, dnorm(x,0,sd.diff),col="green")



#######################################################
##  (5b) Fisher's randomization test continuous response
##  using difference in means divided by naive variance

# Using same data, but different test statistic in randomization test


#using the difference in means divided by naive variance estimate 
# as statistic. 

k <- length(yt)  #treatment group size
n <- length(y)

diff.obs<- (mean(yt)-mean(yc))

hat.var.hat.ace <- (var(yt)/k + var(yc)/(n-k)) 

stat.obs <- diff.obs/sqrt(hat.var.hat.ace)

nsim<-10000
stat.null<-rep(0,nsim)
p.value<-rep(0,nsim)      


for(i in 1:nsim) {
  
  trt.sim<-sample(trt) 
  yt.sim<-y[trt.sim=="C"]
  yc.sim<-y[trt.sim=="T"]
  hat.var.hat.ace.sim <- (var(yt.sim)/k + var(yc.sim)/(n-k)) 
  diff.sim<- (mean(yt.sim)-mean(yc.sim))
  stat.sim <- diff.sim/sqrt(hat.var.hat.ace.sim)
  stat.null[i]<-stat.sim
  p.value[i]<-mean(abs(stat.null[1:i])>=abs(stat.obs))
} 

#has the *estimated* p-value converged?
plot(p.value)
#note that p.value[nsim] will converge to the 
#correct p-value as nsim -> infinity. To decide
#if your nsim is big enough, use the above plot 
#and decide if it has converged

#ok, I think it has converged
hist(stat.null,prob=T,main="diff in means divided by naive sd ace")    # note prob=T makes the height equal to a proportion
abline(v= c( -abs(stat.obs), abs(stat.obs)), lty=2,col="blue")
p.value.rdist.naive.var <-mean( abs(stat.null)>=abs(stat.obs) )
p.value.rdist.naive.var

# compare to a normal distribution

sd.diff <- 	1

x<-seq(-3*sd.diff,3*sd.diff,length=100)
lines(x, dnorm(x,0,sd.diff),col="green")



#######################################################
##  (5c) Fisher's randomization test continuous response
##  using difference in means divided by Aronow variance

# Using same data, but different test statistic in randomization test


#using the difference in means divided by naive variance estimate 
# as statistic. 

k <- length(yt)  #treatment group size
n <- length(y)

diff.obs<- (mean(yt)-mean(yc))


hat.var.hat.ace.aronow <- sharp.var(yt,yc)


stat.obs <- diff.obs/sqrt(hat.var.hat.ace.aronow)

nsim<-10000
stat.null<-rep(0,nsim)
p.value<-rep(0,nsim)      


for(i in 1:nsim) {
  
  trt.sim<-sample(trt) 
  yt.sim<-y[trt.sim=="C"]
  yc.sim<-y[trt.sim=="T"]
  hat.var.hat.ace.aronow.sim <- sharp.var(yt.sim,yc.sim)
  diff.sim<- (mean(yt.sim)-mean(yc.sim))
  stat.sim <- diff.sim/sqrt(hat.var.hat.ace.aronow.sim)
  stat.null[i]<-stat.sim
  p.value[i]<-mean(abs(stat.null[1:i])>=abs(stat.obs))
} 

#has the *estimated* p-value converged?
plot(p.value)
#note that p.value[nsim] will converge to the 
#correct p-value as nsim -> infinity. To decide
#if your nsim is big enough, use the above plot 
#and decide if it has converged

#ok, I think it has converged
hist(stat.null,prob=T,main="diff in means divided by Aronow sd ace")    # note prob=T makes the height equal to a proportion
abline(v= c( -abs(stat.obs), abs(stat.obs)), lty=2,col="blue")
p.value.rdist.aronow.var <-mean( abs(stat.null)>=abs(stat.obs) )
p.value.rdist.aronow.var

# compare to a normal distribution

sd.diff <- 	1

x<-seq(-3*sd.diff,3*sd.diff,length=100)
lines(x, dnorm(x,0,sd.diff),col="green")



