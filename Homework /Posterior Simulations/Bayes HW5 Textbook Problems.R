##Problem 5.6


theta <- seq(-3, 3, .01)
prior <- c(.95, .05)
dens <- prior[1]*dnorm(theta, 1, .5) + prior[2]*dnorm(theta, -1, .5)
plot(theta, dens, ylim=c(0,1.1*max(dens)), type="l", xlab="theta", ylab="", xaxs="i", yaxs="i", yaxt="n", bty="n", cex=2)
mtext("prior density", cex=2, 3)

marg <- dnorm(-.25, c(1,-1), sqrt(c(0.5, 0.5)^2+1/10))
posterior <- (prior *marg)/sum(prior*marg)
 
post <- posterior[1]*dnorm(theta, 1.5/14, sqrt(1/14)) + posterior[2]*dnorm(theta, -6.5/14, sqrt(1/14))
plot(theta, post, ylim=c(0,1.1*max(post)), type="l", xlab="theta", ylab="", xaxs="i", yaxs="i", yaxt="n", bty="n", cex=2)
mtext("posterior density", cex=2, 3)

##PROBLEM 5.11

bikes <- c(16, 9, 10, 13, 19, 20, 18, 17, 35, 55)
other <- c(58, 90, 48, 57, 103, 57, 86, 112, 273, 64)
total <- bikes + other
prop <- bikes/total

bdata <- matrix(c(bikes, total, prop), nrow=10, ncol=3)
mean(prop)
var(prop)

##FIND STARTING ESTIMATES OF ALPHA AND BETA FOR GRID

ab <- ((mean(prop)*(1-mean(prop)))/var(prop))-1
ab*mean(prop)
ab*(1-mean(prop))

2.58*c(1,4)
10.59*c(1,4)

log(2.58/42.36)
log(10.32/10.59)

log(2.58+10.32)
log(10.59+42.36)

log(alpha/beta) grid values -> (-2.79, -0.02)
log(alpha + beta) grid values -> (2.57, 3.96)

OR

alpha grid values -> (2.58, 10.32)
beta grid values ->  (10.59, 42.36)

#RUN GRID TO GET CONTOUR PLOT 

bikes <- c(16, 9, 10, 13, 19, 20, 18, 17, 35, 55)
other <- c(58, 90, 48, 57, 103, 57, 86, 112, 273, 64)
total <- bikes + other
prop <- bikes/total

log.post.tarone <-  function(gamma,delta){
	alpha   <- exp(gamma+delta)/(1+exp(gamma))
	beta    <- exp(delta)/(1+exp(gamma))
	ldens <- sum(lgamma(alpha + bikes) + lgamma(beta + total - bikes) - 
	lgamma(alpha + beta + total)) + length(total)*(lgamma(alpha + beta) - lgamma(alpha) - lgamma(beta))
    ldens - 5/2*log(alpha+beta) + log(alpha) + log(beta)
}

##gamma <- seq(-2.79, -0.02, length=100)
##delta <- seq(2.57, 3.96, length=100)

##CHANGED ABOVE GRID VALUES AFTER FIRST PLOT

gamma <- seq(-2.79, -0.02, length=100)
delta <- seq(1, 4.5, length=100)

contours <- seq(.05, .95, .1)
logdens <- matrix(0, 100, 100)

for (i in 1:100){
	for (j in 1:100){
		logdens[i, j] <- log.post.tarone(gamma[i], delta[j])
	}
}

dens <- exp(logdens-max(logdens))
contour(gamma,delta,dens,levels=contours,xlab="log(alpha/beta)",
	ylab="log(alpha+beta)",la)

##SAMPLE FROM POSTERIOR

pg <- apply(dens, 1, sum)
pg <- pg/sum(pg)
ng <- length(pg)
nd <- ncol(dens)

N <- 1000
post.alpha <- post.beta <- post.theta <- rep(NA, N)

for (k in 1:N){
	i <- sample(1:ng, 1, prob=pg)
# conditional distribution of beta given alpha
	pd.g <- dens[i, ]/sum(dens[i, ])
	j <- sample(1:nd, 1, prob=pd.g)
# No random jitter here; compare with BDA3
	post.alpha[k] <- alpha <- exp(gamma[i]+delta[j])/(1+exp(gamma[i]))
	post.beta[k] <- beta <- exp(delta[j])/(1+exp(gamma[i]))
	}
	
	
##PART C
	
post.theta.mean <- rep(NA, 10)	
for (j in 1:10) {
	for (k in 1:N){
		post.theta[k] <- rbeta(1, post.alpha[k] + bikes[j], post.beta[k] + total[j] - bikes[j])
	}
	post.theta.mean[j] <- mean(post.theta)
}

x <- cbind(bdata, post.theta.mean)
colnames(x) <- c("bikes", "total", "raw proportion", "posterior proportion")

##PART D

N<-1000
overall.post.theta <- rep(NA,N)

for (k in 1:N){
		overall.post.theta[k] <- rbeta(1, 10*post.alpha[k] + sum(bikes), 10*post.beta[k] + sum(total) - sum(bikes))
	}

quantile(overall.post.theta, c(.025, .975))
	
##PART E 

N <- 1000
predictive <- rep(NA, N)

for(k in 1:1000){
	predictive[k] <- rbinom(1, 100, overall.post.theta[k])
}

quantile(predictive, c(.025, .975))

##PART F 

Beta distribution for theta j reasonable 

