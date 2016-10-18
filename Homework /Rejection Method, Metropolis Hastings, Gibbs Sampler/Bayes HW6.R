## Problem 2, Part a

p1 <- rbeta(50000, 7, 5)
p2 <- rbeta(50000, 11, 11)

target <- p1-p2

quantile(target, c(.025, .975))

table(p1 > p2)

##Problem 2, Part b

u <- mean(target)
v <- var(target)

integrand <- function(y) { (1/ (sqrt(2*pi*v))) * exp(- (1 / (2*v))*(y-u)^2)}
integrate(integrand, 0, Inf)

##Problem 3, Part a


beta.rejection.ratio <- function(x, a, b) {
	num <- (x ^ (a-1)) * ((1-x)^(b-1))
	den <- (((1-a) / (2-b-a))^(a-1)) * (1-((1-a)/(2-b-a)))^(b-1)
	num/den
}

p3 <- function(numofgen=1000) {
	
	simulations <- vector()
	
	while (length(simulations) <1000){
		mixed <- runif(1)
		u1 <- runif(1)
		if (mixed < .7) {
			u <- runif(1) 
			ratio <- beta.rejection.ratio(u1, 4, 2)
			if (u < ratio){
				simulations <- append(simulations, u1)
			}	
		}
		if (mixed > .7) {
			u <- runif(1)
			ratio <- beta.rejection.ratio(u1, 1.5, 3)
			if (u < ratio){
				simulations <- append(simulations, u1)
			}	
		}
	}
	simulations
}


sim <- p3()
hist(sim, main="Histogram of Simulated Values Using Rejection Method", xlab="Observations", freq=FALSE)

p3.density <- function(x) {
	0.7 * dbeta(x, 4, 2) + 0.3*dbeta(x, 1.5, 3)
}
grid <- seq(from=0, to=1, length=1000)
dens <- p3.density(grid)

lines(grid, dens)

##Part b

p3.mh <- function(h = function(x){0.7 * dbeta(x, 4, 2) + 0.3*dbeta(x, 1.5, 3)}, numofgen=1000, burn=500) {
	chain <- vector()
	chain[1] <- runif(1)
	MHJumps <- 0
	for (i in 2: (numofgen + burn)) {
		y <- rbeta(1, 1, 1)
		u <- runif(1)
		ratio <- (h(y) * dbeta(chain[i-1], 1, 1)) / (h(chain[i-1])*dbeta(y, 1, 1)) 
		if (u < ratio) {
			chain[i] <- y
			if(i > burn){
			MHJumps <- MHJumps + 1
			}
		}
		else if (u > ratio) {
			chain[i] <- chain[i-1]
		}
	}
	iterations <- seq(1, numofgen, 1)
	gen <- chain[(burn+1) : length(chain)]
	plot(iterations, gen, xlab="iterations", ylab="chain", main="Monte Carlo Chain", type='l')
	cat("\nMH acceptance probability is",100*MHJumps/numofgen,"%\n")
	gen
}

##Part c

p3.ar1.mh <- function(h = function(x){0.7*dbeta(x, 4, 2) + 0.3*dbeta(x, 1.5, 3)}, start = 0.3, sig = 2, b = 0.5, numofgen =1000, burn=500) {
	
	q <- function(x, u, s) {(2*pi*s^2)^(-0.5) * exp( -((x[1] - u*x[2])^2) / s^2)} 
	chain <- vector()	
	chain[1] <- start  
	MHJumps <- 0
	for (i in 2: (numofgen+burn)){
		gen <- 2
		while(gen < 0 | gen > 1){
		gen <- rnorm(1, mean= b*chain[i-1], sd = sig)
		}
		propdens <- q(x = c(chain[i-1], gen), u = b, s = sig)
		ratio = (h(gen)*q(x=c(gen, chain[i-1]), u=b, s=sig)) / (h(chain[i-1]) * propdens)
		if (runif(1) < ratio){
			chain[i] <- gen
			if(i > burn){
			MHJumps <- MHJumps + 1
			}
		}
		else {
			chain[i] <- chain[i-1]
		}
	}
	
	plot(chain, type='l', xlab='iterations', ylab='chain', main='Random Walk Markov Chain')
	cat("\nMH acceptance probability is",100*MHJumps/numofgen,"%\n")
	chain[(burn+1): length(chain)]
}

##Problem 4

##Part b

y <- c(28, 8, -3, 7, -1, 1, 18, 12)
sig <- c(15, 10, 16, 11, 9, 11, 10, 18)

install.packages("pscl")
library("pscl")

p4.gibbs <- function(y, sig, start = c(10, 100), numofgen=1000, burn=500, plot="theta") {
	
	iterations <- numofgen + burn
	theta.chain <- matrix(NA, iterations, 8)
	mean.chain <- vector()
	tau.chain <- vector()
	
	theta.chain[1,] <- y
	mean.chain <- start[1]
	tau.chain <- start[2]
	n <- length(y)
	sig2 <- sig^2
	
	for (i in 2:iterations){
		
		for(j in 1:n) {
			u <- ( y[j]/sig2[j] + mean.chain[i-1]/tau.chain[i-1] ) / (1/sig2[j] + 1/tau.chain[i-1])
			sd <- sqrt( 1 / (1/sig2[j] + 1/tau.chain[i-1]) )
			theta.chain[i, j] <- rnorm(1, u, sd)
		}
		
		theta.mean <- mean(theta.chain[i, ])
		mean.chain <- append(mean.chain, rnorm(1, theta.mean, sqrt(tau.chain[i-1]/n) ) )
		
		tau.chain <- append(tau.chain, rigamma(1, (n-1)/2, .5 * sum( (theta.chain[i,] - mean.chain[i])^2 ) ) )
			
	}
	
	theta.chain.post.burn <- theta.chain[(burn+1):iterations, ]
	mean.chain.post.burn <- mean.chain[(burn+1):iterations]
	tau.chain.post.burn <- tau.chain[(burn+1):iterations]
	
	if (plot == "theta") {
	plot(theta.chain.post.burn[,1], type='l', xlab='iterations', ylab='chain', main='Trace Plot for Theta One')
	return(theta.chain.post.burn)
	}
	
	if (plot == "mean") {
	plot(mean.chain.post.burn, type='l', xlab='iterations', ylab='chain', main='Trace Plot for Mu')
	return(mean.chain.post.burn)
	}
	
	if(plot == "tau"){
	plot(tau.chain.post.burn, type='l', xlab='iterations', ylab='chain', main='Trace Plot for Tau')
	return(tau.chain.post.burn)
	}
	
	if(plot== "difference"){
	diff <- theta.chain.post.burn[,1] - theta.chain.post.burn[,3]
	return(diff)
	}
	
	if(plot == "partf"){
		
		count <- 0 
		for (i in 1: nrow(theta.chain.post.burn)){
			if (theta.chain.post.burn[i, 1] > max(theta.chain.post.burn[i, -1])) {
				count <- count + 1
			}
		}
		
		prob <- count / nrow(theta.chain.post.burn)
		prob
	}
}

##Part c

theta.sim <- p4.gibbs(y, sim, plot="theta")
theta.mean.sim <- mean(theta.sim)

tau.sim <- p4.gibbs(y, sim, plot="tau")
tau.mean.sim <- mean(tau.sim)

mean.sim <- sort(p4.gibbs(y, sig, plot="mean"))
##dens <- dnorm(mean.sim, theta.mean.sim, sqrt(tau.mean.sim/8))

grid <- seq(from=-10, to=30, by=.1)
dens <- dnorm(grid, theta.mean.sim, sqrt(tau.mean.sim/8))
plot(grid, dens, main='Estimated Posterior Density for Mu', type='l')
quantile(mean.sim, c(.025, .975))

##Part d 

sig2 <- sig^2

mu.sim <- p4.gibbs(y, sig, plot="mean")
mu.mean.sim <- mean(mu.sim)

tau.sim <- p4.gibbs(y, sim, plot="tau")
tau.mean.sim <- mean(tau.sim)

u <- ( y[1]/sig2[1] + mu.mean.sim/tau.mean.sim) / (1/sig2[1] + 1/tau.mean.sim)
sd <- sqrt( 1 / (1/sig2[1] + 1/tau.mean.sim) )

grid <- seq(from=-13, to=47, by=.1)
dens <- dnorm(grid, u, sd)

#dens <- dnorm(theta1.sim, u, sd)

plot(grid, dens, main='Estimated Posterior Density for Theta 1', type='l') ##Something is not quite right here

theta1.sim <- p4.gibbs(y, sim, plot="theta") [,1]
theta1.sim <- sort(theta1.sim)
quantile(theta1.sim, c(.025, .975))

##Part e

diff <- sort(p4.gibbs(y, sig, plot="difference"))
quantile(diff, c(.025, .975))

##Part f

p4.gibbs(y, sig, plot="partf")











