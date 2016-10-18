#Problem 5

trunc.pois  <- function(lambda, k) {

s <- vector()

repeat {
		
		random <- rpois(1, lambda)	
		
		if (random <= k) {
			s <- append(s, random)
		}
		
		if (length(s) == 1000) {break}
	
	}
	
	par(mfrow=c(1,2))
	hist(s, main="Histogram of Sample")
	barplot(table(s), main="Bar Plot of Counts")
	cat("The mean is ", mean(s), "and the variance is ", var(s))
	
}



#Problem 6

pairdice <- function() {

	store <- vector()

	repeat {
	
		roll <- sample(1:6, 1) + sample(1:6, 1)
		store <- append(store, roll)
		if(length(table(store)) == 11) {break}
		
	}
	
	plot(table(store), xlab="Sum of Dice Rolls", ylab="Frequency", main="Frequency for Outcome of Two Dice Rolls")
	
}

#Problem 7

#Part i

par(mfrow=c(1,1))

parti <- function(gen=10000) {
	
	gens <- vector()
	dens <- vector()
	u <- runif(gen)
	
	for (j in 1:gen) {
		
		if(u[j]<.5 ) {
			gens[j]  <-  log(2*u[j])/2 
		}
		else if (u[j] >.5){
			gens[j] <- -log(2-2*u[j])/2
		}
	
	}

	for(j in 1: gen){
		
		if(gens[j]<=0){
			dens[j] <- exp(2*gens[j])
		}
		else if (gens[j]>0){
			dens[j] <- exp(-2*gens[j])	
		}
		
	}
	
	plot(gens, dens)
	
}

#Part ii

partii <- function(gen=10000) {
	
	gens <- vector()
	dens <- vector()
	u <- runif(gen)
	
	for (j in 1:gen) {
		sol <- uniroot( function(x, y) {10*x^3-15*x^4 + 6*x^5-y}, c(0,1), y=u[j])
		gens[j] <- sol$root
	}
		
	for(j in 1:gen){
		dens [j] <- 30*(gens[j]^2-2*gens[j]^3+gens[j]^4)
	}
	
	
	plot(gens, dens)
	
}

#Part iii

partiii <- function(gen=10000, a, b) {
	
	gens <- vector()
	dens <- vector()
	u <- runif(gen)
	
	for (j in 1:gen) {
		gens[j] <- (-log(1-u[j])/a)^(1/b)
	}
		
	for(j in 1:gen){
		dens [j] <- a*b*gens[j]^(b-1)*exp(-a*gens[j]^b)
	}
	
	
	plot(gens, dens)
	
}

##Simulation Problem 8

u <- runif(10000000)

#Part i
u1 <- u
u2 <- sqrt(1-u^2)
cov(u1, u2)

#Part ii
u1 <- u^2
u2 <- sqrt(1-u^2)
cov(u1, u2)


#Part iii
u1 <- u
u2 <- exp(u)
cov(u1, u2)

##Simulation Problem 9

install.packages("ggplot2")
library(ggplot2)

prob9 <- function(r, p1, p2, u1, u2, v1, v2){
	
	#Break function
	
	if(p1 + p2 != 1) stop("p1 + p2 does not equal 1")
	
	mixed <- function(x) p1*dnorm(x, u1, sqrt(v1))+p2*dnorm(x, u2, sqrt(v2))
	
	#Set xlim in Plot Function
	lowx <- min(u1, u2)
	if (lowx == u1) {
			lowx <- u1 - 4 * sqrt(v1)
			highx <- u2 + 4 * sqrt(v2)
		}
		
	else if (lowx == u2) {
			lowx <- u2 - 4 * sqrt(v2)
			highx <- u1 + 4 * sqrt(v1)
		}
	
	#Plot 	
	
	plot(mixed, xlim=c(lowx, highx), main="Density Plot of Mixture Distribution using Plot Function")
	
	#Return  r x 1 vector
	rlz <- vector()
	for (i in 1:r){
			rlz[i] <- p1*rnorm(1, u1, sqrt(v1))+p2*dnorm(1, u2, sqrt(v2))
	}
	
	rlz <- t(t(rlz))
	colnames(rlz) < c("realizations") 
	rlz
	
}

#Part iv

#a
rlz <- prob9(1000, .5, .5, 0, 4, 1, 2)
rlz <- data.frame(rlz)
colnames(rlz) <- c("realizations")
ggplot(rlz, aes(x=realizations)) + geom_density()

#b
prob9(1000, .2, .8, 1, 5, 1, 1)
rlz <- data.frame(rlz)
colnames(rlz) <- c("realizations")
ggplot(rlz, aes(x=realizations)) + geom_density()

#c
prob9(1000, .1, .9, 3, 4, .3, .3)
rlz <- data.frame(rlz)
colnames(rlz) <- c("realizations")
ggplot(rlz, aes(x=realizations)) + geom_density()















