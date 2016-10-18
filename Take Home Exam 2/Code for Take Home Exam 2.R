##Andrew Hillard
##Take Home Exam 2

------------------------------------------------------------------------------------------------------------

##Problem 1

acc <- c(24, 25, 31, 31, 22, 21, 26, 20, 16, 22) #Fatal Accidents
death <- c(734, 516, 754, 877, 814, 362, 764, 809, 223, 1066) #Passenger Deaths
rate <- c(.19, .12, .15, .16, .14, .06, .13, .13, .03, .15) #Death Rate - Passenger Deaths per 100 million passenger miles
mi <- death/rate #Passenger miles
year <-1976:1985

##Part I of Take Home Exam

##Plot Code

plot.gamma.poisson <- function(y = acc, a = 1, b = 1, L=20000, alpha = 0.05) {
	
		  	rgamma.prior <- sort(rgamma(1000, a, b))
		  	rgamma.post <-  sort(rgamma(1000, a + sum(y), b + 10))
 	 		densprior <- dgamma(rgamma.prior,a,b)
  			denspost <- dgamma(rgamma.post, a + sum(y), b + 10)
  			m <- mean(rgamma.post)
  			xmax <- max(rgamma.prior, rgamma.post)
  			ymax <- max(densprior,denspost)
  			plot(rgamma.prior,densprior,type="l",lty=2,xlab="theta", ylab="density", main="Solid=Post, Dashed=Prior", xlim=c(0, xmax), ylim=c(0,ymax))
  			points(rgamma.post,denspost,type="l",lty=1)	
  			abline(v=m)
  			
}

##Credible Set Code

cs.gamma.poisson <- function(y = acc, a = 1, b = 1, L=20000, alpha = 0.05) {
	csexact <- qgamma(c(alpha/2, (1-alpha/2)), a + sum(y), b + 10)  
	csexact
}

##HPD Code

hpd.gamma.poisson <- function(y = acc, a = 1, b = 1, L=20000, alpha = 0.05) {
	
  		p <- sort( rgamma(L, a + sum(y), b + 10) )
  		densityatp <- dgamma(p, a + sum(y), b + 10)
  		sorted <- sort(densityatp)
 		percentile <- sorted[floor(L * alpha)]
  		roots <- matrix(0, L, 1)
  		count <- 0
  			for(i in 1:L)
   				 if(abs(densityatp[i] - percentile) <= 0.001) {
     			 ##print(p[i])
     			 count <- count + 1
     			 roots[count] <- p[i]
    		}
  		roots <- roots[1:count]
  		roots <- c(min(roots), max(roots))
  		##plot(p, densityatp, type = "l")
 		##abline(v = roots, h = percentile)
 		return(list(roots, percentile))
 		
}

##Problem 1 Code

p1.gamma.poisson <- function(y = acc, a = 1, b = 1, L=20000, alpha = 0.05) {

	tmp <- y
	plot.gamma.poisson(tmp, a, b, L, alpha)
	cs <- cs.gamma.poisson(tmp, a, b, L, alpha)
	hpd <- hpd.gamma.poisson(tmp, a, b, L, alpha)
	
	cat(paste("\nExact CS =[", cs[1], ",", cs[2],"]"))
	cat(paste("\nHPD = [", hpd[[1]][1], ",", hpd[[1]][2], "]"))
	
}

#Part a, b, c, d

p1.gamma.poisson(acc, 1, 1)
p1.gamma.poisson(acc, 10, 10)
p1.gamma.poisson(acc, 100, 100)

p1.gamma.poisson(death, 1, 1)
p1.gamma.poisson(death, 10, 10)
p1.gamma.poisson(death, 100, 100)

##Problem 1 Part II

##Plot Code

plot.gamma.poisson.m <- function(y = acc, m=mi, a = 1, b = 1, L=20000, alpha = 0.05) {
	
		  	rgamma.prior <- sort(rgamma(1000, a, b))
		  	rgamma.post <-  sort(rgamma(1000, a + sum(y), b + sum(m)))
 	 		densprior <- dgamma(rgamma.prior,a,b)
  			denspost <- dgamma(rgamma.post, a + sum(y), b + sum(m))
  			m <- mean(rgamma.post)
  			xmax.prior <- max(rgamma.prior)
  			ymax.prior <- max(densprior)
  			xmax.post <- max(rgamma.post)
  			ymax.post <- max(denspost)
  			par(mfrow=c(2,1))
  			plot(rgamma.prior,densprior,type="l",lty=2,xlab="theta", ylab="density", main="Prior Distribution", xlim=c(0, xmax.prior), ylim=c(0,ymax.prior))
  			plot(rgamma.post,denspost,type="l",lty=1, main="Posterior Distribution", xlim=c(0, xmax.post), ylim=c(0,ymax.post))
  			abline(v=m)

}

##Credible Set Code

cs.gamma.poisson.m <- function(y = acc, m=mi, a = 1, b = 1, L=20000, alpha = 0.05) {
	csexact <- qgamma(c(alpha/2, (1-alpha/2)), a + sum(y), b + sum(m))  
	csexact
}

##HPD Code

hpd.gamma.poisson.m <- function(y = acc, m=mi, a = 1, b = 1, L=20000, alpha = 0.05) {
	
  		p <- sort( rgamma(L, a + sum(y), b + sum(m)) )
  		densityatp <- dgamma(p, a + sum(y), b + sum(m))
  		sorted <- sort(densityatp)
 		percentile <- sorted[floor(L * alpha)]
  		roots <- matrix(0, L, 1)
  		count <- 0
  			for(i in 1:L)
   				 if(abs(densityatp[i] - percentile) <= 0.001) {
     			 ##print(p[i])
     			 count <- count + 1
     			 roots[count] <- p[i]
    		}
  		roots <- roots[1:count]
  		roots <- c(min(roots), max(roots))
  		##plot(p, densityatp, type = "l")
 		##abline(v = roots, h = percentile)
 		return(list(roots, percentile))
 		
}

##Problem 1 Code

p1.gamma.poisson.m <- function(y = acc, m=mi, a = 1, b = 1, L=20000, alpha = 0.05) {
	
	plot.gamma.poisson.m(y, m, a, b, L, alpha)
	cs <- cs.gamma.poisson.m(y, m, a, b, L, alpha)
	hpd <- hpd.gamma.poisson.m(y, m, a, b, L, alpha)
	
	cat(paste("\nExact CS =[", cs[1], ",", cs[2],"]"))
	##cat(paste("\nHPD = [", hpd[[1]][1], ",", hpd[[1]][2], "]"))
	
}

##Had to rerun HPD multiple times to pick up interval. This is why I put hpd function outside of p1.gamma.poisson.m function. 

p1.gamma.poisson.m(acc, mi, 1, 1)
hpd.gamma.poisson.m(acc, mi, 1, 1)

p1.gamma.poisson.m(acc, mi, 10, 10)
hpd.gamma.poisson.m(acc, mi, 10, 10)

p1.gamma.poisson.m(acc, mi, 100, 100)
hpd.gamma.poisson.m(acc, mi, 100, 100)

------------------------------------------------------------------------------------------------------------

##Problem 2

y <- c(3.7, -6.7, -10.5, -6.1, -17.6, 2.3, -4.5, -7.7, -9.4, -10.4, -10.9, -9.3)

#Part a

p2.marginal.likelihood <- function(obs=y, model.v = 25, prior.u=-5, prior.v=100){
	
	u <- rnorm(50000, prior.u, sqrt(prior.v) )
	model.sd <- sqrt(model.v)
	n <- length(obs)
	py <- vector()
	for (i in 1: length(u)) {
		##py[i] <- sum(log(dnorm(y,u[i], model.sd)))
		py[i] <- prod(dnorm(y, u[i], model.sd)) 
		##py[i] <- exp(sum(log(dnorm(y,u[i], model.sd))))	
	}
	
	mean(py)

}



#Part b

p2.bayes.factor <- function(obs=y, model.v = 25){
	
	##From Notes
	
	n <- length(y)
	ybar <- mean(y)
	z <- (sqrt(n) * abs(ybar--4) ) / sqrt(model.v)
	B <- 1 / exp(.5 * z^2)
	B
	
} 



#Part c - In Notes

p2.posterior.probs <- function(obs=y, model.v = 25, prior.u=-5, prior.v=100) {
	
	n <- length(obs)
	ybar <- mean(y)
	z <- (sqrt(n) * abs(ybar--4) ) / sqrt(model.v)
	p0 <- .4
	p1 <- .6
	num <- exp(.5*(z^2)*(1+model.v / (n*prior.v))^-1)
	den <- (1+ (n*prior.v) / model.v)^(1/2)
	post.0 <- (1+((1-p0)/p0)*(num/den))^-1
	post.0

}



#Part d

p2.cs <- function(obs=y, model.v = 25, prior.u=-5, prior.v=100, alpha=.01){
	
	n <- length(y)
	post.u <- (model.v / (n*prior.v+model.v))*prior.u + ((n*prior.v)/(n*prior.v+model.v))*mean(y)
	post.v <- (prior.v*model.v)/(n*prior.v+model.v)
	
	csexact <- qnorm(c(alpha/2, (1-alpha/2)),  post.u, post.v)  
	csexact
	
}

#Problem 2 Solution Function

p2 <- function(obs=y, model.v = 25, prior.u=-5, prior.v=100, alpha=.01){



	ml <- p2.marginal.likelihood(obs, model.v, prior.u, prior.v)
	B <- p2.bayes.factor(obs, model.v)
	pp <- p2.posterior.probs(obs=y, model.v = 25, prior.u=-5, prior.v=100)
	cs <- p2.cs(obs, model.v, prior.u, prior.v, alpha)

	cat(paste("\nThe Marginal Likelihood for Part D = ", ml))
	cat(paste("\nThe Bayes Factor of H0 vs. H1 = ",  B))
	cat(paste("\nThe Posterior Probability of H0 = ", pp))
	cat(paste("\nExact CS = [", cs[1], ",", cs[2],"]"))
	
}

#Run Problem 2 Function

p2(y)
	
------------------------------------------------------------------------------------------------------------

##Problem 3

##Install Quartz11 to use geoR package - link: http://www.xquartz.org/
install.packages("geoR")
library(geoR)

#Part a

p3.jointsample <- function(obs=x) {
	
	n <- length(obs)
	ybar <- mean(obs)
	s2 <- var(obs)
	sample.post.v <- rinvchisq(100000, n-1, s2)
	sample.post.u <- rnorm(100000, ybar, sample.post.v/n)
	sample.joint <- matrix(NA, nrow=100000, ncol=2)
	sample.joint[,1] <- sample.post.u
	sample.joint[,2] <- sample.post.v
	colnames(sample.joint) <- c("mean", "variance")
	sample.joint
	
}

#Part b

x <- rnorm(30, 10, sqrt(9))

marginal.approximation <- function(obs=x){
	
	tmp <- p3.jointsample(obs)
	sample.post.u <- tmp[,1]
	sample.post.v <- tmp[,2]
	sample.y <- vector()
	for (i in 1:100000) {
		sample.y[i] <- rnorm(1, sample.post.u[i], sqrt(sample.post.v[i]))	
	}
	
	u <- mean(sample.y)
	v <- var(sample.y)
	sd <- sqrt(v)
	cat(paste("\nApproximately a Normal Distribution with Mean =", u, "and Standard Deviation = ", sd))
	
}

marginal.approximation(x)

#Part c

p3.cs <- function(obs=x, model.u=1, model.v=1){
	
	tmp <- obs
	gen <- p3.jointsample(tmp)
	sample.post.u <- gen[,1]
	sample.post.v <- gen[,2]
	
	cs.u <- quantile(sample.post.u, c(.05,.95))
	cs.v <- quantile(sample.post.v, c(.05, .95))
	cat(paste("\nApproximate CS for the Mean = [", cs.u[1], ",", cs.u[2],"]"))
	cat(paste("\nApproximate CS for the Variance = [", cs.v[1], ",", cs.v[2],"]"))
	
}

p3.cs(x)

#Part d

#Part i
gen.reject.student.t <- function(obs = x, numofgen=1) {
	
	xbar <- mean(obs)
	s2 <- var(obs)
	n <- length(obs)
	constant <- 10 
	gens <- matrix(0, numofgen, 1)
  	
  	for(i in 1:numofgen) {
    		y <- rnorm(1, 10, sqrt(10))
    		den <- dnorm(y, 10, sqrt(10))
      		num <-  (1+n * ((y - xbar)^2) / ( (n-1) * s2)) ^ (-n/2)
      		u0 <- 1
      		u <- 1
    			while(u >= u0) {
      					u <- runif(1)
      					y <- rnorm(1, 10, sqrt(10))
      					den <- dnorm(y, 10, sqrt(10))
      					num <-  (1+n * ((y - xbar)^2) / ( (n-1) * s2)) ^ (-n/2)
      					u0 <- num / (constant*den)
    			}
    			
    gens[i] <- y

	}
	
	return(as.vector(gens))
	
}

#Part ii

p3.pd.plot <- function(obs=x, numofgen=1000){
	
	xbar <- mean(obs)
	s2 <- var(obs)
	n <- length(obs)
	u <- sort(gen.reject.student.t(obs, numofgen))
	dens <- vector()
	for(i in 1 : length(u)){
		dens[i] <- (1+n * ((u[i] - xbar)^2) / ( (n-1) * s2)) ^ (-n/2)
	}
	
	xmin <- min(u)
	xmax <- max(u)
	ymin <- min(dens)
	ymax <- max(dens)
	plot(u, dens, xlim=c(xmin,xmax), ylim=c(ymin, ymax), xlab="mean values", ylab="dens", main="Marginal Density")
	
}


par(mfrow=c(1,1))
p3.pd.plot(x)

------------------------------------------------------------------------------------------------------------


##Problem 4 

##Sample using nonfunctional probabilities 

nonfunctional.p <- vector()
for (i in 1:20) {
	nonfunctional.p[i] <- .5 + i/50
}

p4.sample.nonfunctional <- function(numofgen=100000, p=nonfunctional.p){
	
	sample.all <- matrix(NA, numofgen, 20)
	for (i in 1:numofgen) {
		for(j in 1:20){
			sample.all[i, j] <- sample(0:1, 1,  prob=c(1-p[j], p[j]))
		}
	}
	
	sample.all
}


##Or sample using functional probabilities 

nonfunctional.p <- vector()
for (i in 1:20) {
	nonfunctional.p[i] <- .5 + i/50
}

functional.p <- 1-nonfunctional.p

p4.sample.functional <- function(numofgen=100000, p=functional.p){
	
	sample.all <- matrix(NA, numofgen, 20)
	for (i in 1:numofgen) {
		for(j in 1:20){
			sample.all[i, j] <- sample(0:1, 1,  prob=c(1-p[j], p[j]))
		}
	}
	
	sample.all
}


#Part a with non functional probabilities 

sample.pa <- p4.sample.nonfunctional(100000, p=nonfunctional.p)
sim.prob.pa <- table(apply(sample.pa, 1, sum)) / 100000 #Add up probabilities where count is less than or equal to 5

#OR Part a with functional probabilities

sample.func <- p4.sample.functional(100000, p=functional.p)
sim.prob.pa <- table(apply(sample.func, 1, sum)) / 100000 #Add up probabilities where count is greater than or equal to 15

#Part b - Based on Rules of Conditional Probability

sim.prob.pa[4] / sum(sim.prob.pa[1:4])

##Or, through conditional simulation...

p4.sample.conditional <- function(numofgen=100, p=nonfunctional.p) {
	
	sample.all <- matrix(NA, numofgen, 20)
	for (i in 1:numofgen) {
		holder <- rep(1, 20)
		while( sum(holder) > 5) {
				for(j in 1:20){
					holder[j] <- sample(0:1, 1,  prob=c(1-p[j], p[j]))
				}
		}	
		
		sample.all[i, ] <- holder	
			
	}
	
	sample.all
}

##This code will take a long time to run because probability of getting X less than or equal to 5 is very low. 
cond <- p4.sample.conditional(10000, nonfunctional.p)
p5.cond <- table(apply(cond, 1, sum)) 
p5.cond[length(p5.cond)]/sum(p5.cond)

##Alternative code that you can use to test and give you rough estimate. Run chunks at a time. 
tmp1<- p4.sample.conditional(25, nonfunctional.p)
tmp2<- p4.sample.conditional(25, nonfunctional.p)
tmp3<- p4.sample.conditional(25, nonfunctional.p)
tmp4<- p4.sample.conditional(25, nonfunctional.p)
tmp5 <- p4.sample.conditional(25, nonfunctional.p)
tmp6 <- p4.sample.conditional(50, nonfunctional.p)
total <- rbind(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6)
p5.cond <- table(apply(total, 1, sum)) 
p5.cond[length(p5.cond)]/sum(p5.cond)



------------------------------------------------------------------------------------------------------------

##Problem 5

p5 <- function() {
	
	gen <- vector()
	while (sum(gen) < 1) {
		gen <- append(gen, runif(1))
	}
	
	gen
	length(gen)
}
	
	
p5.sim <- function(numofgen) {
	
	sim <- vector()
	for (i in 1:numofgen) {
		sim[i] <- p5()
		
	}
	
	sim
	
}

simulation <- p5.sim(50000)
mean(simulation)
var(simulation)

##The mean is the natural number e

------------------------------------------------------------------------------------------------------------

##Problem 6

##Part b

sample.joint <- function(obs=x, numofgen=10000) {
	
	obs.m <- mean(obs)
	obs.v <- var(obs)
	beta.est <- obs.m/obs.v
	alpha.est <- obs.m*beta.est
	
	obs.no.0 <- obs[obs!=0]
	n1 <- length(obs.no.0)
	
	obs.0 <- obs[obs==0]
	n2 <- length(obs.0)
	
	post.theta <- vector()
	post.lambda <- vector()
	
	for (k in 1:numofgen) {
	
		tmp.theta <- vector()
		tmp.lambda <- vector()
		for (i in 1:n1) {
			tmp.theta[i] <- rbeta(1, 1, 2)
			tmp.lambda[i] <- rgamma(1, obs.no.0[i] + alpha.est, 1+beta.est)
		}
	
		for(j in 1:n2){
			u <- runif(1)
			if (u <.5) {
				lam <- rgamma(1, alpha.est, beta.est)
				tmp.lambda <- append(tmp.lambda, lam)
			} 
			if (u > .5) {
				lam <- rgamma(1, alpha.est, 1+beta.est)
				tmp.lambda <- append(tmp.lambda, lam)
			}
		
			u <- runif(1)
			mix.prob1 <- 1 / (1+exp(-lam))
			if (u < mix.prob1){
				theta <- rbeta(1, 2, 1)
				tmp.theta <- append(tmp.theta, theta)
			}
				if (u > mix.prob1){
				theta <- rbeta(1, 1, 2)
				tmp.theta <- append(tmp.theta, theta)
			}
			
		}
	
		post.theta[k] <- mean(tmp.theta)
		post.lambda[k] <- mean(tmp.lambda)
	
	}
	
	joint <- matrix(NA, numofgen, 2)
	joint[,1] <- post.theta
	joint[,2] <- post.lambda
	joint
	
}

p6.pb <- function(obs=x, numofgen=10000) {
	
	y <- obs
	joint <- sample.joint(y, numofgen=10000)
	
	sample.theta.post <- joint[,1]
	sample.lambda.post <- joint[,2]
	
	theta.post.u <- mean(sample.theta.post)
	lambda.post.u <- mean(sample.lambda.post)
	
	cat(paste("\nThe Posterior Mean of Theta = ", theta.post.u))
	cat(paste("\nThe Posterior Mean of Lambda = ", lambda.post.u))
	
	theta.cs <- quantile(sample.theta.post, c(.05,.95))
	lambda.cs <- quantile(sample.lambda.post, c(.05, .95))
	cat(paste("\nApproximate CS for the Mean of Theta = [", theta.cs[1], ",", theta.cs[2],"]"))
	cat(paste("\nApproximate CS for the Mean of Lambda = [", lambda.cs[1], ",", lambda.cs[2],"]"))

}

gen.y <- p6.pc(numofgen=1000) #Load code in part c to generate these mock observations 

p6.pb(gen.y)

##Part c													

gen.inverse.transform.zip.not.0 <- function(pmf = function(y, theta=.25, lambda=5) { ((1-theta)*(lambda^y)*exp(-lambda)) / factorial(y)}, theta=.25, lambda=5, numofgen = 1, minx = 1){
  	gens <- matrix(0, numofgen, 1)
  	u <- runif(numofgen)
 	prob.0 <- theta +(1-theta)*exp(-lambda) 
	prob.not.0 <- 1 - prob.0
  	for(i in 1:numofgen) {
  	 	cdf <- 0
   	 	count <- minx
   	 	repeat {
      		if(cdf <= u[i] && u[i] < cdf + pmf(count)/prob.not.0)  {
        		gens[i] <- count
        		break
      		}
     		cdf <- cdf + pmf(count)/prob.not.0
      		count <- count + 1
    		}
  	}
  	return(as.vector(gens))
}

p6.pc <- function(numofgen=1, theta=.25, lambda=5){
	
	prob.0 <- theta +(1-theta)*exp(-lambda) 
	gens <- vector()
	for(i in 1:numofgen){
		u <- runif(1, 0, 1)
		if (u <= prob.0) {
			gens[i] <- 0
		}
		if (u > prob.0){
			gens[i] <- gen.inverse.transform.zip.not.0()
		}
	}	
	gens	
}

p6.pc(numofgen=1000)

#Part d

p6.pd <- function(obs=x, numofgen=100000){
	
	obs.m <- mean(obs)
	obs.v <- var(obs)
	beta.est <- obs.m/obs.v
	alpha.est <- obs.m*beta.est
	
	y.sample <- vector()
	for (i in 1:numofgen){
		t <- runif(1)
		l <- rgamma(1, alpha.est, beta.est)
		y.sample[i] <- p6.pc(1, t, l)
		
	}
	
	y.sample
	u <- mean(y.sample)
	v <- var(y.sample)
	cat(paste("\nThe Marginal of y has an approximate ZIP distribution with mean", u, "and variance of", v))
	
}

p6.pd(gen.y)

#Part e

p6.plot.pmf <- function(x=1, theta=.25, lambda=5) {
	(theta+(1-theta)*exp(-lambda))*I(x==0) + ((1-theta)*(lambda^x)*exp(-lambda) / factorial(x) )*I(x>0)
}

p6.pe <- function(obs=x, numofgen=20000){
	
	obs.m <- mean(obs)
	obs.v <- var(obs)
	beta.est <- obs.m/obs.v
	alpha.est <- obs.m*beta.est
	
	tmp1 <- obs
	tmp2 <- numofgen
	joint <- sample.joint(tmp1, tmp2)
	theta <- joint[,1]
	lambda <- joint[,2]
	y.pred.sample <- vector()
	
	for (i in 1:numofgen){
		
		y.pred.sample[i] <- p6.pc(1, theta[i], lambda[i])
		
	}
	
	u <- mean(y.pred.sample)
	v <- var(y.pred.sample)
	cat(paste("\nThe Predictive Marginal Distribution of y has an approximate ZIP distribution with mean", u, "and variance of", v))
	
	y.pred.sample <- sort(y.pred.sample)
	barplot(table(y.pred.sample) / numofgen)
	
}

p6.pe(gen.y, 20000)

------------------------------------------------------------------------------------------------------------

##Problem 7
obs <- c(164, 142, 110, 153, 103, 52, 174, 88, 178, 184, 58, 62, 132, 128)

emp.dist <- function(x, obs) {
	
	n <- length(obs)
	count <- 0
	
	for (i in 1:n){
		if (obs[i] <= x) {
			count <- count + 1
		}
	}
	
	count / n
	
}

ks.stat <- function(obs, min=50, max=200) {
	
	obs.sort <- sort(obs)
	n <- length(obs)
	x.emp <- vector()
	x.f <- vector()
	
	for(i in 1:n){
		x.emp[i] <- emp.dist(obs.sort[i], obs.sort)
		x.f[i] <- punif(obs.sort[i], 50, 200)
	}
	
	dist <- abs(x.emp - x.f)
	D <- max(dist)
	D
	
}

p7.sim.ks <- function(numofgen=1000, min=50, max=200) {
	
	sim.ks <- vector()
	for(i in 1:numofgen){
		x <- runif(14, 50, 200)
		sim.ks[i] <- ks.stat(x, min=50, max=200)
	}
	sim.ks
	
}

#p-value 

p7.pvalue <- function(obs, numofgen=20000, min=50, max=200){
	
	tmp <- obs
	sim.ks <- p7.sim.ks(numofgen=20000, min, max)
	obs.ks <- ks.stat(tmp, min=50, max=200)
	p.value <- length(sim.ks[sim.ks > obs.ks]) / numofgen
	p.value 
	
}

p7.pvalue(obs)















