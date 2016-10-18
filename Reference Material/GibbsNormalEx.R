gibbs.sampler.normal<-function(data = gen.reject.normal(100), ksi = 0, tau = 1, a = 1, b = 1, start = c(0, 1), numofgen = 100)
{
  # we will first generate mu then sigma
  n <- length(data)
  xbar <- mean(data)
  chain <- matrix(0, numofgen, 2)
  chain[1, ] <- start
  a1 <- n/2 + a
  for(i in 2:numofgen) {
    var1 <- (tau^2 * chain[i - 1, 2])/(n * tau^2 + chain[i - 1, 2])
    mean1 <- var1 * (ksi/tau^2 + (n * xbar)/chain[i - 1, 2])
    chain[i, 1] <- rnorm(1, mean1, sd=sqrt(var1))
    b1 <- 1/b + 1/2 * sum((data - chain[i, 1])^2)
    chain[i, 2] <-1/rgamma(1,a1,b1)
    #1/gen.reject.gamma(1, a = a1, b = 1/b1)
  }
  plot(chain[, 1], type = "l", main = "Chain for mu", xlab = "mu")
  plot(chain[, 2], type = "l", main = "Chain for sigma^2", xlab = "sigma^2")
  hist(chain[, 1], main = "Histogram for mu", xlab = "mu")
  hist(chain[, 2], main = "Histogram for sigma^2", xlab = "sigma^2")
  return(chain)
}

runit=function(data=rnorm(10),L=50000,burnin=.1*L,alpha=.1, ksi=0,tau=1,a=2,b=1, start=c(0,1))
{
  x=gibbs.sampler.normal(data=data, ksi=ksi,tau=tau,a=a,b=b, start=start,numofgen = L)
  x=x[(burnin+1):L,]
  postmeans=colSums(x)/L
  #credible sets
  CSmu=quantile(x=x[,1],c(alpha/2,(1-alpha/2)))
  CSsigsq=quantile(x=x[,2],c(alpha/2,(1-alpha/2)))
  cat("\nPosterior means, mu=",postmeans[1],"sigmasq=",postmeans[2],"")
  cat("\n",100*(1-alpha),"% CS for mu: [",CSmu[1],",",CSmu[2],"]")
  cat("\n",100*(1-alpha),"% CS for sigmasq: [",CSsigsq[1],",",CSsigsq[2],"]")
}

#Calls to the function:
# small sample size, but "perfect" hyper-parameters
runit(data=rnorm(10), L=50000,alpha=.05,ksi=0,tau=1,a=2,b=1, start=c(0,1))

# large sample size, but still not good
runit(data=rnorm(100),L=50000,alpha=.05, ksi=0,tau=100,a=1,b=1, start=c(0,1))

# large sample size, with large b, means very small variation in the inverse gamma
runit(data=rnorm(1000), L=50000,alpha=.05, ksi=0,tau=1,a=1,b=100, start=c(0,1))

# (semi) non-informative prior for mu
runit(data=rnorm(1000), L=50000,alpha=.05, ksi=0,tau=1000000000000000,a=1,b=1, start=c(0,1))

# run it with horrible starting values and hyper-parameters, prior dominates the results
runit(data=rnorm(1000),L=50000,alpha=.05,  ksi=-100,tau=.01,a=10,b=1000001, start=c(-10000,100000000))
