# Example: analysis bioassay experiment
# BDA3, p102, 71 experiments, Tarone(1982) data
# plot joint posterior of gamma and delta

tarone  <- matrix(scan("tarone.txt"),byrow=T,ncol=2)
deaths  <- tarone[,1]
rats    <- tarone[,2]
log.post.tarone <-  function(gamma,delta){
	alpha   <- exp(gamma+delta)/(1+exp(gamma))
	beta    <- exp(delta)/(1+exp(gamma))
	ldens <- sum(lgamma(alpha + deaths) + lgamma(beta + rats - deaths) - 
		lgamma(alpha + beta + rats)) + length(rats)*(lgamma(alpha + beta) - 							lgamma(alpha) - lgamma(beta))
    ldens - 5/2*log(alpha+beta) + log(alpha) + log(beta)
}


#postscript ("posterior.tarone.ps")
pdf("posterior.tarone.pdf")
#Plot Figure 5.3(a), page 112
ll <- 50
gamma <- seq(-2.3,-1.3,length=ll)
delta <- seq(1,5,length=ll)
contours<- seq(.05,.95,.1)
logdens <- matrix(0, ll, ll)
for (i in 1:ll){
	for (j in 1:ll){
		logdens[i, j] <- log.post.tarone(gamma[i], delta[j])
	}
}
#  rescale to prevent overflow
dens <- exp(logdens-max(logdens))
contour(gamma,delta,dens,levels=contours,xlab="log(alpha/beta)",
	ylab="log(alpha+beta)",la)

dev.off()

#  Sample from the posterior for lab 71:
pg <- apply(dens, 1, sum)
pg <- pg/sum(pg)
ng <- length(pg)
nd <- ncol(dens)

N <- 1000
post.alpha <- post.beta <- post.theta71 <- rep(NA, N)

for (k in 1:N){
	i <- sample(1:ng, 1, prob=pg)
# conditional distribution of beta given alpha
	pd.g <- dens[i, ]/sum(dens[i, ])
	j <- sample(1:nd, 1, prob=pd.g)
# No random jitter here; compare with BDA3
	post.alpha[k] <- alpha <- exp(gamma[i]+delta[j])/(1+exp(gamma[i]))
	post.beta[k] <- beta <- exp(delta[j])/(1+exp(gamma[i]))
	post.theta71[k] <- rgamma(1, alpha + deaths[71], beta + rats[71] - deaths[71])
	}
#  posterior mean
mean(post.theta71)
# [1] 0.2703153
#  95% credible interval
quantile(post.theta71, c(.025, .975))
#     2.5%     97.5% 
# 0.1048994 0.5207433 

mean(post.alpha)
mean(post.beta)
# [1] 2.357516
# [1] 14.11325
