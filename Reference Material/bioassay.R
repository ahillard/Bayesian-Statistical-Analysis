# R-code for the contour plot in Figure 3.3(a) ---- 
#
x = c(-.86, -.30, -.05, .73)
n = c(5, 5, 5, 5)
y = c(0, 1, 3, 5)

# Likelihood function
      
pl <- function(a, b){
	tmp <- exp(a + b*x)
	p <- tmp/(1 + tmp)
	tmp2 <- y*log(p) + (n - y)*log(1 - p)
	exp(sum(tmp2))
}

alpha = seq(-5,10,0.1) 
beta  = seq(-10,40,0.5)
na    = length(alpha)
nb    = length(beta)
z     = matrix(0,na,nb)

for (i  in (1:na)){
    for (j  in (1:nb)){
        z[i,j] = pl(alpha[i],beta[j])
        }
    }

# Here is how to save a plot as a postscript or pdf file:

#  Figure 3.3(a)
#  Pick one of the next two statements
#postscript(file='graph.ps',horizontal=F,width=7.2,height=10)
#pdf(file='graph1.pdf',width=7.2,height=10)
par(mfrow=c(1,1))                             #  optional
# persp(alpha,beta,z,xlab="x",ylab="y",zlab="z")
contour(alpha,beta,z, nlevels=10)
#dev.off()

#------To get a closer contour plot ----------
alpha2 = seq(0.75,0.90,0.001)
beta2  = seq(7.4,8.0,0.002)
na2    = length(alpha2)
nb2    = length(beta2)
z     = matrix(0,na2,nb2)
for (i  in (1:na2)){
    for (j  in (1:nb2)){
        z[i,j] = pl(alpha2[i],beta2[j])
        }
    }
contour(alpha2,beta2,z,nlevels=25)

#  Sample from the discretized posterior: see BDA3
#  Idea: sample alpha from marginal distribution and sample beta given alpha

# Note that p(i, j) is prop. to z

#  There is a problem with underflow in the edges of z. Set to missing values to zero.

z[is.na(z)] <- 0

# marginal distribution for alpha: renormalize
pa <- apply(z, 1, sum)        # sum across the rows to get marginal pdf
pa <- pa/sum(pa)              # normalize so sum is 1.

set.seed(4029414)    #  Start random number generator in same place for demonstration
post.alpha <- post.beta <- rep(NA, 5000)
for (k in 1:5000){
	i <- sample(1:na, 1, prob=pa)
# conditional distribution of beta given alpha
	pb.a <- z[i, ]/sum(z[i, ])
	j <- sample(1:nb, 1, prob=pb.a)
# add random jitter as in BDA3 so sample looks continuous
	post.alpha[k] <- alpha[i] + runif(1, -.05, .05)
	post.beta[k] <- beta[j] + runif(1, -.25, .25)
}

#  Figure 3.3(b)
plot(post.alpha, post.beta, xlim=c(-5, 10), ylim=c(-10, 50), pch='.', cex=2.5)

# Now look at the posterior distribution of the ld50: 

ld50 <- -post.alpha/post.beta

#  Figure 3.4
hist(ld50, nclass=50, xlim=c(-1, 1))
mean(ld50)
# [1] -0.105453

# Any negative betas?
sum(post.beta < 0)

# approx. Bayes estimates:

mean(post.alpha)
mean(post.beta)
# [1] 1.297378
# [1] 11.52929

#  find posterior mode:

logl <- function(x) -log(pl(x[1], x[2]))

# The optim() function finds the minimizer of a multivariate function. Output includes
# lots of information including the vector where minimum is achieved: $par

theta.hat <- optim(c(1.4, 12), logl)$par
# [1] 0.8464302 7.7462467

############################################################################
#  Normal approximation to the likelihood
X <- cbind(rep(1, 4), x)
tmp <- exp(X%*%theta.hat)
p.hat <- tmp/(1 + tmp)
            # [,1]
# [1,] 0.002585492
# [2,] 0.187093628
# [3,] 0.611511260
# [4,] 0.998675271

#  Set up Fisher information of the sample
In <- matrix(0, 2, 2)
In[1, 1] <- sum(n*p.hat*(1 - p.hat))
In[1, 2] <- In[2, 1] <- sum(n*x*p.hat*(1 - p.hat))
In[2, 2] <- sum(n*x*x*p.hat*(1 - p.hat))

# Take square root of determinant
det.In <- sqrt(det(In))


# Store values of normal posterior approximation for contour plot

z1 <- matrix(0,na,nb)
for (i  in (1:na)){
    for (j  in (1:nb)){
    	theta <- c(alpha[i], beta[j])
        z1[i,j] <- exp(-.5*t(theta - theta.hat)%*%In%*%(theta-theta.hat))*det.In/(2*pi)
        }
    }

#  Figure 4.1(a)
pdf("graph2.pdf", height=5, width=10)
par(mfrow=c(1, 2))
contour(alpha,beta,z, nlevels=12, main='True posterior')
contour(alpha,beta,z1, nlevels=12, main='Normal approximation')
dev.off()

#  Simulate from the normal approximation:

# Covariance matrix for normal approximation
V <- solve(In)
A <- chol(V)         #Cholesky decomposition, so A'A = V

set.seed(5405084)    # Starting value for random numbers so I can reproduce this
N <- 5000
norm.sample <- matrix(0, 2, N)
for (i in 1:N) norm.sample[, i] <- theta.hat + t(A)%*%rnorm(2)

# Figures 3.3(b) and 4.1(b)
pdf("graph3.pdf", height=5, width=10)
par(mfrow=c(1, 2))
plot(post.alpha, post.beta, xlim=c(-5, 10), ylim=c(-10, 50), pch='.', cex=2.5, main='Sample from posterior')
plot(norm.sample[1, ], norm.sample[2, ], xlim=c(-5, 10), ylim=c(-10, 50), pch='.', cex=2.5, main='Sample from norm. approx.')
dev.off()

#  Normal approximation to posterior of ld50:

ld50.norm <- -norm.sample[1, ]/norm.sample[2, ]

#  Extreme outliers:
quantile(ld50.norm)
         # 0%         25%         50%         75%        100% 
# -4.653819e+03 -1.764759e-01 -1.174852e-01 -4.453192e-02  9.165704e+01 

mean(ld50.norm)
# [1] -1.045771

#  Figure 4.2: Histogram of central 95%
q.025 <- quantile(ld50.norm, .025)
q.975 <- quantile(ld50.norm, .975)
ix <- (ld50.norm > q.025) & (ld50.norm < q.975)
hist(ld50.norm[ix], nclass=50, prob=T, xlim=c(-1, 1))

pdf("graph4.pdf", height=5, width=10)
par(mfrow=c(1, 2))
hist(ld50, nclass=50, xlim=c(-1, 1))
hist(ld50.norm[ix], nclass=50, prob=T, xlim=c(-1, 1))
dev.off()

#  Compare the two credible intervals for ld50:

quantile(ld50, c(.025, .975))
      2.5%      97.5% 
# -0.2697082  0.1016277 
 
quantile(ld50.norm, c(.025, .975))
      2.5%      97.5% 
# -0.6632878  0.4663803 


