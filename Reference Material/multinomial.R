# R code for Figure 3.2, BDA3 p. 70

# First create a function for generating random samples from Dirichlet distribution
set.seed(8863155)     # Optional: just to make sure I get the same answer next time
rdirichlet <- function(n,alpha) {
    myMatrix <- matrix(NA,n,length(alpha))
    for (i in 1:length(alpha)) {
      myMatrix[,i] <- rgamma(n,alpha[i],1);
    }
    myMatrix <- myMatrix/apply(myMatrix,1,sum);
    myMatrix;
  }  

# Generate a sample of 100,000 from the posterior

thetas <- rdirichlet(100000, c(728,584,174))

# Setting freq=F generates a density histogram

hist(thetas[,1]-thetas[,2], nclass=50,freq=F)

# Sometimes it's usefult to plot a density estimate

lines(density(thetas[,1]-thetas[,2]))

#  Plot the normal approximation
m <- (728 - 584)/1450
v <- (728*(1450-728) + 584*(1450 - 584) + 2*728*584)/1450^2/1451
s <- sqrt(v)
lines(xx <- seq(0, .2, , 100), dnorm(xx, m, s), lty=2)

# Estimate a 95\% credible interval
post.mc <- thetas[, 1] - thetas[, 2]
quantile(post.mc, c(.025, .975))

      2.5%      97.5% 
0.04953588 0.14422711  

#  Normal approximation

c(m - 1.96*s, m + 1.96*s)

[1] 0.05063317 0.14798752


