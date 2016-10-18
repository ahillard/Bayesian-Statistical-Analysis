#  R code for rejection method for Beta density with beta >= 1
rejbeta <- function(al, bet){
	repeat{
		u1 <- runif(1)
		x <- u1^{1/al}  #  generate a beta(alpha, 1) r.v.
		u2 <- runif(1)
		if (u2 < (1 - x)^{bet - 1}) break
	}
	x
}
#  Generate 10,000 Beta(2, 4) random variables
y <- rep(NA, 10000)
for (i in 1:length(y)) y[i] <- rejbeta(2, 4)
pdf("beta-simulation.pdf", height=6, width=8)
hist(y, freq=FALSE)
zz <- seq(0, 1, , 80)
lines(zz, dbeta(zz, 2, 4))
dev.off()


#  R example of rejection sampling for standard normal
#  Use g(x) = exp(-x), x > 0.

# This function produces a single Z, so no argument needed.
rejnorm <- function(){
# Note use of repeat{} loop: continues until break
	repeat{				
		u1 <- runif(1)
		x <- -log(1 - u1)
		u2 <- runif(1)
		if(u2 < exp(-(x - 1)^2/2)) break
	}
	u3 <- runif(1)
	if(u3 < .5) z = -x
	else z = x
	z
}

y <- rep(NA, 10000)
for (i in 1:length(y)) y[i] <- rejnorm()
mean(y)
var(y)
pdf("std-norm-simulation.pdf", height=6, width=8)
hist(y, freq=FALSE)
zz <- seq(-3, 3, , 100)
lines(zz, dnorm(zz), lty=2)
dev.off()



