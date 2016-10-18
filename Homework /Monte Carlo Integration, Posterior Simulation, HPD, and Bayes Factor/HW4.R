#Problem 3.12 

Scale everything to 1976

x <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
y <- c(24, 25, 31, 31, 22, 21, 26, 20, 16, 22)

crude = glm(y~x,family=poisson(link=log))

#Log Likelihood Function

pl <- function(a, b) {
	sumy <- sum(y)
	sumyx <- sum(y*x)
	summean <- sum(exp(a+b*x))
	exp(a*sumy+b*sumyx-summean)	
}

#Run grid on Posterior

alpha <- seq(from =0, to=4, by = .01)
beta <- seq(from = -.10, to = .05, by = .0001)
na <- length(alpha)
nb <- length(beta)
z <- matrix(0, na, nb)

for (i  in (1:na)){
    for (j  in (1:nb)){
        z[i,j] = pl(alpha[i],beta[j])
        }
    }

par(mfrow=c(1,1))  
contour(z, nlevels=10, xlim=c(.6, 1), ylim=c(0,.8))


#Sample from Joint Density - Part g

pa <- apply(z, 1, sum)
pa <- pa/sum(pa)      
post.alpha <- vector()
post.beta <- vector()


for (k in 1:5000){
	i <- sample(1:na, 1, prob=pa)
# conditional distribution of beta given alpha
	pb.a <- z[i, ]/sum(z[i, ])
	j <- sample(1:nb, 1, prob=pb.a)
# add random jitter as in BDA3 so sample looks continuous
	post.alpha[k] <- alpha[i] + runif(1, -2, 2)
	post.beta[k] <- beta[j] + runif(1, -.075, .075)
}

partg <- exp(post.alpha + post.beta*10)
hist(partg)

# Part h

parth <- vector()
for (i in 1:length(partg)) {
	parth[i] <- rpois(1, partg[i])
}

quantile(parth, c(.025, .975))


#Problem 4.1 
#Part b

newton.raphson.univariate<<-function(f = function(x)
{
  x
}
, fprime = function(x)
{
  1
}
, xo = 1, iterno = 10)
{
  x <- xo
  for(i in 1:iterno) {
    x <- x - f(x)/fprime(x)
    cat(c("iteration ", i, ", current value x=", x, "with f(x)=", f(x), "\n"))
    if(abs(f(x))<.0001)
      break
  }
  return(c(x, f(x), fprime(x)))
}

y <- c(-2, -1, 0, 1.5, 2.5)
fprime <- function(theta, obs) {2*sum((y-theta)/(1+(y-theta)^2))}
fdoubleprime <- function(theta, obs) {2*sum((((y-theta)^2)-1)/(1+(y-theta)^2)^2)}

newton.raphson.univariate(fprime, fdoubleprime, xo=.5, iterno=10)

#Part c Approximate 

(-1*fdoubleprime(-.13764, y))^-1

m <- -.13764
s2 <- 0.7273305
x <- seq(-5,5,0.1)*sqrt(s2) + m
hx <- dnorm(x,m,sqrt(s2))
hx2 <- hx / sum(hx)

plot(x, hx2, xlab="theta values", ylab="posterior values", main="Approximate Normal Distribution", ylim=c(0, .045), xlim=c(-4, 4))

##Part c EXACT

par(mfrow=c(1, 1))

f <- function(theta, obs, n) {prod(1/(1+(obs-theta)^2))}

d <- seq(-5,5,0.1)
vec <- vector()
for (i in 1: length(d)){
	vec[i] <- f(d[i], obs=y, n=5)

}

vec2 <- vec / sum(vec)

par(new=TRUE)
plot(d, vec2, type='l', ylim=c(0, .045), xlim=c(-4, 4), xlab="theta values", ylab="posterior values")

#4.2

x <- c(-.86, -.30, -.05, 0.73)
n <- c(5, 5, 5, 5)
death <- c(0, 1, 3, 5)
a <- 0.8
b <- 7.7
tmp <- a+b*x

info <- matrix(nrow=2, ncol=2)

info[1,1] <- sum((n*exp(tmp))/(1+exp(tmp))^2)
info[2,2] <- sum((n*x^2*exp(tmp))/(1+exp(tmp))^2)
info[1,2] <- sum((n*x*exp(tmp))/(1+exp(tmp))^2)
info[2,1] <- sum((n*x*exp(tmp))/(1+exp(tmp))^2)

solve(info)

#Problem 6

r.duniform <- function(n){
	sample((n-1):(n+1), 1)
}

v6 <- replicate(1000, r.duniform(1000))

mean(v6)


#Problem 7

v7 <- replicate(10000, sample(1:12, 1))
counts <- table(v7)
relative.freq <- counts / 10000


#Problem 8

mc.integral.unif<<-
function(numofgen = 100, p = 1, fun = function(x)
{
  x^0
}
)
{
  unifs <- matrix(0, numofgen, p)
  for(i in 1:numofgen)
    unifs[i, ] <- runif(p)
  mc.sum <- 0
  for(i in 1:numofgen)
    mc.sum <- mc.sum + fun(unifs[i, ])
  return(mc.sum/numofgen)
}

# Part a

fun1 <- function(x){(1/x^2)*exp(-(1-(1/x))/2)}
fun2 <- function(x){(1/x^2)*exp(-((1/x)-1)/2)}

mc.integral.unif(numofgen = 50000, p = 1, fun = fun1)+ mc.integral.unif(numofgen = 50000, p = 1, fun = fun2)

# Part b

fun3 <- function(x){(1/x^2)*((1/x)-1)/(1+((1/x)-1)^2)^2}

mc.integral.unif(numofgen = 50000, p = 1, fun = fun3)

# Part c

fun4 <- function(x){4*exp(4*x-2+(4*x-2)^2)}

mc.integral.unif(numofgen = 50000, p = 1, fun = fun4)

# Part d

fun5 <- function(x){(1/x^2)*exp(-(-(1/x)+1)^2)}
fun6 <- function(x){(1/x^2)*exp(-(-(1/x)+2)^2)}

mc.integral.unif(numofgen = 50000, p = 1, fun = fun5)*mc.integral.unif(numofgen = 50000, p = 1, fun = fun6)

# Part e

Set value of a
fun7 <- function(x, a) {9*((9*x+1)^a)/factorial(9*x)}
mc.integral.unif(numofgen = 50000, p = 1, fun = fun7)

mc.integral.unif.modified <- function(numofgen = 50000, p = 1, a)
{
  	store <- vector()
	unifs <- matrix(0, numofgen, p)
  	for(i in 1:numofgen)
    	unifs[i, ] <- runif(p)
		mc.sum <- 0
		for(i in 1:numofgen){
	    		mc.sum <- mc.sum + 9*((9*unifs[i, ]+1)^a)/factorial(9*unifs[i, ])
		}
		mc.sum/numofgen
	
}

#Example when a=2

mc.integral.unif.modified(50000, 1, 2)

#Problem 9
#Part a

exp.mle <- function(x) {
	length(x)/sum(x)
}

r <- rexp(20, 10)

exp.mle(r)

#Part b

exp.map <- function(x, a, b) {
	(length(x)+a-1)/(b+sum(x))
}

r <- rexp(10, 5)
a <- 2 
b <- 3

exp.map(r, a, b)



#Problem 10

gen.reject.beta<<-
function(numofgen = 1, a = 0.5, b = 0.5)
{
  if(a == 1 && b == 1)
    return("use runif(1)")
  if(2 - a - b == 0) return("oops inappropriate values") #compute the constant c
  t <- (1 - a)/(2 - b - a)
  const <- (t^(a - 1) * (1 - t)^(b - 1) * gamma(a + b))/(gamma(a) * gamma(b))
  cat("\naverage number of iterations to get one value is fixed= ", const, "\n")
  gens <- matrix(0, numofgen, 1)
  for(i in 1:numofgen) {
    #generate two independent uniforms(0,1)
    u1 <- runif(1)
    y <- (u1^(a - 1) * (1 - u1)^(b - 1))/(t^(a - 1) * (1 - t)^(b - 1))
    u <- runif(1)
    while(u >= y) {
      u1 <- runif(1)
      y <- (u1^(a - 1) * (1 - u1)^(b - 1))/(t^(a - 1) * (1 - t)^(b - 1))
      u <- runif(1)
    }
    gens[i] <- u1
  }
  return(as.vector(gens))
}


hpd.beta.binomial <- function(x = 1, a = 1, b = 1, n = 1, numofgen = 100, alpha = 0.05){

  # x=the observe number of successes, a,b the hyper parameters of the Beta prior
  # n=the number of trials, and alpha affects the coverage we want.
 
 pdf.beta <- function(x, a = a, b = b){
    (x^(a - 1) * (1 - x)^(b - 1) * gamma(a + b))/(gamma(a) * gamma(b))
 }
  
	p <- sort(gen.reject.beta(numofgen = numofgen, a = x + a, b = n - x + b))
	densityatp <- pdf.beta(x = p, a = x + a, b = n - x + b)
	sorted <- sort(densityatp)
	percentile <- sorted[floor(numofgen * alpha)]
	roots <- matrix(0, numofgen, 1)
	count <- 0
  for(i in 1:numofgen)
    if(abs(densityatp[i] - percentile) <= 0.001) {
      print(p[i])
      count <- count + 1
      roots[count] <- p[i]
    }
  
roots <- roots[1:count] 
plot(p, densityatp, type = "l")  
abline(v = roots, h = percentile)
return(list(roots, percentile))

}



p10 <- function(x, m, a, b) {
	
	par(mfrow=c(1, 2))
	pgrid <- seq(0,1,length=1000)
	priordensity <- dbeta(pgrid, a, b)
	plot(pgrid, priordensity, ylab="density")
	
	posterior <- vector()
	term1 <- sum(x)+a
	term2 <- b+m*length(x)-sum(x)
	
	posteriordensity <- dbeta(pgrid, term1, term2)

	plot(pgrid, posteriordensity, type="l", ylab="density")

}


n=20 generated data from Binom(5, .45), a = 2, b = 2

testy <- rbinom(20, 5, .45)
hpd.beta.binomial(testy, 2, 2, 5, numofgen=100, alpha=.025)

bayesfactor <- function(x,m,a,b){

	term1 <- sum(x)+a
	term2 <- b+m*length(x)-sum(x)

	posteriornull <- (.5^(term1-1))*(1-po)^(term2-1)
	
	p1 <- mean(rbeta(100000, term1, term2))
	posterioralt <- (p1^(term1-1))*(1-p1)^(term2-1)

	priornull <- (.5^(a-1))*(.5^(b-1))
	prioralt <- (p1^(a-1))*(p1^(b-1))

	bayesfactor <- (posteriornull/posterioralt)/(priornull/prioralt)

	bayesfactor	

}



#Part i

x <- rbinom(20,5,.45)
p10(x, m=5, a=2, b=2)
hpd.beta.binomial(x, a = 2, b = 2, n = 5, numofgen = 100, alpha = 0.05)
bayesfactor(x, m=5, a=2, b=2)

Part ii

x <- rbinom(10, 10, .95)
p10(x,m=10, a=1, b=1)
hpd.beta.binomial(x, a = 1, b = 1, n = 10, numofgen = 100, alpha = 0.05)
bayesfactor(x, m=10, a=1, b=1)

Part iii

x <- rbinom(100, 20, .5)
p10(x, m=20, a=2, b=1)
hpd.beta.binomial(x, a = 2, b = 1, n = 20, numofgen = 100, alpha = 0.05)
bayesfactor(x, m=20, a=2, b=1)

i) n=20 generated data from Binom(5, .45), a = 2, b = 2.
ii) n=10 generated data from Binom(10, .95), a = 1, b = 1.
iii) n=100 generated data from Binom(20, .5), a = 2, b = 1.





















