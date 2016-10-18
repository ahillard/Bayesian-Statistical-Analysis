pmf.binom<<-function(x,par1=10,par2=.5)
{choose(par1,x)*par2^x*(1-par2)^(par1-x)}

pmf.poisson<<-function(x,lambda=5)
  {exp(-lambda)*lambda^x/factorial(x)}

pmf.geometric<<-function(x,p=.3)
  {p*(1-p)^(x-1)}

gen.inverse.transform.discrete<<-function(pmf = function(x, par1 = 1, par2 = 10)
{
  x^0/(par2 - par1+1)
}
, numofgen = 1, minx = 1, maxx = 10)
{
  # notice the first argument is a function itself that takes 3 arguments. The default
  # function is the discrete uniform from par1,par1+1,...,par2-1,par2
  gens <- matrix(0, numofgen, 1)
  u <- runif(numofgen)
  for(i in 1:numofgen) {
    cdf <- 0
    for(x in minx:maxx) {
      if(cdf <= u[i] && u[i] < cdf + pmf(x)) {
        gens[i] <- x
        break
      }
      cdf <- cdf + pmf(x)
    }
  }
  return(as.vector(gens))
}

gen.binom.invtransf<<-function(numofgen = 1,
                               n = 1, p = 0.5)
{
  gens <- matrix(0, numofgen, 1)
  u <- runif(numofgen)
  for(i in 1:numofgen) {
    prob <- (1 - p)^n
    cdf <- 0
    for(x in 0:n) {
      if(cdf <= u[i] && u[i] < cdf + prob) {
        gens[i] <- x
        break
      }
      cdf <- cdf + prob
      prob <- (prob * (n - x) * p)/((x + 1) * (1 - p))
    }
  }
  return(as.vector(gens))
}

gen.poisson.invtransf<<-function(numofgen = 1,
                                 lambda = 5)
{
  gens <- matrix(0, numofgen, 1)
  u <- runif(numofgen)
  for(i in 1:numofgen) {
    prob <- exp( - lambda)
    cdf <- 0
    count <- 0
    repeat {
      if(cdf <= u[i] && u[i] < cdf + prob) {
        gens[i] <- count
        break
      }
      cdf <- cdf + prob
      prob <- (prob * lambda)/(count + 1)
      count <- count + 1
    }
  }
  return(as.vector(gens))
}

gen.inverse.transform.discrete.inf<<-
function(pmf = function(x, par1 = 0.5, par2 = 10)
{
  par1 * (1 - par1)^(x-1)
}
, numofgen = 1, minx = 1)
{
  # maxx is now infinity
  # default function is geometric
  gens <- matrix(0, numofgen, 1)
  u <- runif(numofgen)
  for(i in 1:numofgen) {
    cdf <- 0
    count <- minx
    repeat {
      if(cdf <= u[i] && u[i] < cdf + pmf(count)) {
        gens[i] <- count
        break
      }
      cdf <- cdf + pmf(count)
      count <- count + 1
    }
  }
  return(as.vector(gens))
}

gen.inverse.transform.discrete.finite<<-
function(x=0:1,probs=c(.5,.5), numofgen = 1)
{
  gens <- matrix(0, numofgen, 1)
  u <- runif(numofgen)
  d=length(probs)
  for(i in 1:numofgen) {
    cdf <- 0
    for(j in 1:d){
      if(cdf <= u[i] && u[i]< cdf + probs[j]) {
        gens[i] <- x[j]
        break
      }
      cdf <-cdf + probs[j]
    }
  }
  return(as.vector(gens))
}

gen.reject.discrete<<-
function(targetprobs = c(0.5, 0.5), probs = c(0.5, 0.5), x = c(-1, 1), numofgen = 1)
{
  #compute the constant c, notice division of vectors here; s-plus understands to divide pairwise the coordinates
  const <- max(targetprobs/probs)
  cat("\naverage number of iterations to get one value = ", const, "\n")
  gens <- matrix(0, numofgen, 1)
  for(i in 1:numofgen) {
    y <- gen.inverse.transform.discrete.finite(x = x, probs = probs, numofgen = 1)
    #find the index that this y corresponds to
    for(j in 1:length(x))
      if(x[j] == y) {
        index <- j
        break
      }
    u <- runif(1)
    while(u >= targetprobs[index]/(const * probs[index])) {
      y <- gen.inverse.transform.discrete.finite(x = x, probs = probs, numofgen = 1)
      for(j in 1:length(x))
        if(x[j] == y) {
          index <- j
          break
        }
      u <- runif(1)
    }
    gens[i] <- y
  }
  return(as.vector(gens))
}

gen.multinomial<<-
function(n = 1, probs = c(0.5, 0.5), numofgen = 1)
{
  # this is the dimension of the multinomial vector (the lineraly dependent one)
  k <- length(probs)
  # the generated values are vectors now: each row
  #will correspond to a different simulated value

  gens <- matrix(0, numofgen, k)
  for(i in 1:numofgen) {
    # we need to generate all k-conditional binomials
    gens[i, 1] <- gen.binom.invtransf(numofgen = 1, n = n, p = probs[1])
    if(k > 2) for(j in 2:(k - 1))
      gens[i, j] <- gen.binom.invtransf(numofgen = 1, n = n - sum(gens[i, 1:(j - 1)]), p = probs[j]/(1 - sum(probs[1:
                                                                                                                     (j - 1)]))) #fix the last value xk
    gens[i, k] <- n - sum(gens[i, 1:(k - 1)])
  }
  return(gens)
}
############################################
###############################################
###########################################
gen.inverse.transform.cont<<-
function(cdfinv = function(x){-log(x)}, numofgen = 1)
{
  u <- runif(numofgen)
  gens <-cdfinv(u)
  return(as.vector(gens))
}
gen.inverse.transform.gamma<<-
function(numofgen = 1,a=1,b=1)
{
  gens=matrix(0,numofgen,1)
  for (i in 1:numofgen)
  {
    u <- runif(a)
    gens[i]<--b*sum(log(u))
  }
  return(as.vector(gens))
}

gen.reject.gamma<<-
function(numofgen = 1, a = 1, b = 1)
{
  #default values yield exp(1)
  #compute the constant c
  const <- (a^a * exp( - a + 1))/gamma(a)
  cat("\naverage number of iterations to get one value = ", const, "\n")
  gens <- matrix(0, numofgen, 1)
  for(i in 1:numofgen) {
    #generate exp(1/ab)
    y <- - a * b * log(runif(1))
    u0 <- (y^(a - 1) * exp(( - y * (a - 1))/(a * b)) * exp(a - 1))/(a * b)^(a - 1)
    u <- runif(1)
    while(u >= u0) {
      y <- - a * b * log(runif(1))
      u0 <- (y^(a - 1) * exp(( - y * (a - 1))/(a * b)) * exp(a - 1))/(a * b)^(a - 1)
      u <- runif(1)
    }
    gens[i] <- y
  }
  return(as.vector(gens))
}

pdf.gamma<<-function(x,a=1,b=1)
{x^(a-1)*exp(-x/b)/(b^a*gamma(a))}

gen.reject.normal<<-
function(numofgen = 1, mu = 0, sig= 1)
{
  #compute the constant c
  const <- sqrt(2*exp(1)/pi)
  cat("\naverage number of iterations to get one value is fixed= ", const, "\n")
  gens <- matrix(0, numofgen, 1)
  for(i in 1:numofgen) {
    #generate two independent exp(1)
    y <- - log(runif(2))
    while(y[2] <=.5*(y[1]-1)^2) y <- - log(runif(2))
    u=runif(1)
    if(u<=.5)gens[i] <- -y[1]
    else gens[i] <- y[1]
  }
  return(as.vector(sig*gens+mu))
}

pdf.normal<<-function(x,mu=1,sig=1)
  {(2*pi*sig^2)^(-.5)*exp(-.5*sig^(-2)*(x-mu)^2)}

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
#-------------------------------
#  Composition methods
gen.continuous<<-
function(numofgen = 1, pdf = "normal", par1 = 0, par2 = 1)
{
  # Check to see which r.v. we want and then generate from it
  if(pdf == "normal") gens <- gen.reject.normal(numofgen = numofgen, mu = par1, sig = par2) else if(pdf == "gamma")
    gens <- gen.reject.gamma(numofgen = numofgen, a = par1, b = par2)
  else if(pdf == "beta")
    gens <- gen.reject.beta(numofgen = numofgen, a = par1, b = par2)
  else return("valid p.d.f.'s are normal, gammas, betas")
  return(gens)
}

gen.compos.cont<<-
function(mixprobs = c(0.5, 0.5),
  pdf1 = "normal", par11 = 0, par12 = 1,
  pdf2 = "normal", par21 = 0, par22 = 1,
  numofgen = 1)
{
  gens <- matrix(1, numofgen, 1)
  u <- runif(numofgen)
  for(i in 1:numofgen) {
    if(u[i] < mixprobs[1])
      gens[i] <- gen.continuous(numofgen = 1, pdf = pdf1, par1 = par11, par2 = par12)
    else gens[i] <- gen.continuous(numofgen = 1, pdf = pdf2, par1 = par21, par2 = par22)
  }
  #plotting the generated values
  if(pdf1 == "normal")
    pdfplot1 <- function(x, a = par11, b = par12)
    {
      (2 * pi * b^2)^(-0.5) * exp((-0.5 * (x - a)^2)/b^2)
    }
  else if(pdf1 == "gamma")
    pdfplot1 <- function(x, a = par11, b = par12)
    {
      (x^(a - 1) * exp( - x/b))/(b^a * gamma(a))
    }
  else if(pdf1 == "beta")
    pdfplot1 <- function(x, a = par11, b = par12)
    {
      (x^(a - 1) * (1 - x)^(b - 1) * gamma(a + b))/(gamma(a) * gamma(b))
    }
  if(pdf2 == "normal")
    pdfplot2 <- function(x, a = par21, b = par22)
    {
      (2 * pi * b^2)^(-0.5) * exp((-0.5 * (x - a)^2)/b^2)
    }
  else if(pdf2 == "gamma")
    pdfplot2 <- function(x, a = par21, b = par22)
    {
      (x^(a - 1) * exp( - x/b))/(b^a * gamma(a))
    }
  else if(pdf2 == "beta")
    pdfplot2 <- function(x, a = par21, b = par22)
    {
      (x^(a - 1) * (1 - x)^(b - 1) * gamma(a + b))/(gamma(a) * gamma(b))
    }
  plot(sort(gens), mixprobs[1] * pdfplot1(x=sort(gens),a=par11,b=par12) + mixprobs[2] * pdfplot2(x=sort(gens),a=par21,b=par22), type = "l")
  return(as.vector(gens))
}
#######################################
gen.multinorm<<-
function(numofgen = 1, mu = c(0, 0), sigma = diag(2))
{
  choldecomp <- chol(sigma)
  p=length(mu)
  gens <- matrix(0, numofgen, p)
  for(i in 1:numofgen) {
    w <- gen.reject.normal(numofgen = p)
    gens[i, ] <- mu + t(choldecomp) %*% as.vector(w)
  }
  return(gens)
}

gen.wishart<<-
function(numofgen = 1, p = 1, n = 10, sigma = 1)
{
  #default value is a central chi-square with 1 degree of freedom
  gens <- array(0, c(p, p, numofgen))
  for(i in 1:numofgen) {
    w <- matrix(0, p, p)
    for(j in 1:n) {
      x <- gen.multinorm(numofgen = 1, mu = rep(0, p), sigma = sigma)
      w <- w + t(x) %*% x
    }
    gens[, , i] <- w
  }
  return(gens)
}
##############################################
##############################################
########HPDs and CSs
hpd.beta.binomial<<-function(x = 1, a = 1, b = 1,
    n = 1, numofgen = 100, alpha = 0.05)
{
  # x=the observe number of successes, a,b the hyper parameters of the Beta prior
  # n=the number of trials, and alpha affects the coverage we want.
  pdf.beta <- function(x, a = a, b = b)
  {
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

hpd<<-function(x = gen.inverse.transform.gamma(
  1000, 1, 1), pdf = function(x, a = 1, b = 1)
{
  (x^(a - 1) * exp( - x/b))/(b^a * gamma(a))
}
, tol = 1e-006, alpha = 0.05)
{
  numofgen <- length(x)
  p <- sort(x)
  densityatp <- pdf(p)
  sorted <- sort(densityatp)
  percentile <- sorted[floor(numofgen * alpha)]
  roots <- matrix(0, numofgen, 1)
  count <- 0
  for(i in 1:numofgen)
    if(abs(densityatp[i] - percentile) <= tol) {
      print(p[i])
      count <- count + 1
      roots[count] <- p[i]
    }
  roots <- roots[1:count]
  plot(p, densityatp, type = "l")
  abline(v = roots, h = percentile)
  return(list(roots, percentile))
}
#########################################
#Monte Carlo
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

mc.integral<<-
function(x = 1, fun = function(x)
{
  x^0
}
)
{
  # if we pass a matrix it means we have random vectors
  # then first apply the R^p->R^k, function fun() to each vector (by columns)
  # and if the resulting function is k>1 dimensional apply function mean to each coordinate
  if(is.matrix(x)) {
    z <- apply(x, 1, fun)
    if(is.matrix(z)) {
      w <- apply(z, 1, mean)
    }
    else w <- mean(z)
  }
  else {
    w <- mean(fun(x))
  }
  return(w)
}
