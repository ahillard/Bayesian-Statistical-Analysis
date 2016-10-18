########################HW4
##############problem 6)
#theoretical average should be N
gend1<<-function(N,numofgen=100)
  {
    gens <- matrix(0, numofgen, 1)
    u <- runif(numofgen)
    for(i in 1:numofgen)
    {
        if(u[i]<1/3)

          if(missing(N))
            #use the letter N if not passed
            gens[i] <- "N-1"
          else#the number was passed, use it
            gens[i] <- N-1

        else if(u[i]<2/3)
          if(missing(N))
            #use the letter N if not passed
            gens[i] <- "N"
          else#the number was passed, use it
            gens[i] <- N
        else
          if(missing(N))
            #use the letter N if not passed
            gens[i] <- "N+1"
          else#the number was passed, use it
            gens[i] <- N+1
    }
    return(as.vector(gens))
}
#calls, use whatever N you want
#gend1()

################problem 7)
gendie=function(L=10000)
{
  #use gen.inverse.transform.discrete.finite
  gens=gen.inverse.transform.discrete.finite(
    x=1:12,probs=rep(1/12,12),numofgen=L)
  print(table(gens))
  cat(paste("Expected freq=",L/12))
}
#gendie()


################problem 8)
fun1=function(x){exp(-x/2)}
fun2=function(x){x/(1+x^2)^2}
fun3=function(x){exp(x+x^2)}
fun4=function(x){exp(-x^2)}
fun5=function(x){x/gamma(x)}
#either transform to (0,1) and use
#mc.integral.unif or make the density of a
#rv appear, sample from it, and use mc.integral
#a)this one escapes
#[0,+Inf)
fun11=function(y){1/y^2*(fun1(1/y-1+0))}
#(-Inf,0], this goes off
fun12=function(y){1/y^2*(fun1(-1/y+1+0))}
mc.integral.unif(numofgen=50000, p = 1,
 fun=fun11)+mc.integral.unif(numofgen=50000,
p = 1,fun=fun12)
#answer=Inf
#b)
fun21=function(y){1/y^2*(fun2(1/y-1+0))}
mc.integral.unif(numofgen=50000, p = 1,
                 fun=fun21)
#around 0.5023448
#c)
mc.integral(x=runif(50000,-2,2),fun=fun3)
#around 23.07489
#d) split it up or use biv normal
fun41=function(y){1/y^2*(fun4(-1/y+1+0))}
fun42=function(y){1/y^2*(fun4(-1/y+1+1))}
mc.integral.unif(numofgen=50000, p = 1,
  fun=fun41)*mc.integral.unif(numofgen=50000,
  p = 1,fun=fun42)
#around 1.438625
#e)
mc.integral(x=runif(50000,1,10),fun=fun5)
#around 0.5653258

################problem 9)
#a)
MLEExp=function(x=rexp(20,10))
{#the MLE is simply 1/X-bar
  cat(paste("MLE is=",1/mean(x)))
}
#MLEExp()
#MLE is= 9.4896336854516
#true value is 10, so with n=20 (small)
#you won't get too close to 10
#e.g., try this MLEExp(x=rexp(20000000,10))

#b)theta|x's~Gamma(n+a,(nX-bar+1/b)^-1)
MAPExp=function(x=rexp(10,5),a=2,b=3,
                x0 = 10, iterno = 100)
{
#the MAP is simply (n+a-1)/(n*X-bar+1/b)
#but we also use Newton-Raphson below
  n=length(x)
  xbar=mean(x)
  MAPest=(n+a-1)/(n*xbar+1/b)
  loglik_theta=function(theta)
  {
    (n+a-1)*log(theta)-theta*(n*xbar+1/b)
  }
  f_theta=function(theta)
  {
    (n+a-1)/theta-(n*xbar+1/b)
  }
  f_theta_prime=function(theta)
  {
    -(n+a-1)/theta^2
  }
  ret=newton.raphson.univariate(f =f_theta,
      fprime =f_theta_prime, xo = x0, iterno = iterno)
  print(ret)
  cat(paste("\nMAP from Newton-Raphson is=",ret[1]))
  cat(paste("\nMAP from exact calculations is=",MAPest))
  #it helps to see how bad it performs for n=10
  #if you plot the posterior
  gens=sort(rgamma(100000,
          shape=n+a,scale=(n*xbar+1/b)^-1))
  densgamma=dgamma(gens,shape=n+a,scale=(n*xbar+1/b)^-1)
  plot(gens,densgamma,
       xlab="x",ylab="pdf",main="",type="l")
  title("Gamma pdf")
  abline(v = c(ret[1],MAPest), h = max(densgamma))
#  cat(paste("\nshape=",n+a,", scale=",n*xbar+1/b))
}
#MAPExp()
#if you start it at a bad x0, N-R won't work
#i.e., won't converge, easy to crush it
#try different bad x0
#MAPExp(x0=12,iterno=100)

################problem 10)
#Xi~Bi(m,p),p~Beta(a,b)
#p|x's~Beta(a+sum(x),m*n-sum(x)+b)
BetaBinomProb1=function(x=rbinom(20,5,.45),
  m=5,n=20,a=2,b=2,alpha=0.05,L=20000,p0=.5)
{
  #a)plot prior and posterior
  windows()
  gridx=seq(0,1,length=1000)
  densprior=dbeta(gridx,a,b)
  denspost=dbeta(gridx,a+sum(x),m*n-sum(x)+b)
  ymax=max(densprior,denspost)
  plot(gridx,densprior,type="l",lty=2,xlab="p",
    ylab="density",main="Solid=Post, Dashed=Prior",
    ylim=c(0,ymax))
  points(gridx,denspost,type="l",lty=1)
  #b)either approximate them or get the true values
  postgens=rbeta(L,a+sum(x),m*n-sum(x)+b)
  CSapprox=quantile(postgens,c(alpha/2,
    (1-alpha/2)))
  #same as this
  #CSapprox=c(postgens[floor(L*alpha/2)],
  # postgens[floor(L*(1-alpha/2))])

  CSexact=qbeta(c(alpha/2,
    (1-alpha/2)),a+sum(x),
    m*n-sum(x)+b) # the true values
  cat(paste("\nApproximate CS=[",CSapprox[1],
            ",",CSapprox[2],"]"))
  cat(paste("\nExact CS=[",CSexact[1],
            ",",CSexact[2],"]"))
  abline(v=CSexact)
  #c) for the HPD use the routine hpd.beta.binomial
  cat("\n")
#  hpd.beta.binomial(x = sum(x), a = a, b = b,
#    n = m*n, numofgen = L,alpha = alpha)
  windows()
  p=sort(postgens)
  densityatp=dbeta(p,a+sum(x),m*n-sum(x)+b)
  sorteddens=sort(densityatp)
  percentile=sorteddens[floor(L * alpha)]
  roots=matrix(0, L, 1)
  count=0
  for(i in 1:L)
    if(abs(densityatp[i] - percentile) <= 0.001) {
      print(p[i])
      count=count + 1
      roots[count]=p[i]
    }
  roots=roots[1:count]
  plot(p, densityatp, type = "l")
  abline(v = roots, h = percentile)
  #d) Bayes factor, Ho:p=p0 vs H1:p~=p0
  #from our notes: B=f(x|p0)/m1(x)
  #compute the density first
  fx_given_p0=1
  for(i in 1:n)
    fx_given_p0=fx_given_p0*dbinom(x[i],m,p0)
  #now grab the marginal, first integrate
  # f(x's|p)*pi(p) and get m1_atxs of this form
  m1_atxs=prod(choose(m,x))*
    beta(a+sum(x),m*n-sum(x)+b)/beta(a,b)
  Bfactor=fx_given_p0/m1_atxs
  cat(paste("\nBayes factor=",Bfactor))

}
#e)calls, might have to rerun to pick up
#both roots
#i) BetaBinomProb1() #default values
#Bayes Factor could be > or < 1, and
#it depends on thedata you generated
#since data comes from p=p0=.45, close to .5
#ii) BetaBinomProb1(x=rbinom(10,10,.95),m=10,n=10,a=1,b=1,alpha=0.05,L=20000)
#Bayes Factor should be small since .95 is away from .5
#iii) BetaBinomProb1(x=rbinom(100,20,.5),m=20,n=100,a=2,b=1,alpha=0.05,L=20000)
#Bayes Factor should be huge since data comes
#exactly from p=p0=.5
