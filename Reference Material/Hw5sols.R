########################HW5
##############problem 5)
gen_trunc_Pois<<-function(numofgen=1000,lambda=5,k=10)
{
  #first get the probs corresponding to this
  #truncated Poisson then generate from the
  #corresponding discrete rv
  x=0:k#support of the truncated Poisson
  truncPoisprobs=dpois(x,lambda=lambda)/ppois(k,lambda=lambda)
  gens=gen.inverse.transform.discrete.finite(
    x=x,probs=truncPoisprobs,numofgen=numofgen)
  hist(gens,main="Truncated Poisson")
  #get the frequency table
  tab=table(gens)
#approximate the mean of the truncated Poisson
  approx_mean=mean(gens)
  cat("\nApproximate mean is",approx_mean)
}
#calls, use whatever lambda you want
#gen_trunc_Pois(numofgen=1000,lambda=5,k=10)
#gen_trunc_Pois(numofgen=1000,lambda=8,k=10)

################problem 6)
dice_ex=function()
{
  #use gen.inverse.transform.discrete.finite
  counts=rep(0,11)#11 possible outcomes,2:12
  iter=1
  notdone=TRUE
  while(notdone)
  #if one of them is still 0 continue the loop
  {
    #generate the outcome of the toss
    die1=sample(x=1:6,size=1,prob=rep(1/6,6)) #same function as gen.inverse.transform.discrete.finite
    die2=sample(x=1:6,size=1,prob=rep(1/6,6)) #same function as gen.inverse.transform.discrete.finite
    roll=die1+die2#a number from 2-12
 #   cat("\n",die1,die2,roll,"\n")
    #see to which count this belongs to
    for(i in 1:11)
      if(i+1==roll)
      {
        counts[i]=counts[i]+1
        break
      }
    if(min(counts)>0)
      notdone=FALSE
    cat("\nIteration",iter,"\n")
    iter=iter+1
  }
  table(counts)
  print(counts)
}
#calls
#dice_ex()
################problem 7)
gen_func1_prob7=function(numofgen=10000)
{
  #use rejection depending on x>0 or x<=0
  #this is a mixture:
  #f(x)=.5*2exp(2x)I(x<=0)+.5*2exp(-2x)I(x>0)
  #so first choose the component, and then just
  #generate an exp(2), this is a Laplace distr
  f<-function(x)
  {
    exp(2*x)*(x<=0)+exp(-2*x)*(x>0)
  }
  gens=rep(0,numofgen)
  for(i in 1:numofgen)
  {
    u=runif(1)
    if(u<.5)#this is - an exp(2) rv
      gens[i]=-rexp(1,2)
    else#this branch is just exp(2)
      gens[i]=rexp(1,2)
  }
  plot(sort(gens),f(x=sort(gens)),main="Laplace pdf",ylab="f(x)",xlab="x",type = "l")
}
#gen_func1_prob7()
gen_func2_prob7=function(numofgen=10000)
{
  #f(x)=30*(x^2-2*x^3+x^4),0<=x<=1
  #use rejection algorithm
  #use as proposal the obvious Uniform(0,1)
  #with density g(x)=I(x is in [0,1])
  #you can either maximize f(x)/g(x)
  #to get c=max{f(x)/g(x)}=max{f(x)}
  #but who has time, we recognize first that
  #since 0<=x<=1, f(x)<=30*(1^2-2*0^3+1^4)=60
  #so that choosing c=60 we satisfy the requirement
  #for the rejection algorithm, it will just take
  #us longer to accept a value
  f<-function(x)
  {
    30*(x^2-2*x^3+x^4)
  }
  c=60
  gens=rep(0,numofgen)
  for(i in 1:numofgen)
  {
    u=runif(1)
    Y=runif(1)
    #g(y)=1,c=60
    while(u>=f(Y)/c)
    {
      u=runif(1)
      Y=runif(1)
    }
    gens[i]=Y
  }
  plot(sort(gens),f(x=sort(gens)),ylab="f(x)",xlab="x",type = "l")
#we have a fixed domain here, we could do this too
#  x=seq(0,1,length=100)
#  plot(x,f(x),type="l")
}
#gen_func2_prob7()
gen_func3_prob7=function(L=10000,a = .1, b = 2)
{
  #use inverse transform
  #y=F(x)=1-exp(-a*x^b)->
  #x=F^(-1)(y)=[(-1/a)*(log(1-y))]^(1/b)
  #f(x)=dF(x)/dx=a*b*x^(b-1)*exp(-a*x^b)
  u=runif(L)
  gens=((-1/a)*(log(1-u)))^(1/b)
  pdfWeibull<-function(x,a=a,b=b)
  {
    a*b*x^(b-1)*exp(-a*x^b)
  }
  plot(sort(gens),pdfWeibull(x=sort(gens),a=a,b=b),
       main="Weibull pdf",ylab="pdf",xlab="x",type = "l")
}
#gen_func3_prob7()
################problem 8)
unif_MC_ex=function(L=50000)
{
  u=runif(L)
  cov1=mean(u*sqrt(1-u))-mean(u)*mean(sqrt(1-u))
  cat("\nApproximate cov1 is",cov1)
  cov2=mean(u^2*sqrt(1-u))-mean(u^2)*mean(sqrt(1-u))
  cat("\nApproximate cov2 is",cov2)
  cov3=mean(u*exp(u))-mean(u)*mean(exp(u))
  cat("\nApproximate cov3 is",cov3)
}
#unif_MC_ex()
################problem 9)
pdf.normal<<-function(x,mu=1,sigsq=1)
{(2*pi*sigsq)^(-.5)*exp(-.5*sigsq^(-1)*(x-mu)^2)}
normmix2=function(r=1000,theta=c(.5,.5, 0, 4, 1, 2))
{
  p1=theta[1]
  p2=theta[2]
  if(p1+p2!=1)
    stop("Component probabilities do not sum up to 1.")
  mu1=theta[3]
  mu2=theta[4]
  sig1sq=theta[5]
  sig2sq=theta[6]
  mixpdf=function(x,p1=p1,p2=p2,mu1=mu1,mu2=mu2,
                  sig1sq=sig1sq,sig2sq=sig2sq)
  {
    p1*pdf.normal(x,mu=mu1,sig=sig1sq)+
    p2*pdf.normal(x,mu=mu2,sig=sig2sq)
  }
  gens=rep(0,r)
  for(i in 1:r)
  {
    u=runif(1)
    if(u<p1)#this is the first normal rv
      gens[i]=rnorm(1,mean = mu1, sd = sqrt(sig1sq))
    else#this is the second normal rv
      gens[i]=rnorm(1,mean = mu2, sd = sqrt(sig2sq))
  }
  x=sort(gens)
  pdf=mixpdf(x=x,
             p1=p1,p2=p2,mu1=mu1,mu2=mu2,
             sig1sq=sig1sq,sig2sq=sig2sq)
  plot(x,pdf,
    main="Univariate mixture of 2 normal components",ylab="f(x)",xlab="x",type = "l")
  library(ggplot2)
  print(ggplot2::qplot(x, pdf, geom = "path")+
          ggplot2::theme_classic() +
          ggplot2::theme(panel.border =
    ggplot2::element_rect(fill = NA, size = 1)))
  return(as.vector(gens))
}
#normmix2()
#normmix2(theta=c(.2,.8, 1, 5, 1, 1))
#normmix2(theta=c(.1,.9, 3, 4, .3, .3))
