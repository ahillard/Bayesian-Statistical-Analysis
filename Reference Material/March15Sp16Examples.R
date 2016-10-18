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

# solve f(x)=x^2-1

newton.raphson.univariate(f = function(x){x^2-1}, fprime = function(x){2*x}, xo = 1, iterno = 10)
newton.raphson.univariate(f = function(x){x^2-1}, fprime = function(x){2*x}, xo = 11, iterno = 10)
newton.raphson.univariate(f = function(x){x^2-1}, fprime = function(x){2*x}, xo = -11, iterno = 10)

# solve f(x)=ln(x)-1

newton.raphson.univariate(f = function(x){log(x)-1}, fprime = function(x){1/x}, xo = 1, iterno = 10)

# MLE of x, Y~Exp(x), Y is the data

pdf1 <- function(x, y=10 ){x* exp(-y*x)}
pdf1.prime <- function(x, y=10 ){ exp(-y*x)-x*y*exp(-y*x)}
pdf1.doubleprime <- function(x, y=10 ){ -y*exp(-y*x)-y*exp(-y*x)+x*y^2*exp(-y*x)}

newton.raphson.univariate(f = pdf1.prime, fprime = pdf1.doubleprime, xo = .1, iterno = 10)
newton.raphson.univariate(f = pdf1.prime, fprime = pdf1.doubleprime, xo = 1, iterno =10)

# Take a look at the function pdf1.prime
x=seq(0,1,length=300)
plot(x,pdf1.prime(x),type="l")
abline(v=.1,h=0)
