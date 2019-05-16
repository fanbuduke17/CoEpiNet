# 05/08/2019
# Truncated Exponential Distribution

# 1. a function to generate samples from Exp(rate) truncated at (a,b)
sam_texp <- function(n, rate, a=0, b=Inf){
  quans = runif(n)
  sams = -(log(exp(-rate*a)-quans*(exp(-rate*a)-exp(-rate*b))))/rate
  return(sams)
}

# # test it against a function found online
# par(mfrow=c(1,2))
# N = 10000; gam = 0.01
# a.vals = c(0:9)
# b.vals = c(5, 8, 10, 25, 36, 100, 300, 1000, 1e4, Inf)
# for(ix in 1:length(a.vals)){
#   a = a.vals[ix]; b = b.vals[ix]
#   hist(TruncatedDistributions::rtexp(N, gam, a, b),
#        main = paste("rate =", gam,"; Truncated within [",a, ",",b,"]"))
#   hist(sam_texp(N, gam, a, b),
#        main = paste("rate =", gam,"; Truncated within [",a, ",",b,"]"))
# }
# # seems to be doing the right thing

# 2. a function to evaluate density of truncated exponential
dtexp <- function(x, rate, a=0, b=Inf, log=F){
  if(any(x < a | x > b)){
    stop("Truncated Exponential: input values not within the bounds!\n")
  }
  dx = rate * exp(-rate * x)/(exp(-rate * a) - exp(-rate * b))
  if(length(x)==1){
    res = dx
  }else{
    res = prod(dx)
  }
  if(log){
    return(log(res))
  }else{
    return(res)
  }
}

# # try it out
# x = sam_texp(20,0.5,1,7)
# log(prod(TruncatedDistributions::dtexp(x,0.5,1,7)))
# dtexp(x,0.5,1,7,log = T)
# # seems to be working!!
