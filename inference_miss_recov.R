# 05/16/2019
# Inference on partially observed data
# --Missing (some) recovery times 

# preparations: loading stuff
setwd("~/Documents/EpiNet/")

source("./inference_util.R")

# example data
# dats = readRDS("~/Documents/EpiNet_coupled_1.rds")
miss1 = readRDS("~/Documents/miss_recov.rds")

# function to do Bayesian inference on data w/ recovery times missingness
# output/plot results every `output.sams` recorded samples
# priors: data frame of vars `count` & `avg`
# assume "quarantine" case!
infer_miss_recov <- function(dats, priors, init.params = NULL,
                             verbose = T, plot = T, output.sams = 100, 
                             samples = 1000, burn = 100, thin = 1,
                             impute = "filter", model="SIR", 
                             true.params = c(.03, .15, .005, 0.001, .005, .05, 0.1, .05),
                             timing = T, seed=42){
  if(timing){ time.st = Sys.time()}
  
  set.seed(seed)
  
  # preparations
  G0 = dats$G0; I0 = dats$I0; events = dats$events
  reports = dats$report; report.times = dats$report.times
  
  if(plot){ par(mfrow=c(2,2)) }
  
  # get time intervals to operate on, and those who need exact recovery times imputed
  MR = get_miss_recov(reports, report.times, events)
  recov.persons = unlist(MR$recover)
  intervals = MR$intervals
  
  # get neighborhood info for all infection cases at their infection times
  nei_infec = get_nei_infection(G0, events, reports, report.times)
  
  # initialize param values
  if(is.null(init.params)){
    # draw from priors
    if(nrow(priors)==4){
      a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
      avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
      b.pr = a.pr/avg.pr
    }else{
      a.pr = priors$count
      b.pr = a.pr/priors$avg
    }
    params.cur = rgamma(8, shape = a.pr, rate = b.pr)
  }else{
    params.cur = init.params
  }
  
  # parameter values storage
  params = matrix(ncol=8, nrow=samples)
  vars = c("beta","gamma",
           "alpha.SS","alpha.SI","alpha.II",
           "omega.SS","omega.SI","omega.II")
  colnames(params) = vars
  
  # run iterations
  S = samples * thin + burn
  for(it in 1:S){
    
    if(verbose){ cat("\nIteration",it,"..\n") }
    
    # (1) propose recovery times
    gam.cur = params.cur[2]
    times = NULL
    for(ix in 1:length(MR$recover)){
      recovs = MR$recover[[ix]]
      lb = intervals$lb[ix]; ub = intervals$ub[ix]
      if(impute == "filter"){
        imputed = propose_recov_filter(lb, ub, recovers = recovs, events, nei_infec, gam = gam.cur)
      }else{
        imputed = propose_recov_rej(lb, ub, recovers = recovs, events, nei_infec, gam = gam.cur)
      }
      times = c(times, imputed)
    }
    recover.dat = data.frame(time = times, per1 = recov.persons)
    
    #print(recover.dat)

    
    if(verbose){ cat("Recovery times imputation done.\n") }
    
    # (2) compute event counts and summations
    PA = parse_augment(G0, I0, events, recover.dat)
    event.counts = PA$event.counts
    big.sums = PA$big.sums
    
    # (3) sample params
    a.post = a.pr + event.counts
    b.post = b.pr + big.sums
    
    params.cur = rgamma(8, shape = a.post, rate = b.post)
    
    if(verbose){ 
      cat("Parameter values sampled:",params.cur,"\n")
    }
    
    ## record this sample after burn-in and thinning
    if((it - burn) %% thin == 0){
      s = (it - burn)/thin
      params[s,] = params.cur
      
      ### make plots periodically
      if(plot & s %% output.sams == 0){
        for(ix in 1:8){
          v = vars[ix]
          sams = params[1:s,ix]
          if(!is.null(true.params)){
            truth = true.params[ix]
            xl = range(sams, truth)
            hist(sams, main=paste("Posterior samples for",v),
                 xlab = v, xlim=xl, col="lightblue")
            abline(v=truth,col="red")
          }else{
            xl = range(sams)
            hist(sams, main=paste("Posterior samples for",v),
                 xlab = v, xlim=xl, col="lightblue")
          }
        }
      }
    }
    
  }
  
  # traceplot
  for(ix in 1:8){
    v = vars[ix]
    sams = params[,ix]
    if(!is.null(true.params)){
      truth = true.params[ix]
      xl = range(sams, truth)
      plot(sams~c(1:samples), main=paste("Traceplot for",v), 
           xlab = "sample", ylab = v, type="l")
      abline(h=truth,col="red")
    }else{
      xl = range(sams)
      plot(sams~c(1:samples), main=paste("Traceplot for",v), 
           xlab = "sample", ylab = v, type="l")
    }
  }
  
  # timing it
  if(timing){ 
    time.en = Sys.time(); 
    cat("\n\n")
    print(time.en - time.st)
    #cat("\n\nTime spent (in seconds):",time.en - time.st,"\n")
  }
  
  return(params)
}


# try it out
pr = data.frame(count = rep(1,4), avg = c(0.05, 0.1, 0.005, 0.05))
## 1. use "filter"
pdf("./miss1_filter.pdf")
res.fil = infer_miss_recov(miss1, priors = pr, output.sams = 100, 
                           samples = 1000, burn = 100)
# time: 13.75 mins
dev.off()
## 2. use "reject"
pdf("./miss1_reject.pdf")
res.rej = infer_miss_recov(miss1, priors = pr, output.sams = 100, 
                           samples = 1000, burn = 100,
                           impute = "rej")
# time: 11.27 mins
dev.off()

