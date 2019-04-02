# 03/04/2019
# Simulator for continuous-time co-evolution of social networks and epidemics

# 03/09/2019
# Run on server

# 03/18/2019
# Debugged


library(Matrix)
library(igraph)
library(ggplot2)

stochastic_coevolve <- function(N=200, tmax=100, 
                                bet=0.05, gam=0.2, init.infec = 1, model = "SIS",
                                alpha.r=0.005, alpha.d=0.05, quarantine=F,
                                init.net = NULL, init.p = 0.1,
                                verbose = F, plot = F){
  ## alpha.* = c(SSrate, SIrate, IIrate) ##
  # set.seed(seed)
  
  # check edge rates validity
  if(quarantine & (length(alpha.r) < 3 | length(alpha.d) < 3)){
    stop("Under quarantine, need to specify 3 values for edge connection/disconnection rates respectively!")
  }
  
  # output basic info
  cat("N =", N, " Model:", model,"\n",
      "beta =", bet, " gama =", gam, "\n",
      "init.infec =", init.infec, " init.p =", init.p, "\n",
      "alpha.r =", alpha.r, "\n", "alpha.d =", alpha.d, "\n",
      "quarantine?", quarantine,"\n")
  
  # summarization stats to be recorded
  time = NULL
  event.type = NULL
  preval = NULL
  dens = NULL
  
  
  # initialize ER(p) graph if the initial netwok is not provided
  if(is.null(init.net)){
    init.net = sample_gnp(n = N, p = init.p, directed = F)
  }
  # convert to adjacency matrix
  adjmat = get.adjacency(init.net, type=c("both"), names=F)
  adjmat = as(adjmat,"dgTMatrix")
  
  
  ## check whether adjmat is symmetric and nnzero() returns an even number
  if(!isSymmetric(adjmat) | length(adjmat@i)%%2 == 1){
    stop(paste("The initial adjmat isn't symmetric, where nonzero =", 
               nnzero(adjmat),"\n"))
  }
  
  # initialize epidemic
  infected = sample(N, init.infec)
  epid = rep(0,N); epid[infected] = 1
  cat("Initially infected: ",infected,"\n")
  
  # initialize data objects to keep track of
  susceptible = which(epid==0)
  N.I = init.infec; N.S = N - N.I
  t.cur = 0
  if(quarantine){
    N.c = numeric(3); N.d = numeric(3)
    adjmat.ss = adjmat[susceptible,susceptible]
    N.c[1] = nnzero(adjmat.ss)/2
    N.d[1] = N.S * (N.S - 1)/2 - N.c[1]
    adjmat.si = adjmat[susceptible,infected]
    N.c[2] = nnzero(adjmat.si)
    N.d[2] = N.S * N.I - N.c[2]
    adjmat.ii = adjmat[infected,infected]
    N.c[3] = nnzero(adjmat.ii)/2
    N.d[3] = N.I * (N.I - 1)/2 - N.c[3]
  
    N.SI = N.c[2]
  }else{
    N.c = nnzero(adjmat)/2
    N.d = N * (N-1)/2 - N.c
    
    adjmat.si = adjmat[susceptible,infected]
    N.SI = nnzero(adjmat.si)
  }
  
  
  # iterations
  while (t.cur < tmax) {
    # NOTES #
    # a. check to make sure N.S also > 0
    # b. handle edge cases (N.S = 1 and/or N.I = 1)
    
    # # check whether all N's are non-negative
    # Ns = c(N.SI, N.I, N.d, N.c)
    # if(any(Ns < 0) | all(Ns == 0)){
    #   msg = paste("Weird stuff happened! All the Ns =", paste(Ns,collapse = ", "), "\n")
    #   stop(msg)
    # }
    
    # 1. calculate sum of risks and sample the next event time
    Lambda = bet * N.SI + N.I * gam + sum(alpha.r * N.d) + sum(alpha.d * N.c)
    t.next = t.cur + rexp(1, rate = Lambda)
    
    # 2. determine the event type
    ratios = c(bet * N.SI, N.I * gam, alpha.r * N.d, alpha.d * N.c)
    # 2.1 deal with network change
    if(quarantine){
      z = sample(8, 1, prob = ratios)
      ## NEED TO FILL OUT z==3:8 ##
      
      if(z==3){
        # reconnect a random S-S pair
        num.dis = (N.S - 1) - rowSums(adjmat.ss)
        from.ind = sample(N.S,1,prob = num.dis); from = susceptible[from.ind]
        to.cands = which(adjmat.ss[from.ind,]==0)
        to.cands = to.cands[to.cands!=from.ind]
        if(length(to.cands) == 1){
          to.ind = to.cands
        }else{
          to.ind = sample(to.cands,1)
        }
        to = susceptible[to.ind]
        adjmat.ss[as.integer(from.ind),as.integer(to.ind)] = 1
        adjmat.ss[as.integer(to.ind),as.integer(from.ind)] = 1
        adjmat[as.integer(from),as.integer(to)] = 1
        adjmat[as.integer(to),as.integer(from)] = 1
        N.c[1] = N.c[1] + 1; N.d[1] = N.d[1] - 1
        
        msg = paste0("reconnect an S-S pair, (",from,",",to,")")
      }
      if(z==4){
        # reconnect a random S-I pair
        if(N.S == 1 | N.I == 1){
          if(N.S == 1){
            if(N.I == 1){
              # only one possible S-I pair
              from = susceptible; to = infected
            }else{
              # N.S == 1 but N.I != 1
              from = susceptible
              to.ind = sample(which(adjmat.si==0),1); to = infected[to.ind]
            }
          }else{
            # N.S != 1 but N.I == 1
            num.dis = N.I - adjmat.si
            from.ind = sample(N.S,1,prob = num.dis); from = susceptible[from.ind]
            to = infected
          }
        }else{
          # N.S != 1 and N.I != 1
          num.dis = N.I - rowSums(adjmat.si)
          from.ind = sample(N.S,1,prob = num.dis); from = susceptible[from.ind]
          to.ind = sample(which(adjmat.si[from.ind,]==0),1); to = infected[to.ind]
        }
        adjmat[as.integer(from),as.integer(to)] = 1
        adjmat[as.integer(to),as.integer(from)] = 1
        adjmat.si = adjmat[susceptible,infected]
        
        N.c[2] = N.c[2] + 1; N.d[2] = N.d[2] - 1
        N.SI = N.c[2]
        
        msg = paste0("reconnect an S-I pair, (",from,",",to,").")
      }
      if(z==5){
        # reconnect a random I-I pair
        num.dis = (N.I - 1) - rowSums(adjmat.ii)
        from.ind = sample(N.I,1,prob = num.dis); from = infected[from.ind]
        to.cands = which(adjmat.ii[from.ind,]==0)
        to.cands = to.cands[to.cands!=from.ind]
        if(length(to.cands) == 1){
          to.ind = to.cands
        }else{
          to.ind = sample(to.cands,1)
        }
        to = infected[to.ind]
        adjmat.ii[as.integer(from.ind),as.integer(to.ind)] = 1
        adjmat.ii[as.integer(to.ind),as.integer(from.ind)] = 1
        adjmat[as.integer(from),as.integer(to)] = 1
        adjmat[as.integer(to),as.integer(from)] = 1
        N.c[3] = N.c[3] + 1; N.d[3] = N.d[3] - 1
        
        msg = paste0("reconnect an I-I pair, (",from,",",to,").")
      }
      if(z==6){
        # disconnect a random S-S pair
        dis.ind = sample(N.c[1], 1) 
        from.ind = adjmat.ss@i[dis.ind] + 1; to.ind = adjmat.ss@j[dis.ind] + 1
        adjmat.ss[as.integer(from.ind),as.integer(to.ind)] = 0
        adjmat.ss[as.integer(to.ind),as.integer(from.ind)] = 0
        
        from = susceptible[from.ind]; to = susceptible[to.ind]
        adjmat[as.integer(from),as.integer(to)] = 0
        adjmat[as.integer(to),as.integer(from)] = 0
        
        N.c[1] = N.c[1] - 1; N.d[1] = N.d[1] + 1
        
        # ## sometimes from/to is just empty... check it out
        # if(length(from)!=1){
        #   stop("from.ind = ", from.ind, 
        #        " but susceptible length is ", length(susceptible),
        #        " N.c[1] =", N.c[1], " dis.ind =", dis.ind, 
        #        " while adjmat.ss@i has length ", length(adjmat.ss@i))
        # }
        # ##
        
        msg = paste0("disconnect an S-S pair, (",from,",",to,").")
      }
      if(z==7){
        # disconnect a random S-I pair
        dis.ind = sample(N.c[2], 1)
        
        if(N.S == 1 | N.I == 1){
          if(N.S == 1){
            if(N.I == 1){
              # only one possible S-I pair
              from = susceptible; to = infected
            }else{
              # N.S == 1 but N.I != 1
              from = susceptible
              to.ind = sample(which(adjmat.si==1),1); to = infected[to.ind]
            }
          }else{
            # N.S != 1 but N.I == 1
            num.con = adjmat.si
            from.ind = sample(N.S,1,prob = num.con); from = susceptible[from.ind]
            to = infected
          }
          adjmat[as.integer(from),as.integer(to)] = 0
          adjmat[as.integer(to),as.integer(from)] = 0
          adjmat.si = adjmat[susceptible,infected]
        }else{
          # N.S != 1 and N.I != 1
          dis.ind = sample(N.c[2], 1)
          from = adjmat.si@i[dis.ind] + 1; to = adjmat.si@j[dis.ind] + 1
          adjmat.si[as.integer(from),as.integer(to)] = 0
          
          from = susceptible[from]; to = infected[to]
          adjmat[as.integer(from),as.integer(to)] = 0
          adjmat[as.integer(to),as.integer(from)] = 0
        }
        
        N.c[2] = N.c[2] - 1; N.d[2] = N.d[2] + 1
        N.SI = N.c[2]
        
        msg = paste0("reconnect an S-I pair, (",from,",",to,").")
      }
      if(z==8){
        # disconnect a random I-I pair
        dis.ind = sample(N.c[3], 1) 
        from = adjmat.ii@i[dis.ind] + 1; to = adjmat.ii@j[dis.ind] + 1
        adjmat.ii[as.integer(from),as.integer(to)] = 0
        adjmat.ii[as.integer(to),as.integer(from)] = 0
        
        from = infected[from]; to = infected[to]
        adjmat[as.integer(from),as.integer(to)] = 0
        adjmat[as.integer(to),as.integer(from)] = 0
        
        N.c[3] = N.c[3] - 1; N.d[3] = N.d[3] + 1
        
        msg = paste0("disconnect an I-I pair, (",from,",",to,").")
      }
    }else{
      z = sample(4, 1, prob = ratios)
      if(z==3){
        # reconnect a random pair
        num.dis = (N - 1) - rowSums(adjmat)
        from = sample(N,1,prob = num.dis)
        to.cands = which(adjmat[from,]==0)
        to.cands = to.cands[to.cands!=from]
        if(length(to.cands) == 1){
          to = to.cands
        }else{
          to= sample(to.cands,1)
        }
        adjmat[as.integer(from), as.integer(to)] = 1
        adjmat[as.integer(to),as.integer(from)] = 1
        
        N.c = N.c + 1; N.d = N.d - 1
        
        # UPDATE: if any of S or I is empty now
        if(is.null(susceptible) | is.null(infected)){
          adjmat.si = NULL; N.SI = 0
        }else{
          adjmat.si = adjmat[susceptible,infected]
          N.SI = nnzero(adjmat.si)
        }
        
        msg = paste0("reconnect a pair, (",from,",",to,").")
      }
      if(z==4){
        # disconnect a random pair
        dis.ind = sample(N.c, 1)
        from = adjmat@i[dis.ind] + 1; to = adjmat@j[dis.ind] + 1
        adjmat[from,to] = 0; adjmat[to,from] = 0
        N.c = N.c - 1; N.d = N.d + 1
        
        # UPDATE: if any of S or I is empty now
        if(is.null(susceptible) | is.null(infected)){
          adjmat.si = NULL; N.SI = 0
        }else{
          adjmat.si = adjmat[susceptible,infected]
          N.SI = nnzero(adjmat.si)
        }
        
        msg = paste0("disconnect a pair, (",from,",",to,").")
      }
    }
    #2.2 deal with epidemic progress
    if(z==2){
      # recover an infected person
      rec.ind = sample(N.I,1)
      recover = infected[rec.ind]
      infected = infected[-rec.ind]
      N.I = N.I - 1 
      if(model == "SIS"){
        susceptible = c(susceptible,recover)
        epid[recover] = 0
        N.S = N.S + 1
      }else{
        epid[recover] = -1
      }
      
      msg = paste0("individual ", recover,
                   " is recovered. Disease prevalence = ",
                   N.I/N, ".")
    }
    if(z==1){
      # select a susceptible person to infect
      if(N.S == 1){
        # only one person to infect
        new.infec = susceptible
        susceptible = NULL
      }else{
        # multiple susceptible inviduals to choose from
        if(N.I == 1){
          # only one infected person
          num.nei.i = adjmat.si
        }else{
          num.nei.i = rowSums(adjmat.si)
        }
        infec.ind = sample(N.S,1,prob = num.nei.i)
        new.infec = susceptible[infec.ind]
        susceptible = susceptible[-infec.ind]
      }
      epid[new.infec] = 1
      N.S = N.S - 1
      infected = c(infected, new.infec); N.I = N.I + 1
      
      msg = paste0("individual ", new.infec,
                   " is infected. Disease prevalence = ",
                   N.I/N, ".")
    }
    
    # A sanity check: if N.I==0, stop the simulation
    if(N.I==0){
      if(verbose){
        cat(paste0("At time ", t.next, ", "),msg,"\n")
      }
      # still have to record the event
      time = c(time, t.next)
      event.type = c(event.type, z)
      preval = c(preval, N.I)
      dens = c(dens, nnzero(adjmat))
      
      cat("Disease extinct at time", t.next, ", simulation stops.\n")
      break
    }
    
    
    # # check whether adjmat is still symmetric--things have gone wrong somewhere
    # if(!isSymmetric(adjmat) | length(adjmat@i)%%2 == 1){
    #   stop(paste("The adjmat isn't symmetric any more, where nonzero =", 
    #              nnzero(adjmat), "and that's after z =", z,
    #              "from.ind =", from.ind, "to.ind = ", to.ind, 
    #              "to.cands =", paste(to.cands,collapse = ","),
    #              "from =", from, "to =", to,"\n"))
    # }
    
    
    # 2.2.* make updates about S-I, S-S, I-I connections if z==1 or 2
    # while handling edge cases of N.S <= 1 and/or N.I == 1
    if(z %in% c(1,2)){
      if(N.S == 0){
        adjmat.si = NULL
        N.SI = 0
      }else{
        adjmat.si = adjmat[susceptible,infected]
        N.SI = nnzero(adjmat.si)
      }
      if(quarantine){
        if(N.S <= 1){
          adjmat.ss = 0
          N.c[1] = 0
          N.d[1] = 0
        }else{
          adjmat.ss = adjmat[susceptible,susceptible]
          N.c[1] = nnzero(adjmat.ss)/2
          N.d[1] = N.S * (N.S - 1)/2 - N.c[1]
        }
        N.c[2] = N.SI
        N.d[2] = N.S * N.I - N.c[2]
        if(N.I == 1){
          adjmat.ii = 0
          N.c[3] = 0
          N.d[3] = 0
        }else{
          adjmat.ii = adjmat[infected,infected]
          N.c[3] = nnzero(adjmat.ii)/2
          N.d[3] = N.I * (N.I - 1)/2 - N.c[3]
        }
      }
    }
    
    # 3. report the event and summarize and set t.cur
    
    t.cur = t.next
    
    # record the event
    time = c(time, t.cur)
    event.type = c(event.type, z)
    preval = c(preval, N.I)
    dens = c(dens, nnzero(adjmat))
    
    if(verbose){
      cat(paste0("At time ", t.cur, ", "),msg,"\n")
    }
    
  }
  
  # return results/summary
  cat("Simulation finished. At the end: \n", 
      "Disease prevalence =", N.I/N, "\n",
      "Network density = ", nnzero(adjmat)/(N*(N-1)),"\n")
  event.log = data.frame(time = time, event = event.type,
                         preval = preval/N, dens = dens/(N*(N-1)))
  
  if(plot){
    subt = paste(model," beta =",bet," gamma =", gam, 
                 " alpha.r =", alpha.r, " alpha.d =", alpha.d)
    plot(preval ~ time, data = event.log, 
         xlab="Time", ylab="Disease prevalence", 
         main = "Disease Prevalence vs. Time", sub = subt,
         ylim = c(0,1), type="l",lwd=2)
    plot(dens ~ time, data = event.log,
         xlab="Time", ylab="Edge density", 
         main = "Network Edge Density vs. Time", sub = subt,
         type="l",lwd=2)
  }
  
  return(event.log)
}



#### the function to repeat the simulation multiple times ####

rep_stochastic_coevolve <- function(N=200, tmax=100, 
                                    bet=0.05, gam=0.2, init.infec = 1, model = "SIS",
                                    alpha.r=0.005, alpha.d=0.05, quarantine=F,
                                    init.net = NULL, init.p = 0.1,
                                    seed = 42, plot = T, n.sim = 10){
  # repeat `stochastic_coevolve` for n.sim times
  set.seed(seed)
  
  dats = NULL
  for(re in 1:n.sim){
    cat("\nSimulation",re,":\n")
    res = stochastic_coevolve(N = N, tmax = tmax, bet = bet, gam = gam,
                              init.infec = init.infec, model = model,
                              alpha.r = alpha.r, alpha.d = alpha.d, 
                              quarantine = quarantine, 
                              init.net = init.net, init.p = init.p)
    res$REP = re
    dats = rbind(dats, res)
  }
  
  dats$REP = as.factor(dats$REP)
  dats$event = as.factor(dats$event)
  
  # plots
  if(plot){
    cap = paste(model," beta =",bet," gamma =", gam, 
                " alpha.r =", paste(alpha.r,collapse = ","), 
                " alpha.d =", paste(alpha.d,collapse = ","))
    
    if(quarantine){
      br = as.character(c(1:8))
      event.code = c("infection", "recovery", 
                     "S-S recon","S-I recon", "I-I recon",
                     "S-S discon", "S-I discon", "I-I discon")
    }else{
      br = as.character(c(1:4))
      event.code = c("infection", "recovery", "reconnection", "disconnection")
    }
    
    # 1. prevalence
    p.preval = ggplot(data = dats, aes(x=time, y=preval,group = REP))+
      geom_line(size=0.2, color="darkred") + ylim(0,1) + 
      labs(y = "disease prevalence", caption = cap) +
      theme_bw()
    # 2. network density
    p.dens = ggplot(data = dats, aes(x=time, y=dens,group = REP))+
      geom_line(size=0.2) +  
      labs(y = "edge density", caption = cap) +
      theme_bw()
    # 3. event frequency
    p.event = ggplot(data = dats, aes(x=event)) + 
      geom_bar(fill="deepskyblue",aes(y=..count../sum(..count..)), width = 0.2) +
      scale_x_discrete(breaks = br, labels = event.code) +
      labs(x="event type", y="proportion",caption = cap) +
      ggtitle("Average frequencies of event types") + 
      theme_bw()
    
    print(p.preval)
    print(p.dens)
    print(p.event)
  }
  
  return(dats)
}




#### RUN CODES BELOW ON SERVER ####

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
ss = as.numeric(slurm_arrayid)
set.seed(ss)

outdir = "cont_test_SIS/"


## 1. no quarantine; experiment to observe the "baseline"

# ix = ss %% 10 + 1
# 
# beta.vals= seq(from = 0.01, to = 0.1, by = .01)
# BE = beta.vals[ix]
# 
# if(ss <= 10){
#   pdf(paste0(outdir,"SIS_noQua_",BE,".pdf"))
# 
#   dats = rep_stochastic_coevolve(bet = BE, init.infec = 10,
#                                  seed = sample(1000,1), n.sim = 10)
# 
#   saveRDS(dats,file=paste0(outdir,"SIS_noQua_",BE,".rds"))
#   dev.off()
# }else{
#   pdf(paste0(outdir,"SIR_noQua_",BE,".pdf"))
# 
#   dats = rep_stochastic_coevolve(model = "SIR", bet = BE, init.infec = 10,
#                                  seed = sample(1000,1), n.sim = 10)
# 
#   saveRDS(dats,file=paste0(outdir,"SIR_noQua_",BE,".rds"))
#   dev.off()
# }


## 2. with quarantine: fix beta=0.09, try varying alpha.*.SI from .000 to .009

# a change: try running for 1000 time units

# another change: try out huge d.vals

BE = .09

r.vals= seq(from = 0, to = .009, by = .001)

ix = ss %% length(r.vals) + 1
r.rate = r.vals[ix]

#d.vals = seq(from = 0, to = .09, by = .01)
d.vals = seq(from = 0, to = .9, by = .1)

if(ss <= length(r.vals)){
  pdf(paste0(outdir,"SIS_Qua_",r.rate,".pdf"))
  for(iy in 1:length(d.vals)){
    d.rate = d.vals[iy]
    dats = rep_stochastic_coevolve(bet = BE, init.infec = 10,
                                   alpha.r = c(.005, r.rate, .005),
                                   alpha.d = c(.05, d.rate, .05),
                                   quarantine = T,
                                   seed = sample(1000,1), n.sim = 10,
                                   tmax = 1000)
    saveRDS(dats,file=paste0(outdir,"SIS_Qua_",r.rate,"_",d.rate,".rds"))
  }
  dev.off()
}else{
  pdf(paste0(outdir,"SIR_Qua_",r.rate,".pdf"))
  for(iy in 1:length(d.vals)){
    d.rate = d.vals[iy]
    dats = rep_stochastic_coevolve(model = "SIR", bet = BE, init.infec = 10,
                                   alpha.r = c(.005, r.rate, .005),
                                   alpha.d = c(.05, d.rate, .05),
                                   quarantine = T,
                                   seed = sample(1000,1), n.sim = 10,
                                   tmax = 1000)
    saveRDS(dats,file=paste0(outdir,"SIR_Qua_",r.rate,"_",d.rate,".rds"))
  }
  dev.off()
}




#### try it out locally -- outdated ####
# stochastic_coevolve(tmax = 10, verbose = F)
# stochastic_coevolve(tmax = 10, model = "SIR")
# 
# stochastic_coevolve(tmax = 10, 
#                     alpha.r = c(0.008, 0.008, 0.004), 
#                     alpha.d = c(0.03, 0.03, 0.06),
#                     quarantine = T, init.infec = 20)
# 
# stochastic_coevolve(tmax = 8.6, model = "SIR",
#                     alpha.r = c(0.008, 0.008, 0.004), 
#                     alpha.d = c(0.03, 0.03, 0.06),
#                     quarantine = T, init.infec = 50)
