# 05/06/2019
# Simulation & Inference
# work with COMPLETE data


library(Matrix)
library(igraph)
library(ggplot2)

# Basic function: simulate full event sequence and do MLE & Bayesian

# do inference every "infer.interval" events
# if "Bayes"=T, specify priors in "priors" as a dataframe

# 05/09/2019
# add "demo.sleep" to pause the plotting for better presentation
stochastic_coevolve_infer <- 
  function(N=50, tmax=100, 
           bet=0.05, gam=0.2, init.infec = 1, model = "SIS",
           alpha.r=0.005, alpha.d=0.05, quarantine=F,
           init.net = NULL, init.p = 0.1,
           verbose = F, plot = F,
           infer.interval = 100, 
           Bayes = F, priors = NULL, samples = 500,
           infer.verbose = T, infer.plot = T, 
           demo.sleep = F, seed = 42){
    ## alpha.* = c(SSrate, SIrate, IIrate) ##
    set.seed(seed)
    
    # check edge rates validity
    if(quarantine & (length(alpha.r) < 3 | length(alpha.d) < 3)){
      stop("Under quarantine, need to specify 3 values for edge connection/disconnection rates respectively!")
    }
    
    # output basic info
    cat("N =", N, " Model:", model,"\n",
        "beta =", bet, " gamma =", gam, "\n",
        "init.infec =", init.infec, " init.p =", init.p, "\n",
        "alpha.r =", alpha.r, "\n", "alpha.d =", alpha.d, "\n",
        "quarantine?", quarantine,"\n",
        "Bayesian inference?", Bayes, "\n",
        "inference interval =", infer.interval, "\n")
    
    # summarization stats to be recorded
    time = NULL
    event.type = NULL
    preval = NULL
    dens = NULL
    # add involved individual labels 
    per1 = NULL
    per2 = NULL
    
    # more stats to record for inference
    if(quarantine){
      event.counts = numeric(8)
      big.sums = numeric(8)
    }else{
      event.counts = numeric(4)
      big.sums = numeric(4)
    }
    tot.event = 0
    
    if(Bayes){
      if(quarantine & nrow(priors)==4){
        a.pr = c(priors$count[1:2],rep(priors$count[3:4],each=3))
        avg.pr = c(priors$avg[1:2],rep(priors$avg[3:4],each=3))
        b.pr = a.pr/avg.pr
      }else{
        a.pr = priors$count
        b.pr = a.pr/priors$avg
      }
      if(quarantine){
        vars = c("beta","gamma",
                 "alpha.r.SS","alpha.r.SI","alpha.r.II",
                 "alpha.d.SS","alpha.d.SI","alpha.d.II")
      }else{
        vars = c("beta","gamma","alpha.r","alpha.d")
      }
    }else{
      MLE.tab = NULL
    }
    params = c(bet, gam, alpha.r, alpha.d)
    
    
    # initialize ER(p) graph if the initial netwok is not provided
    if(is.null(init.net)){
      init.net = sample_gnp(n = N, p = init.p, directed = F)
    }
    # convert to adjacency matrix
    adjmat = get.adjacency(init.net, type=c("both"), names=F)
    adjmat = as(adjmat,"dgTMatrix")
    # save the initial network structure (adj.mat)
    init.adj = adjmat
    
    
    ## check whether adjmat is symmetric and nnzero() returns an even number
    if(!isSymmetric(adjmat) | length(adjmat@i)%%2 == 1){
      stop(paste("The initial adjmat isn't symmetric, where nonzero =", 
                 nnzero(adjmat),"\n"))
    }
    
    # initialize epidemic
    infected = sample(N, init.infec)
    epid = rep(0,N); epid[infected] = 1
    cat("Initially infected: ",infected,"\n")
    # save initial cases
    init.infected = infected
    
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
        # still update inference statistics
        del.t = t.next - t.cur
        z = as.integer(z)
        event.counts[z] = event.counts[z] + 1
        tot.event = tot.event + 1
        
        big.sums = big.sums + c(N.SI, N.I, N.d, N.c) * del.t
        
        # still have to record the event
        time = c(time, t.next)
        event.type = c(event.type, z)
        preval = c(preval, N.I)
        dens = c(dens, nnzero(adjmat))
        # here: a recovery must be the latest event
        per1 = c(per1, recover)
        per2 = c(per2, NA)
        
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
      
      # 3. update inference statistics and do inference occasionally
      # 3.1 update stats
      del.t = t.next - t.cur
      z = as.integer(z)
      event.counts[z] = event.counts[z] + 1
      tot.event = tot.event + 1
      
      big.sums = big.sums + c(N.SI, N.I, N.d, N.c) * del.t
      
      # 3.2 inference if ...
      if(tot.event %% infer.interval == 0){
        if(Bayes){
          # Bayesian inference: Gamma sampling
          a.post = a.pr + event.counts
          b.post = b.pr + big.sums
          cat("With", tot.event,"events, posterior means and variances:\n",
              a.post/b.post, "\n",
              a.post/(b.post^2), "\n")
          if(infer.plot){
            par(mfrow=c(2,2))
            for(ix in 1:length(params)){
              v = vars[ix]; p = params[ix]
              sams = rgamma(samples, shape=a.post[ix], rate = b.post[ix])
              xl = range(sams, p)
              hist(sams, main=paste("Posterior samples for",vars[ix]),
                   xlab = vars[ix], xlim=xl, col="lightblue")
              abline(v=p,col="red")
              if(demo.sleep & ix%%4 == 0){
                Sys.sleep(0.2)
              }
            }
          }
        }else{
          # MLE
          MLEs = event.counts/big.sums
          if(infer.verbose){
            cat("With",tot.event, "events, MLEs:","\n", MLEs, "\n")
          }
          MLE.tab = rbind(MLE.tab, c(tot.event, MLEs))
        }
      }
      
      
      
      # 4. report the event and summarize and set t.cur
      
      t.cur = t.next
      
      # record the event
      time = c(time, t.cur)
      event.type = c(event.type, z)
      preval = c(preval, N.I)
      dens = c(dens, nnzero(adjmat))
      # record labels of individual(s) involved
      if (z >= 3){
        per1 = c(per1, from)
        per2 = c(per2, to)
      }else if (z==2){
        per1 = c(per1, recover)
        per2 = c(per2, NA)
      }else{
        per1 = c(per1, new.infec)
        per2 = c(per2, NA)
      }
      
      
      if(verbose){
        cat(paste0("At time ", t.cur, ", "),msg,"\n")
      }
      
    }
    
    # return results/summary
    cat("Simulation finished. At the end: \n", 
        "Disease prevalence =", N.I/N, "\n",
        "Network density = ", nnzero(adjmat)/(N*(N-1)),"\n")
    event.log = data.frame(time = time, event = event.type,
                           per1 = per1, per2 = per2,
                           preval = preval/N, dens = dens/(N*(N-1)))
    
    # inference at the end
    if(Bayes){
      # Bayesian inference: Gamma sampling
      a.post = a.pr + event.counts
      b.post = b.pr + big.sums
      cat("Posterior means and variances with", tot.event, "events in total:\n",
          a.post/b.post, "\n",
          a.post/(b.post^2), "\n")
      cat("Event type counts:\n", event.counts,"\n")
      if(infer.plot){
        par(mfrow=c(2,2))
        for(ix in 1:length(params)){
          v = vars[ix]; p = params[ix]
          sams = rgamma(samples, shape=a.post[ix], rate = b.post[ix])
          xl = range(sams, p)
          hist(sams, main=paste("Posterior samples for",vars[ix]),
               xlab = vars[ix], xlim=xl, col="lightblue")
          abline(v=p,col="red")
          if(demo.sleep & ix%%4 == 0){
            Sys.sleep(0.2)
          }
        }
      }
    }else{
      # MLE
      MLEs = event.counts/big.sums
      cat("MLEs with", tot.event, "events in total:","\n", MLEs, "\n")
      cat("Event type counts:\n", event.counts,"\n")
      
      MLE.tab = rbind(MLE.tab, c(tot.event, MLEs))
      MLE.tab = as.data.frame(MLE.tab)
      if(quarantine){
        names(MLE.tab) = c("event.num","beta","gamma",
                           "alpha.r.SS","alpha.r.SI","alpha.r.II",
                           "alpha.d.SS","alpha.d.SI","alpha.d.II")
      }else{
        names(MLE.tab) = c("event.num","beta","gamma","alpha.r","alpha.d")
      }
    }
    
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
    
    if(infer.plot){
      #params = c(bet, gam, alpha.r, alpha.d)
      if(Bayes){
        par(mfrow=c(1,1))
      }else{
        vars = names(MLE.tab)[-1]
        e.num = MLE.tab$event.num
        if(quarantine){
          par(mfrow=c(2,2))
        }else{
          par(mfrow=c(2,2))
        }
        for(ix in 1:length(vars)){
          v = vars[ix]
          s.m = mean(MLE.tab[[v]], na.rm=T); 
          s.sd = sd(MLE.tab[[v]], na.rm=T)
          yl = range(s.m - 3*s.sd, s.m + 3*s.sd, params[ix])
          if(all(is.na(MLE.tab[[v]]))){
            plot.new()
          }else{
            plot(MLE.tab[[v]]~e.num, xlab="# events", ylab=v,
                 main = paste("MLE for", v, "v.s. # events"), 
                 type="l", ylim=yl)
          }
          abline(h=params[ix], col="red")
        }
        par(mfrow=c(1,1))
      }
    }
    
    # return event data
    return(list(G0 = init.adj, I0 = init.infected,
                events = event.log))
  }





#######################
####### trials: #######
#######################
# 1. MLEs
EL = stochastic_coevolve_infer(N=100, model="SIR", seed=53) 
# disease-free in the end (~33days), fairly accurate
# but perhaps there are way too many events (>1600)
EL = stochastic_coevolve_infer(N=100, bet=0.03, gam=0.15, model="SIR", 
                               init.p = .07, alpha.r = .003,
                               seed=83) 
# disease-free after ~49 days, reasonably accurate
# event types: 25 26 733 747 

EL = stochastic_coevolve_infer(bet=0.1, quarantine=T,seed=83,
                               alpha.r = c(.005, 0.001, .005),
                               alpha.d = c(.05, 0.1, .05))
# still disease-free, not terrible
EL = stochastic_coevolve_infer(N=100, bet=0.1, gam = 0.1, quarantine=T,
                               alpha.r = c(.005, 0.001, .005),
                               alpha.d = c(.05, 0.1, .05),
                               seed=83)
# the counts are crazy... so not sure this is a good reference
EL = stochastic_coevolve_infer(N=100, bet=0.03, gam = 0.15, model = "SIR",
                               quarantine=T,
                               alpha.r = c(.005, 0.001, .005),
                               alpha.d = c(.05, 0.1, .05),
                               seed=83)
# only need a few infections/recoveries to get beta and gamma
# but can get zero edge events if temporary N.I/N.S tiny
# event counts: 26 27 470 17 4 479 88 8

# 2. Bayesian: w/ Gamma priors
pr = data.frame(count = rep(1,4), avg = c(0.05, 0.1, 0.005, 0.05))
EL = stochastic_coevolve_infer(seed=53, Bayes = T, priors = pr) 
# looks pretty good
#############
#############
EL = stochastic_coevolve_infer(N=100, bet=0.03, gam=0.15, model="SIR", 
                               init.p = .07,alpha.r = .003,
                               seed=71,demo.sleep = T,
                               Bayes = T, priors = pr, samples = 1000)
# results are a bit different than the MLE
# because of the posterior sampling
# ~ 1000 events after ~33 days (disease-free)
# event types: 37 38 468 536
#############
#############
EL = stochastic_coevolve_infer(N =100, bet=0.03, gam=0.15, 
                               quarantine=T,seed=83,
                               alpha.r = c(.005, 0.001, .005),
                               alpha.d = c(.05, 0.1, .05),
                               Bayes = T, priors = pr)
# somehow alpha.d.SI not very good (under-estimated), 
# but there are 133 events???
##########
##########
EL = stochastic_coevolve_infer(N=100, bet=0.03, gam = 0.15, model = "SIR",
                               quarantine=T,
                               alpha.r = c(.005, 0.001, .005),
                               alpha.d = c(.05, 0.1, .05),
                               seed=83,
                               Bayes = T, priors = pr,
                               demo.sleep = T)
# alpha.r.II has event count 0, but luckily the prior is good (!)
# still: only need very few events to get accurate estimation
# ~ 1600 events after ~44 days
# event types: 15 16 757 7 0 795 60 2
##########
##########
