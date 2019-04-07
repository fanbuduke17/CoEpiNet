# 04/06/2019

#### TOY EXAMPLE ####
## dynamically plot network structure

library(Matrix)
library(igraph)

stochastic_coevolve_d <- function(N=200, tmax=100, 
                                  bet=0.05, gam=0.2, init.infec = 1, model = "SIS",
                                  alpha.r=0.005, alpha.d=0.05, quarantine=F,
                                  init.net = NULL, init.p = 0.1,
                                  verbose = F, plot = F, 
                                  d.int=5, seed = 43){
  
  # "d.int": dynamic plot interval (if NULL: no plots)
  set.seed(seed)
  
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
  
  
  ## if plot dynamic networks, set up edgelist to store dynamic network
  # if(!is.null(d.int)){
  #   edgelist = data.frame(from=adjmat@i+1, to=adjmat@j+1, time=0)
  # }
  
  
  ## check whether adjmat is symmetric and nnzero() returns an even number
  if(!isSymmetric(adjmat) | length(adjmat@i)%%2 == 1){
    stop(paste("The initial adjmat isn't symmetric, where nonzero =", 
               nnzero(adjmat),"\n"))
  }
  
  # initialize epidemic
  infected = sample(N, init.infec)
  epid = rep(0,N); epid[infected] = 1
  cat("Initially infected: ",infected,"\n")
  
  
  ## if plot dynamic networks, set up margin
  ## plot time 0
  ## and keep track of plotting times
  if(!is.null(d.int)){
    par(mar=c(1,1,1,2))
    title = paste("Time:",0)
    g.cur = graph_from_adjacency_matrix(adjmat, mode = "undirected")
    lo = layout_in_circle(g.cur)
    vcol = rep("floralwhite",N)
    vcol[epid==1] = "firebrick1"
    vcol[epid==-1] = "deepskyblue2"
    V(g.cur)$color = vcol
    plot(g.cur,vertex.label.color="black", vertex.size=18,
         vertex.label.cex=1, vertex.label.font=2, edge.width=2,
         layout=lo, main=title)
    # legend("bottomright", legend = c("Susceptible","Infected","Recovered"), pch=21,
    #        col="black", 
    #        pt.bg=c("floralwhite","firebrick1","deepskyblue2"), 
    #        pt.cex=1.5, cex=.8, bty="n", ncol=1)
    plotted = 0
  }
  
  
  
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
    
    
    ## plot dynamic networks with approx. d.int interval
    if(!is.null(d.int) & t.cur-plotted >= d.int){
      # edgelist.cur = data.frame(from=adjmat@i+1, to=adjmat@j+1, time=t.int)
      # edgelist = rbind(edgelist, edgelist.cur)
      t.int = round(t.cur)
      title = paste("Time:",t.int)
      g.cur = graph_from_adjacency_matrix(adjmat, mode = "undirected")
      vcol = rep("floralwhite",N)
      vcol[epid==1] = "firebrick1"
      vcol[epid==-1] = "deepskyblue2"
      V(g.cur)$color = vcol
      plot(g.cur,vertex.label.color="black", vertex.size=18,
           vertex.label.cex=1, vertex.label.font=2, edge.width=2,
           layout=lo, main=title)
      plotted = t.cur
    }
    
  }
  
  # return results/summary
  cat("Simulation finished. At the end: \n", 
      "Disease prevalence =", N.I/N, "\n",
      "Network density = ", nnzero(adjmat)/(N*(N-1)),"\n")
  event.log = data.frame(time = time, event = event.type,
                         preval = preval/N, dens = dens/(N*(N-1)))
  
  # summary plot
  if(plot){
    par(mar=c(5,4,4,2)+0.1)
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
  
  ## Final state dynamic plot
  if(!is.null(d.int)){
    # edgelist.cur = data.frame(from=adjmat@i+1, to=adjmat@j+1, time=t.int)
    # edgelist = rbind(edgelist, edgelist.cur)
    title = "At the end"
    g.cur = graph_from_adjacency_matrix(adjmat, mode = "undirected")
    vcol = rep("floralwhite",N)
    vcol[epid==1] = "firebrick1"
    vcol[epid==-1] = "deepskyblue2"
    V(g.cur)$color = vcol
    plot(g.cur,vertex.label.color="black", vertex.size=18,
         vertex.label.cex=1, vertex.label.font=2, edge.width=2,
         layout=lo, main=title)
  }
  
  return(event.log)
}


#### Codes to produce toy example plots ####
# Showcase the toy example
# 1. No quarantine
par(mfrow=c(2,4))
d1 = stochastic_coevolve_d(N=20, alpha.r = .001, alpha.d = .01, 
                           bet = .3, gam = .1, model = "SIR",
                           seed = 43, d.int = 5, verbose = T)
plot.new()
legend(.01,.9, legend = c("Susceptible","Infected","Recovered"), pch=21,
       col="black", 
       pt.bg=c("floralwhite","firebrick1","deepskyblue2"), 
       pt.cex=2.5, cex=2, bty="n", ncol=1)


# 2. Simple quarantine: 
# alpha.d.SI large, alpha.r.SI=0, the rest is same
par(mfrow=c(1,4), mar=c(1,1,2,2))
d2 = stochastic_coevolve_d(N=20, bet = .3, gam = .1, model = "SIR",
                           alpha.r = c(.001,0,.001), alpha.d = c(.01, .2, .01),
                           quarantine = T, seed = 43, 
                           d.int = 2, verbose = T)