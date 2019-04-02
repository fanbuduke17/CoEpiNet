## 03/20/2019
# post processing the results from continuous-time simulations
# make plots of peak prevalence (SIR) and equilibrium prevalence (SIS)
# for various patameter settings


library(dplyr)
library(stringr)
library(ggplot2)
library(scales)

rootdir = "~/Documents/"


##### a function to summarize results from a collection of repeated experiments ####
# different with model (SIR, SIS) and with parameter settings #

#### A COMPROMISE FOR SIS: 
#### calculate the average prevalence after reaching the max preval
#### as a substitute for the equilibrium prevalence

#### Adjustment for SIR: 
#### take the average of the top 20 prevalence values as the max preval

get_preval <- function(path="cont_test/", model = "SIR", par1 = NULL,
                       plot = TRUE, errorbar = TRUE, parname = "beta",
                       equimode = "max"){
  dir = paste0(rootdir, path)
  setwd(dir)
  
  pattern = paste0(model,"\\_.*",par1,"\\_.*rds")
  dats = list.files(pattern = pattern)
  summ = NULL
  
  if(length(dats)==0){
    stop("Queried files not found!")
  }
  
  # read all the files and summarize
  for(dat in dats){
    # print info (just to make sure)
    cat("Processing file",dat,"......")
    # extract the last param value
    s = str_sub(dat,end = -5)
    s.l = unlist(str_split(s,"\\_"))
    v = as.numeric(s.l[length(s.l)])
    # read the file
    dat = readRDS(dat)
    if(model == "SIS"){
      if(equimode == "max"){
        # use the values after the peak point to take average
        equi.res = dat %>% group_by(REP) %>% 
          summarize(equi.prev = mean(preval[which.max(preval)[1]:n()])) %>%
          summarise(avg.equi = mean(equi.prev), sd.equi = sd(equi.prev))
      }else{
        # use the values in the last 20% observations to take average
        tmax = max(dat$time)
        equi.res = dat %>% group_by(REP) %>% 
          filter(time >= 0.8*max(time)) %>%
          summarise(equi.prev = mean(preval)) %>%
          summarise(avg.equi = mean(equi.prev), sd.equi = sd(equi.prev))
      }
    }else{
      equi.res = dat %>% group_by(REP) %>% 
        summarize(equi.prev = mean(preval[order(preval,decreasing = T)[1:20]])) %>%
        summarise(avg.equi = mean(equi.prev), sd.equi = sd(equi.prev))
    }
    equi.res$val = v
    summ = rbind(summ,equi.res)
    cat("done.\n")
  }
  
  # make a plot
  if(plot){
    if(is.null(par1)){
      cap = paste0("model: ",model)
    }else{
      cap = paste0("model: ",model,", fixed parameter=",par1)
    }
    p = ggplot(data = summ, aes(x=val, y=avg.equi))+
      geom_line(color="coral3") +
      scale_x_continuous(breaks=sort(summ$val), trans = log_trans()) +
      scale_y_continuous(limits = c(0,1),trans = exp_trans())
    if(errorbar){
      p = p +
        geom_errorbar(aes(ymin=avg.equi-sd.equi, ymax=avg.equi+sd.equi), color="deepskyblue2")
    }
    p = p + 
      labs(y = ifelse(model=="SIS","equilibrium prevalence","peak prevalence"), 
           x = parname, caption = cap) +
      theme_bw()
    print(p)
  }
  return(summ)
}


## get plots for noQua results
get_preval(model = "SIR")
get_preval(model = "SIS",equimode = 20)

## try out a couple of Qua results
get_preval(path="cont_test_SIR", model = "SIR", par1=0.001, 
           parname = "alpha.d.SI")
get_preval(path="cont_test_SIS", model = "SIS", par1=0.001, 
           parname = "alpha.d.SI")



##### a function to summarize multiple collections of repeated experiments #####
multi_get_preval <- function(path="cont_test_SIS",model="SIS",
                             par.list = seq(from=0, to=0.009, by=0.001),
                             groupby = "par", equimode = "max"){
  # groupby=="par": group by alpha.r.SI (0~0.009)
  # otherwise: group by alpha.d.SI (0~0.09)
  
  # put data together
  res.all = NULL
  for(par in par.list){
    res = get_preval(path, model, par, FALSE, FALSE, NULL, equimode)
    res$par1 = par
    res.all = rbind(res.all,res)
  }
  
  # plot
  if(groupby=="par"){
    p = ggplot(data=res.all, aes(y=avg.equi,x=val,color=as.factor(par1))) + 
      geom_line() +
      scale_x_continuous(breaks=sort(unique(res.all$val))) +
      labs(caption=paste0("model: ",model), x="alpha.d.SI", 
           y=ifelse(model=="SIS","equilibrium prevalence","peak prevalence"),
           color = "alpha.r.SI") +
      theme_bw()
  }else{
    p = ggplot(data=res.all, aes(y=avg.equi,x=par1,color=as.factor(val))) + 
      geom_line() +
      scale_x_continuous(breaks=par.list) +
      labs(caption=paste0("model: ",model), x="alpha.r.SI", 
           y=ifelse(model=="SIS","equilibrium prevalence","peak prevalence"),
           color = "alpha.d.SI") +
      theme_bw()
  }
  
  # if SIS: re-scale the axes
  if(model == "SIS"){
    p = p +
      scale_y_continuous(trans = exp_trans())
  }
  
  print(p)
  
  return(res.all)
}


# examine the results for SIR and SIS
res.all.SIR = multi_get_preval("cont_test_SIR","SIR")
res.all.SIR2 = multi_get_preval("cont_test_SIR","SIR",groupby = "val")
res.all.SIS = multi_get_preval()
res.all.SIS2 = multi_get_preval(groupby = "val")

# new way to calculate average prevalence at equilibrium
res.all.SIS = multi_get_preval(equimode = 20)
res.all.SIS = multi_get_preval(groupby = "val",equimode = 20)



#### DON'T RUN CODES BELOW ####

# try smoothing for better visual

dat.sel = dat %>% filter(REP==1) %>% dplyr::select(time, preval)

lo = loess(preval ~ time, dat.sel, span=0.3)

plot(preval~time, type="l", data=dat.sel)
lines(predict(lo)~time, col="red", data=dat.sel)
