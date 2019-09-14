# 06/02/2019
# Move all `sim_inference` running commands to this script
# So that things are less messy


setwd("~/Documents/EpiNet/")
source("./sim_inference.R")
pr = data.frame(count = rep(1,4), avg = c(0.05, 0.1, 0.005, 0.05))


#################################
####### single experiment #######
#################################
# 1. MLEs
EL = stochastic_coevolve_infer(N=100, model="SIR", seed=53) 
# disease-free in the end (~33days), fairly accurate
# but perhaps there are way too many events (>1600)

# set.seed(53)
# EL2 = stochastic_coevolve_infer2(N=100, model="SIR") 

EL = stochastic_coevolve_infer(N=100, bet=0.03, gam=0.15, model="SIR", 
                               init.p = .07, alpha.r = .003,
                               seed=83) 
# disease-free after ~49 days, reasonably accurate
# event types: 25 26 733 747 

# 06/02/2019
# generate a decoupled case for
# a coupled inference
set.seed(71)
EL = stochastic_coevolve_infer2(N=100, bet=0.03, gam=0.12, model="SIR", 
                                alpha.r = .003, alpha.d = .07) 
saveRDS(EL, "~/Documents/EpiNet_decoup_2.rds")


EL = stochastic_coevolve_infer(bet=0.1, quarantine=T,seed=83,
                               alpha.r = c(.005, 0.001, .005),
                               alpha.d = c(.05, 0.1, .05))
# still disease-free, not terrible
EL = stochastic_coevolve_infer(N=100, bet=0.1, gam = 0.1, quarantine=T,
                               alpha.r = c(.005, 0.001, .005),
                               alpha.d = c(.05, 0.1, .05),
                               seed=83)
# the counts are crazy... so not sure this is a good reference
set.seed(83)
EL = stochastic_coevolve_infer(N=100, bet=0.03, gam = 0.15, model = "SIR",
                               quarantine=T,
                               alpha.r = c(.005, 0.001, .005),
                               alpha.d = c(.05, 0.1, .05))
# only need a few infections/recoveries to get beta and gamma
# but can get zero edge events if temporary N.I/N.S tiny
# event counts: 26 27 470 17 4 479 88 8

## try the updated function:
set.seed(71)
EL2 = stochastic_coevolve_infer2(N=100, bet=0.1, gam = 0.15, model = "SIR",
                                 quarantine=T,
                                 alpha.r = c(.005, 0.001, .005),
                                 alpha.d = c(.05, 0.1, .05), return.infer = F)
# the updated function seems to be working
# save the dataset: looks interesting (almost everybody got sick once)
# Event type counts: 93 94 601 47 46 536 254 39
# saveRDS(EL2, "~/Documents/EpiNet_coupled_2.rds")

set.seed(73)
EL2 = stochastic_coevolve_infer2(N=100, bet=0.05, gam = 0.15, model = "SIR",
                                 quarantine=T,
                                 alpha.r = c(.005, 0.001, .005),
                                 alpha.d = c(.05, 0.1, .05), return.infer = F)
# Event type counts: 60 61 1151 38 9 1049 208 15 
# saveRDS(EL2, "~/Documents/EpiNet_coupled_3.rds")

set.seed(83)
EL2 = stochastic_coevolve_infer2(N=100, bet=0.03, gam = 0.12, model = "SIR",
                                 quarantine=T,
                                 alpha.r = c(.005, 0.001, .005),
                                 alpha.d = c(.05, 0.1, .05), return.infer = T)
# a somewhat unfriendly example
# Event type counts: 25 26 665 10 1 590 119 9 
# saveRDS(EL2, "~/Documents/EpiNet_coupled_baddie.rds")



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




## try the updated function too:
set.seed(73)
EL3 = stochastic_coevolve_infer2(N=100, bet=0.05, gam = 0.15, model = "SIR",
                                 quarantine=T,
                                 alpha.r = c(.005, 0.001, .005),
                                 alpha.d = c(.05, 0.1, .05), 
                                 Bayes = T, priors = pr,
                                 return.infer = F, samples = 1000)
# after repeating the above command for a few times:
# Event type counts: 53 54 904 30 9 887 176 10 
# saveRDS(EL3, "~/Documents/EpiNet_coupled_4.rds")

# in some other try:
# Event type counts: 59 60 781 32 12 737 209 18 
# big.sums: 
# 1248.997 419.1147 161640 33239.59 2748.679 14592.05 2005.444 375.3079 
# saveRDS(EL3, "~/Documents/EpiNet_coupled_7.rds")

# go back to the old function
set.seed(73)
EL5 = stochastic_coevolve_infer(N=100, bet=0.03, gam = 0.12, model = "SIR",
                                quarantine=T,
                                alpha.r = c(.005, 0.001, .005),
                                alpha.d = c(.05, 0.1, .05), 
                                Bayes = T, priors = pr,
                                return.infer = F, samples = 1000)
# first with bet = .04 (?)
# Event type counts: 68 69 338 28 24 342 137 22 
# saveRDS(EL5, "~/Documents/EpiNet_coupled_11.rds")
# Event type counts: 52 53 700 22 6 658 137 9 
# saveRDS(EL5, "~/Documents/EpiNet_coupled_12.rds")
# then with bet = .03
# Event type counts: 26 27 771 15 1 814 89 5 
# saveRDS(EL5, "~/Documents/EpiNet_coupled_13.rds")


set.seed(1003)
EL3 = stochastic_coevolve_infer2(N=100, bet=0.04, gam = 0.15, model = "SIR",
                                 quarantine=T,
                                 alpha.r = c(.005, 0.001, .005),
                                 alpha.d = c(.05, 0.1, .05), 
                                 Bayes = T, priors = pr,
                                 return.infer = F, samples = 1000)
# Event type counts: 35 36 1095 15 2 1035 108 4 
# saveRDS(EL3, "~/Documents/EpiNet_coupled_5.rds")

EL3 = stochastic_coevolve_infer2(N=100, bet=0.03, gam = 0.12, model = "SIR",
                                 quarantine=T,
                                 alpha.r = c(.005, 0.001, .005),
                                 alpha.d = c(.05, 0.1, .05), 
                                 Bayes = T, priors = pr,
                                 return.infer = F, samples = 1000)
# Event type counts: 25 26 747 16 4 777 108 8 
# saveRDS(EL3, "~/Documents/EpiNet_coupled_6.rds")


# 3. Try "hubnet": one person connected to all, all others ER
# 06/02/2019

hubnet = readRDS("~/Documents/hubnet100.rds")

set.seed(71)
pdf("./hubnet100_plots.pdf", width=9, height=6)
EL_hub1 = stochastic_coevolve_infer2(N=100, bet=0.03, gam = 0.12, model = "SIR",
                                     quarantine=T,
                                     alpha.r = c(.005, 0.001, .005),
                                     alpha.d = c(.05, 0.1, .05), return.infer = F)
# Event type counts:
#   21 22 1287 19 0 1297 91 6 

EL_hub2 = stochastic_coevolve_infer2(N=100, bet=0.03, gam = 0.12, 
                                     model = "SIR", quarantine=T,
                                     alpha.r = c(.005, 0.001, .005),
                                     alpha.d = c(.05, 0.1, .05), 
                                     Bayes = T, priors = pr, infer.interval = 1000,
                                     return.infer = F, samples = 1000,
                                     init.net = hubnet)
# Event type counts:
#   44 45 1136 34 9 1097 170 15 
# save this one
saveRDS(EL_hub2, "~/Documents/EpiNet_coupled_hubnet2.rds")

# save this one (or rather, someone along the way)
# saveRDS(EL_hub2, "~/Documents/EpiNet_coupled_hubnet1.rds")

dev.off()



# 4. Try with larger social networks
# N = 500, for example


# MLEs
set.seed(73)
pdf("./N500_MLEs.pdf",height = 6, width = 9)
EL_big1 = stochastic_coevolve_infer2(N=500, bet=0.03, gam = 0.12, model = "SIR",
                                 quarantine=T,
                                 alpha.r = c(.005, 0.001, .005),
                                 alpha.d = c(.05, 0.1, .05), return.infer = T)
dev.off()
# save this one 
# MLEs with 50549 events in total:(which is good estimation)
#   0.03006962 0.1232847 
#   0.004999269 0.001043554 0.004855036 0.0499428 0.0994639 0.05106905 
# Event type counts: (pretty much everyone got sick once)
#   497 498 20770 1210 1745 17227 6855 1747 
# saveRDS(EL_big1,"~/Documents/EpiNet_coupled_big_1.rds")


# Bayesian
pdf("./N500_Bayes_2.pdf",height = 6, width = 9)
EL_big2 = stochastic_coevolve_infer2(N=500, bet=0.03, gam = 0.12, model = "SIR",
                                 quarantine=T,
                                 alpha.r = c(.005, 0.001, .005),
                                 alpha.d = c(.05, 0.1, .05), 
                                 Bayes = T, priors = pr,
                                 return.infer = F, samples = 2000,
                                 infer.interval = 10000)
dev.off()
# Posterior means and variances with 53466 events in total:
#   0.02883917 0.1237284 0.00500658 0.0009646837 0.005125877 0.04902018 0.09998641 0.05088475 
# 1.673436e-06 3.074038e-05 1.121113e-09 8.368837e-10 1.431859e-08 1.293245e-07 1.448042e-06 1.533012e-06 
# Event type counts:
#   496 497 22357 1111 1834 18580 6903 1688 


#################################
##### repeated experiments ######
#################################


res1 = rep_stochastic_coevolve_infer(N=100, bet=0.03, gam = 0.12, model = "SIR",
                                     quarantine=T,
                                     alpha.r = c(.005, 0.001, .005),
                                     alpha.d = c(.05, 0.1, .05),
                                     seed=73, return.infer=T,
                                     Bayesian = T, priors = pr, 
                                     demo.sleep = F, rep= 4,
                                     event.thres = 1000,
                                     pdfpath = "./ex1_all",
                                     savedat = F)

res2 = rep_stochastic_coevolve_infer(N=100, bet=0.03, gam = 0.12, model = "SIR",
                                     quarantine=T,
                                     alpha.r = c(.005, 0.001, .005),
                                     alpha.d = c(.05, 0.1, .05),
                                     seed=83, return.infer=T,
                                     Bayesian = T, priors = pr, 
                                     demo.sleep = F, rep= 5,
                                     event.thres = 1000,
                                     pdfpath = "./ex2",
                                     savedat = T, fname = "./exdat2")

# try MLEs
res3 = rep_stochastic_coevolve_infer(N=100, bet=0.03, gam = 0.12, model = "SIR",
                                     quarantine=T,
                                     alpha.r = c(.005, 0.001, .005),
                                     alpha.d = c(.05, 0.1, .05),
                                     seed=83, return.infer=T,
                                     Bayesian = F, priors = pr, 
                                     demo.sleep = F, rep= 5,
                                     event.thres = 1000,
                                     pdfpath = "./ex3",
                                     savedat = T, fname = "./ex3dat")

# Note ex3dat_5:
# Event type ounts: 47 48 621 35 13 573 189 17 
# MLE plots (w/ CIs) saved

# 06/02/2019
# try the `hubnet`
# Bayesian
res_hub = rep_stochastic_coevolve_infer(N=100, bet=0.03, gam = 0.12, model = "SIR",
                                     quarantine=T, init.net = hubnet,
                                     alpha.r = c(.005, 0.001, .005),
                                     alpha.d = c(.05, 0.1, .05),
                                     seed=83, return.infer=T,
                                     Bayesian = T, priors = pr, 
                                     demo.sleep = F, rep= 4,
                                     event.thres = 1000,
                                     pdfpath = "./hubnet2",
                                     savedat = F)

# MLE
res_hub2 = rep_stochastic_coevolve_infer(N=100, bet=0.03, gam = 0.12, model = "SIR",
                                     quarantine=T,
                                     alpha.r = c(.005, 0.001, .005),
                                     alpha.d = c(.05, 0.1, .05),
                                     seed=83, return.infer=T,
                                     Bayesian = F, priors = pr, 
                                     demo.sleep = F, rep= 5,
                                     event.thres = 1000,
                                     pdfpath = "./hubnet_MLEs",
                                     savedat = T, fname = "./hubnet2")
# Note that "hubnet2_5" seems quite okay
# 1543 events in total (end at 33.46)
# Event type counts:
#   47 48 621 35 13 573 189 17
