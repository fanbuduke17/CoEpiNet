
# CoEpiNet
Codes for _Likelihood-based Inference for Partially Observed Epidemics on Dynamic Networks_

## Things included
1. simulation codes (both SIR- and SIS-type) epidemic processes on adaptive networks
2. inference codes on complete data
3. inference codes on partially observed data w/ missing recovery times
4. codes for generating visualizations of toy examples
5. example datasets from simulations

## Things **not** included
1. real data (proprietary)
2. codes for inference on real data

## How to run simulations

### Simulate complete event data and estimate parameters

The function `stochastic_coevolve_infer2` in `sim_inference.R` simulates one realization of a temporal network epidemic process (it can be coupled or decoupled) and carried out maximum likelihood estimation or Bayesian estimation for the parameters.

Another function `rep_stochastic_coevolve_infer` defined in the same file does the "simulate+infer" procedure repeatedly.

Some examples of running simulations and complete data inference are included in `run_sim_inference.R`:

 - Simulations with N=100 people, starting with an ER(100, 0.1) network
 - Simulations with N=100 people, starting with a "hubnet" (one person is connected to all, and the rest form an ER(N-1, 0.1) network)
 - Simulations with N=500 people, starting with an ER(500, 0.1) network


### Inference from partially observed epidemic events

The function

The run times of the Chewbacca algorithm and the rejection sampler for data augmentation are compared using the following codes

```{r}
source("./inference_util.R)
bp_res =
bench::press(
  ix = c(1:length(recovers)),
  {
    bench::mark(
      length(propose_recov_rej(lb=intervals$lb[ix],ub = intervals$ub[ix], recovers = recovers[[ix]],
                                  events = miss1$events, nei_infec = nei_infec_miss)),
      length(propose_recov_filter(lb=intervals$lb[ix],ub = intervals$ub[ix], recovers = recovers[[ix]],
                           events = miss1$events, nei_infec = nei_infec_miss))
    )
  }
 )
}
```
