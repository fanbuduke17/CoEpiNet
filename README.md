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
2. codes for running inference on real data

## How to run simulations

### Simulate complete event data and estimate parameters

The function `stochastic_coevolve_infer2` in `sim_inference.R` simulates one realization of a temporal network epidemic process (it can be coupled or decoupled) and carried out maximum likelihood estimation or Bayesian estimation for the parameters $\Theta = \{\beta, \gamma, \alpha_{SS},\alpha_{SI},\alpha_{II}, \omega_{SS},\omega_{SI},\omega_{II}\}$.

Another function `rep_stochastic_coevolve_infer` defined in the same file does the "simulate+infer" procedure repeatedly.

Some examples of running simulations are included in `run_sim_inference.R`.

```{r}
siyy <- function(x){
  sum(x^2)
}
```
