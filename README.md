 Monte Carlo EM 
 ====

# Introduction

* This project described a Monte Carlo EM (**MCEM**) method to derive Maximum Likelihood Estimates (**MLE**) of the log-likelihood function. 

* In the E-step, perform *K = 500* Gibbs sampling incorporated with a Metropolis-Hastings step, and drop the first *100* as a burn-in procedure.

* Read article: *[Maximum Likelihood Algorithms for Generalized Linear Mixed Models (McCulloch 1997)](www.jstor.org/stable/2291460)*

# Procedures
## Flow Chart

## Results

|Variables      | True Value | Initial Value| Converged Value 
|------------|------------|------------|------------|
| $\alpha$      | right-aligned |------------|------------|
| col 2 is      | centered      |------------|------------|
| zebra stripes | are neat      | ------------|------------|


## ![beta](https://github.com/Gaochenyin/MCEM/blob/master/3.2.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/300)
## Metropolis-Hasting Step
由基本的贝叶斯思想可得<a href="https://www.codecogs.com/eqnedit.php?latex=F(Z|\Omega,Y)&space;=&space;F(Z|\Omega)L(\Omega|Y)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?F(Z|\Omega,Y)&space;=&space;F(Z|\Omega)L(\Omega|Y)" title="F(Z|\Omega,Y) = F(Z|\Omega)L(\Omega|Y)" /></a>
## Gibbs Sampler
