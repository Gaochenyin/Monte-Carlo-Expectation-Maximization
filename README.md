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
| $\beta_1$      | 1 |0|0.9953680|
| $\beta_2$     | 1     |0|1.4076125|
| $\sigma_1$ | 2      | 1| 1.387342|
| $\sigma_2$ | 10      | 5|9.132040|
| $\pi_1$ | 0.6    | 0.8| 0.480500|


##
<div style="float:left;border:solid 1px 000;margin:2px;"><img src="https://github.com/Gaochenyin/MCEM/blob/master/3.2.png"  width="400" ></div>

<div style="float:left;border:solid 1px 000;margin:2px;"><img src="https://github.com/Gaochenyin/MCEM/blob/master/3.3.png" width="400"></div>

<div style="float:left;border:solid 1px 000;margin:2px;"><img src="https://github.com/Gaochenyin/MCEM/blob/master/3.4.png" width="400"></div>

## Metropolis-Hasting Step
由基本的贝叶斯思想可得<a href="https://www.codecogs.com/eqnedit.php?latex=F(Z|\Omega,Y)&space;=&space;F(Z|\Omega)L(\Omega|Y)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?F(Z|\Omega,Y)&space;=&space;F(Z|\Omega)L(\Omega|Y)" title="F(Z|\Omega,Y) = F(Z|\Omega)L(\Omega|Y)" /></a>
## Gibbs Sampler
