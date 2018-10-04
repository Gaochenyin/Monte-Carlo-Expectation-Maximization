 Monte Carlo EM 
 ====

# Introduction

* This project described a Monte Carlo EM (**MCEM**) method to derive Maximum Likelihood Estimates (**MLE**) of the log-likelihood function. 

* In the E-step, perform *K = 500* Gibbs sampling incorporated with a Metropolis-Hastings step, and drop the first *100* as a burn-in procedure.

* Read article: *[Maximum Likelihood Algorithms for Generalized Linear Mixed Models (McCulloch 1997)](www.jstor.org/stable/2291460)*

# Procedures

## Flow Chart


## Results

$\frac{\hat{\beta_1^{t+1}}}{2}$
### Values

1. My convergence is pretty good, all parameters are **converged** in less than *50* steps, which cost about 1 minute.

2. Besides, my project also contain different simulation with different initial value, which also obtain similar result. However, EM alogorithm is highly rely on **random numbers**, the final evaluation of these results is essential.

|Variables  | True Value | Initial Value| Converged Value 
|------------|------------|------------|------------|
| $\beta_1$      | 1 |0|0.9953680|
| $\beta_2$     | 1     |0|1.4076125|
| $\sigma_1$ | 2      | 1| 1.387342|
| $\sigma_2$ | 10      | 5|9.132040|
| $\pi_1$ | 0.6    | 0.8| 0.480500|

<div style="float:left;border:solid 1px 000;margin:2px;"><img src="https://github.com/Gaochenyin/MCEM/blob/master/beta.png"  width="400" ></div>

<div style="float:left;border:solid 1px 000;margin:2px;"><img src="https://github.com/Gaochenyin/MCEM/blob/master/sigma.png" width="400"></div>

<div style="float:left;border:solid 1px 000;margin:2px;"><img src="https://github.com/Gaochenyin/MCEM/blob/master/pi.png" width="400"></div>

### Evaluation

1. MCEM could obtain a fair results based on the intialization mentioned before. Because EM Algorithm is sensitive to the initial values of parameters, we choose two fixed initialization later

2. The MSE of $\beta_2$ and $\sigma_2$ by MCEM are much bigger than other parameters. This may be the result of the difference of the magnitudes.

|N| $\beta_1$ |$\beta_2$|$\sigma_1$|$\sigma_2$|$\pi_1$|
|--------|--------|--------|--------|--------|--------|
|100|0.01931882| 0.04455292|0.06768088| 0.5722849 |0.02300298
|200|0.01833736| 0.04397230| 0.06758311| 0.5673952 |0.02303765
|300|0.01820673| 0.04591263| 0.06918503| 0.6156710 |0.02475426
|400|0.01826401| 0.04716113| 0.06819809| 0.6093218 |0.02446871
|500|0.01859774| 0.04598016| 0.06920732| 0.6096665 |0.02488787
|600|0.01844890| 0.04572807| 0.06838667| 0.6032845 |0.02462379
|700|0.01873417| 0.04545865| 0.06950417| 0.6325660 |0.02525724
|800|0.01852523| 0.04516546| 0.06837866| 0.6270182 |0.02496853
|900|0.01822611| 0.04475583| 0.06697709| 0.6238291 |0.02465581
|1000|0.01824120| 0.04517658| 0.06683521| 0.6233873| 0.02449014

<div style="float:left;border:solid 1px 000;margin:2px;"><img src="https://github.com/Gaochenyin/MCEM/blob/master/MSE.png"  width="600" ></div>

