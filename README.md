 Monte Carlo EM 
 ====

# Introduction

* This project described a Monte Carlo EM (**MCEM**) method to derive Maximum Likelihood Estimates (**MLE**) of the log-likelihood function. 

* In the E-step, perform *K = 500* Gibbs sampling incorporated with a Metropolis-Hastings step, and drop the first *100* as a burn-in procedure.

* Read article: *[Maximum Likelihood Algorithms for Generalized Linear Mixed Models (McCulloch 1997)](www.jstor.org/stable/2291460)*

# Model Notation

In this project, we consider a clustering problem. Suppose we have observed n observations, each observation is a binary process, i.e. the response $Y_{ij}=0 or 1$,$i=1,\cdots,n$,$j=1,\cdots,T$. Here n is the number of subjects and T is the length of observation. In general, T might vary across subjects, time points may also be different. In this project, however, we simply assume that all subjects have common time length and time points. We also assume that these subjects belong to two clusters. For each cluster, the conditional
expectation of response variable is
$$P_{ij}=E(Y_{ij}|U_i=1,X_{1,ij},Z_{1,i})=g^{-1}(\beta_1X_{1,ij}+Z_{1,i})$$
$$P_{ij}=E(Y_{ij}|U_i=2,X_{2,ij},Z_{2,i})=g^{-1}(\beta_2X_{2,ij}+Z_{2,i})$$
where $U$ is cluster membership, $X_{c,ij}$ and $Z_{c,i} (c = 1,2)$ are fixed and random effects, respectively. The link
function $g^{−1}(x)=\frac{exp(x)}{1+exp(x)}$ is given. In a typical clustering problem, $U$ is usually unknown, and hence we
treat $U$ as another random effect.

For random effects, we assume that $Z_{c,i}\sim N(0,\sigma_c^2)$ and $P(U = 1)=\pi$ (then $\pi_2=1-\pi_1$ ). Then the
parameter to be estimated is $\Omega=(\beta_1,\beta_2,\sigma_1,\sigma_2,\pi_1)$. Treating random effects as missing data, one can write
the complete data likelihood function as

$$L(\Omega|Y_{ij},U_i,Z_{U_i,i})=\prod_{i=1}^n\prod_{c=1}^2(\pi_cf_c(Z_{c,i})\[\prod_{j=1}^Tf_c(Y_{ij}|Z_{c,i})\])^{w_{i,c}}$$

where $f_c(Z_{c,i})$ is the density function of Normal distribution, $f_c(Y_ij|Z_{c,i})=P_{ij}^{Y_ij}(1−P_{ij})^{1−Y_{ij}}$. $ω_{ic}$ is the dummy variable of $U_i$

# Procedures

## Flow Chart

<div style="float:left;border:solid 1px 000;margin:2px;"><img src="https://github.com/Gaochenyin/MCEM/blob/master/flow.png"  width="400" ></div>

## Results
* EM Algorithm is sensitive to the initial values of parameters. We choose two fixed initialization.
* For each parameter, we calculate the changing rate of it and let $\epsilon^t$ be as following

$$\epsilon^{t}=max(\frac{\hat{\beta}_1^{t+1}-\hat{\beta}_1^{t}}{\hat{\beta}_1^{t}+\delta},\frac{\hat{\beta}_2^{t+1}-\hat{\beta}_2^{t}}{\hat{\beta}_2^{t}+\delta},\frac{\hat{\sigma}_1^{t+1}-\hat{\sigma}_1^{t}}{\hat{\sigma}_1^{t}+\delta},\frac{\hat{\sigma}_2^{t+1}-\hat{\sigma}_2^{t}}{\hat{\sigma}_2^{t}+\delta},\frac{\hat{\pi}_1^{t+1}-\hat{\pi}_1^{t}}{\hat{\pi}_1^{t}+\delta})$$

where $\delta>0$ is to assure that the denominator is positive. Setting the threshold $\epsilon_0$, if $\epsilon^t<\epsilon_0$ then we will consider the simulation converges.

* In this project, we choose $\delta=10^{-12},\epsilon_0=2.5*10^{-2}$
### Values

* Our convergences are pretty good, all parameters are **converged** in less than *50* steps, which cost about 1 minute.

* Besides, our project also contain different simulation with different initial value, which also obtain similar result. However, EM alogorithm is highly rely on **random numbers**, the final evaluation of these results is essential.

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

* MCEM could obtain a fair results based on the intialization mentioned before. 

* The MSE of $\beta_2$ and $\sigma_2$ by MCEM are much bigger than other parameters. This may be the result of the difference of the magnitudes.

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

