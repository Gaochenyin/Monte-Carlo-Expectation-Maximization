 Monte Carlo Expectation Maximization Algorithm
 ====

# Introduction

* This project described a Monte Carlo EM (**MCEM**) method to derive Maximum Likelihood Estimates (**MLE**) of the log-likelihood function. 

* In the E-step, perform *K = 500* Gibbs sampling incorporated with a Metropolis-Hastings step, and drop the first *100* as a burn-in procedure.

* Read article: *[Maximum Likelihood Algorithms for Generalized Linear Mixed Models (McCulloch 1997)](www.jstor.org/stable/2291460)*

# Model Notation

In this project, we consider a clustering problem. Suppose we have observed n observations, each observation is a binary process, i.e. the response <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;Y_{ij}=0~or~1,i=1,\cdots,n,j=1,\cdots,T" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;Y_{ij}=0~or~1,i=1,\cdots,n,j=1,\cdots,T" title="Y_{ij}=0~or~1,i=1,\cdots,n,j=1,\cdots,T" /></a>. Here n is the number of subjects and T is the length of observation. In general, T might vary across subjects, time points may also be different. In this project, however, we simply assume that all subjects have common time length and time points. We also assume that these subjects belong to two clusters. For each cluster, the conditional
expectation of response variable is


<a href="https://www.codecogs.com/eqnedit.php?latex=P_{ij}=E(Y_{ij}|U_i=1,X_{1,ij},Z_{1,i})=g^{-1}(\beta_1X_{1,ij}&plus;Z_{1,i})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?P_{ij}=E(Y_{ij}|U_i=1,X_{1,ij},Z_{1,i})=g^{-1}(\beta_1X_{1,ij}&plus;Z_{1,i})" title="P_{ij}=E(Y_{ij}|U_i=1,X_{1,ij},Z_{1,i})=g^{-1}(\beta_1X_{1,ij}+Z_{1,i})" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=P_{ij}=E(Y_{ij}|U_i=2,X_{2,ij},Z_{2,i})=g^{-1}(\beta_2X_{2,ij}&plus;Z_{2,i})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?P_{ij}=E(Y_{ij}|U_i=2,X_{2,ij},Z_{2,i})=g^{-1}(\beta_2X_{2,ij}&plus;Z_{2,i})" title="P_{ij}=E(Y_{ij}|U_i=2,X_{2,ij},Z_{2,i})=g^{-1}(\beta_2X_{2,ij}+Z_{2,i})" /></a>


where <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;U" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;U" title="U" /></a> is cluster membership, <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;X_{c,ij}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;X_{c,ij}" title="X_{c,ij}" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;Z_{c,i}~(c&space;=&space;1,2)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;Z_{c,i}~(c&space;=&space;1,2)" title="Z_{c,i}~(c = 1,2)" /></a> are fixed and random effects, respectively. The link
function <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;g^{-1}(x)=\frac{exp(x)}{1&plus;exp(x)}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;g^{-1}(x)=\frac{exp(x)}{1&plus;exp(x)}" title="g^{-1}(x)=\frac{exp(x)}{1+exp(x)}" /></a> is given. In a typical clustering problem, <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;U" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;U" title="U" /></a> is usually unknown, and hence we
treat <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;U" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;U" title="U" /></a> as another random effect.

For random effects, we assume that <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;Z_{c,i}\sim&space;N(0,\sigma_c^2)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;Z_{c,i}\sim&space;N(0,\sigma_c^2)" title="Z_{c,i}\sim N(0,\sigma_c^2)" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;P(U&space;=&space;1)=\pi" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;P(U&space;=&space;1)=\pi" title="P(U = 1)=\pi" /></a> (then <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\pi_2=1-\pi_1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\pi_2=1-\pi_1" title="\pi_2=1-\pi_1" /></a> ). Then the
parameter to be estimated is <a href="https://www.codecogs.com/eqnedit.php?latex=\inline&space;\Omega=(\beta_1,\beta_2,\sigma_1,\sigma_2,\pi_1)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\inline&space;\Omega=(\beta_1,\beta_2,\sigma_1,\sigma_2,\pi_1)" title="\Omega=(\beta_1,\beta_2,\sigma_1,\sigma_2,\pi_1)" /></a>. Treating random effects as missing data, one can write
the complete data likelihood function as

<a href="https://www.codecogs.com/eqnedit.php?latex=L(\Omega|Y_{ij},U_i,Z_{U_i,i})=\prod_{i=1}^n\prod_{c=1}^2(\pi_cf_c(Z_{c,i})\[\prod_{j=1}^Tf_c(Y_{ij}|Z_{c,i})\])^{w_{i,c}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?L(\Omega|Y_{ij},U_i,Z_{U_i,i})=\prod_{i=1}^n\prod_{c=1}^2(\pi_cf_c(Z_{c,i})\[\prod_{j=1}^Tf_c(Y_{ij}|Z_{c,i})\])^{w_{i,c}}" title="L(\Omega|Y_{ij},U_i,Z_{U_i,i})=\prod_{i=1}^n\prod_{c=1}^2(\pi_cf_c(Z_{c,i})\[\prod_{j=1}^Tf_c(Y_{ij}|Z_{c,i})\])^{w_{i,c}}" /></a>

where <a href="https://www.codecogs.com/eqnedit.php?latex=f_c(Z_{c,i})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f_c(Z_{c,i})" title="f_c(Z_{c,i})" /></a> is the density function of Normal distribution, <a href="https://www.codecogs.com/eqnedit.php?latex=f_c(Y_ij|Z_{c,i})=P_{ij}^{Y_ij}(1−P_{ij})^{1−Y_{ij}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f_c(Y_ij|Z_{c,i})=P_{ij}^{Y_ij}(1−P_{ij})^{1−Y_{ij}}" title="f_c(Y_ij|Z_{c,i})=P_{ij}^{Y_ij}(1−P_{ij})^{1−Y_{ij}}" /></a>. <a href="https://www.codecogs.com/eqnedit.php?latex=ω_{ic}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?ω_{ic}" title="ω_{ic}" /></a> is the dummy variable of <a href="https://www.codecogs.com/eqnedit.php?latex=U_i" target="_blank"><img src="https://latex.codecogs.com/gif.latex?U_i" title="U_i" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=w_{ic}=\begin{cases}&1~~,if~subject~i~belongs~to~cluster~c\\&space;&0~~,otherwise&space;\end{cases}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?w_{ic}=\begin{cases}&1~~,if~subject~i~belongs~to~cluster~c\\&space;&0~~,otherwise&space;\end{cases}" title="w_{ic}=\begin{cases}&1~~,if~subject~i~belongs~to~cluster~c\\ &0~~,otherwise \end{cases}" /></a>

the random effects <a href="https://www.codecogs.com/eqnedit.php?latex=U" target="_blank"><img src="https://latex.codecogs.com/gif.latex?U" title="U" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=Z" target="_blank"><img src="https://latex.codecogs.com/gif.latex?Z" title="Z" /></a> are called the [latent varaibles](https://en.wikipedia.org/wiki/Expectation%E2%80%93maximization_algorithm#Description) and <a href="https://www.codecogs.com/eqnedit.php?latex=(Y,U,Z)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?(Y,U,Z)" title="(Y,U,Z)" /></a> is called
complete data. The distribution of <a href="https://www.codecogs.com/eqnedit.php?latex=U" target="_blank"><img src="https://latex.codecogs.com/gif.latex?U" title="U" /></a> depends on <a href="https://www.codecogs.com/eqnedit.php?latex=\pi_1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\pi_1" title="\pi_1" /></a> and the distribution of <a href="https://www.codecogs.com/eqnedit.php?latex=Z" target="_blank"><img src="https://latex.codecogs.com/gif.latex?Z" title="Z" /></a> depends on <a href="https://www.codecogs.com/eqnedit.php?latex=U" target="_blank"><img src="https://latex.codecogs.com/gif.latex?U" title="U" /></a>, <a href="https://www.codecogs.com/eqnedit.php?latex=\sigma_1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sigma_1" title="\sigma_1" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=\sigma_2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sigma_2" title="\sigma_2" /></a>.

# Procedures

## Flow Chart

<div style="float:left;border:solid 1px 000;margin:2px;"><img src="https://github.com/Gaochenyin/MCEM/blob/master/flow_chart.png"  width="800" ></div>

## Estimation (Netwon-Raphson algorithm)

Partial differentiate the complete log-likelihood of the parameters and set derivatives to 0, we get the maximum likelihood estimators 

<a href="https://www.codecogs.com/eqnedit.php?latex=\hat{\pi}_{MLE,1}=\frac{1}{n}\sum_{i=1}^n\mathbb{I}_{\{U_i=1\}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\hat{\pi}_{MLE,1}=\frac{1}{n}\sum_{i=1}^n\mathbb{I}_{\{U_i=1\}}" title="\hat{\pi}_{MLE,1}=\frac{1}{n}\sum_{i=1}^n\mathbb{I}_{\{U_i=1\}}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\hat{\sigma}_{MLE,c}=\sqrt{\frac{\sum_{i=1}^n\mathbb{I}_{\{U_i=c\}}Z_{c,i}^2}{\sum_{i=1}^n\mathbb{I}_{\{U_i=c\}}}}" target="_blank"><img src="https://latex.codecogs.com/png.latex?\hat{\sigma}_{MLE,c}=\sqrt{\frac{\sum_{i=1}^n\mathbb{I}_{\{U_i=c\}}Z_{c,i}^2}{\sum_{i=1}^n\mathbb{I}_{\{U_i=c\}}}}" title="\hat{\sigma}_{MLE,c}=\sqrt{\frac{\sum_{i=1}^n\mathbb{I}_{\{U_i=c\}}Z_{c,i}^2}{\sum_{i=1}^n\mathbb{I}_{\{U_i=c\}}}}" /></a>

However the MLE of <a href="https://www.codecogs.com/eqnedit.php?latex=\beta_c" target="_blank"><img src="https://latex.codecogs.com/png.latex?\beta_c" title="\beta_c" /></a> is hard to obtain as a close form. We can use Netwon-Raphson algorithm to iteratively approximate it.

<a href="https://www.codecogs.com/eqnedit.php?latex=\beta_c^{(t&plus;1)}&space;=\beta^{(t)}_c&plus;\frac{\sum_{i=1}^n\mathbb{I}_{\{U_i=c\}}\sum_{j=1}^T(Y_{ij}X_{ij}-\frac{X_{c,ij}e^{(\beta_c^tX_{c,ij}&plus;Z_{c,i})}}{1&plus;e^{\beta_c^tX_{c,ij}&plus;Z_{c,i}}})}{\sum_{i=1}^n\mathbb{I}_{\{U_i=c\}}\sum_{j=1}^T\frac{X_{c,ij}^2e^{\beta_c^tX_{c,ij}&plus;Z_{c,i}}}{(1&plus;e^{\beta_c^tX_{c,ij}&plus;Z_{c,i}})^2}}" target="_blank"><img src="https://latex.codecogs.com/png.latex?\beta_c^{(t&plus;1)}&space;=\beta^{(t)}_c&plus;\frac{\sum_{i=1}^n\mathbb{I}_{\{U_i=c\}}\sum_{j=1}^T(Y_{ij}X_{ij}-\frac{X_{c,ij}e^{(\beta_c^tX_{c,ij}&plus;Z_{c,i})}}{1&plus;e^{\beta_c^tX_{c,ij}&plus;Z_{c,i}}})}{\sum_{i=1}^n\mathbb{I}_{\{U_i=c\}}\sum_{j=1}^T\frac{X_{c,ij}^2e^{\beta_c^tX_{c,ij}&plus;Z_{c,i}}}{(1&plus;e^{\beta_c^tX_{c,ij}&plus;Z_{c,i}})^2}}" title="\beta_c^{(t+1)} =\beta^{(t)}_c+\frac{\sum_{i=1}^n\mathbb{I}_{\{U_i=c\}}\sum_{j=1}^T(Y_{ij}X_{ij}-\frac{X_{c,ij}e^{(\beta_c^tX_{c,ij}+Z_{c,i})}}{1+e^{\beta_c^tX_{c,ij}+Z_{c,i}}})}{\sum_{i=1}^n\mathbb{I}_{\{U_i=c\}}\sum_{j=1}^T\frac{X_{c,ij}^2e^{\beta_c^tX_{c,ij}+Z_{c,i}}}{(1+e^{\beta_c^tX_{c,ij}+Z_{c,i}})^2}}" /></a>

## Sampling (Markov chain Monte Carlo)

Since it is hard to calculate the prior distribution of <a href="https://www.codecogs.com/eqnedit.php?latex=Y_i" target="_blank"><img src="https://latex.codecogs.com/png.latex?Y_i" title="Y_i" /></a> , it is difficult to sample directly from the multivariate distribution <a href="https://www.codecogs.com/eqnedit.php?latex=f(U_i,Z_i|Y_{ij},\Omega)" target="_blank"><img src="https://latex.codecogs.com/png.latex?f(U_i,Z_i|Y_{ij},\Omega)" title="f(U_i,Z_i|Y_{ij},\Omega)" /></a>. We can use Gibbs Sampling, a Markov chain Monte Carlo (MCMC) algorithm to obtain a sequence of observations which are approximated from the multivariate distribution.

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{aligned}f(U_i|Z_{U_i},\mathbf{Y_i})&space;&=\frac{f(U_i,Z_i,\mathbf{Y_i}|\Omega)}{f(Z_{i},\mathbf{Y_i}|\Omega)}\\&space;&=\frac{\pi_{U_i}f_{U_i}(Z_{i}|\sigma_1,\sigma_2)\prod\limits_{j=1}^Tf_{U_i}(Y_{ij}|Z_{i},\Omega)}{\sum\limits_{c=1}^2\pi_cf_c(Z_{i}|\sigma_1,\sigma_2)\prod\limits_{j=1}^Tf_{c}(Y_{ij}|Z_{i},\Omega)}&space;\end{aligned}" target="_blank"><img src="https://latex.codecogs.com/png.latex?\begin{aligned}f(U_i|Z_{U_i},\mathbf{Y_i})&space;&=\frac{f(U_i,Z_i,\mathbf{Y_i}|\Omega)}{f(Z_{i},\mathbf{Y_i}|\Omega)}\\&space;&=\frac{\pi_{U_i}f_{U_i}(Z_{i}|\sigma_1,\sigma_2)\prod\limits_{j=1}^Tf_{U_i}(Y_{ij}|Z_{i},\Omega)}{\sum\limits_{c=1}^2\pi_cf_c(Z_{i}|\sigma_1,\sigma_2)\prod\limits_{j=1}^Tf_{c}(Y_{ij}|Z_{i},\Omega)}&space;\end{aligned}" title="\begin{aligned}f(U_i|Z_{U_i},\mathbf{Y_i}) &=\frac{f(U_i,Z_i,\mathbf{Y_i}|\Omega)}{f(Z_{i},\mathbf{Y_i}|\Omega)}\\ &=\frac{\pi_{U_i}f_{U_i}(Z_{i}|\sigma_1,\sigma_2)\prod\limits_{j=1}^Tf_{U_i}(Y_{ij}|Z_{i},\Omega)}{\sum\limits_{c=1}^2\pi_cf_c(Z_{i}|\sigma_1,\sigma_2)\prod\limits_{j=1}^Tf_{c}(Y_{ij}|Z_{i},\Omega)} \end{aligned}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{aligned}f(Z_i|U_{i},\mathbf{Y_i})&space;&=\frac{f(U_i,Z_i,\mathbf{Y_i}|\Omega)}{f(Z_{i},\mathbf{Y_i}|\Omega)}\\&space;&=\frac{\pi_{U_i}f_{U_i}(Z_{i}|\sigma_1,\sigma_2)\prod\limits_{j=1}^Tf_{U_i}(Y_{ij}|Z_{i},\Omega)}{\int_{\mathbb{R}}\pi_{U_i}f_{U_i}(z|\sigma_1,\sigma_2)\prod\limits_{j=1}^Tf_{U_i}(Y_{ij}|z,\Omega)dz}&space;\end{aligned}" target="_blank"><img src="https://latex.codecogs.com/png.latex?\begin{aligned}f(Z_i|U_{i},\mathbf{Y_i})&space;&=\frac{f(U_i,Z_i,\mathbf{Y_i}|\Omega)}{f(Z_{i},\mathbf{Y_i}|\Omega)}\\&space;&=\frac{\pi_{U_i}f_{U_i}(Z_{i}|\sigma_1,\sigma_2)\prod\limits_{j=1}^Tf_{U_i}(Y_{ij}|Z_{i},\Omega)}{\int_{\mathbb{R}}\pi_{U_i}f_{U_i}(z|\sigma_1,\sigma_2)\prod\limits_{j=1}^Tf_{U_i}(Y_{ij}|z,\Omega)dz}&space;\end{aligned}" title="\begin{aligned}f(Z_i|U_{i},\mathbf{Y_i}) &=\frac{f(U_i,Z_i,\mathbf{Y_i}|\Omega)}{f(Z_{i},\mathbf{Y_i}|\Omega)}\\ &=\frac{\pi_{U_i}f_{U_i}(Z_{i}|\sigma_1,\sigma_2)\prod\limits_{j=1}^Tf_{U_i}(Y_{ij}|Z_{i},\Omega)}{\int_{\mathbb{R}}\pi_{U_i}f_{U_i}(z|\sigma_1,\sigma_2)\prod\limits_{j=1}^Tf_{U_i}(Y_{ij}|z,\Omega)dz} \end{aligned}" /></a>

Now we just assume that both of  can be attained. Then, we suppose that <a href="https://www.codecogs.com/eqnedit.php?latex=(U_{(k),i},Z_{(k),U_{(k)},i})" target="_blank"><img src="https://latex.codecogs.com/png.latex?(U_{(k),i},Z_{(k),U_{(k)},i})" title="(U_{(k),i},Z_{(k),U_{(k)},i})" /></a> is the <a href="https://www.codecogs.com/eqnedit.php?latex=i" target="_blank"><img src="https://latex.codecogs.com/png.latex?i" title="i" /></a>th component of the <a href="https://www.codecogs.com/eqnedit.php?latex=k" target="_blank"><img src="https://latex.codecogs.com/png.latex?k" title="k" /></a>th sample, we want to draw the <a href="https://www.codecogs.com/eqnedit.php?latex=i" target="_blank"><img src="https://latex.codecogs.com/png.latex?i" title="i" /></a>th component of the <a href="https://www.codecogs.com/eqnedit.php?latex=(k&plus;1)" target="_blank"><img src="https://latex.codecogs.com/png.latex?(k&plus;1)" title="(k+1)" /></a>th sample from these two conditional distributions. We draw

<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{aligned}&space;&U_{(k&plus;1),i}\sim&space;f_{U_i|Z_{U_i},i,\mathbf{Y}_i}(u|Z_{U_{(k)},i},\mathbf{Y}_i,\Omega)\\&space;&Z_{(k&plus;1),U_{(k&plus;1)},i}\sim&space;f_{Z_{U_i}|U_i,\mathbf{Y}_i}(z|U_{(k&plus;1)},\mathbf{Y}_i,\Omega)&space;\end{aligned}" target="_blank"><img src="https://latex.codecogs.com/png.latex?\begin{aligned}&space;&U_{(k&plus;1),i}\sim&space;f_{U_i|Z_{U_i},i,\mathbf{Y}_i}(u|Z_{U_{(k)},i},\mathbf{Y}_i,\Omega)\\&space;&Z_{(k&plus;1),U_{(k&plus;1)},i}\sim&space;f_{Z_{U_i}|U_i,\mathbf{Y}_i}(z|U_{(k&plus;1)},\mathbf{Y}_i,\Omega)&space;\end{aligned}" title="\begin{aligned} &U_{(k+1),i}\sim f_{U_i|Z_{U_i},i,\mathbf{Y}_i}(u|Z_{U_{(k)},i},\mathbf{Y}_i,\Omega)\\ &Z_{(k+1),U_{(k+1)},i}\sim f_{Z_{U_i}|U_i,\mathbf{Y}_i}(z|U_{(k+1)},\mathbf{Y}_i,\Omega) \end{aligned}" /></a>

Then we use the use Monte Carlo Integrating to approximate the expectation of the log-likelihood. 
## Results


* EM Algorithm is sensitive to the initial values of parameters. We choose two fixed initialization.
* For each parameter, we calculate the changing rate of it and let <a href="https://www.codecogs.com/eqnedit.php?latex=\epsilon^t" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\epsilon^t" title="\epsilon^t" /></a> be as following

<a href="https://www.codecogs.com/eqnedit.php?latex=\epsilon^{t}=\max\{\frac{\hat{\beta}_1^{t&plus;1}-\hat{\beta}_1^{t}}{\hat{\beta}_1^{t}&plus;\delta},\frac{\hat{\beta}_2^{t&plus;1}-\hat{\beta}_2^{t}}{\hat{\beta}_2^{t}&plus;\delta},\frac{\hat{\sigma}_1^{t&plus;1}-\hat{\sigma}_1^{t}}{\hat{\sigma}_1^{t}&plus;\delta},\frac{\hat{\sigma}_2^{t&plus;1}-\hat{\sigma}_2^{t}}{\hat{\sigma}_2^{t}&plus;\delta},\frac{\hat{\pi}_1^{t&plus;1}-\hat{\pi}_1^{t}}{\hat{\pi}_1^{t}&plus;\delta}\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\epsilon^{t}=\max\{\frac{\hat{\beta}_1^{t&plus;1}-\hat{\beta}_1^{t}}{\hat{\beta}_1^{t}&plus;\delta},\frac{\hat{\beta}_2^{t&plus;1}-\hat{\beta}_2^{t}}{\hat{\beta}_2^{t}&plus;\delta},\frac{\hat{\sigma}_1^{t&plus;1}-\hat{\sigma}_1^{t}}{\hat{\sigma}_1^{t}&plus;\delta},\frac{\hat{\sigma}_2^{t&plus;1}-\hat{\sigma}_2^{t}}{\hat{\sigma}_2^{t}&plus;\delta},\frac{\hat{\pi}_1^{t&plus;1}-\hat{\pi}_1^{t}}{\hat{\pi}_1^{t}&plus;\delta}\}" title="\epsilon^{t}=\max\{\frac{\hat{\beta}_1^{t+1}-\hat{\beta}_1^{t}}{\hat{\beta}_1^{t}+\delta},\frac{\hat{\beta}_2^{t+1}-\hat{\beta}_2^{t}}{\hat{\beta}_2^{t}+\delta},\frac{\hat{\sigma}_1^{t+1}-\hat{\sigma}_1^{t}}{\hat{\sigma}_1^{t}+\delta},\frac{\hat{\sigma}_2^{t+1}-\hat{\sigma}_2^{t}}{\hat{\sigma}_2^{t}+\delta},\frac{\hat{\pi}_1^{t+1}-\hat{\pi}_1^{t}}{\hat{\pi}_1^{t}+\delta}\}" /></a>

where <a href="https://www.codecogs.com/eqnedit.php?latex=\delta>0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\delta>0" title="\delta>0" /></a> is to assure that the denominator is positive. Setting the threshold <a href="https://www.codecogs.com/eqnedit.php?latex=\epsilon_0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\epsilon_0" title="\epsilon_0" /></a>, if <a href="https://www.codecogs.com/eqnedit.php?latex=\epsilon^t<\epsilon_0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\epsilon^t<\epsilon_0" title="\epsilon^t<\epsilon_0" /></a> then we will consider the simulation converges.

* In this project, we choose <a href="https://www.codecogs.com/eqnedit.php?latex=\delta=10^{-12},\epsilon_0=2.5*10^{-2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\delta=10^{-12},\epsilon_0=2.5*10^{-2}" title="\delta=10^{-12},\epsilon_0=2.5*10^{-2}" /></a>
### Values

* Our convergences are pretty good, all parameters are **converged** in less than *50* steps, which cost about 1 minute.

* Besides, our project also contain different simulation with different initial value, which also obtain similar result. However, EM alogorithm is highly rely on **random numbers**, the final evaluation of these results is essential.

|Variables  | True Value | Initial Value| Converged Value 
|------------|------------|------------|------------|
| <a href="https://www.codecogs.com/eqnedit.php?latex=\beta_1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\beta_1" title="\beta_1" /></a>      | 1 |0|0.9953680|
| <a href="https://www.codecogs.com/eqnedit.php?latex=\beta_2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\beta_2" title="\beta_2" /></a>     | 1     |0|1.4076125|
| <a href="https://www.codecogs.com/eqnedit.php?latex=\sigma_1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sigma_1" title="\sigma_1" /></a> | 2      | 1| 1.387342|
| <a href="https://www.codecogs.com/eqnedit.php?latex=\sigma_2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sigma_2" title="\sigma_2" /></a> | 10      | 5|9.132040|
| <a href="https://www.codecogs.com/eqnedit.php?latex=\pi_1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\pi_1" title="\pi_1" /></a> | 0.6    | 0.8| 0.480500|

<div style="float:left;border:solid 1px 000;margin:2px;"><img src="https://github.com/Gaochenyin/MCEM/blob/master/beta.png"  width="600" ></div>
<div style="float:left;border:solid 1px 000;margin:2px;"><img src="https://github.com/Gaochenyin/MCEM/blob/master/sigma.png" width="600"></div>
<div style="float:left;border:solid 1px 000;margin:2px;"><img src="https://github.com/Gaochenyin/MCEM/blob/master/pi.png" width="600"></div>

### Evaluation

* We conducted different number of simulations:<a href="https://www.codecogs.com/eqnedit.php?latex=100,200,\cdots,1000" target="_blank"><img src="https://latex.codecogs.com/gif.latex?100,200,\cdots,1000" title="100,200,\cdots,1000" /></a> and evaluate the corresponding *MSE*. From the result, we concluded that MCEM could obtain a fair results based on the intialization mentioned before. 

* The MSE of <a href="https://www.codecogs.com/eqnedit.php?latex=\beta_2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\beta_2" title="\beta_2" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=\sigma_2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sigma_2" title="\sigma_2" /></a> by MCEM are much bigger than other parameters. This may be the result of the difference of the magnitudes.

|N| <a href="https://www.codecogs.com/eqnedit.php?latex=\beta_1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\beta_1" title="\beta_1" /></a>  |<a href="https://www.codecogs.com/eqnedit.php?latex=\beta_2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\beta_2" title="\beta_2" /></a> |<a href="https://www.codecogs.com/eqnedit.php?latex=\sigma_1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sigma_1" title="\sigma_1" /></a>|<a href="https://www.codecogs.com/eqnedit.php?latex=\sigma_2" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sigma_2" title="\sigma_2" /></a>|<a href="https://www.codecogs.com/eqnedit.php?latex=\pi_1" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\pi_1" title="\pi_1" /></a>|
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

