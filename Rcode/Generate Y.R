set.seed(123456)
n <- 100#simulations
t <- 10
beta1 <- 1
beta2 <- 1
p1 <- 0.6
sigma1 <-2
sigma2 <- 10
x <- rnorm(t)
y_prob <- matrix(numeric(n*t),ncol=n)
y <- matrix(numeric(n*t),ncol=n)
U <- numeric(n)

for (i in 1:n)
{
  U[i] <- rbinom(1,1,prob = p1)
  if(U[i]==1)
  {
    z1 <- rnorm(1,sd=sigma1)
    y_prob[,i] <- exp(beta1*x+z1)/(1+exp(beta1*x+z1))
  }
  else
  {
    z2 <- rnorm(1,sd=sigma2)
    y_prob[,i] <- exp(beta2*x+z2)/(1+exp(beta2*x+z2))
  }
  y[,i] <- rbinom(10,1,prob = y_prob[,i])
}