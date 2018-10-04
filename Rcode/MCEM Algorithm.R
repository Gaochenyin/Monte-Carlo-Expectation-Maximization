
MCEM <- function(seed,N,Time,beta1_0,beta2_0,sigma1_0,sigma2_0,p1_0,iterstep){
  set.seed(123456)
  n <- N
  t <- Time
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
  set.seed(12)#seed)
  allstep <- iterstep
  beta1_t <- array(NA,dim=allstep)
  beta2_t <- array(NA,dim=allstep)
  beta1_t_ac <- array(NA,dim=allstep)
  beta2_t_ac <- array(NA,dim=allstep)
  sigma1_t <- array(NA,dim=allstep)
  sigma2_t <- array(NA,dim=allstep)
  sigma1_t_ac <- array(NA,dim=allstep)
  sigma2_t_ac <- array(NA,dim=allstep)
  p1_t <- array(NA,dim=allstep)
  p1_t_ac <- array(NA,dim=allstep)
  beta1_t[1] <- beta1_0
  beta2_t[1] <- beta2_0
  beta1_t_ac[1] <- 0
  beta2_t_ac[1] <- 0
  sigma1_t[1] <-sigma1_0
  sigma2_t[1] <-sigma2_0
  sigma1_t_ac[1] <- 1
  sigma2_t_ac[1] <- 5
  p1_t[1] <-p1_0
  p1_t_ac[1] <- 0.8
  step_t <- 1
  tol1 <- array(1,dim = allstep)
  tol2 <- array(1,dim = allstep)
  while(step_t<allstep&tol1[step_t]>10^{-10})
  {
    UZlist <- list()
    #generate random variable
    ##对每个人进行抽样U和Z
    for (member in 1:n) {
      N <- 500
      if(step_t<40)N <- 300
      if(step_t<20)N <- 150
      burnin <- 100
      X <- matrix(0,N,2)
      X[1,] <- c(rbinom(1,1,prob = p1_t[step_t]),rnorm(1,0,sigma1_t[step_t]))
      ypresent <- y[,member]
      for (i in 2:N)
      {
        ##generate U
        
        Zlast <- X[i-1,2]
        prob1 <- prod(exp(ypresent*(beta1_t[step_t]*x+Zlast))/(1+exp(beta1_t[step_t]*x+Zlast)))*p1_t[step_t]/sigma1_t[step_t]*exp(-Zlast^2/(2*sigma1_t[step_t]^2))
        prob2 <- prod(exp(ypresent*(beta2_t[step_t]*x+Zlast))/(1+exp(beta2_t[step_t]*x+Zlast)))*(1-p1_t[step_t])/sigma2_t[step_t]*exp(-Zlast^2/(2*sigma2_t[step_t]^2))
        X[i,1] <- rbinom(1,1,prob = prob1/(prob1+prob2))
        
        
        ##generate Z
        
        Ulast <- X[i,1]
        if(Ulast==1)
        {
          Z_star <- rnorm(1,0,sigma1_t[step_t])
          prod_star <- prod(exp(ypresent*(beta1_t[step_t]*x+Z_star))/(1+exp(beta1_t[step_t]*x+Z_star)))
          prod_last <- prod(exp(ypresent*(beta1_t[step_t]*x+Zlast))/(1+exp(beta1_t[step_t]*x+Zlast)))
          
        }
        else
        {
          Z_star <- rnorm(1,0,sigma2_t[step_t])
          prod_star <- prod(exp(ypresent*(beta2_t[step_t]*x+Z_star))/(1+exp(beta2_t[step_t]*x+Z_star)))
          prod_last <- prod(exp(ypresent*(beta2_t[step_t]*x+Zlast))/(1+exp(beta2_t[step_t]*x+Zlast)))
          
        }
        alpha <- min(1,prod_star/prod_last) 
        U <- runif(1)
        if(U<alpha)
          X[i,2] <- Z_star
        else
          X[i,2] <- Zlast
        
        
      }
      X <- X[burnin:N,] 
      UZlist[[member]] <- X
    }
    ##进行M-step
    step_t <- step_t+1
    Zsquare_1 <- numeric(N-100)
    Zsquare_0 <- numeric(N-100)
    U_mean <- numeric(N-100)
    beta1k <- numeric(N-100)
    beta2k <- numeric(N-100)
    sigma1k <- numeric(N-100)
    sigma2k <- numeric(N-100)
    p1k <- numeric(N-100)
    for(number in 1:(N-100))
    {
      presentdata <- sapply(UZlist, function(x){x[number,]})
      presentdata_U <- presentdata[1,]
      comp1 <- numeric(n)
      comp2 <- numeric(n)
      for(m in 1:length(presentdata_U))
      {
        if(presentdata_U[m]==1)
        {
          comp1[m] <- sum(y[,m]*x-x*exp((beta1_t[step_t-1]*x+presentdata[2,m]))/(1+exp(beta1_t[step_t-1]*x+presentdata[2,m])))
          comp2[m] <- sum(exp(beta1_t[step_t-1]*x+presentdata[2,m])*x^2/(1+exp(beta1_t[step_t-1]*x+presentdata[2,m]))^2)
        }
        else
        {
          comp1[m] <- sum(y[,m]*x-x*exp((beta2_t[step_t-1]*x+presentdata[2,m]))/(1+exp(beta2_t[step_t-1]*x+presentdata[2,m])))
          comp2[m] <- sum(exp(beta2_t[step_t-1]*x+presentdata[2,m])*x^2/(1+exp(beta2_t[step_t-1]*x+presentdata[2,m]))^2)
        }
      }
      #newthon
      beta1k[number] <- sum(comp1[presentdata_U==1])/sum(comp2[presentdata_U==1])+beta1_t[step_t-1]
      beta2k[number] <- sum(comp1[presentdata_U==0])/sum(comp2[presentdata_U==0])+beta2_t[step_t-1]
      presentdata_1 <- matrix(presentdata[,presentdata[1,]==1],nrow=2)
      presentdata_0 <- matrix(presentdata[,presentdata[1,]==0],nrow=2)
      sigma1k[number] <- sum(sigma1_t_ac[step_t-1]^3-presentdata_1[2,]^2*sigma1_t_ac[step_t-1])/sum(-3*presentdata_1[2,]^2+sigma1_t_ac[step_t-1]^2)+sigma1_t_ac[step_t-1]
      sigma2k[number] <- sum(sigma2_t_ac[step_t-1]^3-presentdata_0[2,]^2*sigma2_t_ac[step_t-1])/sum(-3*presentdata_0[2,]^2+sigma2_t_ac[step_t-1]^2)+sigma2_t_ac[step_t-1]
      p1k[number] <- sum((1-presentdata[1,])/(1-p1_t_ac[step_t-1])-(presentdata[1,])/(p1_t_ac[step_t-1]))/sum(-presentdata[1,]/p1_t_ac[step_t-1]^2-(1-presentdata[1,])/(1-p1_t_ac[step_t-1])^2)+p1_t_ac[step_t-1]
      
      ##ordinary
      Zsquare_1[number] <- mean(presentdata_1[2,]^2)
      Zsquare_0[number] <- mean(presentdata_0[2,]^2)
      U_mean[number] <- mean(presentdata[1,])
    }
    beta1_t[step_t] <- mean(beta1k)
    beta2_t[step_t] <- mean(beta2k)
    
    ##locuis ac
    deltabeta1_t <- beta1_t[step_t]-beta1_t_ac[step_t-1]
    deltabeta2_t <- beta2_t[step_t]-beta2_t_ac[step_t-1]
    
    beta1k_ac <- numeric(N-100)
    beta2k_ac <- numeric(N-100)
    for(number in 1:(N-100))
    {
      presentdata <- sapply(UZlist, function(x){x[number,]})
      presentdata_U <- presentdata[1,]
      comp2_ac <- matrix(NA,nrow = 100,ncol = 10)
      for(m in 1:length(presentdata_U))
      {
        if(presentdata_U[m]==1)
        {
          comp2_ac[m,] <- exp(beta1_t[step_t]*x+presentdata[2,m])*x^2/(1+exp(beta1_t[step_t]*x+presentdata[2,m]))^2
        }
        else
        {
          comp2_ac[m,] <- exp(beta2_t[step_t]*x+presentdata[2,m])*x^2/(1+exp(beta2_t[step_t]*x+presentdata[2,m]))^2
        }
      }
      beta1k_ac[number] <- (sum(comp2_ac[presentdata_U==1,])-var(as.numeric(comp2_ac[presentdata_U==1,])))*deltabeta1_t/(sum(comp2[presentdata_U==1]))+beta1_t[step_t-1]
      beta2k_ac[number] <- (sum(comp2_ac[presentdata_U==0,])-var(as.numeric(comp2_ac[presentdata_U==0,])))*deltabeta2_t/(sum(comp2[presentdata_U==0]))+beta2_t[step_t-1]
    }
    sigma1_t[step_t] <- sqrt(mean(Zsquare_1))
    sigma2_t[step_t] <- sqrt(mean(Zsquare_0))
    p1_t[step_t] <- mean(U_mean)
    ##ac
    beta1_t_ac[step_t] <- mean(beta1k_ac)
    beta2_t_ac[step_t] <- mean(beta2k_ac)
    sigma1_t_ac[step_t] <- mean(sigma1k)
    sigma2_t_ac[step_t] <- mean(sigma2k)
    p1_t_ac[step_t] <- mean(p1k)
    ##tol
    now1 <- c(beta1_t[step_t],beta2_t[step_t],sigma1_t[step_t],sigma2_t[step_t],p1_t[step_t])
    last1 <- c(beta1_t[step_t-1],beta2_t[step_t-1],sigma1_t[step_t-1],sigma2_t[step_t-1],p1_t[step_t-1])
    tol1[step_t] <- mean((now1-last1)/now1)^2
    now2 <- c(beta1_t[step_t],beta2_t[step_t],sigma1_t_ac[step_t],sigma2_t_ac[step_t],p1_t_ac[step_t])
    last2 <- c(beta1_t[step_t-1],beta2_t[step_t-1],sigma1_t_ac[step_t-1],sigma2_t_ac[step_t-1],p1_t_ac[step_t-1])
    tol2[step_t] <- mean((now2-last2)/now2)^2
    print(step_t)
  }
  return(data.frame(beta1_t,beta2_t,sigma1_t,sigma2_t,p1_t,beta1_t_ac,beta2_t_ac,sigma1_t_ac,sigma2_t_ac,p1_t_ac))
}
