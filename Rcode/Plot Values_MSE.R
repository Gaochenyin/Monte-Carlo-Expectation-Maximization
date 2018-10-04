##算法1
plot(beta1_t,pch=3,xlab = '',ylab = '',type = 'b',xlim = c(0,step_t))
mtext(expression(beta[1]),side=2,line=1.5,cex=1)
par(new=T)
plot(beta2_t,axes=F,xlab = '',ylab = '',type = 'b',col='blue',pch=3,xlim = c(0,step_t))
axis(4,col='blue',las=1)
mtext(expression(beta[2]),side=4,line=-1.5,cex=1)
legend('bottomright',legend = c('beta1','beta2'),pch = 3,col = c('black','blue'))

plot(sigma1_t,pch=3,xlab = '',ylab = '',type = 'b',xlim = c(0,step_t))
mtext(expression(sigma[1]),side=2,line=1.5,cex=1)
par(new=T)
plot(sigma2_t,axes=F,xlab = '',ylab = '',type = 'b',col='blue',pch=3,xlim = c(0,step_t))
axis(4,col='blue',las=1)
mtext(expression(sigma[2]),side=4,line=-1.5,cex=1)
legend('bottomright',legend = c('sigma1','sigma2'),pch = 3,col = c('black','blue'))

plot(p1_t,pch=3,xlab = '',ylab = '',type = 'b',xlim = c(0,step_t))
mtext(expression(pi[1]),side=2,line=1.5,cex=1)
legend('topright',legend = c('pi1'),pch = 3,col = 'black')



##加速前后对比
par(cex.main=2)
plot(beta2_t,pch=3,xlab = '迭代次数',ylab = '',type = 'b',xlim = c(0,50),main = 'beta Locuis加速')
par(new=T)
plot(beta2_t_ac,axes=F,xlab = '',ylab = '',type = 'b',col='blue',pch=3,xlim = c(0,50))
legend('bottomright',legend = c('加速前','加速后'),pch = 3,col = c('black','blue'))

plot(diff(beta2_t)[1:145]/beta2_t[2:146],ylab = '变化率',type='b',pch=3,ylim=c(-0.5,1),xlab='迭代次数',main='beta 变化率')
lines(diff(beta2_t_ac)[1:145]/beta2_t_ac[2:146],type='b',pch=3,col='blue')
legend('topright',legend = c('加速前','加速后'),pch = 3,col = c('black','blue'))

##进行1000次模拟
coef100_1000 <- list()
iterstep <- 50
for (i in 1:1000) {
  coef100_1000[[i]] <- MCEM(i,iterstep)
  print(paste('完成',i,seq=''))
}


##模拟100次
par(cex.main=2)
plot(apply(matrix(sapply(coef100_1000,function(x){x[,i]}),nrow = 50)[,1:100],1,mean),xlab = '迭代次数',type = 'b',ylab = '',main = mainnames[i])
lines(apply(matrix(sapply(coef100_1000,function(x){x[,i+3]}),nrow = 50)[,1:100],1,mean),type = 'b',ylab = '',pch=2,col=2)
legend('topright',legend = c('算法1','算法2'),pch=1:2,col=1:2)
#100次模拟效果图
matplot((matrix(sapply(coef100_1000,function(x){x[,8]}),nrow = 50)[,1:100]),xlab='迭代次数',ylab = '',type='b',main = '100次模拟的beta2',pch=rep(1:10,10),col=1:100)

##参数MSE
plot(apply(matrix(sapply(coef100_1000,function(x){x[,1]}),nrow = 50)[,1:100],1,sd),xlab='迭代次数',ylab = '',type='b',col=1,pch=1,ylim = c(0,.1))
for(i in c(2,3,5))
  lines(apply(matrix(sapply(coef100_1000,function(x){x[,i]}),nrow = 50)[,1:100],1,sd),xlab='迭代次数',ylab = '',type='b',col=i*10,pch=i)
par(new=T)
plot(apply(matrix(sapply(coef100_1000,function(x){x[,4]}),nrow = 50)[,1:100],1,sd),axes=F,xlab='迭代次数',ylab = '',type='b',col=4*10,pch=4)
axis(4,col=4,las=1)
legend('topleft',legend = c('beta1','beta2','sigma1','sigma2','pi1'),pch =1:5,col=1:5 )


##1000次MSE
par(cex.main=2)
plot(apply(matrix(sapply(coef100_1000,function(x){x[,5]}),nrow = 50)[,1:100],1,sd),main='多次模拟pi1的MSE',xlab = '迭代次数',type = 'b',pch=100/100,col=100/100,ylab = '',ylim = c(0,0.05))
for(i in seq(200,1000,by=100))
{
  lines(apply(matrix(sapply(coef100_1000,function(x){x[,5]}),nrow = 50)[,1:i],1,sd),type = 'b',pch=i/100,col=i/100)
}
legend('bottomright',pch = 1:10,col = 1:10,legend =paste(seq(100,1000,by =100),'次',sep=''),cex=0.8 )
MSE200_1000 <- matrix(NA,nrow = 9,ncol = 5)
for (i in 1:5) {
  for (j in seq(200,1000,by=100)) {
    MSE200_1000[j/100-1,i] <- apply(matrix(sapply(coef100_1000,function(x){x[,i]}),nrow = 50)[,1:j],1,sd)[50]
  }
}
Time <- seq(200,1000,by=100)
MSE200_1000 <- data.frame(Time,MSE200_1000)
colnames(MSE200_1000) <- c('Time','beta1','beta2','sigma1','sigma2','pi1')
par(cex.main=2)
matplot(MSE200_1000[,2:6],pch=1:5,col=c('green','blue','purple','red','black'),main='MSE',ylab = '',xlab = '模拟次数(100)',type = 'b',cex = 2)
legend('right',legend = c('beta1','beta2','sigma1','sigma2','pi1'),pch=1:5,col = c('green','blue','purple','red','black'))



##构造经验分布直方图
coef_dis <- matrix(NA,nrow = 1000,ncol = 5)
for (i in 1:1000) {
  coef_dis[i,] <- unlist(coef100_1000[[i]][50,1:5])
  print(i)
}
coef_dis <- data.frame(coef_dis)

par(cex.main=2)
par(mfrow=c(2,3))
for (i in 1:5){
  hist(coef_dis[,i],col='steelblue',border='white',main = paste(mainnames[i],'直方图'),xlab='')
}

##poorman 100
#poorman <- list()
for(number in 1:100)
{
  set.seed(number)
  b1_0 <- coef_dis[sample(c(1:1000),1),1]
  b2_0 <- coef_dis[sample(c(1:1000),1),2]
  s1_0 <- coef_dis[sample(c(1:1000),1),3]
  s2_0 <- coef_dis[sample(c(1:1000),1),4]
  p1_0 <- coef_dis[sample(c(1:1000),1),5]
  poorman[[number]] <- MCEM(1015,100,10,b1_0,b2_0,s1_0,s2_0,p1_0,10)
  print(paste('已完成',number))
}


mean(coef_dis[,1])
mean(sapply(poorman, function(x){x[,1]})[10,])

par(cex.main=2)
par(mfrow=c(2,3))
for (i in 1:5) {
  matplot((sapply(poorman,function(x){x[,i]})),type='l',xlab = '迭代次数',ylab='',main = paste(mainnames[i]))
  #lines(h=quantile(sapply(poorman,function(x){x[,i]})[10,],c(0.025,0.975)))
}
quantile(sapply(poorman,function(x){x[,5]})[10,],c(0.5))
sample(c(1:1000),1)


b_1 <- quantile(sapply(poorman,function(x){x[,1]})[10,],c(0.5))
b_2 <- quantile(sapply(poorman,function(x){x[,2]})[10,],c(0.5))
s_1 <- quantile(sapply(poorman,function(x){x[,3]})[10,],c(0.5))
s_2 <- quantile(sapply(poorman,function(x){x[,4]})[10,],c(0.5))
p_1 <- quantile(sapply(poorman,function(x){x[,5]})[10,],c(0.5))
