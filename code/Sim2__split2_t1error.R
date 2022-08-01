
set.seed(20)

reps = 1000

K     = 100 #number of stages in trial
S     = 2
T     = 3

nk_.ge   = 80
Nx.ge    = c(501, 1279, 324) + c(1034, 2604, 683)
Px.ge    = Nx.ge/sum(Nx.ge)
nk_.ge   = round(nk_.ge*Px.ge) #patients per stage
Py.ge   = matrix(c(0.140, 0.262, 0.414, 0.140, 0.262, 0.414), nrow=3)
#Py    = c(0.25,0.22) #probability of death on each treatment
#S     = length(Py[1,]) #Number of treatments
#T     = length(Px) #Number of care

alpha.ge       = 1 #uninformative beta
beta.ge        = 1 #uninformative beta
Y.ge           = array(dim=c(K,S,T))
postmean.ge    = array(dim=c(K,S,T))
mle.ge         = array(dim=c(K,S,T))
PostAlpha.ge   = array(dim=c(K,S,T))
PostBeta.ge    = array(dim=c(K,S,T))

nk.ge          = array(dim=c(K,S,T))
for (h in 1:length(nk_.ge)){
  nk.ge[,,h] = round(nk_.ge[h]/S)
}

rand_ratio.ge = rep(0.5, T)
rand_ratios.ge = matrix(nrow=K, ncol=T)

UPPER_LIMIT <- 0.9
LOWER_LIMIT <- 0.1

nka.ge = array(dim = c(K, S, T,reps))
Ya.ge = array(dim = c(K, S, T,reps))
postmean_a.ge = array(dim = c(K, S, T,reps))
mle_a.ge = array(dim=c(K, S,T, reps))

thetas.ge = array(dim = c(reps,K,T))
#y = vector(length=3)

conf_sup.ge <- matrix(nrow=K/10, ncol=reps)

# Stage 1
for (n in 1:reps){

rand_ratio.ge = rep(0.5, T)
  
for (h in 1:T){
for(j in 1:S){
y.ge = rbinom(1, round(nk.ge[1,j,h]), 0.25)
Y.ge[1,j,h]         = y.ge
PostAlpha.ge[1,j,h] = alpha.ge + Y.ge[1,j,h]
PostBeta.ge[1,j,h]  = beta.ge  + nk.ge[1,j,h] -Y.ge[1,j,h]
postmean.ge[1,j,h]  = PostAlpha.ge[1,j,h]/(PostAlpha.ge[1,j,h]+PostBeta.ge[1,j,h])
mle.ge[1,j,h]       = sum(Y.ge[1,j,h])/(nk.ge[1,j,h])
}
}

for(i in 2:K){

for (h in 1:T){
    
nk.ge[i, 1,h] = rbinom(1,nk_.ge[h], 1-rand_ratio.ge[h])
nk.ge[i, 2,h] = rbinom(1,nk_.ge[h], rand_ratio.ge[h])

for(j in 1:S){
  y.ge = rbinom(1, round(nk.ge[i,j,h]), 0.25)

Y.ge[i,j,h]         = y.ge
PostAlpha.ge[i,j,h] = alpha.ge + sum(Y.ge[1:i,j,h], na.rm=T)
PostBeta.ge[i,j,h]  = beta.ge + sum(nk.ge[1:i,j,h]-Y.ge[1:i,j,h], na.rm=T)
postmean.ge[i,j,h]  = PostAlpha.ge[i,j,h]/(PostAlpha.ge[i,j,h]+PostBeta.ge[i,j,h])
mle.ge[i,j,h]       = sum(Y.ge[1:i,j,h])/sum(nk.ge[1:i,j,h])


}

if (i > 28){
    
if (i %% 7 == 0){
X1 <- rbeta(100, PostAlpha.ge[i-28, 1,h], PostBeta.ge[i-28, 1,h])
X2 <- rbeta(100, PostAlpha.ge[i-28, 2,h], PostBeta.ge[i-28, 2,h])
theta <- length(X1[X1>=X2])/100
thetas.ge[n,i,h] <- theta
k = i/(K)
rand_ratio.ge[h] <- theta^k / (theta^k + (1-theta)^k)
rand_ratios.ge[i,h] = rand_ratio.ge[h]

#making sure samples don't go to 0
if (rand_ratio.ge[h]>UPPER_LIMIT){
  rand_ratio.ge[h] = UPPER_LIMIT
}
if (rand_ratio.ge[h]<LOWER_LIMIT){
  rand_ratio.ge[h] = LOWER_LIMIT
}

}
}

}
if (i %% 10 == 0){
  pat_df <- data.frame(t=rep(1, sum(nk.ge[1:i,,])), o=rep(0, sum(nk.ge[1:i,,])))
  pat_df[1:sum(nk.ge[1:i,1,]), 't'] = 0
  pat_df[1:(sum(Y.ge[1:i,1,])), 'o'] = 1
  pat_df[(sum(nk.ge[1:i,1,])+1):(sum(nk.ge[1:i,1,])+sum(Y.ge[1:i,2,])), 'o'] = 1
  
  fit = glm(o~t,data=pat_df, family=binomial)
  conf_sup.ge[i/10,n] = summary(fit)$coef['t','Pr(>|z|)']
}


}

nka.ge[,,,n] = nk.ge
Ya.ge[,,,n] = Y.ge
postmean_a.ge[,,,n] = postmean.ge
mle_a.ge[,,,n] = mle.ge
}


t1e_metric.ge = vector(length=10)
for (i in 1:10){
  t1e_metric.ge[i] = mean(conf_sup.ge[i,]<0.05)

}

t1e_metric.ge

plot(t1e_metric.ge*100, lwd=4, type='b', col='blue', ylim=c(0,15),
     axes=FALSE, xlab='Trial progress(%)', ylab='Family-wise error rate')
axis(1, seq(1,10), labels = seq(12.5, 125,12.5))
axis(2, seq(0, 15, 2.5), labels = seq(0, 15, 2.5))




#mean(apply(conf_sup.ge[8,,]<0.05, 2, any))  #0.201
#mean(apply(conf_sup.ge[8,,]<(0.05/3), 2, any))  #0.082 #bonferroni correction










nk_avg.ge = array(dim=c(K,S,T))
Y_avg.ge  = array(dim=c(K,S,T))

for (h in 1:T){
for (j in 1:S){
  nk_avg.ge[,j,h] = apply(nka.ge[,j,h,], 1, mean)
  Y_avg.ge[,j,h]  = apply(Ya.ge[,j,h,], 1, mean, na.rm=T)
}
}


for (h in 1:T){
print(sum(nk_avg.ge[1:80,1,h]))
print(sum(nk_avg.ge[1:80,2,h]))
}

sum(nk_avg.ge)

for (h in 1:T){
  print(sum(Y_avg.ge[1:80,1,h]))
  print(sum(Y_avg.ge[1:80,2,h]))
}


#mortality at n=6450
sum(nk_avg.ge[1:80,,])
sum(Y_avg.ge[1:80,,])



for (h in 1:T){
  conf_sup_avg = apply(conf_sup[,h,],1, mean)
  plot(conf_sup_avg)
  lines(conf_sup_avg)
}


h=1

c_vector = c('red', 'blue', 'brown')

power_avg = apply(conf_sup.h[,h,]<0.05, 1, mean)
power_avg = c(0, power_avg)
plot(power_avg, xlab='Trial Progress (%)', ylab = 'Power', ylim= c(0, 1), xaxt = 'n')

for (h in 1:T){
  power_avg = apply(conf_sup.ge[,h,]<0.05, 1, mean)
  power_avg = c(0, power_avg)
  points(power_avg, lwd=2)
  lines(power_avg,lty='dashed', col=c_vector[h], lwd=3)
  
  power_avg = apply(conf_sup.h[,h,]<0.05, 1, mean)
  power_avg = c(0, power_avg)
  points(power_avg, lwd=2)
  lines(power_avg, col=c_vector[h], lwd=3)
  
  
  
}

#legend('bottomright', legend = c('REMAP-CAP', 'Tuning'), fill=c('red', 'blue'))

legend("bottomright",legend = c('subgroup iii', 'subgroup ii', 'subgroup i'),col = c_vector[c(3,2,1)], lwd=2)
legend('bottom', legend=c(expression('T'[s]), expression('RMC'[s])), lty = c('solid','dashed'))

axis(1, seq(1,11,2), labels = seq(0, 125,25))

abline(h=0.8, lty=4)
abline(h=0.9, lty=4)






c_vector = c('red', 'blue', 'brown')
plot(nk_avg.ge[c(2, 1:10*10),2,1]/nk_.ge[1]*100, ylim = c(0,100), col='red', xlab= 'Trial progress (%)', 
     ylab='Proportion of patients allocated dexamethasone(%)', type = 'b', axes=FALSE)
axis(1, seq(1,11,2), labels = seq(0, 125,25))
axis(2, seq(0, 100, 25), labels = seq(0, 100, 25))

for (h in 1:T){
points(nk_avg.ge[c(2, 1:10*10),2,h]/nk_.ge[h]*100, col=c_vector[h], type = 'b', lwd=3)
points(nk_avg.h[c(2, 1:10*10),2,h]/nk_.h[h]*100, col=c_vector[h], type = 'b', lty = 'dashed', lwd=3)
}


legend("right",legend = c('subgroup iii', 'subgroup ii','subgroup i'),col = c_vector[c(3,2,1)], lwd=4)
legend('bottomleft', legend = c(expression('T'[s]), expression('RMC'[s])), lwd =3, lty = c('solid','dashed'))
abline(v=4.4, lty='dashed')


sum(Ya)/reps
#Y_avg = matrix(nrow=K, ncol=S)
#for (j in 1:S){
#  Y_avg[,j] = apply(Ya[,j,], 1, mean)
#}

treatment_superior = matrix(nrow=K, ncol=T)
for (h in 1:T){
for (i in 1:K){
  treatment_superior[i,h] = mean(postmean_a[i,1,h,] >= postmean_a[i,2,h,])
}
}

for (h in 1:T){
plot(treatment_superior[,h])
}

for (h in 1:T){
thetas_avg = apply(thetas.ge[,,h], 2, mean)
plot(thetas_avg)
}

for (h in 2:T){
above_theta = thetas.ge > 0.95
above_theta_avg = apply(above_theta, 2, mean)
plot(above_theta_avg)
}
#plot(postmean[,1],ylim=range(postmean),type="l",lwd=4,
#     ylab="Posterior mean",xlab="Interim Analysis")
#lines(postmean[,2],type="l",lwd=4,col="red")
#legend("bottomright",c("Treatment 1","Treatment 2"),col=c("black","red"),lwd=4,bty="n")

#lines(mle[,1],ylim=range(mle),type="l",lwd=4,lty=2)
#lines(mle[,2],type="l",lwd=4,col="red",lty=2)

mle_avg = matrix(nrow=T, ncol=S)
mle_sd = matrix(nrow=T, ncol=S)

for (h in 1:T){
  for (j in 1:S){
    mle_avg[h,j] = mean(mle_a[100,j,h,])
    mle_sd[h,j] = sd(mle_a[100,j,h,])
  }
}  


bias = mle_avg - Py
c_vector =c('blue', 'red','yellow')
for (h in 1:T){
  points(c(9.5+(h-1)*18,10.5+(h-1)*18), bias[h,], col=c_vector)
  arrows(c(9.5+(h-1)*18,10.5+(h-1)*18), bias[h,]-2*mle_sd[h,], c(9.5+(h-1)*18,10.5+(h-1)*18), 
         bias[h,]+2*mle_sd[h,], code=3, angle=90, lty=8, col=c_vector)
}


#treatment effect bias at n=6400

mle_t_avg = vector(length=3)
mle_t_sd = vector(length=3)

treatment_effect = Py.ge[,1] - Py.ge[,2]

for (h in 1:T){
  mle_t_avg[h] = mean(((mle_a.ge[80,1, h,] - mle_a.ge[80,2,h,])- treatment_effect[h])/treatment_effect[h])
  mle_t_sd[h] = sd(((mle_a.ge[80,1, h,] - mle_a.ge[80,2,h,])- treatment_effect[h])/treatment_effect[h])
}


sample_size = round(6400*Px.ge)


upper = mle_t_avg + 1.96*mle_t_sd/sqrt(sample_size)
mid = mle_t_avg
lower = mle_t_avg - 1.96*mle_t_sd/sqrt(sample_size)


for (h in 1:T){
  points(1.5+(h-1)*4, mid[h], col='blue', lwd=6)
  arrows(c(1.5+(h-1)*4), lower[h], c(1.5+(h-1)*4),
         upper[h], code=3, angle=90, lty=8, col='blue', lwd=2)
}

