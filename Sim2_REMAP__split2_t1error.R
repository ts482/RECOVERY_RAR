
set.seed(20)

reps = 1000

K     = 100 #number of stages in trial
S     = 2
T     = 3

nk_.he   = 80 #patients per stage
Nx.he    = c(501, 1279, 324) + c(1034, 2604, 683)
Px.he    = Nx.he/sum(Nx.he)
nk_.he   = round(nk_.he*Px.he)
Py.he   = matrix(c(0.140, 0.262, 0.414, 0.140, 0.262, 0.414), nrow=3)
#Py    = c(0.25,0.22) #probability of death on each treatment
#S     = length(Py[1,]) #Number of treatments
#T     = length(Px) #Number of care

pre_ratio.he <- vector(length=S)

alpha.he       = 1 #uninformative beta
beta.he        = 1 #uninformative beta
Y.he           = array(dim=c(K,S,T))
postmean.he    = array(dim=c(K,S,T))
mle.he         = array(dim=c(K,S,T))
PostAlpha.he   = array(dim=c(K,S,T))
PostBeta.he    = array(dim=c(K,S,T))

nk.he          = array(dim=c(K,S,T))
for (h in 1:length(nk_.he)){
  nk.he[,,h] = round(nk_.he[h]/S)
}
rand_ratio.he = matrix(0.5, nrow=S, ncol=T)
rand_ratios.he = array(dim=c(K,S,T))

UPPER_LIMIT <- 0.9
LOWER_LIMIT <- 0.1

nka.he = array(dim = c(K, S,T, reps))
Ya.he = array(dim = c(K, S, T,reps))
postmean_a.he = array(dim = c(K, S, T,reps))
mle_a.he = array(dim=c(K, S,T, reps))

X.he = array(dim=c(100,S, T))

prob_optimal.he = vector(length=S)
prob_optimals.he = array(dim=c(reps, S, K, T))

#y = vector(length=3)

conf_sup.he <- array(dim=c(K/10, T,reps))
conf_sup_dir.he <- array(dim=c(K/10, T,reps))


# Stage 1
for (n in 1:reps){
  rand_ratio.he = matrix(0.5, nrow=S, ncol=T)
  for (h in 1:T){
    for(j in 1:S){
      y.he = rbinom(1, round(nk.he[1,j,h]), Py.he[h,j])
      Y.he[1,j,h]         = y.he
      PostAlpha.he[1,j,h] = alpha.he + Y.he[1,j,h]
      PostBeta.he[1,j,h]  = beta.he  + nk.he[1,j,h] -Y.he[1,j,h]
      postmean.he[1,j,h]  = PostAlpha.he[1,j,h]/(PostAlpha.he[1,j,h]+PostBeta.he[1,j,h])
      mle.he[1,j,h]       = sum(Y.he[1,j,h])/(nk.he[1,j,h])
    }
  }
  
for(i in 2:K){

for (h in 1:T){
for(j in 1:S){
  nk.he[i, j,h] = rbinom(1,nk_.he[h], rand_ratio.he[j,h])
  y.he = rbinom(1, round(nk.he[i,j,h]), Py.he[h,j])

Y.he[i,j,h]         = y.he
PostAlpha.he[i,j,h] = alpha.he + sum(Y.he[1:i,j,h], na.rm=T)
PostBeta.he[i,j,h]  = beta.he + sum(nk.he[1:i,j,h]-Y.he[1:i,j,h], na.rm=T)
postmean.he[i,j,h]  = PostAlpha.he[i,j,h]/(PostAlpha.he[i,j,h]+PostBeta.he[i,j,h])
mle.he[i,j,h]       = sum(Y.he[1:i,j,h])/sum(nk.he[1:i,j,h])

}
if (i > 28){
    
if (i %% 7 == 0){

for (j in 1:S){
  X.he[,j, h] = rbeta(100, PostAlpha.he[i-28, j,h], PostBeta.he[i-28, j,h])
}
optimum = apply(X.he[,,h], 1, which.min)
for (j in 1:S){
  prob_optimal.he[j] = mean(optimum == j)
  pre_ratio.he[j] = sqrt(prob_optimal.he[j]/(sum(nk.he[1:(i-28),j,h])+1))
}
prob_optimals.he[n,,i,h] = prob_optimal.he

rand_ratio.he[,h] = pre_ratio.he/sum(pre_ratio.he)
rand_ratios.he[i,,h] = rand_ratio.he[,h]

#making sure samples don't go to 0
for (j in 1:S){
if (rand_ratio.he[j]>1-(S-1)*LOWER_LIMIT){
  rand_ratio.he[j] = UPPER_LIMIT
}
if (rand_ratio.he[j]<LOWER_LIMIT){
  rand_ratio.he[j] = LOWER_LIMIT
}
}
}
}
}
  
if (i %% 10 == 0){
  for (h in 1:T){
  #pat_df <- data.frame(t=rep(1, i*nk_.he[h]), o=rep(0, i*nk_.he[h]), 
  #                     c=rep(as.integer(h==2), i*nk_.he[h]), v=rep(as.integer(h==3), i*nk_[h]))
  pat_df <- data.frame(t=rep(1, i*nk_.he[h]), o=rep(0, i*nk_.he[h]))
  pat_df[1:sum(nk.he[1:i,1,h]), 't'] = 0
  pat_df[1:(sum(Y.he[1:i,1,h])), 'o'] = 1
  pat_df[(sum(nk.he[1:i,1,h])+1):(sum(nk.he[1:i,1,h])+sum(Y.he[1:i,2,h])), 'o'] = 1
  
  
  fit = glm(o~t,data=pat_df, family=binomial)
  conf_sup.he[i/10,h,n] = summary(fit)$coef['t','Pr(>|z|)']
  conf_sup_dir.he[i/10,h,n] = summary(fit)$coef['t','Estimate']
  
  #if (h==1){
  #  master_df = pat_df
  #} else {
  #  master_df <- rbind(master_df, pat_df)
  #}
  
  }
  #fit = glm(o~t*c+t*v, data = master_df, family=binomial)
  #conf_sup2[i/10,,n] = summary(fit)$coef[c('t','t:c', 't:v'), 'Pr(>|z|)']
}
}

nka.he[,,,n] = nk.he
Ya.he[,,,n] = Y.he
postmean_a.he[,,,n] = postmean.he
mle_a.he[,,,n] = mle.he
}

sum(Ya.he)/reps

mean(apply(conf_sup.he[8,,]<0.05, 2, any))  #0.213
mean(apply(conf_sup.he[8,,]<(0.05/3), 2, any))  #0.089 #bonferroni correction



nk_avg.he = array(dim=c(K,S,T))
Y_avg.he = array(dim=c(K,S,T))
for (h in 1:T){
  for (j in 1:S){
    nk_avg.he[,j,h] = apply(nka.he[,j,h,], 1, mean)
    Y_avg.he[,j,h] = apply(Ya.he[,j,h,], 1, mean, na.rm=T)
  }
}

for (h in 1:T){
  print(sum(nk_avg.he[1:80,1,h]))
  print(sum(nk_avg.he[1:80,2,h]))
}

for (h in 1:T){
  print(sum(Y_avg.he[,1,h]))
  print(sum(Y_avg.he[,2,h]))
}


#mortality at n=6450
sum(nk_avg.he[1:80,,])
sum(Y_avg.he[1:80,,])





#plot(nk_avg[,1], ylim = c(0,50), col='red')
#points(nk_avg[,2], col='blue')
#legend("bottomleft",legend = c('no treatment', 'dexamethasone'),col = c("red", "blue"), lwd=4)

for (h in 1:T){
  plot(nk_avg[,1,h], ylim = c(0,nk_[h]), col='red')
  points(nk_avg[,2,h], col='blue')
  legend("bottomleft",legend = c('no treatment', 'dexamethasone'),col = c("red", "blue"), lwd=4)
}

sum(Ya)/reps


treatment_superior = matrix(nrow=K, ncol=T)
for (h in 1:T){
  for (i in 1:K){
    treatment_superior[i,h] = mean(postmean_a[i,1,h,] >= postmean_a[i,2,h,])
  }
}

for (h in 1:T){
  plot(treatment_superior[,h])
}


thetas_avg = apply(prob_optimals[,2,], 2, mean)
plot(thetas_avg)
above_theta = prob_optimals[,2,] > 0.95
above_theta_avg = apply(above_theta, 2, mean)
plot(above_theta_avg)

for (h in 1:T){
  thetas_avg = apply(prob_optimals[,2,,h], 2, mean)
  plot(thetas_avg)
}

for (h in 1:T){
  above_theta = prob_optimals[,2,,h] > 0.95
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
c_vector =c('blue', 'red','brown')
for (h in 1:T){
  points(c(13.5+(h-1)*18,14.5+(h-1)*18), bias[h,], col=c_vector)
  arrows(c(13.5+(h-1)*18,14.5+(h-1)*18), bias[h,]-2*mle_sd[h,], c(13.5+(h-1)*18,14.5+(h-1)*18), 
         bias[h,]+2*mle_sd[h,], code=3, angle=90, lty=8, col=c_vector)
}


tick_intervals = c(2,6,10,14)
axis(1, labels = rep(c('FeR', 'FuR', 'Ts', 'RMCs'), 3), 
     at = c(tick_intervals, tick_intervals+18, tick_intervals+36))
abline(v = c(17, 35))
legend("bottomleft",legend = c('no treatment', 'dexamethasone'),col = c("blue", "red"), lwd=4)
text(x=9, y=0.12, labels='subgroup i')
text(x=27, y=0.12, labels='subgroup ii')
text(x=43, y=0.12, labels='subgroup iii')


#treatment effect bias at n=6400
treatment_effect = Py.g[,1] - Py.g[,2]

for (h in 1:T){
  mle_t_avg[h] = mean(((mle_a.he[80,1, h,] - mle_a.he[80,2,h,])- treatment_effect[h])/treatment_effect[h])
  mle_t_sd[h] = sd(((mle_a.he[80,1, h,] - mle_a.he[80,2,h,])- treatment_effect[h])/treatment_effect[h])
}

sample_size = round(6400*Px.he)


upper = mle_t_avg + 1.96*mle_t_sd/sqrt(sample_size)
mid = mle_t_avg
lower = mle_t_avg - 1.96*mle_t_sd/sqrt(sample_size)

h=1
plot(2.5+(h-1)*4, mid[h],xlim = c(0,16), ylim = c(-0.1, 0.1),
     col='red', xaxt='n', ylab='Relative Bias', xlab='Patient group', lwd=6)
for (h in 1:T){
  points(2.5+(h-1)*4, mid[h], col='red', lwd=6)
  arrows(2.5+(h-1)*4, lower[h], 2.5+(h-1)*4, 
         upper[h], code=3, angle=90, lty=8, col='red', lwd=2)
}

