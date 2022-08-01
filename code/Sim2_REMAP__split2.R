
set.seed(20)

reps = 1000

K     = 100 #number of stages in trial
S     = 2
T     = 3

nk_.h   = 80 #patients per stage
Nx.h    = c(501, 1279, 324) + c(1034, 2604, 683)
Px.h    = Nx.h/sum(Nx.h)
nk_.h   = round(nk_.h*Px.h)
Py.h   = matrix(c(0.140, 0.262, 0.414, 0.178, 0.233, 0.293), nrow=3)
#Py    = c(0.25,0.22) #probability of death on each treatment
#S     = length(Py[1,]) #Number of treatments
#T     = length(Px) #Number of care

pre_ratio.h <- vector(length=S)

alpha.h       = 1 #uninformative beta
beta.h        = 1 #uninformative beta
Y.h           = array(dim=c(K,S,T))
postmean.h    = array(dim=c(K,S,T))
mle.h         = array(dim=c(K,S,T))
PostAlpha.h   = array(dim=c(K,S,T))
PostBeta.h    = array(dim=c(K,S,T))

nk.h          = array(dim=c(K,S,T))
for (h in 1:length(nk_.h)){
  nk.h[,,h] = round(nk_.h[h]/S)
}
rand_ratio.h = matrix(0.5, nrow=S, ncol=T)
rand_ratios.h = array(dim=c(K,S,T))

UPPER_LIMIT <- 0.9
LOWER_LIMIT <- 0.1

nka.h = array(dim = c(K, S,T, reps))
Ya.h = array(dim = c(K, S, T,reps))
postmean_a.h = array(dim = c(K, S, T,reps))
mle_a.h = array(dim=c(K, S,T, reps))

X.h = array(dim=c(100,S, T))

prob_optimal.h = vector(length=S)
prob_optimals.h = array(dim=c(reps, S, K, T))

#y = vector(length=3)

conf_sup.h <- array(dim=c(K/10, T,reps))
#conf_sup2 <- array(dim=c(K/10, T, reps))
# Stage 1
for (n in 1:reps){
  rand_ratio.h = matrix(0.5, nrow=S, ncol=T)
  for (h in 1:T){
    for(j in 1:S){
      y.h = rbinom(1, round(nk.h[1,j,h]), Py.h[h,j])
      Y.h[1,j,h]         = y.h
      PostAlpha.h[1,j,h] = alpha.h + Y.h[1,j,h]
      PostBeta.h[1,j,h]  = beta.h  + nk.h[1,j,h] -Y.h[1,j,h]
      postmean.h[1,j,h]  = PostAlpha.h[1,j,h]/(PostAlpha.h[1,j,h]+PostBeta.h[1,j,h])
      mle.h[1,j,h]       = sum(Y.h[1,j,h])/(nk.h[1,j,h])
    }
  }
  
for(i in 2:K){

for (h in 1:T){
for(j in 1:S){
  nk.h[i, j,h] = rbinom(1,nk_.h[h], rand_ratio.h[j,h])
  y.h = rbinom(1, round(nk.h[i,j,h]), Py.h[h,j])

Y.h[i,j,h]         = y.h
PostAlpha.h[i,j,h] = alpha.h + sum(Y.h[1:i,j,h], na.rm=T)
PostBeta.h[i,j,h]  = beta.h + sum(nk.h[1:i,j,h]-Y.h[1:i,j,h], na.rm=T)
postmean.h[i,j,h]  = PostAlpha.h[i,j,h]/(PostAlpha.h[i,j,h]+PostBeta.h[i,j,h])
mle.h[i,j,h]       = sum(Y.h[1:i,j,h])/sum(nk.h[1:i,j,h])

}
if (i > 28){
    
if (i %% 7 == 0){

for (j in 1:S){
  X.h[,j, h] = rbeta(100, PostAlpha.h[i-28, j,h], PostBeta.h[i-28, j,h])
}
optimum = apply(X.h[,,h], 1, which.min)
for (j in 1:S){
  prob_optimal.h[j] = mean(optimum == j)
  pre_ratio.h[j] = sqrt(prob_optimal.h[j]/(sum(nk.h[1:(i-28),j,h])+1))
}
prob_optimals.h[n,,i,h] = prob_optimal.h

rand_ratio.h[,h] = pre_ratio.h/sum(pre_ratio.h)
rand_ratios.h[i,,h] = rand_ratio.h[,h]

#making sure samples don't go to 0
for (j in 1:S){
if (rand_ratio.h[j]>1-(S-1)*LOWER_LIMIT){
  rand_ratio.h[j] = UPPER_LIMIT
}
if (rand_ratio.h[j]<LOWER_LIMIT){
  rand_ratio.h[j] = LOWER_LIMIT
}
}
}
}
}
  
if (i %% 10 == 0){
  for (h in 1:T){
  #pat_df <- data.frame(t=rep(1, i*nk_.h[h]), o=rep(0, i*nk_.h[h]), 
  #                     c=rep(as.integer(h==2), i*nk_.h[h]), v=rep(as.integer(h==3), i*nk_[h]))
  pat_df <- data.frame(t=rep(1, i*nk_.h[h]), o=rep(0, i*nk_.h[h]))
  pat_df[1:sum(nk.h[1:i,1,h]), 't'] = 0
  pat_df[1:(sum(Y.h[1:i,1,h])), 'o'] = 1
  pat_df[(sum(nk.h[1:i,1,h])+1):(sum(nk.h[1:i,1,h])+sum(Y.h[1:i,2,h])), 'o'] = 1
  
  
  fit = glm(o~t,data=pat_df, family=binomial)
  conf_sup.h[i/10,h,n] = summary(fit)$coef['t','Pr(>|z|)']
  
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

nka.h[,,,n] = nk.h
Ya.h[,,,n] = Y.h
postmean_a.h[,,,n] = postmean.h
mle_a.h[,,,n] = mle.h
}

sum(Ya.h)/reps

nk_avg.h = array(dim=c(K,S,T))
Y_avg.h = array(dim=c(K,S,T))
nk_sd.h = array(dim=c(K,S,T))
Y_sd.h = array(dim=c(K,S,T))
for (h in 1:T){
  for (j in 1:S){
    nk_avg.h[,j,h] = apply(nka.h[,j,h,], 1, mean)
    Y_avg.h[,j,h] = apply(Ya.h[,j,h,], 1, mean, na.rm=T)
    nk_sd.h[,j,h] = apply(nka.h[,j,h,], 1, sd)
    Y_sd.h[,j,h] = apply(Ya.h[,j,h,], 1, sd, na.rm=T)
  }
}

for (h in 1:T){
  print(sum(nk_avg.h[1:80,1,h]))
  print(sum(nk_avg.h[1:80,2,h]))
}

for (h in 1:T){
  print(sum(Y_avg.h[,1,h]))
  print(sum(Y_avg.h[,2,h]))
}


#mortality at n=6450
sum(nk_avg.h[1:80,,])
sum(Y_avg.h[1:80,,])
sd(apply(Ya.h[1:80, 1,1, ] + Ya.h[1:80, 2,1, ] +
           Ya.h[1:80, 1,2, ] + Ya.h[1:80, 2,2, ] +
           Ya.h[1:80, 1,3, ] + Ya.h[1:80, 2,3, ], 2, sum))/(reps ** 0.5)




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






#mean squared error
true_t_effect = Py.h[,1] - Py.h[,2]
mse.h = apply((mle_a.h[43,1,,] - mle_a.h[43,2,,] - true_t_effect)**2, 1, mean)
mse.h







#treatment effect bias at n=6400
treatment_effect = Py.g[,1] - Py.g[,2]

for (h in 1:T){
  mle_t_avg[h] = mean(((mle_a.h[80,1, h,] - mle_a.h[80,2,h,])- treatment_effect[h])/treatment_effect[h])
  mle_t_sd[h] = sd(((mle_a.h[80,1, h,] - mle_a.h[80,2,h,])- treatment_effect[h])/treatment_effect[h])
}


upper = mle_t_avg + 1.96*mle_t_sd/sqrt(1000)
mid = mle_t_avg
lower = mle_t_avg - 1.96*mle_t_sd/sqrt(1000)

h=1
plot(2.5+(h-1)*4, mid[h],xlim = c(0,12), ylim = c(-0.1, 0.1), #xlim 16 if plotting Tf and RMCf
     col='red', xaxt='n', ylab='Relative Bias', xlab='Patient group', lwd=6)
for (h in 1:T){
  points(2.5+(h-1)*4, mid[h], col='red', lwd=6)
  arrows(2.5+(h-1)*4, lower[h], 2.5+(h-1)*4, 
         upper[h], code=3, angle=90, lty=8, col='red', lwd=2)
}
