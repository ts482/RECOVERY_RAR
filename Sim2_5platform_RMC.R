
set.seed(20)

reps = 1000

K     = 100 #number of stages in trial
Si    = 4 #number of treatments
T     = 3 #number of subgroups

#3 treatment groups at ratio of 2:1:1:1
#therefore randomisation aim will be to have treatment groups at 1.5, control group at 1

nk_.m  = 120 #patients per stage so ar days =80, sample size = 9600 (as opposed to 9602)


rand_ratio.m  = c(0.4, 0.2, 0.2, 0.2)

Nx.m1   = c(501, 1279, 324) + c(1034, 2604, 683) #dex group, control group as reported for dex   #no ox, ox, vent
Nx.m2  = c(362, 938, 261) + c(425, 1131, 60)       #hydroxychloroquine, lopinavir
N.control_group = c(1034, 2604, 683)
N.treatment_groups = c(362, 938, 261) + c(425, 1131, 60) + c(501, 1279, 324)
Nx.m = Nx.m1 + Nx.m2
Px.m   = Nx.m/sum(Nx.m)
Py.m1   = matrix(c(0.140, 0.262, 0.414, 0.178, 0.233, 0.293), nrow=3)
Py.m2  = matrix(c(0.160, 0.270, 0.421, 0.167, 0.247, 0.4), nrow = 3)         #hydroxychloroquine, lopinavir
Py.m = matrix(c(Py.m1, Py.m2), nrow=3)
#Py    = c(0.25,0.22) #probability of death on each treatment
#S     = length(Py[1,]) #Number of treatments
#T     = length(Px) #Number of care


alpha.m       = 1 #uninformative beta
beta.m        = 1 #uninformative beta
Y.m           = matrix(nrow=K,ncol=Si)
postmean.m    = matrix(nrow=K,ncol=Si)
mle.m         = matrix(nrow=K,ncol=Si)
PostAlpha.m   = matrix(nrow=K,ncol=Si)
PostBeta.m    = matrix(nrow=K,ncol=Si)

nk.m          = matrix(round(nk_.m*rand_ratio.m), nrow=K,ncol=Si, byrow=T)
rand_ratios.m = matrix(nrow=K, ncol=Si)

UPPER_LIMIT <- 0.5
LOWER_LIMIT <- 0.05

nka.m = array(dim = c(K, Si, reps))
Ya.m = array(dim = c(K, Si, reps))
postmean_a.m = array(dim = c(K, Si, reps))
mle_a.m = array(dim=c(K, Si, reps))

y.m = vector(length=3)

conf_sup.m <- array(dim = c(Si-1,K/10, reps))

prob_optimal.m = vector(length=Si-1)
prob_optimals.m = array(dim=c(reps, Si-1, K))
pre_ratio.m = vector(length=Si-1)


X.m = matrix(nrow=100, ncol=Si)
# Stage 1
for (n in 1:reps){

rand_ratio.m = c(0.4, 0.2, 0.2, 0.2)

for(j in 1:Si){
for (h in 1:T){
  y.m[h] = rbinom(1, round(nk.m[1,j]* Px.m[h]), Py.m[h,j])
}
Y.m[1,j]         = sum(y.m)
PostAlpha.m[1,j] = alpha.m + Y.m[1,j]
PostBeta.m[1,j]  = beta.m  + nk.m[1,j] -Y.m[1,j]
postmean.m[1,j]  = PostAlpha.m[1,j]/(PostAlpha.m[1,j]+PostBeta.m[1,j])
mle.m[1,j]       = sum(Y.m[1,j])/(nk.m[1,j])
}


for(i in 2:K){
  
for(j in 1:Si){
  nk.m[i, j] = round(nk_.m * rand_ratio.m[j])
for (h in 1:T){
  y.m[h] = rbinom(1, round(nk.m[i,j]* Px.m[h]), Py.m[h,j])
}
Y.m[i,j]         = sum(y.m)
PostAlpha.m[i,j] = alpha.m + sum(Y.m[1:i,j], na.rm=T)
PostBeta.m[i,j]  = beta.m + sum(nk.m[1:i,j]-Y.m[1:i,j], na.rm=T)
postmean.m[i,j]  = PostAlpha.m[i,j]/(PostAlpha.m[i,j]+PostBeta.m[i,j])
mle.m[i,j]       = sum(Y.m[1:i,j])/sum(nk.m[1:i,j])


}

if (i > 28){

if (i %% 7 == 0){
  
  for (j in 1:Si){
    X.m[,j] = rbeta(100, PostAlpha.m[i-28, j], PostBeta.m[i-28, j])
  }
optimum = apply(X.m, 1, which.min)
k = i/(K)

for (j in 2:Si){
  prob_optimal.m[j-1] = mean(optimum == j)
  pre_ratio.m[j-1] = sqrt(prob_optimal.m[j-1]/(sum(nk.m[1:(i-28),j])+1))
}

for (j in 2:Si){
  rand_ratio.m[j] = pre_ratio.m[j-1]*0.6/sum(pre_ratio.m)
}
while (min(rand_ratio.m)<0.05){
  maxi = which.max(rand_ratio.m[2:Si]) + 1
  mini = which.min(rand_ratio.m[2:Si]) + 1
  diff = 0.05 - rand_ratio.m[mini]
  rand_ratio.m[maxi] = rand_ratio.m[maxi] - diff
  rand_ratio.m[mini] = rand_ratio.m[mini] + diff
}
}
rand_ratios.m[i,] = rand_ratio.m

}
if (i %% 10 == 0){
  #calculating power for dex
  pat_df <- data.frame(t=rep(1, sum(nk.m[1:i, c(1,2)])), 
                       o=rep(0, sum(nk.m[1:i, c(1,2)])))
  pat_df[1:sum(nk.m[1:i,1]), 't'] = 0
  pat_df[1:(sum(Y.m[1:i,1])), 'o'] = 1
  pat_df[(sum(nk.m[1:i,1])+1):(sum(nk.m[1:i,1])+sum(Y.m[1:i,2])+1), 'o'] = 1
  
  fit = glm(o~t,data=pat_df, family=binomial)
  conf_sup.m[1,i/10,n] = summary(fit)$coef['t','Pr(>|z|)']
  
  
  #calculating power for hydroxy
  pat_df <- data.frame(t=rep(1, sum(nk.m[1:i, c(1,3)])), 
                       o=rep(0, sum(nk.m[1:i, c(1,3)])))
  pat_df[1:sum(nk.m[1:i,1]), 't'] = 0
  pat_df[1:(sum(Y.m[1:i,1])), 'o'] = 1
  pat_df[(sum(nk.m[1:i,1])+1):(sum(nk.m[1:i,1])+sum(Y.m[1:i,3])+1), 'o'] = 1
  
  fit = glm(o~t,data=pat_df, family=binomial)
  conf_sup.m[2,i/10,n] = summary(fit)$coef['t','Pr(>|z|)']
  
  
  #calculating power for lopinavir
  pat_df <- data.frame(t=rep(1, sum(nk.m[1:i, c(1,4)])), 
                       o=rep(0, sum(nk.m[1:i, c(1,4)])))
  pat_df[1:sum(nk.m[1:i,1]), 't'] = 0
  pat_df[1:(sum(Y.m[1:i,1])), 'o'] = 1
  pat_df[(sum(nk.m[1:i,1])+1):(sum(nk.m[1:i,1])+sum(Y.m[1:i,4])+1), 'o'] = 1
  
  fit = glm(o~t,data=pat_df, family=binomial)
  conf_sup.m[3,i/10,n] = summary(fit)$coef['t','Pr(>|z|)']
}
}


nka.m[,,n] = nk.m
Ya.m[,,n] = Y.m
postmean_a.m[,,n] = postmean.m
mle_a.m[,,n] = mle.m

}




Y_avg.m = matrix(nrow=K, ncol=Si)
nk_avg.m = matrix(nrow=K, ncol=Si)
Y_sd.m = matrix(nrow=K, ncol=Si)
nk_sd.m = matrix(nrow=K, ncol=Si)
for (j in 1:Si){
  nk_avg.m[,j] = apply(nka.m[,j,], 1, mean)
  Y_avg.m[,j] = apply(Ya.m[,j,], 1, mean, na.rm=T)
  nk_sd.m[,j] = apply(nka.m[,j,], 1, sd)
  Y_sd.m[,j] = apply(Ya.m[,j,], 1, sd, na.rm=T)
}


lines(nk_avg.m[c(1,seq(10, 100, 10)),2]/nk_.m*100, col='blue', type='b', lty='dashed', lwd=3)
arrows(5:11 + 0.05, (nk_avg.m[seq(40,100,10),2] - nk_sd.m[seq(40,100,10),2])/nk_.i*100,
       y1 = (nk_avg.m[seq(40,100,10),2] + nk_sd.m[seq(40,100,10),2])/nk_.i*100,
       code=3, angle=90, lty=9, col='blue', lwd=2)
lines(nk_avg.m[c(1,seq(10, 100, 10)),3]/nk_.m*100, col='yellow', type='b', lty='dashed', lwd=3)
arrows(5:11 + 0.05, (nk_avg.m[seq(40,100,10),3] - nk_sd.m[seq(40,100,10),3])/nk_.i*100,
       y1 = (nk_avg.m[seq(40,100,10),3] + nk_sd.m[seq(40,100,10),3])/nk_.i*100,
       code=3, angle=90, lty=9, col='yellow', lwd=2)
lines(nk_avg.m[c(1,seq(10, 100, 10)),4]/nk_.m*100, col='green', type='b', lty='dashed', lwd=3)
arrows(5:11 + 0.05, (nk_avg.m[seq(40,100,10),4] - nk_sd.m[seq(40,100,10),4])/nk_.i*100,
       y1 = (nk_avg.m[seq(40,100,10),4] + nk_sd.m[seq(40,100,10),4])/nk_.i*100,
       code=3, angle=90, lty=9, col='green', lwd=2)

legend('topright', legend = c(expression('T'[f]), expression('RMC'[f])), lwd =3, lty = c('solid','dashed'))
legend('top', legend=c(expression('T'[f]*' SD'), expression('RMC'[f]*' SD')), lty = c(8,9))


#lines(nk_avg.m[c(1,1:10*10),2]/nk_.m*100, lwd=4, type='b', col='brown')


for (j in 1:Si){
  print(sum(nk_avg.m[1:80,j]))
}


for (j in 1:Si){
print(sum(Y_avg.m[1:80,j]))
}

for (j in 1:Si){
  print(sum(Y_avg.m[,j])/sum(nk_avg.m[,j]))
}


#deaths at n=9600
sum(nk_avg.m[1:80,])
sum(Y_avg.m[1:80,])
sd(apply(Ya.m[1:80, 1, ] + Ya.m[1:80, 2, ] + 
           Ya.m[1:80, 3, ] + Ya.m[1:80, 4, ], 2, sum))/(reps ** 0.5)


sum(Y_avg.m[1:80,])/sum(nk_avg.m[1:80,])



#power calculation and plotting

treatment_power_avg = matrix(nrow=11, ncol=3)
treatment_power_err = matrix(nrow=10, ncol=3)
for (j in 1:3){
  treatment_power_avg[,j] = c(0, apply(conf_sup.m[j,,]<0.05, 1, mean))
  treatment_power_err[,j] = apply(conf_sup.m[j,,]<0.05, 1, sd)/(reps**0.5)
  
}




points(treatment_power_avg[,1], col='red', type='b', lwd=2)
#points(treatment_power_avg[,2], col='red', type='b', lwd=2, lty='dashed')
points(treatment_power_avg[,3], col='red', type='b', lwd=2, lty='dashed')

arrows(2:11, treatment_power_avg[2:11,1] - treatment_power_err[,1],
       y1 = treatment_power_avg[2:11,1] + treatment_power_err[,1],
       code=3, angle=90, lty=8, col='red', lwd=2.5)
arrows(2:11, treatment_power_avg[2:11,3] - treatment_power_err[,3],
       y1 = treatment_power_avg[2:11,3] + treatment_power_err[,3],
       code=3, angle=90, lty=9, col='red', lwd=2.5)




plot(nk_avg[,1], ylim = c(0,60), col='red', xlab='Trial completion (%)', ylab='Patients allocated')
points(nk_avg[,2], col='blue')
legend("bottomleft",legend = c('no treatment', 'dexamethasone'),col = c("red", "blue"), lwd=4)

sum(Ya)/reps
#Y_avg = matrix(nrow=K, ncol=Si)
#for (j in 1:Si){
#  Y_avg[,j] = apply(Ya[,j,], 1, mean)
#}

treatment_superior = vector(length=K)
for (i in 1:K){
  treatment_superior[i] = mean(postmean_a[i,1,] >= postmean_a[i,2,])
}

plot(treatment_superior)

thetas_avg = apply(thetas, 2, mean)
plot(thetas_avg)
above_theta = thetas > 0.95
above_theta_avg = apply(above_theta, 2, mean)
plot(above_theta_avg)
#plot(postmean[,1],ylim=range(postmean),type="l",lwd=4,
#     ylab="Posterior mean",xlab="Interim Analysis")
#lines(postmean[,2],type="l",lwd=4,col="red")
#legend("bottomright",c("Treatment 1","Treatment 2"),col=c("black","red"),lwd=4,bty="n")

#lines(mle[,1],ylim=range(mle),type="l",lwd=4,lty=2)
#lines(mle[,2],type="l",lwd=4,col="red",lty=2)


#trial progress bias

mle_avg = vector(length= 10) #vector(length=Si)
mle_sd = vector(length=10)#vector(length=Si)

treatment_effect = 0.257 - 0.229

for (i in 1:10){
  mle_avg[i] = mean(((mle_a.m[i*10,1,] - mle_a.m[i*10,2,])-treatment_effect)/treatment_effect)
  mle_sd[i] = sd(((mle_a.m[i*10,1,] - mle_a.m[i*10,2,])-treatment_effect)/treatment_effect)
}



upper = mle_avg + 1.96*mle_sd/sqrt(6450)
mid = mle_avg
lower = mle_avg - 1.96*mle_sd/sqrt(6450)

points(1:10-0.2, mid, ylim=c(-0.06,0.06), col= 'blue', lwd=6)
arrows(1:10-0.2, lower, 1:10-0.2, upper, code=3, angle=90, lty=8, col='blue', lwd=3)

legend('bottomright', legend = c(expression('T'[f]), expression('RMC'[f])),col = c('blue','red'), 
       lty = c('dashed', 'dashed'), lwd=3)




# end of trial treatment effect bias


treatment_effect = 0.257 - 0.229
mle_t_avg = mean(((mle_a.m[80,1,] - mle_a.m[80,2,])-treatment_effect)/treatment_effect)
mle_t_sd = sd(((mle_a.m[80,1,] - mle_a.m[80,2,])-treatment_effect)/treatment_effect)


upper = mle_t_avg + 1.96*mle_t_sd/sqrt(6450)
mid = mle_t_avg
lower = mle_t_avg - 1.96*mle_t_sd/sqrt(6450)

points(13.5, mid, col = 'blue', lwd=6)
arrows(13.5, lower, 13.5, upper, code=3, angle=90, lty=8, col='blue', lwd=2)





#mean squared error

arm_estimates.m = apply(Px.m * Py.m, 2, sum)

treatment_effects = arm_estimates.m[1] - arm_estimates.m[2:4]

treatment_effects

mse.m = vector(length=3)

for (j in 1:3){
  mse.m[j] = mean((mle_a.m[80,1,] - mle_a.m[80,j+1,]- treatment_effects[j])**2)
}
mse.m






#end of trial bias

mle_t_avg = vector(length=3)
mle_t_sd = vector(length=3)


arm_estimates.m = apply(Px.m * Py.m, 2, sum)

treatment_effects = arm_estimates.m[1] - arm_estimates.m[2:4]

treatment_effects


for (j in 1:3){
  mle_t_avg[j] = mean(((mle_a.m[80,1,] - mle_a.m[80,j+1,])- treatment_effects[j])/treatment_effects[j])
  mle_t_sd[j] = sd(((mle_a.m[80,1,] - mle_a.m[80,j+1,])- treatment_effects[j])/treatment_effects[j])
}



upper = mle_t_avg + 1.96*mle_t_sd/sqrt(1000)
mid = mle_t_avg
lower = mle_t_avg - 1.96*mle_t_sd/sqrt(1000)

for (j in 1:3){
  points(1.125+(j-1), mid[j], col='red', lwd=6)
  arrows(1.125+c(j-1), lower[j], 1.125+c(j-1),
         upper[j], code=3, angle=90, lty=8, col='red', lwd=2)
}

axis(1, labels = c('dexamethasone','hydroxychloroquine','lopinavir'), at = c(1,2,3))
legend('bottomright', legend = c(expression('T'[f]), expression('RMC'[f])),col = c('blue','red'), lty = c('dashed','dashed'))

abline(h=0, lty='dashed')
