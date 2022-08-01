
set.seed(20)

reps = 1000

K     = 100 #number of stages in trial
Si    = 4 #number of treatments
T     = 3 #number of subgroups

#3 treatment groups at ratio of 2:1:1:1
#therefore randomisation aim will be to have treatment groups at 1.5, control group at 1

nk_.ke  = 120 #patients per stage so ar days =80, sample size = 9600 (as opposed to 9602)
rand_ratio.ke  = c(0.4, 0.2, 0.2, 0.2)
pre_ratio.ke <- vector(length=Si)


Nx.ke1   = c(501, 1279, 324) + c(1034, 2604, 683) #dex group, control group as reported for dex   #no ox, ox, vent
Nx.ke2  = c(362, 938, 261) + c(425, 1131, 60)       #hydroxychloroquine, lopinavir
N.control_group = c(1034, 2604, 683)
N.treatment_groups = c(362, 938, 261) + c(425, 1131, 60) + c(501, 1279, 324)
Nx.ke = Nx.ke1 + Nx.ke2
Px.ke   = Nx.ke/sum(Nx.ke)
Py.ke1   = matrix(c(0.140, 0.262, 0.414, 0.178, 0.233, 0.293), nrow=3)
Py.ke2  = matrix(c(0.160, 0.270, 0.421, 0.167, 0.247, 0.4), nrow = 3)         #hydroxychloroquine, lopinavir
Py.ke = matrix(c(Py.ke1, Py.ke2), nrow=3)
#Py    = c(0.25,0.22) #probability of death on each treatment
#S     = length(Py[1,]) #Number of treatments
#T     = length(Px) #Number of care


alpha.ke       = 1 #uninformative beta
beta.ke        = 1 #uninformative beta
Y.ke           = matrix(nrow=K,ncol=Si)
postmean.ke    = matrix(nrow=K,ncol=Si)
mle.ke         = matrix(nrow=K,ncol=Si)
PostAlpha.ke   = matrix(nrow=K,ncol=Si)
PostBeta.ke    = matrix(nrow=K,ncol=Si)

nk.ke          = matrix(round(nk_.ke*rand_ratio.ke), nrow=K,ncol=Si, byrow=T)
rand_ratios.ke = matrix(nrow=K, ncol=Si)

UPPER_LIMIT <- 0.5
LOWER_LIMIT <- 0.05

nka.ke = array(dim = c(K, Si, reps))
Ya.ke = array(dim = c(K, Si, reps))
postmean_a.ke = array(dim = c(K, Si, reps))
mle_a.ke = array(dim=c(K, Si, reps))

y.ke = vector(length=3)

conf_sup.ke <- array(dim = c(Si-1,K/10, reps))

prob_optimal.ke = vector(length=Si-1)
prob_optimals.ke = array(dim=c(reps, Si-1, K))
pre_ratio.ke = vector(length=Si-1)


X.ke = matrix(nrow=100, ncol=Si)
# Stage 1
for (n in 1:reps){

rand_ratio.ke = c(0.4, 0.2, 0.2, 0.2)

for(j in 1:Si){
  nk.ke[1, j] = rbinom(1,nk_.ke, rand_ratio.ke[j])
#for (h in 1:T){
 # y.ke[h] = rbinom(1, rbinom(1,nk.ke[1,j], Px.ke[h]), Py.ke[h,j])
#}
Y.ke[1,j]         = rbinom(1,nk.ke[1,j], 0.25)
PostAlpha.ke[1,j] = alpha.ke + Y.ke[1,j]
PostBeta.ke[1,j]  = beta.ke  + nk.ke[1,j] -Y.ke[1,j]
postmean.ke[1,j]  = PostAlpha.ke[1,j]/(PostAlpha.ke[1,j]+PostBeta.ke[1,j])
mle.ke[1,j]       = sum(Y.ke[1,j])/(nk.ke[1,j])
}


for(i in 2:K){
  
for(j in 1:Si){
  nk.ke[i, j] = rbinom(1,nk_.ke, rand_ratio.ke[j])
#for (h in 1:T){
#  y.ke[h] = rbinom(1, rbinom(1,nk.ke[i,j], Px.ke[h]), Py.ke[h,j])
#}
Y.ke[i,j]         = rbinom(1,nk.ke[i,j], 0.25)
PostAlpha.ke[i,j] = alpha.ke + sum(Y.ke[1:i,j], na.rm=T)
PostBeta.ke[i,j]  = beta.ke + sum(nk.ke[1:i,j]-Y.ke[1:i,j], na.rm=T)
postmean.ke[i,j]  = PostAlpha.ke[i,j]/(PostAlpha.ke[i,j]+PostBeta.ke[i,j])
mle.ke[i,j]       = sum(Y.ke[1:i,j])/sum(nk.ke[1:i,j])


}


if (i %% 10 == 0){
  #calculating power for dex
  pat_df <- data.frame(t=rep(1, sum(nk.ke[1:i, c(1,2)])), 
                       o=rep(0, sum(nk.ke[1:i, c(1,2)])))
  pat_df[1:sum(nk.ke[1:i,1]), 't'] = 0
  pat_df[1:(sum(Y.ke[1:i,1])), 'o'] = 1
  pat_df[(sum(nk.ke[1:i,1])+1):(sum(nk.ke[1:i,1])+sum(Y.ke[1:i,2])+1), 'o'] = 1
  
  fit = glm(o~t,data=pat_df, family=binomial)
  conf_sup.ke[1,i/10,n] = summary(fit)$coef['t','Pr(>|z|)']
  
  
  #calculating power for hydroxy
  pat_df <- data.frame(t=rep(1, sum(nk.ke[1:i, c(1,3)])), 
                       o=rep(0, sum(nk.ke[1:i, c(1,3)])))
  pat_df[1:sum(nk.ke[1:i,1]), 't'] = 0
  pat_df[1:(sum(Y.ke[1:i,1])), 'o'] = 1
  pat_df[(sum(nk.ke[1:i,1])+1):(sum(nk.ke[1:i,1])+sum(Y.ke[1:i,3])+1), 'o'] = 1
  
  fit = glm(o~t,data=pat_df, family=binomial)
  conf_sup.ke[2,i/10,n] = summary(fit)$coef['t','Pr(>|z|)']
  
  
  #calculating power for lopinavir
  pat_df <- data.frame(t=rep(1, sum(nk.ke[1:i, c(1,4)])), 
                       o=rep(0, sum(nk.ke[1:i, c(1,4)])))
  pat_df[1:sum(nk.ke[1:i,1]), 't'] = 0
  pat_df[1:(sum(Y.ke[1:i,1])), 'o'] = 1
  pat_df[(sum(nk.ke[1:i,1])+1):(sum(nk.ke[1:i,1])+sum(Y.ke[1:i,4])+1), 'o'] = 1
  
  fit = glm(o~t,data=pat_df, family=binomial)
  conf_sup.ke[3,i/10,n] = summary(fit)$coef['t','Pr(>|z|)']
}
}


nka.ke[,,n] = nk.ke
Ya.ke[,,n] = Y.ke
postmean_a.ke[,,n] = postmean.ke
mle_a.ke[,,n] = mle.ke

}


#type 1 error



t1e_metric.ke = vector(length=10)
for (i in 1:10){
  t1e_metric.ke[i] = mean(apply(conf_sup.ke[,i,]<(0.05/3), 2, any))
}

t1e_metric.ke[8]

lines(t1e_metric.ke*100, lwd=4, type='b', col='green')






Y_avg.ke = matrix(nrow=K, ncol=Si)
nk_avg.ke = matrix(nrow=K, ncol=Si)
Y_sd.ke = vector(length=Si)
nk_sd.ke = vector(length=Si)
for (j in 1:Si){
  nk_avg.ke[,j] = apply(nka.ke[,j,], 1, mean)
  Y_avg.ke[,j] = apply(Ya.ke[,j,], 1, mean, na.rm=T)
  nk_sd.ke[j] = sd(nka.ke[,j,])
  Y_sd.ke[j] = sd(Ya.ke[,j,])
}


#number of patients randomised to dexamethasone
sum(nk_avg.ke[,2])

#lines(nk_avg.ke[c(1,1:10*10),2]/nk_.ke*100, lwd=4, type='b', col='brown')


## sanity check: comparing death rates in the same treatment in different trials
apply(Y_avg.i[1:80,], 2, sum)/apply(nk_avg.i[1:80,], 2, sum)
apply(Y_avg.m[1:80,], 2, sum)/apply(nk_avg.m[1:80,], 2, sum)
apply(Y_avg.l[1:80,], 2, sum)/apply(nk_avg.l[1:80,], 2, sum)
apply(Y_avg.ke[1:80,], 2, sum)/apply(nk_avg.ke[1:80,], 2, sum)



for (j in 1:Si){
  print(sum(nk_avg.ke[,j]))
}


for (j in 1:Si){
print(sum(Y_avg.ke[,j]))
}

for (j in 1:Si){
  print(sum(Y_avg.ke[,j])/sum(nk_avg.ke[,j]))
}


#deaths at n=6450
sum(nk_avg.ke[1:80,])
sum(Y_avg.ke[1:80,])
sd(apply(Ya.ke[1:80, 1, ] + Ya.ke[1:80, 2, ] + 
           Ya.ke[1:80, 3, ] + Ya.ke[1:80, 4, ], 2, sum))/(reps ** 0.5)

sum(Y_avg.ke[1:80,])/sum(nk_avg.ke[1:80,])


#power calculation and plotting

treatment_power_avg = matrix(nrow=11, ncol=3)
treatment_power_err = matrix(nrow=10, ncol=3)
for (j in 1:3){
  treatment_power_avg[,j] = c(0, apply(conf_sup.ke[j,,]<0.05, 1, mean))
  treatment_power_err[,j] = apply(conf_sup.ke[j,,]<0.05, 1, sd)/(reps**0.5)
}




points(treatment_power_avg[,1], col='green', type='b', lwd=2)
#points(treatment_power_avg[,2], col='green', type='b', lwd=2, lty='dashed')
points(treatment_power_avg[,3], col='green', type='b', lwd=2, lty='dashed')

arrows(2:11, treatment_power_avg[2:11,1] - treatment_power_err[,1],
       y1 = treatment_power_avg[2:11,1] + treatment_power_err[,1],
       code=3, angle=90, lty=8, col='green', lwd=2.5)
arrows(2:11, treatment_power_avg[2:11,3] - treatment_power_err[,3],
       y1 = treatment_power_avg[2:11,3] + treatment_power_err[,3],
       code=3, angle=90, lty=9, col='green', lwd=2.5)


abline(h=0.8, lty=4)
abline(h=0.9, lty=4)

legend("topleft",legend = c(expression('T'[f]),expression('RMC'[f]), 'FeR','FuR'),
       col = c('blue','red','brown','green'), lwd=2)
legend('left', legend=c('dexamethasone', 'lopinavir'), 
       lty = c('solid','dashed'))
legend('top', legend=c('dexamethasone SE', 'lopinavir SE'), lty = c(8,9))



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
  mle_avg[i] = mean(((mle_a.ke[i*10,1,] - mle_a.ke[i*10,2,])-treatment_effect)/treatment_effect)
  mle_sd[i] = sd(((mle_a.ke[i*10,1,] - mle_a.ke[i*10,2,])-treatment_effect)/treatment_effect)
}



upper = mle_avg + 1.96*mle_sd/sqrt(6450)
mid = mle_avg
lower = mle_avg - 1.96*mle_sd/sqrt(6450)

points(1:10-0.2, mid, ylim=c(-0.06,0.06), col= 'blue', lwd=6)
arrows(1:10-0.2, lower, 1:10-0.2, upper, code=3, angle=90, lty=8, col='blue', lwd=3)

legend('bottomright', legend = c(expression('T'[f]), expression('RMC'[f])),col = c('blue','red'), 
       lty = c('dashed', 'dashed'), lwd=3)



#mean squared error


arm_estimates.ke = apply(Px.ke * Py.ke, 2, sum)

treatment_effects = arm_estimates.ke[1] - arm_estimates.ke[2:4]

treatment_effects

mse.ke = vector(length=3)

for (j in 1:3){
  mse.ke[j] = mean((mle_a.ke[80,1,] - mle_a.ke[80,j+1,]- treatment_effects[j])**2)
}

mse.ke*10^4








# end of trial treatment effect bias


treatment_effect = 0.257 - 0.229
mle_t_avg = mean(((mle_a.ke[80,1,] - mle_a.ke[80,2,])-treatment_effect)/treatment_effect)
mle_t_sd = sd(((mle_a.ke[80,1,] - mle_a.ke[80,2,])-treatment_effect)/treatment_effect)


upper = mle_t_avg + 1.96*mle_t_sd/sqrt(6450)
mid = mle_t_avg
lower = mle_t_avg - 1.96*mle_t_sd/sqrt(6450)

points(13.5, mid, col = 'blue', lwd=6)
arrows(13.5, lower, 13.5, upper, code=3, angle=90, lty=8, col='blue', lwd=2)

#end of trial bias




mle_avg = vector(length=Si)
mle_sd = vector(length=Si)
for (j in 1:Si){
  mle_avg[j] = mean(mle_a[100,j,])
  mle_sd[j] = sd(mle_a[100,j,])
}

bias = mle_avg - c(0.257, 0.229)
mle_sd
c_vector =c('blue', 'red')
points(c(9.5, 10.5), bias, col=c_vector)
arrows(c(9.5,10.5), bias-2*mle_sd, c(9.5,10.5), bias+2*mle_sd, code=3, angle=90, lty=8, col=c_vector)

axis(1, labels = c('FeR', 'FuR',expression('T'[f]), expression('RMC'[f])), at = c(2,6,10,14))
legend('bottomright', legend = c('Standard care', 'Dexanethasone'),col = c('blue','red'), lty = c('dashed','dashed'))
