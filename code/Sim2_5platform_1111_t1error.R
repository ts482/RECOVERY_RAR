
set.seed(20)

reps = 1000

K     = 100 #number of stages in trial
Si    = 4 #number of treatments
T     = 3 #number of subgroups

#3 treatment groups at ratio of 2:1:1:1
#therefore randomisation aim will be to have treatment groups at 1.5, control group at 1

nk_.le  = 120 #patients per stage so ar days =80, sample size = 9600 (as opposed to 9602)
rand_ratio.le  = c(0.25, 0.25, 0.25, 0.25)
pre_ratio.le <- vector(length=Si)


Nx.le1   = c(501, 1279, 324) + c(1034, 2604, 683) #dex group, control group as reported for dex   #no ox, ox, vent
Nx.le2  = c(362, 938, 261) + c(425, 1131, 60)       #hydroxychloroquine, lopinavir
N.control_group = c(1034, 2604, 683)
N.treatment_groups = c(362, 938, 261) + c(425, 1131, 60) + c(501, 1279, 324)
Nx.le = Nx.le1 + Nx.le2
Px.le   = Nx.le/sum(Nx.le)
Py.le1   = matrix(c(0.140, 0.262, 0.414, 0.178, 0.233, 0.293), nrow=3)
Py.le2  = matrix(c(0.160, 0.270, 0.421, 0.167, 0.247, 0.4), nrow = 3)         #hydroxychloroquine, lopinavir
Py.le = matrix(c(Py.le1, Py.le2), nrow=3)
#Py    = c(0.25,0.22) #probability of death on each treatment
#S     = length(Py[1,]) #Number of treatments
#T     = length(Px) #Number of care


alpha.le       = 1 #uninformative beta
beta.le        = 1 #uninformative beta
Y.le           = matrix(nrow=K,ncol=Si)
postmean.le    = matrix(nrow=K,ncol=Si)
mle.le         = matrix(nrow=K,ncol=Si)
PostAlpha.le   = matrix(nrow=K,ncol=Si)
PostBeta.le    = matrix(nrow=K,ncol=Si)

nk.le          = matrix(round(nk_.le*rand_ratio.le), nrow=K,ncol=Si, byrow=T)
rand_ratios.le = matrix(nrow=K, ncol=Si)

UPPER_LIMIT <- 0.5
LOWER_LIMIT <- 0.05

nka.le = array(dim = c(K, Si, reps))
Ya.le = array(dim = c(K, Si, reps))
postmean_a.le = array(dim = c(K, Si, reps))
mle_a.le = array(dim=c(K, Si, reps))

y.le = vector(length=3)

conf_sup.le <- array(dim = c(Si-1,K/10, reps))

prob_optimal.le = vector(length=Si-1)
prob_optimals.le = array(dim=c(reps, Si-1, K))
pre_ratio.le = vector(length=Si-1)


X.le = matrix(nrow=100, ncol=Si)
# Stage 1
for (n in 1:reps){

rand_ratio.le = c(0.25, 0.25, 0.25, 0.25)

for(j in 1:Si){
  nk.le[1, j] = rbinom(1,nk_.le, rand_ratio.le[j])
#for (h in 1:T){
#  y.le[h] = rbinom(1, rbinom(1,nk.le[1,j], Px.le[h]), Py.le[h,j])
#}
Y.le[1,j]         = rbinom(1,nk.le[1,j],0.25)
PostAlpha.le[1,j] = alpha.le + Y.le[1,j]
PostBeta.le[1,j]  = beta.le  + nk.le[1,j] -Y.le[1,j]
postmean.le[1,j]  = PostAlpha.le[1,j]/(PostAlpha.le[1,j]+PostBeta.le[1,j])
mle.le[1,j]       = sum(Y.le[1,j])/(nk.le[1,j])
}


for(i in 2:K){
  
for(j in 1:Si){
  nk.le[i, j] = rbinom(1,nk_.le, rand_ratio.le[j])
#for (h in 1:T){
#  y.le[h] = rbinom(1, rbinom(1,nk.le[i,j], Px.le[h]), Py.le[h,j])
#}
Y.le[i,j]         = rbinom(1,nk.le[i,j],0.25)
PostAlpha.le[i,j] = alpha.le + sum(Y.le[1:i,j], na.rm=T)
PostBeta.le[i,j]  = beta.le + sum(nk.le[1:i,j]-Y.le[1:i,j], na.rm=T)
postmean.le[i,j]  = PostAlpha.le[i,j]/(PostAlpha.le[i,j]+PostBeta.le[i,j])
mle.le[i,j]       = sum(Y.le[1:i,j])/sum(nk.le[1:i,j])


}


if (i %% 10 == 0){
  #calculating power for dex
  pat_df <- data.frame(t=rep(1, sum(nk.le[1:i, c(1,2)])), 
                       o=rep(0, sum(nk.le[1:i, c(1,2)])))
  pat_df[1:sum(nk.le[1:i,1]), 't'] = 0
  pat_df[1:(sum(Y.le[1:i,1])), 'o'] = 1
  pat_df[(sum(nk.le[1:i,1])+1):(sum(nk.le[1:i,1])+sum(Y.le[1:i,2])+1), 'o'] = 1
  
  fit = glm(o~t,data=pat_df, family=binomial)
  conf_sup.le[1,i/10,n] = summary(fit)$coef['t','Pr(>|z|)']
  
  
  #calculating power for hydroxy
  pat_df <- data.frame(t=rep(1, sum(nk.le[1:i, c(1,3)])), 
                       o=rep(0, sum(nk.le[1:i, c(1,3)])))
  pat_df[1:sum(nk.le[1:i,1]), 't'] = 0
  pat_df[1:(sum(Y.le[1:i,1])), 'o'] = 1
  pat_df[(sum(nk.le[1:i,1])+1):(sum(nk.le[1:i,1])+sum(Y.le[1:i,3])+1), 'o'] = 1
  
  fit = glm(o~t,data=pat_df, family=binomial)
  conf_sup.le[2,i/10,n] = summary(fit)$coef['t','Pr(>|z|)']
  
  
  #calculating power for lopinavir
  pat_df <- data.frame(t=rep(1, sum(nk.le[1:i, c(1,4)])), 
                       o=rep(0, sum(nk.le[1:i, c(1,4)])))
  pat_df[1:sum(nk.le[1:i,1]), 't'] = 0
  pat_df[1:(sum(Y.le[1:i,1])), 'o'] = 1
  pat_df[(sum(nk.le[1:i,1])+1):(sum(nk.le[1:i,1])+sum(Y.le[1:i,4])+1), 'o'] = 1
  
  fit = glm(o~t,data=pat_df, family=binomial)
  conf_sup.le[3,i/10,n] = summary(fit)$coef['t','Pr(>|z|)']
}
}


nka.le[,,n] = nk.le
Ya.le[,,n] = Y.le
postmean_a.le[,,n] = postmean.le
mle_a.le[,,n] = mle.le

}



#type 1 error

t1e_metric.le = vector(length=10)
for (i in 1:10){
  t1e_metric.le[i] = mean(apply(conf_sup.le[,i,]<(0.05/3), 2, any))
}

t1e_metric.le[8]


lines(t1e_metric.le*100, lwd=4, type='b', col='brown')





Y_avg.le = matrix(nrow=K, ncol=Si)
nk_avg.le = matrix(nrow=K, ncol=Si)
Y_sd.le = vector(length=Si)
nk_sd.le = vector(length=Si)
for (j in 1:Si){
  nk_avg.le[,j] = apply(nka.le[,j,], 1, mean)
  Y_avg.le[,j] = apply(Ya.le[,j,], 1, mean, na.rm=T)
  nk_sd.le[j] = sd(nka.le[,j,])
  Y_sd.le[j] = sd(Ya.le[,j,])
}


#number of patients randomised to dexamethasone
sum(nk_avg.le[,2])

#lines(nk_avg.le[c(1,1:10*10),2]/nk_.le*100, lwd=4, type='b', col='brown')


for (j in 1:Si){
  print(sum(nk_avg.le[,j]))
}


for (j in 1:Si){
print(sum(Y_avg.le[,j]))
}

for (j in 1:Si){
  print(sum(Y_avg.le[,j])/sum(nk_avg.le[,j]))
}


#deaths at n=6450
sum(nk_avg.le[1:80,])
sum(Y_avg.le[1:80,])
sd(apply(Ya.le[1:80, 1, ] + Ya.le[1:80, 2, ] + 
           Ya.le[1:80, 3, ] + Ya.le[1:80, 4, ], 2, sum))/(reps ** 0.5)





sum(Y_avg.le[1:80,])/sum(nk_avg.le[1:80,])


#power calculation and plotting

treatment_power_avg = matrix(nrow=11, ncol=3)
treatment_power_err = matrix(nrow=10, ncol=3)
for (j in 1:3){
  treatment_power_avg[,j] = c(0, apply(conf_sup.le[j,,]<0.05, 1, mean))
  treatment_power_err[,j] = apply(conf_sup.le[j,,]<0.05, 1, sd)/(reps**0.5)
}




points(treatment_power_avg[,1], col='brown', type='b', lwd=2)
#points(treatment_power_avg[,2], col='yellow', type='b', lwd=2, lty='dashed')
points(treatment_power_avg[,3], col='brown', type='b', lwd=2, lty='dashed')

arrows(2:11, treatment_power_avg[2:11,1] - treatment_power_err[,1],
       y1 = treatment_power_avg[2:11,1] + treatment_power_err[,1],
       code=3, angle=90, lty=8, col='brown', lwd=2.5)
arrows(2:11, treatment_power_avg[2:11,3] - treatment_power_err[,3],
       y1 = treatment_power_avg[2:11,3] + treatment_power_err[,3],
       code=3, angle=90, lty=9, col='brown', lwd=2.5)


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
  mle_avg[i] = mean(((mle_a.le[i*10,1,] - mle_a.le[i*10,2,])-treatment_effect)/treatment_effect)
  mle_sd[i] = sd(((mle_a.le[i*10,1,] - mle_a.le[i*10,2,])-treatment_effect)/treatment_effect)
}



upper = mle_avg + 1.96*mle_sd/sqrt(6450)
mid = mle_avg
lower = mle_avg - 1.96*mle_sd/sqrt(6450)

points(1:10-0.2, mid, ylim=c(-0.06,0.06), col= 'blue', lwd=6)
arrows(1:10-0.2, lower, 1:10-0.2, upper, code=3, angle=90, lty=8, col='blue', lwd=3)

legend('bottomright', legend = c(expression('T'[f]), expression('RMC'[f])),col = c('blue','red'), 
       lty = c('dashed', 'dashed'), lwd=3)





#mean squared error


arm_estimates.le = apply(Px.le * Py.le, 2, sum)

treatment_effects = arm_estimates.le[1] - arm_estimates.le[2:4]

treatment_effects

mse.le = vector(length=3)

for (j in 1:3){
  mse.le[j] = mean((mle_a.le[80,1,] - mle_a.le[80,j+1,]- treatment_effects[j])**2)
}

mse.le*10^4





# end of trial treatment effect bias


treatment_effect = 0.257 - 0.229
mle_t_avg = mean(((mle_a.le[80,1,] - mle_a.le[80,2,])-treatment_effect)/treatment_effect)
mle_t_sd = sd(((mle_a.le[80,1,] - mle_a.le[80,2,])-treatment_effect)/treatment_effect)


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
