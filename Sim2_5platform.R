
set.seed(20)

reps = 1000

K     = 100 #number of stages in trial
Si    = 4 #number of treatments
T     = 3 #number of subgroups

#3 treatment groups at ratio of 2:1:1:1
#therefore randomisation aim will be to have treatment groups at 1.5, control group at 1

nk_.i  = 120 #patients per stage so ar days =80, sample size = 9600 (as opposed to 9602)


rand_ratio.i  = c(0.4, 0.2, 0.2, 0.2)


Nx.i1   = c(501, 1279, 324) + c(1034, 2604, 683) #dex group, control group as reported for dex   #no ox, ox, vent
Nx.i2  = c(362, 938, 261) + c(425, 1131, 60)       #hydroxychloroquine, lopinavir
N.control_group = c(1034, 2604, 683)
N.treatment_groups = c(362, 938, 261) + c(425, 1131, 60) + c(501, 1279, 324)
Nx.i = Nx.i1 + Nx.i2
Px.i   = Nx.i/sum(Nx.i)
Py.i1   = matrix(c(0.140, 0.262, 0.414, 0.178, 0.233, 0.293), nrow=3)
Py.i2  = matrix(c(0.160, 0.270, 0.421, 0.167, 0.247, 0.4), nrow = 3)         #hydroxychloroquine, lopinavir
Py.i = matrix(c(Py.i1, Py.i2), nrow=3)
#Py    = c(0.25,0.22) #probability of death on each treatment
#S     = length(Py[1,]) #Number of treatments
#T     = length(Px) #Number of care


alpha.i       = 1 #uninformative beta
beta.i        = 1 #uninformative beta
Y.i           = matrix(nrow=K,ncol=Si)
postmean.i    = matrix(nrow=K,ncol=Si)
mle.i         = matrix(nrow=K,ncol=Si)
PostAlpha.i   = matrix(nrow=K,ncol=Si)
PostBeta.i    = matrix(nrow=K,ncol=Si)

nk.i          = matrix(round(nk_.i*rand_ratio.i), nrow=K,ncol=Si, byrow=T)
rand_ratios.i = matrix(nrow=K, ncol=Si)

UPPER_LIMIT <- 0.5
LOWER_LIMIT <- 0.05

nka.i = array(dim = c(K, Si, reps))
Ya.i = array(dim = c(K, Si, reps))
postmean_a.i = array(dim = c(K, Si, reps))
mle_a.i = array(dim=c(K, Si, reps))

y.i = vector(length=3)

conf_sup.i <- array(dim = c(Si-1,K/10, reps))

prob_optimal.i = vector(length=Si-1)
prob_optimals.i = array(dim=c(reps, Si-1, K))
pre_ratio.i = vector(length=Si-1)


X.i = matrix(nrow=100, ncol=Si)
# Stage 1
for (n in 1:reps){

rand_ratio.i = c(0.4, 0.2, 0.2, 0.2)

for(j in 1:Si){
for (h in 1:T){
  y.i[h] = rbinom(1, round(nk.i[1,j]* Px.i[h]), Py.i[h,j])
}
Y.i[1,j]         = sum(y.i)
PostAlpha.i[1,j] = alpha.i + Y.i[1,j]
PostBeta.i[1,j]  = beta.i  + nk.i[1,j] -Y.i[1,j]
postmean.i[1,j]  = PostAlpha.i[1,j]/(PostAlpha.i[1,j]+PostBeta.i[1,j])
mle.i[1,j]       = sum(Y.i[1,j])/(nk.i[1,j])
}


for(i in 2:K){
  
for(j in 1:Si){
  nk.i[i, j] = round(nk_.i * rand_ratio.i[j])
for (h in 1:T){
  y.i[h] = rbinom(1, round(nk.i[i,j]* Px.i[h]), Py.i[h,j])
}
Y.i[i,j]         = sum(y.i)
PostAlpha.i[i,j] = alpha.i + sum(Y.i[1:i,j], na.rm=T)
PostBeta.i[i,j]  = beta.i + sum(nk.i[1:i,j]-Y.i[1:i,j], na.rm=T)
postmean.i[i,j]  = PostAlpha.i[i,j]/(PostAlpha.i[i,j]+PostBeta.i[i,j])
mle.i[i,j]       = sum(Y.i[1:i,j])/sum(nk.i[1:i,j])


}

if (i > 28){

if (i %% 7 == 0){
  
  for (j in 1:Si){
    X.i[,j] = rbeta(100, PostAlpha.i[i-28, j], PostBeta.i[i-28, j])
  }
optimum = apply(X.i, 1, which.min)
k = i/(K)
for (j in 2:Si){
  prob_optimal.i[j-1] = mean(optimum == j)  # multiplied by 0.6 as 0.4 is for control group
  pre_ratio.i[j-1] = prob_optimal.i[j-1]^k / (prob_optimal.i[j-1]^k + (1-prob_optimal.i[j-1])^k)
}

for (j in 2:Si){
  rand_ratio.i[j] = pre_ratio.i[j-1]*0.6/sum(pre_ratio.i)
}
while (min(rand_ratio.i)<0.05){
  maxi = which.max(rand_ratio.i[2:Si]) + 1
  mini = which.min(rand_ratio.i[2:Si]) + 1
  diff = 0.05 - rand_ratio.i[mini]
  rand_ratio.i[maxi] = rand_ratio.i[maxi] - diff
  rand_ratio.i[mini] = rand_ratio.i[mini] + diff
}

}
rand_ratios.i[i,] = rand_ratio.i

}
if (i %% 10 == 0){
  #calculating power for dex
  pat_df <- data.frame(t=rep(1, sum(nk.i[1:i, c(1,2)])), 
                       o=rep(0, sum(nk.i[1:i, c(1,2)])))
  pat_df[1:sum(nk.i[1:i,1]), 't'] = 0
  pat_df[1:(sum(Y.i[1:i,1])), 'o'] = 1
  pat_df[(sum(nk.i[1:i,1])+1):(sum(nk.i[1:i,1])+sum(Y.i[1:i,2])+1), 'o'] = 1
  
  fit = glm(o~t,data=pat_df, family=binomial)
  conf_sup.i[1,i/10,n] = summary(fit)$coef['t','Pr(>|z|)']
  
  
  #calculating power for hydroxy
  pat_df <- data.frame(t=rep(1, sum(nk.i[1:i, c(1,3)])), 
                       o=rep(0, sum(nk.i[1:i, c(1,3)])))
  pat_df[1:sum(nk.i[1:i,1]), 't'] = 0
  pat_df[1:(sum(Y.i[1:i,1])), 'o'] = 1
  pat_df[(sum(nk.i[1:i,1])+1):(sum(nk.i[1:i,1])+sum(Y.i[1:i,3])+1), 'o'] = 1
  
  fit = glm(o~t,data=pat_df, family=binomial)
  conf_sup.i[2,i/10,n] = summary(fit)$coef['t','Pr(>|z|)']

  
  #calculating power for lopinavir
  pat_df <- data.frame(t=rep(1, sum(nk.i[1:i, c(1,4)])), 
                       o=rep(0, sum(nk.i[1:i, c(1,4)])))
  pat_df[1:sum(nk.i[1:i,1]), 't'] = 0
  pat_df[1:(sum(Y.i[1:i,1])), 'o'] = 1
  pat_df[(sum(nk.i[1:i,1])+1):(sum(nk.i[1:i,1])+sum(Y.i[1:i,4])+1), 'o'] = 1
  
  fit = glm(o~t,data=pat_df, family=binomial)
  conf_sup.i[3,i/10,n] = summary(fit)$coef['t','Pr(>|z|)']
}
}


nka.i[,,n] = nk.i
Ya.i[,,n] = Y.i
postmean_a.i[,,n] = postmean.i
mle_a.i[,,n] = mle.i

}



Y_avg.i = matrix(nrow=K, ncol=Si)
nk_avg.i = matrix(nrow=K, ncol=Si)
Y_sd.i = matrix(nrow=K, ncol=Si)
nk_sd.i = matrix(nrow=K, ncol=Si)
for (j in 1:Si){
  nk_avg.i[,j] = apply(nka.i[,j,], 1, mean)
  Y_avg.i[,j] = apply(Ya.i[,j,], 1, mean, na.rm=T)
  nk_sd.i[,j] = apply(nka.i[,j,], 1, sd)
  Y_sd.i[,j] = apply(Ya.i[,j,], 1, sd, na.rm=T)
}


#number of people allocated to dexamethasone
sum(nk_avg.i[,2]) # 4101.2


apply(nk_avg.i[1:80,],2, sum)
apply(nk_avg.m[1:80,],2, sum)



plot(nk_avg.i[c(1,seq(10, 100, 10)),1]/nk_.i*100, xlim = c(1,11), ylim = c(0,100), lwd=3, col='red', type = 'b',
     xlab='Trial progress(%)', ylab='Proportion of patients receiving each treatment(%)', axes=FALSE)
lines(nk_avg.i[c(1,seq(10, 100, 10)),2]/nk_.i*100, col='blue', type='b', lwd=3)
arrows(5:11 - 0.05, (nk_avg.i[seq(40,100,10),2] - nk_sd.i[seq(40,100,10),2])/nk_.i*100,
       y1 = (nk_avg.i[seq(40,100,10),2] + nk_sd.i[seq(40,100,10),2])/nk_.i*100,
       code=3, angle=90, lty=8, col='blue', lwd=2)
lines(nk_avg.i[c(1,seq(10, 100, 10)),3]/nk_.i*100, col='yellow', type='b', lwd=3)
arrows(5:11 - 0.05, (nk_avg.i[seq(40,100,10),3] - nk_sd.i[seq(40,100,10),3])/nk_.i*100,
       y1 = (nk_avg.i[seq(40,100,10),3] + nk_sd.i[seq(40,100,10),3])/nk_.i*100,
       code=3, angle=90, lty=8, col='yellow', lwd=2)
lines(nk_avg.i[c(1,seq(10, 100, 10)),4]/nk_.i*100, col='green', type='b', lwd=3)
arrows(5:11 - 0.05, (nk_avg.i[seq(40,100,10),4] - nk_sd.i[seq(40,100,10),4])/nk_.i*100,
       y1 = (nk_avg.i[seq(40,100,10),4] + nk_sd.i[seq(40,100,10),4])/nk_.i*100,
       code=3, angle=90, lty=8, col='green', lwd=2)

legend("topleft",legend = c('standard care', 'dexamethasone', 'hydroxychloroquine', 'lopinavir'),
       col = c("red", "blue", 'yellow','green'), lwd=4)
axis(1, seq(1,11,2), labels = seq(0, 125,25))
axis(2, seq(0, 100, 25), labels = seq(0, 100, 25))

abline(v=4.4, lty='dashed')


for (j in 1:Si){
  print(sum(nk_avg.i[1:80,j]))
}


for (j in 1:Si){
print(sum(Y_avg.i[,j]))
}

for (j in 1:Si){
  print(sum(Y_avg.i[,j])/sum(nk_avg.i[,j]))
}



#deaths at n=9600
sum(nk_avg.i[1:80,])
sum(Y_avg.i[1:80,])
sd(apply(Ya.i[1:80, 1, ] + Ya.i[1:80, 2, ] + 
           Ya.i[1:80, 3, ] + Ya.i[1:80, 4, ], 2, sum))/(reps ** 0.5)


sum(Y_avg.i[1:80,])/sum(nk_avg.i[1:80,])




#power calculation and plotting

treatment_power_avg = matrix(nrow=11, ncol=3)
treatment_power_err = matrix(nrow=10, ncol=3)
for (j in 1:3){
treatment_power_avg[,j] = c(0, apply(conf_sup.i[j,,]<0.05, 1, mean))
treatment_power_err[,j] = apply(conf_sup.i[j,,]<0.05, 1, sd)/(reps**0.5)
}


plot(treatment_power_avg[,1], ylim = c(0,1.1), col = 'blue', axes=F,
     type = 'b', lwd=2, xlab='trial progress (%)',ylab='Power of study')
axis(1, at = seq(1,11,2), labels = seq(0, 125,25))
axis(2, at = seq(0,1,0.2))

#points(treatment_power_avg[,2], col='blue', type='b', lwd=2, lty='dashed')
points(treatment_power_avg[,3], col='blue', type='b', lwd=2, lty='dashed')
arrows(2:11, treatment_power_avg[2:11,1] - treatment_power_err[,1],
       y1 = treatment_power_avg[2:11,1] + treatment_power_err[,1],
       code=3, angle=90, lty=8, col='blue', lwd=2.5)
arrows(2:11, treatment_power_avg[2:11,3] - treatment_power_err[,3],
       y1 = treatment_power_avg[2:11,3] + treatment_power_err[,3],
       code=3, angle=90, lty=9, col='blue', lwd=2.5)


#points(power_avg, lwd=4)
#lines(power_avg, col = 'brown',lwd=4)
#power_avg



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
  mle_avg[i] = mean(((mle_a.i[i*10,1,] - mle_a.i[i*10,2,])-treatment_effect)/treatment_effect)
  mle_sd[i] = sd(((mle_a.i[i*10,1,] - mle_a.i[i*10,2,])-treatment_effect)/treatment_effect)
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
mle_t_avg = mean(((mle_a.i[80,1,] - mle_a.i[80,2,])-treatment_effect)/treatment_effect)
mle_t_sd = sd(((mle_a.i[80,1,] - mle_a.i[80,2,])-treatment_effect)/treatment_effect)


upper = mle_t_avg + 1.96*mle_t_sd/sqrt(6450)
mid = mle_t_avg
lower = mle_t_avg - 1.96*mle_t_sd/sqrt(6450)

points(13.5, mid, col = 'blue', lwd=6)
arrows(13.5, lower, 13.5, upper, code=3, angle=90, lty=8, col='blue', lwd=2)


#mean squared error

arm_estimates.i = apply(Px.i * Py.i, 2, sum)

treatment_effects = arm_estimates.i[1] - arm_estimates.i[2:4]

treatment_effects

mse.i = vector(length=3)

for (j in 1:3){
  mse.i[j] = mean((mle_a.i[80,1,] - mle_a.i[80,j+1,]- treatment_effects[j])**2)
}
mse.i




#treatment effect bias at n=9600

mle_t_avg = vector(length=3)
mle_t_sd = vector(length=3)


arm_estimates.i = apply(Px.i * Py.i, 2, sum)

treatment_effects = arm_estimates.i[1] - arm_estimates.i[2:4]

treatment_effects

for (j in 1:3){
  mle_t_avg[j] = mean(((mle_a.i[80,1,] - mle_a.i[80,j+1,])- treatment_effects[j])/treatment_effects[j])
  mle_t_sd[j] = sd(((mle_a.i[80,1,] - mle_a.i[80,j+1,])- treatment_effects[j])/treatment_effects[j])
}



upper = mle_t_avg + 1.96*mle_t_sd/sqrt(1000)
mid = mle_t_avg
lower = mle_t_avg - 1.96*mle_t_sd/sqrt(1000)

j=1
plot(0.875+(j-1), mid[j],
     col='blue', xaxt='n',xlim = c(0.5,3.5), ylim=c(-1.5, 1.5),
     ylab='Relative Bias', xlab='Treatment arm', lwd=6)
for (j in 1:3){
  points(0.875+(j-1), mid[j], col='blue', lwd=6)
  arrows(0.875+c(j-1), lower[j], 0.875+c(j-1),
         upper[j], code=3, angle=90, lty=8, col='blue', lwd=2)
}

