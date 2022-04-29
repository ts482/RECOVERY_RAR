
set.seed(20)

reps = 1000

K     = 100 #number of stages in trial
Si    = 4 #number of treatments
T     = 3 #number of subgroups

#3 treatment groups at ratio of 2:1:1:1
#therefore randomisation aim will be to have treatment groups at 1.5, control group at 1

nk_.je  = 120 #patients per stage so ar days =80, sample size = 9600 (as opposed to 9602)
rand_ratio.je  = c(0.4, 0.2, 0.2, 0.2)
pre_ratio.je <- vector(length=Si)


Nx.je1   = c(501, 1279, 324) + c(1034, 2604, 683) #dex group, control group as reported for dex   #no ox, ox, vent
Nx.je2  = c(362, 938, 261) + c(425, 1131, 60)       #hydroxychloroquine, lopinavir
N.control_group = c(1034, 2604, 683)
N.treatment_groups = c(362, 938, 261) + c(425, 1131, 60) + c(501, 1279, 324)
Nx.je = Nx.je1 + Nx.je2
Px.je   = Nx.je/sum(Nx.je)
Py.je = matrix(rep(c(0.140, 0.262, 0.414), 4), nrow=3)
#Py    = c(0.25,0.22) #probability of death on each treatment
#S     = length(Py[1,]) #Number of treatments
#T     = length(Px) #Number of care


alpha.je       = 1 #uninformative beta
beta.je        = 1 #uninformative beta
Y.je           = matrix(nrow=K,ncol=Si)
postmean.je    = matrix(nrow=K,ncol=Si)
mle.je         = matrix(nrow=K,ncol=Si)
PostAlpha.je   = matrix(nrow=K,ncol=Si)
PostBeta.je    = matrix(nrow=K,ncol=Si)

nk.je          = matrix(round(nk_.je*rand_ratio.je), nrow=K,ncol=Si, byrow=T)
rand_ratios.je = matrix(nrow=K, ncol=Si)

UPPER_LIMIT <- 0.5
LOWER_LIMIT <- 0.05

nka.je = array(dim = c(K, Si, reps))
Ya.je = array(dim = c(K, Si, reps))
postmean_a.je = array(dim = c(K, Si, reps))
mle_a.je = array(dim=c(K, Si, reps))

y.je = vector(length=3)

conf_sup.je <- array(dim = c(Si-1,K/10, reps))
conf_sup_dir.je <- array(dim = c(Si-1,K/10, reps))

prob_optimal.je = vector(length=Si-1)
prob_optimals.je = array(dim=c(reps, Si-1, K))
pre_ratio.je = vector(length=Si-1)


X.je = matrix(nrow=100, ncol=Si)
# Stage 1
for (n in 1:reps){

rand_ratio.je = c(0.4, 0.2, 0.2, 0.2)

for(j in 1:Si){
for (h in 1:T){
  y.je[h] = rbinom(1, round(nk.je[1,j]* Px.je[h]), Py.je[h,j])
}
Y.je[1,j]         = sum(y.je)
PostAlpha.je[1,j] = alpha.je + Y.je[1,j]
PostBeta.je[1,j]  = beta.je  + nk.je[1,j] -Y.je[1,j]
postmean.je[1,j]  = PostAlpha.je[1,j]/(PostAlpha.je[1,j]+PostBeta.je[1,j])
mle.je[1,j]       = sum(Y.je[1,j])/(nk.je[1,j])
}


for(i in 2:K){
  
for(j in 1:Si){
  nk.je[i, j] = round(nk_.je * rand_ratio.je[j])
for (h in 1:T){
  y.je[h] = rbinom(1, round(nk.je[i,j]* Px.je[h]), Py.je[h,j])
}
Y.je[i,j]         = sum(y.je)
PostAlpha.je[i,j] = alpha.je + sum(Y.je[1:i,j], na.rm=T)
PostBeta.je[i,j]  = beta.je + sum(nk.je[1:i,j]-Y.je[1:i,j], na.rm=T)
postmean.je[i,j]  = PostAlpha.je[i,j]/(PostAlpha.je[i,j]+PostBeta.je[i,j])
mle.je[i,j]       = sum(Y.je[1:i,j])/sum(nk.je[1:i,j])


}

if (i > 28){

  if (i %% 7 == 0){
    
    for (j in 1:Si){
      X.je[,j] = rbeta(100, PostAlpha.je[i-28, j], PostBeta.je[i-28, j])
    }
    optimum = apply(X.je, 1, which.min)
    k = i/(K)
    for (j in 2:Si){
      prob_optimal.je[j-1] = mean(optimum == j)  # multiplied by 0.6 as 0.4 is for control group
      pre_ratio.je[j-1] = prob_optimal.je[j-1]^k / (prob_optimal.je[j-1]^k + (1-prob_optimal.je[j-1])^k)
    }
    if (sum(pre_ratio.je) == 0){
      pre_ratio.je[1:3] = 1/3
    }
    for (j in 2:Si){
      rand_ratio.je[j] = pre_ratio.je[j-1]*0.6/sum(pre_ratio.je)
    }
    while (min(rand_ratio.je)<0.05){
      maxi = which.max(rand_ratio.je[2:Si]) + 1
      mini = which.min(rand_ratio.je[2:Si]) + 1
      diff = 0.05 - rand_ratio.je[mini]
      rand_ratio.je[maxi] = rand_ratio.je[maxi] - diff
      rand_ratio.je[mini] = rand_ratio.je[mini] + diff
    }
    
  }
  rand_ratios.je[i,] = rand_ratio.je

}
if (i %% 10 == 0){
  #calculating power for dex
  pat_df <- data.frame(t=rep(1, sum(nk.je[1:i, c(1,2)])), 
                       o=rep(0, sum(nk.je[1:i, c(1,2)])))
  pat_df[1:sum(nk.je[1:i,1]), 't'] = 0
  pat_df[1:(sum(Y.je[1:i,1])), 'o'] = 1
  pat_df[(sum(nk.je[1:i,1])+1):(sum(nk.je[1:i,1])+sum(Y.je[1:i,2])+1), 'o'] = 1
  
  fit = glm(o~t,data=pat_df, family=binomial)
  conf_sup.je[1,i/10,n] = summary(fit)$coef['t','Pr(>|z|)']
  conf_sup_dir.je[1,i/10,n] = summary(fit)$coef['t','Estimate']
  
  
  #calculating power for hydroxy
  pat_df <- data.frame(t=rep(1, sum(nk.je[1:i, c(1,3)])), 
                       o=rep(0, sum(nk.je[1:i, c(1,3)])))
  pat_df[1:sum(nk.je[1:i,1]), 't'] = 0
  pat_df[1:(sum(Y.je[1:i,1])), 'o'] = 1
  pat_df[(sum(nk.je[1:i,1])+1):(sum(nk.je[1:i,1])+sum(Y.je[1:i,3])+1), 'o'] = 1
  
  fit = glm(o~t,data=pat_df, family=binomial)
  conf_sup.je[2,i/10,n] = summary(fit)$coef['t','Pr(>|z|)']
  conf_sup_dir.je[2,i/10,n] = summary(fit)$coef['t','Estimate']
  
  
  #calculating power for lopinavir
  pat_df <- data.frame(t=rep(1, sum(nk.je[1:i, c(1,4)])), 
                       o=rep(0, sum(nk.je[1:i, c(1,4)])))
  pat_df[1:sum(nk.je[1:i,1]), 't'] = 0
  pat_df[1:(sum(Y.je[1:i,1])), 'o'] = 1
  pat_df[(sum(nk.je[1:i,1])+1):(sum(nk.je[1:i,1])+sum(Y.je[1:i,4])+1), 'o'] = 1
  
  fit = glm(o~t,data=pat_df, family=binomial)
  conf_sup.je[3,i/10,n] = summary(fit)$coef['t','Pr(>|z|)']
  conf_sup_dir.je[3,i/10,n] = summary(fit)$coef['t','Estimate']
  
}
}


nka.je[,,n] = nk.je
Ya.je[,,n] = Y.je
postmean_a.je[,,n] = postmean.je
mle_a.je[,,n] = mle.je

}


mean(apply(conf_sup.je[,8,]<0.05, 2, any)) #0.186
mean(apply(conf_sup.je[,8,]<(0.05/3), 2, any)) #0.077


Y_avg.je = matrix(nrow=K, ncol=Si)
nk_avg.je = matrix(nrow=K, ncol=Si)
Y_sd.je = vector(length=Si)
nk_sd.je = vector(length=Si)
for (j in 1:Si){
  nk_avg.je[,j] = apply(nka.je[,j,], 1, mean)
  Y_avg.je[,j] = apply(Ya.je[,j,], 1, mean, na.rm=T)
  nk_sd.je[j] = sd(nka.je[,j,])
  Y_sd.je[j] = sd(Ya.je[,j,])
}


lines(nk_avg.je[c(1,1:10*10),2]/nk_.je*100, lwd=4, type='b', col='brown')



for (j in 1:Si){
  print(sum(nk_avg.je[,j]))
}


for (j in 1:Si){
print(sum(Y_avg.je[,j]))
}

for (j in 1:Si){
  print(sum(Y_avg.je[,j])/sum(nk_avg.je[,j]))
}




#deaths at n=6450
sum(nk_avg.je[1:80,])
sum(Y_avg.je[1:80,])


dex_conf = apply(conf_sup.je[1,,], 1, mean)
hydroxy_conf = apply(conf_sup.je[2,,], 1, mean)
lopinavir_conf = apply(conf_sup.je[3,,], 1, mean)

dex_conf
hydroxy_conf
lopinavir_conf

power_avg = apply(conf_sup.je<0.05, 1, mean)
power_avg = c(0, power_avg)
points(power_avg, lwd=4)
lines(power_avg, col = 'brown',lwd=4)
power_avg

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
  mle_avg[i] = mean(((mle_a.je[i*10,1,] - mle_a.je[i*10,2,])-treatment_effect)/treatment_effect)
  mle_sd[i] = sd(((mle_a.je[i*10,1,] - mle_a.je[i*10,2,])-treatment_effect)/treatment_effect)
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
mle_t_avg = mean(((mle_a.je[80,1,] - mle_a.je[80,2,])-treatment_effect)/treatment_effect)
mle_t_sd = sd(((mle_a.je[80,1,] - mle_a.je[80,2,])-treatment_effect)/treatment_effect)


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
