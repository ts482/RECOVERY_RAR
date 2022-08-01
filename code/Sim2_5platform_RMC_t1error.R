
set.seed(20)

reps = 1000

K     = 100 #number of stages in trial
Si    = 4 #number of treatments
T     = 3 #number of subgroups

#3 treatment groups at ratio of 2:1:1:1
#therefore randomisation aim will be to have treatment groups at 1.5, control group at 1

nk_.me  = 120 #patients per stage so ar days =80, sample size = 9600 (as opposed to 9602)

rand_ratio.me  = c(0.4, 0.2, 0.2, 0.2)

Nx.me1   = c(501, 1279, 324) + c(1034, 2604, 683) #dex group, control group as reported for dex   #no ox, ox, vent
Nx.me2  = c(362, 938, 261) + c(425, 1131, 60)       #hydroxychloroquine, lopinavir
N.control_group = c(1034, 2604, 683)
N.treatment_groups = c(362, 938, 261) + c(425, 1131, 60) + c(501, 1279, 324)
Nx.me = Nx.me1 + Nx.me2
Px.me   = Nx.me/sum(Nx.me)
Py.me1   = matrix(c(0.140, 0.262, 0.414, 0.178, 0.233, 0.293), nrow=3)
Py.me2  = matrix(c(0.160, 0.270, 0.421, 0.167, 0.247, 0.4), nrow = 3)         #hydroxychloroquine, lopinavir
Py.me = matrix(rep(c(0.140, 0.262, 0.414), 4), nrow=3)
#Py    = c(0.25,0.22) #probability of death on each treatment
#S     = length(Py[1,]) #Number of treatments
#T     = length(Px) #Number of care


alpha.me       = 1 #uninformative beta
beta.me        = 1 #uninformative beta
Y.me           = matrix(nrow=K,ncol=Si)
postmean.me    = matrix(nrow=K,ncol=Si)
mle.me         = matrix(nrow=K,ncol=Si)
PostAlpha.me   = matrix(nrow=K,ncol=Si)
PostBeta.me    = matrix(nrow=K,ncol=Si)

nk.me          = matrix(round(nk_.me*rand_ratio.me), nrow=K,ncol=Si, byrow=T)
rand_ratios.me = matrix(nrow=K, ncol=Si)

UPPER_LIMIT <- 0.5
LOWER_LIMIT <- 0.05

nka.me = array(dim = c(K, Si, reps))
Ya.me = array(dim = c(K, Si, reps))
postmean_a.me = array(dim = c(K, Si, reps))
mle_a.me = array(dim=c(K, Si, reps))

y.me = vector(length=3)

conf_sup.me <- array(dim = c(Si-1,K/10, reps))
conf_sup_dir.me <- array(dim = c(Si-1,K/10, reps))
t_metric.me <- array(dim = c(Si-1,K/10, reps))


prob_optimal.me = vector(length=Si-1)
prob_optimals.me = array(dim=c(reps, Si-1, K))
pre_ratio.me = vector(length=Si-1)


X.me = matrix(nrow=100, ncol=Si)
# Stage 1
for (n in 1:reps){

rand_ratio.me = c(0.4, 0.2, 0.2, 0.2)

for(j in 1:Si){
#for (h in 1:T){
#  y.me[h] = rbinom(1, rbinom(1,nk.me[1,j], Px.me[h]), Py.me[h,j])
#}
Y.me[1,j]         = rbinom(1, round(nk_.me * rand_ratio.me[j]), 0.25)   #sum(y.me)
PostAlpha.me[1,j] = alpha.me + Y.me[1,j]
PostBeta.me[1,j]  = beta.me  + nk.me[1,j] -Y.me[1,j]
postmean.me[1,j]  = PostAlpha.me[1,j]/(PostAlpha.me[1,j]+PostBeta.me[1,j])
mle.me[1,j]       = sum(Y.me[1,j])/(nk.me[1,j])
}


for(i in 2:K){
  
for(j in 1:Si){
  nk.me[i, j] = rbinom(1,nk_.me, rand_ratio.me[j])
#for (h in 1:T){
#  y.me[h] = rbinom(1, rbinom(1, nk.me[i,j], Px.me[h]), Py.me[h,j])
#}
Y.me[i,j]         = rbinom(1, nk.me[i,j], 0.25)
PostAlpha.me[i,j] = alpha.me + sum(Y.me[1:i,j], na.rm=T)
PostBeta.me[i,j]  = beta.me + sum(nk.me[1:i,j]-Y.me[1:i,j], na.rm=T)
postmean.me[i,j]  = PostAlpha.me[i,j]/(PostAlpha.me[i,j]+PostBeta.me[i,j])
mle.me[i,j]       = sum(Y.me[1:i,j])/sum(nk.me[1:i,j])


}

if (i > 28){

  if (i %% 7 == 0){
    
    for (j in 1:Si){
      X.me[,j] = rbeta(100, PostAlpha.me[i-28, j], PostBeta.me[i-28, j])
    }
    optimum = apply(X.me, 1, which.min)
    k = i/(K)
    
    if (all(optimum == 1)){
      pre_ratio.me = c(1/3, 1/3, 1/3)
    } else {
    for (j in 2:Si){
      prob_optimal.me[j-1] = mean(optimum == j)
      pre_ratio.me[j-1] = sqrt(prob_optimal.me[j-1]/(sum(nk.me[1:(i-28),j])+1))
    }
    }
    
    for (j in 2:Si){
      rand_ratio.me[j] = pre_ratio.me[j-1]*0.6/sum(pre_ratio.me)
    }
    while (min(rand_ratio.me)<0.05){
      maxi = which.max(rand_ratio.me[2:Si]) + 1
      mini = which.min(rand_ratio.me[2:Si]) + 1
      diff = 0.05 - rand_ratio.me[mini]
      rand_ratio.me[maxi] = rand_ratio.me[maxi] - diff
      rand_ratio.me[mini] = rand_ratio.me[mini] + diff
    }
  }
  rand_ratios.me[i,] = rand_ratio.me

}
if (i %% 10 == 0){
  #calculating power for dex
  pat_df <- data.frame(t=rep(1, sum(nk.me[1:i, c(1,2)])), 
                       o=rep(0, sum(nk.me[1:i, c(1,2)])))
  pat_df[1:sum(nk.me[1:i,1]), 't'] = 0
  pat_df[1:(sum(Y.me[1:i,1])), 'o'] = 1
  pat_df[(sum(nk.me[1:i,1])+1):(sum(nk.me[1:i,1])+sum(Y.me[1:i,2])+1), 'o'] = 1
  
#  fit = glm(o~t,data=pat_df, family=binomial)
#  conf_sup.me[1,i/10,n] = summary(fit)$coef['t','Pr(>|z|)']
#  conf_sup_dir.me[1,i/10,n] = summary(fit)$coef['t','Estimate']
  t_metric.me[1,i/10,n] = t.test(o~t, data = pat_df)$p.value
  
  
  #calculating power for hydroxy
  pat_df <- data.frame(t=rep(1, sum(nk.me[1:i, c(1,3)])), 
                       o=rep(0, sum(nk.me[1:i, c(1,3)])))
  pat_df[1:sum(nk.me[1:i,1]), 't'] = 0
  pat_df[1:(sum(Y.me[1:i,1])), 'o'] = 1
  pat_df[(sum(nk.me[1:i,1])+1):(sum(nk.me[1:i,1])+sum(Y.me[1:i,3])+1), 'o'] = 1
  
#  fit = glm(o~t,data=pat_df, family=binomial)
#  conf_sup.me[2,i/10,n] = summary(fit)$coef['t','Pr(>|z|)']
#  conf_sup_dir.me[2,i/10,n] = summary(fit)$coef['t','Estimate']
  t_metric.me[2,i/10,n] = t.test(o~t, data = pat_df)$p.value
  
  
  #calculating power for lopinavir
  pat_df <- data.frame(t=rep(1, sum(nk.me[1:i, c(1,4)])), 
                       o=rep(0, sum(nk.me[1:i, c(1,4)])))
  pat_df[1:sum(nk.me[1:i,1]), 't'] = 0
  pat_df[1:(sum(Y.me[1:i,1])), 'o'] = 1
  pat_df[(sum(nk.me[1:i,1])+1):(sum(nk.me[1:i,1])+sum(Y.me[1:i,4])+1), 'o'] = 1
  
#  fit = glm(o~t,data=pat_df, family=binomial)
#  conf_sup.me[3,i/10,n] = summary(fit)$coef['t','Pr(>|z|)']
#  conf_sup_dir.me[3,i/10,n] = summary(fit)$coef['t','Estimate']
  t_metric.me[3,i/10,n] = t.test(o~t, data = pat_df)$p.value
}
}


nka.me[,,n] = nk.me
Ya.me[,,n] = Y.me
postmean_a.me[,,n] = postmean.me
mle_a.me[,,n] = mle.me

}



t1e_metric.me = vector(length=10)
for (i in 1:10){
  t1e_metric.me[i] = mean(apply(t_metric.me[,i,]<(0.05/3), 2, any))
}

lines(t1e_metric.me*100, lwd=4, type='b', col='red')

#legend('topright',legend=c(expression('T'[f]), expression('RMC'[f])), col=c('blue','red'), lwd=4)

legend('topright',legend=c(expression('T'[f]), expression('RMC'[f]), 'FeR','FuR'),
       col=c('blue','red','brown','green'), lwd=4)




Y_avg.me = matrix(nrow=K, ncol=Si)
nk_avg.me = matrix(nrow=K, ncol=Si)
Y_sd.me = vector(length=Si)
nk_sd.me = vector(length=Si)
for (j in 1:Si){
  nk_avg.me[,j] = apply(nka.me[,j,], 1, mean)
  Y_avg.me[,j] = apply(Ya.me[,j,], 1, mean, na.rm=T)
  nk_sd.me[j] = sd(nka.me[,j,])
  Y_sd.me[j] = sd(Ya.me[,j,])
}



## TYPE 1 error


type_1_err.me = vector(length=10)
type_1_err_adj.me = vector(length=10)
t1e.me = vector(length=10)

for (i in 1:10){
  type_1_err.me[i] = mean(apply(conf_sup.me[,i,]<0.05, 2, any)) #0.186
  type_1_err_adj.me[i] =mean(apply(conf_sup.me[,i,]<(0.05/3), 2, any)) #0.077
  t1e.me[i] = mean(apply(t_metric.me[,i,]<0.05/3, 2, any)) #0.186
}

type_1_err_adj.me
t1e.me





#lines(nk_avg.me[c(1,1:10*10),2]/nk_.me*100, lwd=4, type='b', col='brown')


for (j in 1:Si){
  print(sum(nk_avg.me[,j]))
}


for (j in 1:Si){
print(sum(Y_avg.me[,j]))
}

for (j in 1:Si){
  print(sum(Y_avg.me[,j])/sum(nk_avg.me[,j]))
}


#deaths at n=6450
sum(nk_avg.me[1:80,])
sum(Y_avg.me[1:80,])


sum(Y_avg.me[1:80,])/sum(nk_avg.me[1:80,])


dex_conf = apply(conf_sup.me[1,,], 1, mean)
hydroxy_conf = apply(conf_sup.me[2,,], 1, mean)
lopinavir_conf = apply(conf_sup.me[3,,], 1, mean)

dex_conf
hydroxy_conf
lopinavir_conf

power_avg = apply(conf_sup.me<0.05, 1, mean)
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
  mle_avg[i] = mean(((mle_a.me[i*10,1,] - mle_a.me[i*10,2,])-treatment_effect)/treatment_effect)
  mle_sd[i] = sd(((mle_a.me[i*10,1,] - mle_a.me[i*10,2,])-treatment_effect)/treatment_effect)
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
mle_t_avg = mean(((mle_a.me[80,1,] - mle_a.me[80,2,])-treatment_effect)/treatment_effect)
mle_t_sd = sd(((mle_a.me[80,1,] - mle_a.me[80,2,])-treatment_effect)/treatment_effect)


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
