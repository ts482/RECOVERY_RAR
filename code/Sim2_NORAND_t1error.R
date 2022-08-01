

set.seed(20)

reps = 1000
K       = 100 #number of stages in trial
S       = 2
T       = 3

nk_.ae   = 80 #patients per stage
Nx.ae    = c(501, 1279, 324) + c(1034, 2604, 683)
Px.ae    = Nx.ae/sum(Nx.ae)
Py.ae   = matrix(c(0.140, 0.262, 0.414, 0.178, 0.233, 0.293), nrow=3)
#Py    = c(0.25,0.22) #probability of death on each treatment
#S       = length(Py.ae[1,]) #Number of treatments
#T       = length(Px.ae) #Number of care

rand_ratio.ae <- rep(0.5, S)
pre_ratio.ae <- vector(length=2)

alpha.ae       = 1 #uninformative beta
beta.ae        = 1 #uninformative beta
Y.ae           = matrix(nrow=K,ncol=S)
postmean.ae    = matrix(nrow=K,ncol=S)
mle.ae         = matrix(nrow=K,ncol=S)
PostAlpha.ae   = matrix(nrow=K,ncol=S)
PostBeta.ae    = matrix(nrow=K,ncol=S)

nk.ae          = matrix(round(nk_.ae*rand_ratio.ae), nrow=K,ncol=S)
rand_ratios.ae = matrix(nrow=K, ncol=S)

UPPER_LIMIT <- 0.9
LOWER_LIMIT <- 0.1

nka.ae = array(dim = c(K, S, reps))
Ya.ae = array(dim = c(K, S, reps))
postmean_a.ae = array(dim = c(K, S, reps))
mle_a.ae = array(dim=c(K, S, reps))

X.ae = matrix(nrow=100, ncol=S)

y.ae = vector(length=3)

prob_optimal.ae = vector(length=S)
prob_optimals.ae = array(dim=c(reps, S, K))

conf_sup.ae <- matrix(nrow=K/10, ncol=reps)
# Stage 1
for (n in 1:reps){
for(j in 1:S){
for (h in 1:T){
  y.ae[h] = rbinom(1, round(nk.ae[1,j]* Px.ae[h]), 0.25)
}
Y.ae[1,j]         = sum(y.ae)
PostAlpha.ae[1,j] = alpha.ae + Y.ae[1,j]
PostBeta.ae[1,j]  = beta.ae  + nk.ae[1,j] -Y.ae[1,j]
postmean.ae[1,j]  = PostAlpha.ae[1,j]/(PostAlpha.ae[1,j]+PostBeta.ae[1,j])
mle.ae[1,j]       = sum(Y.ae[1,j])/(nk.ae[1,j])
}

for(i in 2:K){

  n_stand = rbinom(1, nk_.ae, prob = rand_ratio.ae)
  nk.ae[i, ] = c(n_stand, nk_.ae-n_stand)
for(j in 1:S){
  
for (h in 1:T){
  y.ae[h] = rbinom(1, round(nk.ae[i,j]* Px.ae[h]), 0.25)
}
Y.ae[i,j]         = sum(y.ae)
PostAlpha.ae[i,j] = alpha.ae + sum(Y.ae[1:i,j], na.rm=T)
PostBeta.ae[i,j]  = beta.ae + sum(nk.ae[1:i,j]-Y.ae[1:i,j], na.rm=T)
postmean.ae[i,j]  = PostAlpha.ae[i,j]/(PostAlpha.ae[i,j]+PostBeta.ae[i,j])
mle.ae[i,j]       = sum(Y.ae[1:i,j])/sum(nk.ae[1:i,j])

#X[,j] = rbeta(100, PostAlpha[i, j], PostBeta[i, j])
}
#optimum = apply(X, 1, which.min)
#for (j in 1:S){
#  prob_optimal[j] = mean(optimum == j)
#}


#prob_optimals[n,,i] = prob_optimal

if (i %% 10 == 0){
pat_df <- data.frame(t=rep(1, i*nk_.ae), o=rep(0, i*nk_.ae))
pat_df[1:sum(nk.ae[1:i,1]), 't'] = 0
pat_df[1:(sum(Y.ae[1:i,1])), 'o'] = 1
pat_df[(sum(nk.ae[1:i,1])+1):(sum(nk.ae[1:i,1])+sum(Y.ae[1:i,2])+1), 'o'] = 1

fit = glm(o~t,data=pat_df, family=binomial)
conf_sup.ae[i/10,n] = summary(fit)$coef['t','Pr(>|z|)']
}
}

  
nka.ae[,,n] = nk.ae
Ya.ae[,,n] = Y.ae
postmean_a.ae[,,n] = postmean.ae
mle_a.ae[,,n] = mle.ae
}

sum(Ya.ae)/reps

## TYPE 1 error



t1e_metric.ae = vector(length=10)

for (i in 1:10){
  t1e_metric.ae[i] = mean(conf_sup.ae[i,]<=0.05) #0.05 exactly!
}

t1e_metric.ae[8]

lines(t1e_metric.ae*100, lwd=4, type='b', col='brown')






Y_avg.ae = matrix(nrow=K, ncol=S)
nk_avg.ae = matrix(nrow=K, ncol=S)
Y_sd.ae = matrix(nrow=K, ncol=S)
nk_sd.ae = matrix(nrow=K, ncol=S)
for (j in 1:S){
  nk_avg.ae[,j] = apply(nka.ae[,j,], 1, mean)
  Y_avg.ae[,j] = apply(Ya.ae[,j,], 1, mean, na.rm=T)
  nk_sd.ae[,j] = apply(nka.ae[,j,], 1, sd)
  Y_sd.ae[,j] = apply(nka.ae[,j,], 1, sd)
}


sum(nk_avg.ae[,1])
sum(nk_avg.ae[,2])

sum(Y_avg.ae[,1])
sum(Y_avg.ae[,2])

#deaths at n=6450
sum(Y_avg.ae[1:80,])
sd(apply(Ya.ae[1:80, 1, ] + Ya.ae[1:80, 2, ], 2, sum))/(reps ** 0.5)




power_avg = apply(conf_sup.ae<0.05, 1, mean)
power_avg = c(0, power_avg)
plot(power_avg, ylim=c(0,1), lwd=4, xaxt = 'n', xlim = c(1,11),
     xlab = 'Trial Progress(%)', ylab= 'Power of study')
axis(1, seq(1,11,2), labels = seq(0, 125,25))
lines(power_avg, col='blue',lwd=4)
power_err = apply(conf_sup.ae<0.05, 1, sd)/(reps ** 0.5)
arrows(2:11, power_avg[2:11] - power_err,
       y1 = power_avg[2:11] + power_err,
       code=3, angle=90, lty=8, col='blue', lwd=2.5)



plot(nk_avg.ae[c(2,seq(10, 100, 10)),2]/nk_.ae*100, xlim = c(1,11), ylim = c(0,100), lwd=4, col='red', type = 'b',
     xlab='Trial progress(%)', ylab='Proportion of patients receiving dexamethasone(%)', axes=FALSE)
axis(1, seq(1,11,2), labels = seq(0, 125,25))
axis(2, seq(0, 100, 25), labels = seq(0, 100, 25))
arrows(1:11, (nk_avg.ae[c(2,seq(10, 100, 10)),2]- nk_sd.ae[c(2,seq(10, 100, 10)),2])/nk_.ae*100,
       y1 = (nk_avg.ae[c(2,seq(10, 100, 10)),2] + nk_sd.ae[c(2,seq(10, 100, 10)),2])/nk_.ae*100,
       code=3, angle=90, lty=8, col='red', lwd=2)
#points(nk_avg[,2], col='blue')
#legend("bottomleft",legend = c('no treatment', 'dexamethasone'),col = c("red", "blue"), lwd=4)

sum(Ya)/reps

treatment_superior = vector(length=K)
for (i in 1:K){
  treatment_superior[i] = mean(postmean_a[i,1,] >= postmean_a[i,2,])
}


plot(treatment_superior)

points(1:10*10,conf_sup_avg)
plot(1-thetas_avg)
par(mfrow=c(1,1))
thetas_avg = apply(prob_optimals[,2,], 2, mean)
points(1-thetas_avg)
above_theta = prob_optimals[,2,] > 0.95
above_theta_avg = apply(above_theta, 2, mean)
plot(above_theta_avg, col='blue')

legend('bottomright',legend = c('Tuning method', 'REMAP-CAP', 'fixed randomisation'), 
       fill=c('red','blue','black'))

tail(above_theta_avg)
#plot(postmean[,1],ylim=range(postmean),type="l",lwd=4,
#     ylab="Posterior mean",xlab="Interim Analysis")
#lines(postmean[,2],type="l",lwd=4,col="red")
#legend("bottomright",c("Treatment 1","Treatment 2"),col=c("black","red"),lwd=4,bty="n")

#lines(mle[,1],ylim=range(mle),type="l",lwd=4,lty=2)
#lines(mle[,2],type="l",lwd=4,col="red",lty=2)


mle_avg = vector(length= 10) #vector(length=S)
mle_sd = vector(length=10)#vector(length=S)
for (i in 1:10){
  mle_avg[i] = mean(mle_a.ae[i*10,1,] - mle_a.ae[i*10,2,])
  mle_sd[i] = sd(mle_a.ae[i*10,1,] - mle_a.ae[i*10,2,])
}


bias = mle_avg - (0.257 - 0.229)


plot(1:10-0.2, bias, ylim=c(-0.05,0.05), xaxt = 'n', xlab = 'Trial Progress(%)', 
     ylab= 'Power of study', col= 'blue')
axis(1, at = 1:10, labels = seq(10,100, 10))
arrows(1:10-0.2, bias-2*mle_sd, 1:10-0.2, bias+2*mle_sd, code=3, angle=90, lty=8, col='blue')
abline(h=0, lty='dashed')

c_vector =c('blue', 'red')

# end of trial treatment effect bias



mle_t_avg = mean(mle_a.ae[43,1,] - mle_a.ae[43,2,])
mle_t_sd = sd(mle_a.ae[43,1,] - mle_a.ae[43,2,])

bias_t = (mle_t_avg - (0.257 - 0.229))

upper = (bias_t + 2*mle_t_sd)/(0.257 - 0.229)
mid = bias_t/(0.257 - 0.229)
lower = (bias_t - 2*mle_t_sd)/(0.257 - 0.229)

plot(2, mid ,xlim = c(0,16), ylim = c(lower-1, upper+1), xaxt='n', 
     ylab='Bias (% of treatment effect)', xlab='Algorithm')
arrows(2, lower, 2, upper, code=3, angle=90, lty=8, )



#mean squared error
mort_rates = apply(Px.ae * Py.ae,2,sum)
true_t_effect = mort_rates[1] - mort_rates[2]
mse.ae = mean((mle_a.ae[80,1,] - mle_a.ae[80,2,] - true_t_effect)**2)
mse.ae






# end of trial bias

mle_avg = vector(length=S)
mle_sd = vector(length=S)
for (j in 1:S){
  mle_avg[j] = mean(mle_a.ae[100,j,])
  mle_sd[j] = sd(mle_a.ae[100,j,])
}

bias = mle_avg - c(0.257, 0.229)
c_vector =c('blue', 'red')

plot(c(1.5, 2.5), bias,xlim = c(0,16), ylim=c(-0.04, 0.04), col=c_vector, xaxt='n', ylab='Bias', xlab='Algorithm')
arrows(c(1.5,2.5), bias-2*mle_sd, c(1.5,2.5), bias+2*mle_sd, code=3, angle=90, lty=8, col=c_vector)
