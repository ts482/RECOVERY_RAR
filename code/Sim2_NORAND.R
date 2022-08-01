

set.seed(20)

reps = 1000
K       = 100 #number of stages in trial
S       = 2
T       = 3

nk_.a   = 80 #patients per stage
Nx.a    = c(501, 1279, 324) + c(1034, 2604, 683)
Px.a    = Nx.a/sum(Nx.a)
Py.a   = matrix(c(0.140, 0.262, 0.414, 0.178, 0.233, 0.293), nrow=3)
#Py    = c(0.25,0.22) #probability of death on each treatment
#S       = length(Py.a[1,]) #Number of treatments
#T       = length(Px.a) #Number of care

rand_ratio.a <- rep(0.5, S)
pre_ratio.a <- vector(length=2)

alpha.a       = 1 #uninformative beta
beta.a        = 1 #uninformative beta
Y.a           = matrix(nrow=K,ncol=S)
postmean.a    = matrix(nrow=K,ncol=S)
mle.a         = matrix(nrow=K,ncol=S)
PostAlpha.a   = matrix(nrow=K,ncol=S)
PostBeta.a    = matrix(nrow=K,ncol=S)

nk.a          = matrix(round(nk_.a*rand_ratio.a), nrow=K,ncol=S)
rand_ratios.a = matrix(nrow=K, ncol=S)

UPPER_LIMIT <- 0.9
LOWER_LIMIT <- 0.1

nka.a = array(dim = c(K, S, reps))
Ya.a = array(dim = c(K, S, reps))
postmean_a.a = array(dim = c(K, S, reps))
mle_a.a = array(dim=c(K, S, reps))

X.a = matrix(nrow=100, ncol=S)

y.a = vector(length=3)

prob_optimal.a = vector(length=S)
prob_optimals.a = array(dim=c(reps, S, K))

conf_sup.a <- matrix(nrow=K/10, ncol=reps)
# Stage 1
for (n in 1:reps){
for(j in 1:S){
for (h in 1:T){
  y.a[h] = rbinom(1, round(nk.a[1,j]* Px.a[h]), Py.a[h,j])
}
Y.a[1,j]         = sum(y.a)
PostAlpha.a[1,j] = alpha.a + Y.a[1,j]
PostBeta.a[1,j]  = beta.a  + nk.a[1,j] -Y.a[1,j]
postmean.a[1,j]  = PostAlpha.a[1,j]/(PostAlpha.a[1,j]+PostBeta.a[1,j])
mle.a[1,j]       = sum(Y.a[1,j])/(nk.a[1,j])
}

for(i in 2:K){

  n_stand = rbinom(1, nk_.a, prob = rand_ratio.a)
  nk.a[i, ] = c(n_stand, nk_.a-n_stand)
for(j in 1:S){
  
for (h in 1:T){
  y.a[h] = rbinom(1, round(nk.a[i,j]* Px.a[h]), Py.a[h,j])
}
Y.a[i,j]         = sum(y.a)
PostAlpha.a[i,j] = alpha.a + sum(Y.a[1:i,j], na.rm=T)
PostBeta.a[i,j]  = beta.a + sum(nk.a[1:i,j]-Y.a[1:i,j], na.rm=T)
postmean.a[i,j]  = PostAlpha.a[i,j]/(PostAlpha.a[i,j]+PostBeta.a[i,j])
mle.a[i,j]       = sum(Y.a[1:i,j])/sum(nk.a[1:i,j])

#X[,j] = rbeta(100, PostAlpha[i, j], PostBeta[i, j])
}
#optimum = apply(X, 1, which.min)
#for (j in 1:S){
#  prob_optimal[j] = mean(optimum == j)
#}


#prob_optimals[n,,i] = prob_optimal

if (i %% 10 == 0){
pat_df <- data.frame(t=rep(1, i*nk_.a), o=rep(0, i*nk_.a))
pat_df[1:sum(nk.a[1:i,1]), 't'] = 0
pat_df[1:(sum(Y.a[1:i,1])), 'o'] = 1
pat_df[(sum(nk.a[1:i,1])+1):(sum(nk.a[1:i,1])+sum(Y.a[1:i,2])+1), 'o'] = 1

fit = glm(o~t,data=pat_df, family=binomial)
conf_sup.a[i/10,n] = summary(fit)$coef['t','Pr(>|z|)']
}
}

  
nka.a[,,n] = nk.a
Ya.a[,,n] = Y.a
postmean_a.a[,,n] = postmean.a
mle_a.a[,,n] = mle.a
}

sum(Ya.a)/reps

Y_avg.a = matrix(nrow=K, ncol=S)
nk_avg.a = matrix(nrow=K, ncol=S)
Y_sd.a = matrix(nrow=K, ncol=S)
nk_sd.a = matrix(nrow=K, ncol=S)
for (j in 1:S){
  nk_avg.a[,j] = apply(nka.a[,j,], 1, mean)
  Y_avg.a[,j] = apply(Ya.a[,j,], 1, mean, na.rm=T)
  nk_sd.a[,j] = apply(nka.a[,j,], 1, sd)
  Y_sd.a[,j] = apply(nka.a[,j,], 1, sd)
}


sum(nk_avg.a[,1])
sum(nk_avg.a[,2])

sum(Y_avg.a[,1])
sum(Y_avg.a[,2])

#deaths at n=6450
sum(Y_avg.a[1:80,])
sd(apply(Ya.a[1:80, 1, ] + Ya.a[1:80, 2, ], 2, sum))/(reps ** 0.5)




power_avg = apply(conf_sup.a<0.05, 1, mean)
power_avg = c(0, power_avg)
plot(power_avg, ylim=c(0,1), lwd=4, xaxt = 'n', xlim = c(1,11),
     xlab = 'Trial Progress(%)', ylab= 'Power of study')
axis(1, seq(1,11,2), labels = seq(0, 125,25))
lines(power_avg, col='blue',lwd=4)
power_err = apply(conf_sup.a<0.05, 1, sd)/(reps ** 0.5)
arrows(2:11, power_avg[2:11] - power_err,
       y1 = power_avg[2:11] + power_err,
       code=3, angle=90, lty=8, col='blue', lwd=2.5)



plot(nk_avg.a[c(2,seq(10, 100, 10)),2]/nk_.a*100, xlim = c(1,11), ylim = c(0,100), lwd=4, col='red', type = 'b',
     xlab='Trial progress(%)', ylab='Proportion of patients receiving dexamethasone(%)', axes=FALSE)
axis(1, seq(1,11,2), labels = seq(0, 125,25))
axis(2, seq(0, 100, 25), labels = seq(0, 100, 25))
arrows(1:11, (nk_avg.a[c(2,seq(10, 100, 10)),2]- nk_sd.a[c(2,seq(10, 100, 10)),2])/nk_.a*100,
       y1 = (nk_avg.a[c(2,seq(10, 100, 10)),2] + nk_sd.a[c(2,seq(10, 100, 10)),2])/nk_.a*100,
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
  mle_avg[i] = mean(mle_a.a[i*10,1,] - mle_a.a[i*10,2,])
  mle_sd[i] = sd(mle_a.a[i*10,1,] - mle_a.a[i*10,2,])
}


bias = mle_avg - (0.257 - 0.229)


plot(1:10-0.2, bias, ylim=c(-0.05,0.05), xaxt = 'n', xlab = 'Trial Progress(%)', 
     ylab= 'Power of study', col= 'blue')
axis(1, at = 1:10, labels = seq(10,100, 10))
arrows(1:10-0.2, bias-2*mle_sd, 1:10-0.2, bias+2*mle_sd, code=3, angle=90, lty=8, col='blue')
abline(h=0, lty='dashed')

c_vector =c('blue', 'red')

# end of trial treatment effect bias



mle_t_avg = mean(mle_a.a[43,1,] - mle_a.a[43,2,])
mle_t_sd = sd(mle_a.a[43,1,] - mle_a.a[43,2,])

bias_t = (mle_t_avg - (0.257 - 0.229))

upper = (bias_t + 2*mle_t_sd)/(0.257 - 0.229)
mid = bias_t/(0.257 - 0.229)
lower = (bias_t - 2*mle_t_sd)/(0.257 - 0.229)

plot(2, mid ,xlim = c(0,16), ylim = c(lower-1, upper+1), xaxt='n', 
     ylab='Bias (% of treatment effect)', xlab='Algorithm')
arrows(2, lower, 2, upper, code=3, angle=90, lty=8, )



#mean squared error
mort_rates = apply(Px.a * Py.a,2,sum)
true_t_effect = mort_rates[1] - mort_rates[2]
mse.a = mean((mle_a.a[80,1,] - mle_a.a[80,2,] - true_t_effect)**2)
mse.a






# end of trial bias

mle_avg = vector(length=S)
mle_sd = vector(length=S)
for (j in 1:S){
  mle_avg[j] = mean(mle_a.a[100,j,])
  mle_sd[j] = sd(mle_a.a[100,j,])
}

bias = mle_avg - c(0.257, 0.229)
c_vector =c('blue', 'red')

plot(c(1.5, 2.5), bias,xlim = c(0,16), ylim=c(-0.04, 0.04), col=c_vector, xaxt='n', ylab='Bias', xlab='Algorithm')
arrows(c(1.5,2.5), bias-2*mle_sd, c(1.5,2.5), bias+2*mle_sd, code=3, angle=90, lty=8, col=c_vector)
