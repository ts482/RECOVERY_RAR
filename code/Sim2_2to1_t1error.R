
set.seed(20)

reps = 1000

K     = 100 #number of stages in trial
S     = 2
T     = 3

nk_.be   = 80 #patients per stage
Nx.be    = c(501, 1279, 324) + c(1034, 2604, 683)
Px.be    = Nx.be/sum(Nx.be)
Py.be   = matrix(c(0.140, 0.262, 0.414, 0.178, 0.233, 0.293), nrow=3)
#Py    = c(0.257,0.229) #probability of death on each treatment
#S     = length(Py[1,]) #Number of treatments
#T     = length(Px) #Number of care


rand_ratio.be <- c(2/3, 1/3)
pre_ratio.be <- vector(length=2)

alpha.be       = 1 #uninformative beta
beta.be        = 1 #uninformative beta
Y.be           = matrix(nrow=K,ncol=S)
postmean.be    = matrix(nrow=K,ncol=S)
mle.be         = matrix(nrow=K,ncol=S)
PostAlpha.be   = matrix(nrow=K,ncol=S)
PostBeta.be    = matrix(nrow=K,ncol=S)

nk.be          = matrix(nrow=K,ncol=S)
for (j in 1:S){
  nk.be[,j] = rand_ratio.be[j]*nk_.be
}
rand_ratios.be = matrix(nrow=K, ncol=S)

UPPER_LIMIT <- 0.9
LOWER_LIMIT <- 0.1

nka.be = array(dim = c(K, S, reps))
Ya.be = array(dim = c(K, S, reps))
postmean_a.be = array(dim = c(K, S, reps))
mle_a.be = array(dim=c(K, S, reps))

X.be = matrix(nrow=100, ncol=S)

y.be = vector(length=3)

prob_optimal.be = vector(length=S)
prob_optimals.be = array(dim=c(reps, S, K))


conf_sup.be <- matrix(nrow=K/10, ncol=reps)
# Stage 1
for (n in 1:reps){
for(j in 1:S){
for (h in 1:T){
  y.be[h] = rbinom(1, round(nk.be[1,j]* Px.be[h]), 0.25)
}
Y.be[1,j]         = sum(y.be)
PostAlpha.be[1,j] = alpha.be + Y.be[1,j]
PostBeta.be[1,j]  = beta.be  + nk.be[1,j] -Y.be[1,j]
postmean.be[1,j]  = PostAlpha.be[1,j]/(PostAlpha.be[1,j]+PostBeta.be[1,j])
mle.be[1,j]       = sum(Y.be[1,j])/(nk.be[1,j])
}

for(i in 2:K){

n_stand = rbinom(1, nk_.be, prob = c(rand_ratio.be[1], rand_ratio.be[2]))
nk.be[i, ] = c(n_stand, nk_.be-n_stand)
  
for(j in 1:S){
  
for (h in 1:T){
  y.be[h] = rbinom(1, round(nk.be[i,j]* Px.be[h]), 0.25)
}
Y.be[i,j]         = sum(y.be)
PostAlpha.be[i,j] = alpha.be + sum(Y.be[1:i,j], na.rm=T)
PostBeta.be[i,j]  = beta.be + sum(nk.be[1:i,j]-Y.be[1:i,j], na.rm=T)
postmean.be[i,j]  = PostAlpha.be[i,j]/(PostAlpha.be[i,j]+PostBeta.be[i,j])
mle.be[i,j]       = sum(Y.be[1:i,j])/sum(nk.be[1:i,j])

#X[,j] = rbeta(100, PostAlpha[i, j], PostBeta[i, j])
}
#optimum = apply(X, 1, which.min)
#for (j in 1:S){
#  prob_optimal[j] = mean(optimum == j)
#}
#prob_optimals[n,,i] = prob_optimal

if (i %% 10 == 0){
  pat_df <- data.frame(t=rep(1, i*nk_.be), o=rep(0, i*nk_.be))
  pat_df[1:sum(nk.be[1:i,1]), 't'] = 0
  pat_df[1:(sum(Y.be[1:i,1])), 'o'] = 1
  pat_df[(sum(nk.be[1:i,1])+1):(sum(nk.be[1:i,1])+sum(Y.be[1:i,2])+1), 'o'] = 1
  
  fit = glm(o~t,data=pat_df, family=binomial)
  conf_sup.be[i/10,n] = summary(fit)$coef['t','Pr(>|z|)']
}
}

  
nka.be[,,n] = nk.be
Ya.be[,,n] = Y.be
postmean_a.be[,,n] = postmean.be
mle_a.be[,,n] = mle.be
}



sum(Ya.be)/reps




## TYPE 1 ERROR

t1e_metric.be = vector(length=10)

for (i in 1:10){
  t1e_metric.be[i] = mean(conf_sup.be[i,]<0.05) #0.05 exactly!
}

t1e_metric.be[8]

lines(t1e_metric.be*100, lwd=4, type='b', col='green')














Y_avg.be = matrix(nrow=K, ncol=S)
nk_avg.be = matrix(nrow=K, ncol=S)
Y_sd.be = matrix(nrow=K, ncol=S)
nk_sd.be = matrix(nrow=K, ncol=S)
for (j in 1:S){
  nk_avg.be[,j] = apply(nka.be[,j,], 1, mean)
  Y_avg.be[,j] = apply(Ya.be[,j,], 1, mean, na.rm=T)
  nk_sd.be[,j] = apply(nka.be[,j,], 1, sd)
  Y_sd.be[,j] = apply(Ya.be[,j,], 1, sd, na.rm=T)
}
  
  
sum(nk_avg.be[,1])
sum(nk_avg.be[,2])

sum(Y_avg.be[,1])
sum(Y_avg.be[,2])


#deaths at n=6450
sum(Y_avg.be[1:80,])
sd(apply(Ya.be[1:80, 1, ] + Ya.be[1:80, 2, ], 2, sum))/(reps ** 0.5)



lines(nk_avg.be[c(2,1:10*10),2]/nk_.be*100, lwd=4, type='b', col='blue')
arrows(1:11, (nk_avg.be[c(2,seq(10, 100, 10)),2]- nk_sd.be[c(2,seq(10, 100, 10)),2])/nk_.be*100,
       y1 = (nk_avg.be[c(2,seq(10, 100, 10)),2] + nk_sd.be[c(2,seq(10, 100, 10)),2])/nk_.be*100,
       code=3, angle=90, lty=8, col='blue', lwd=2)
#points(nk_avg[,2]/nk_*100, col='blue')
#legend("bottomleft",legend = c('no treatment', 'dexamethasone'),col = c("red", "blue"), lwd=4)

sum(Ya)/reps




conf_sup_avg = apply(conf_sup,1, mean)
plot(conf_sup_avg)
lines(conf_sup_avg)


power_avg = apply(conf_sup.be<0.05, 1, mean)
power_avg = c(0, power_avg)
points(power_avg,lwd=4)
lines(power_avg, col='red',lwd=4)
power_err = apply(conf_sup.be<0.05, 1, sd)/(reps ** 0.5)
arrows(2:11, power_avg[2:11] - power_err,
       y1 = power_avg[2:11] + power_err,
       code=3, angle=90, lty=8, col='red', lwd=2.5)


treatment_superior = vector(length=K)
for (i in 1:K){
  treatment_superior[i] = mean(postmean_a[i,1,] >= postmean_a[i,2,])
}


plot(treatment_superior)

thetas_avg = apply(prob_optimals[,2,], 2, mean)
plot(thetas_avg)
above_theta = prob_optimals[,2,] > 0.95
above_theta_avg = apply(above_theta, 2, mean)
points(above_theta_avg, col='red')

legend('bottomright', legend = c('1:1 randomisation', '2:1 randomisation'), fill=c('blue', 'red'))



#legend('bottomright',legend = c('Tuning method', 'REMAP-CAP', 'fixed randomisation'), 
#       fill=c('red','blue','black'))

tail(above_theta_avg)
#plot(postmean[,1],ylim=range(postmean),type="l",lwd=4,
#     ylab="Posterior mean",xlab="Interim Analysis")
#lines(postmean[,2],type="l",lwd=4,col="red")
#legend("bottomright",c("Treatment 1","Treatment 2"),col=c("black","red"),lwd=4,bty="n")

#lines(mle[,1],ylim=range(mle),type="l",lwd=4,lty=2)
#lines(mle[,2],type="l",lwd=4,col="red",lty=2)


# end of trial treatment effect bias



mle_t_avg = mean(mle_a.be[43,1,] - mle_a.be[43,2,])
mle_t_sd = sd(mle_a.be[43,1,] - mle_a.be[43,2,])

bias_t = (mle_t_avg - (0.257 - 0.229))

upper = (bias_t + 2*mle_t_sd)/(0.257 - 0.229)
mid = bias_t/(0.257 - 0.229)
lower = (bias_t - 2*mle_t_sd)/(0.257 - 0.229)

points(6, mid ,xlim = c(0,16), ylim = c(lower-1, upper+1), xaxt='n', ylab='Bias', xlab='Algorithm')
arrows(6, lower, 6, upper, code=3, angle=90, lty=8, )





#mean squared error
mort_rates = apply(Px.be * Py.be,2,sum)
true_t_effect = mort_rates[1] - mort_rates[2]
mse.be = mean((mle_a.be[80,1,] - mle_a.be[80,2,] - true_t_effect)**2)
mse.be








#end of trial bias
mle_avg = vector(length=S)
mle_sd = vector(length=S)
for (j in 1:S){
  mle_avg[j] = mean(mle_a.be[100,j,])
  mle_sd[j] = sd(mle_a.be[100,j,])
}

bias = mle_avg - c(0.257, 0.229)

c_vector =c('blue', 'red')
points(c(5.5, 6.5), bias, col=c_vector)
arrows(c(5.5,6.5), bias-2*mle_sd, c(5.5,6.5), bias+2*mle_sd, code=3, angle=90, lty=8, col=c_vector)
