
set.seed(20)

reps = 1000

K     = 100 #number of stages in trial
S     = 2
T     = 3

nk_.b   = 150 #patients per stage
Nx.b    = c(501, 1279, 324) + c(1034, 2604, 683)
Px.b    = Nx.b/sum(Nx.b)
Py.b   = matrix(c(0.140, 0.262, 0.414, 0.178, 0.233, 0.293), nrow=3)
#Py    = c(0.257,0.229) #probability of death on each treatment
#S     = length(Py[1,]) #Number of treatments
#T     = length(Px) #Number of care


rand_ratio.b <- c(2/3, 1/3)
pre_ratio.b <- vector(length=2)

alpha.b       = 1 #uninformative beta
beta.b        = 1 #uninformative beta
Y.b           = matrix(nrow=K,ncol=S)
postmean.b    = matrix(nrow=K,ncol=S)
mle.b         = matrix(nrow=K,ncol=S)
PostAlpha.b   = matrix(nrow=K,ncol=S)
PostBeta.b    = matrix(nrow=K,ncol=S)

nk.b          = matrix(nrow=K,ncol=S)
for (j in 1:S){
  nk.b[,j] = rand_ratio.b[j]*nk_.b
}
rand_ratios.b = matrix(nrow=K, ncol=S)

UPPER_LIMIT <- 0.9
LOWER_LIMIT <- 0.1

nka.b = array(dim = c(K, S, reps))
Ya.b = array(dim = c(K, S, reps))
postmean_a.b = array(dim = c(K, S, reps))
mle_a.b = array(dim=c(K, S, reps))

X.b = matrix(nrow=100, ncol=S)

y.b = vector(length=3)

prob_optimal.b = vector(length=S)
prob_optimals.b = array(dim=c(reps, S, K))


conf_sup.b <- matrix(nrow=K/10, ncol=reps)
# Stage 1
for (n in 1:reps){
for(j in 1:S){
for (h in 1:T){
  y.b[h] = rbinom(1, round(nk.b[1,j]* Px.b[h]), Py.b[h,j])
}
Y.b[1,j]         = sum(y.b)
PostAlpha.b[1,j] = alpha.b + Y.b[1,j]
PostBeta.b[1,j]  = beta.b  + nk.b[1,j] -Y.b[1,j]
postmean.b[1,j]  = PostAlpha.b[1,j]/(PostAlpha.b[1,j]+PostBeta.b[1,j])
mle.b[1,j]       = sum(Y.b[1,j])/(nk.b[1,j])
}

for(i in 2:K){

n_stand = rbinom(1, nk_.b, prob = c(rand_ratio.b[1], rand_ratio.b[2]))
nk.b[i, ] = c(n_stand, nk_-n_stand)
  
for(j in 1:S){
  
for (h in 1:T){
  y.b[h] = rbinom(1, round(nk.b[i,j]* Px.b[h]), Py.b[h,j])
}
Y.b[i,j]         = sum(y.b)
PostAlpha.b[i,j] = alpha.b + sum(Y.b[1:i,j], na.rm=T)
PostBeta.b[i,j]  = beta.b + sum(nk.b[1:i,j]-Y.b[1:i,j], na.rm=T)
postmean.b[i,j]  = PostAlpha.b[i,j]/(PostAlpha.b[i,j]+PostBeta.b[i,j])
mle.b[i,j]       = sum(Y.b[1:i,j])/sum(nk.b[1:i,j])

#X[,j] = rbeta(100, PostAlpha[i, j], PostBeta[i, j])
}
#optimum = apply(X, 1, which.min)
#for (j in 1:S){
#  prob_optimal[j] = mean(optimum == j)
#}
#prob_optimals[n,,i] = prob_optimal

if (i %% 10 == 0){
  pat_df <- data.frame(t=rep(1, i*nk_.b), o=rep(0, i*nk_.b))
  pat_df[1:sum(nk.b[1:i,1]), 't'] = 0
  pat_df[1:(sum(Y.b[1:i,1])), 'o'] = 1
  pat_df[(sum(nk.b[1:i,1])+1):(sum(nk.b[1:i,1])+sum(Y.b[1:i,2])+1), 'o'] = 1
  
  fit = glm(o~t,data=pat_df, family=binomial)
  conf_sup.b[i/10,n] = summary(fit)$coef['t','Pr(>|z|)']
}
}

  
nka.b[,,n] = nk.b
Ya.b[,,n] = Y.b
postmean_a.b[,,n] = postmean.b
mle_a.b[,,n] = mle.b
}

sum(Ya.b)/reps

Y_avg.b = matrix(nrow=K, ncol=S)
nk_avg.b = matrix(nrow=K, ncol=S)
Y_sd.b = vector(length=S)
nk_sd.b = vector(length=S)
for (j in 1:S){
  nk_avg.b[,j] = apply(nka.b[,j,], 1, mean)
  Y_avg.b[,j] = apply(Ya.b[,j,], 1, mean, na.rm=T)
  nk_sd.b[j] = sd(nka.b[,j,])
  Y_sd.b[j] = sd(Ya.b[,j,])
}
nk_sd.b

sum(nk_avg.b[,1])
sum(nk_avg.b[,2])

sum(Y_avg.b[,1])
sum(Y_avg.b[,2])


#deaths at n=6450
sum(Y_avg.b[1:43,])


lines(nk_avg.b[seq(10, 100, 10),2]/nk_.b*100, lwd=4, type='b', col='blue')
#points(nk_avg[,2]/nk_*100, col='blue')
#legend("bottomleft",legend = c('no treatment', 'dexamethasone'),col = c("red", "blue"), lwd=4)

sum(Ya)/reps




conf_sup_avg = apply(conf_sup,1, mean)
plot(conf_sup_avg)
lines(conf_sup_avg)


power_avg = apply(conf_sup.b<0.05, 1, mean)
points(power_avg,lwd=4)
lines(power_avg, col='red',lwd=4)
power_avg

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



mle_t_avg = mean(mle_a.b[43,1,] - mle_a.b[43,2,])
mle_t_sd = sd(mle_a.b[43,1,] - mle_a.b[43,2,])

bias_t = (mle_t_avg - (0.257 - 0.229))

upper = (bias_t + 2*mle_t_sd)/(0.257 - 0.229)
mid = bias_t/(0.257 - 0.229)
lower = (bias_t - 2*mle_t_sd)/(0.257 - 0.229)

points(6, mid ,xlim = c(0,16), ylim = c(lower-1, upper+1), xaxt='n', ylab='Bias', xlab='Algorithm')
arrows(6, lower, 6, upper, code=3, angle=90, lty=8, )


#end of trial bias
mle_avg = vector(length=S)
mle_sd = vector(length=S)
for (j in 1:S){
  mle_avg[j] = mean(mle_a.b[100,j,])
  mle_sd[j] = sd(mle_a.b[100,j,])
}

bias = mle_avg - c(0.257, 0.229)

c_vector =c('blue', 'red')
points(c(5.5, 6.5), bias, col=c_vector)
arrows(c(5.5,6.5), bias-2*mle_sd, c(5.5,6.5), bias+2*mle_sd, code=3, angle=90, lty=8, col=c_vector)
