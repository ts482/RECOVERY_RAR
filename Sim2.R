
set.seed(20)

reps = 1000

K     = 100 #number of stages in trial
S     = 2
T     = 3


nk_.c  = 150 #patients per stage
rand_ratio.c  = 0.5
Nx.c   = c(501, 1279, 324) + c(1034, 2604, 683)
Px.c   = Nx.c/sum(Nx.c)
Py.c   = matrix(c(0.140, 0.262, 0.414, 0.178, 0.233, 0.293), nrow=3)
#Py    = c(0.25,0.22) #probability of death on each treatment
#S     = length(Py[1,]) #Number of treatments
#T     = length(Px) #Number of care

alpha.c       = 1 #uninformative beta
beta.c        = 1 #uninformative beta
Y.c           = matrix(nrow=K,ncol=S)
postmean.c    = matrix(nrow=K,ncol=S)
mle.c         = matrix(nrow=K,ncol=S)
PostAlpha.c   = matrix(nrow=K,ncol=S)
PostBeta.c    = matrix(nrow=K,ncol=S)

nk.c          = matrix(round(nk_.c*rand_ratio.c), nrow=K,ncol=S)
rand_ratios.c = vector(length=K)

UPPER_LIMIT <- 0.9
LOWER_LIMIT <- 0.1

nka.c = array(dim = c(K, S, reps))
Ya.c = array(dim = c(K, S, reps))
postmean_a.c = array(dim = c(K, S, reps))
mle_a.c = array(dim=c(K, S, reps))

thetas.c = matrix(nrow=reps, ncol=K)

y.c = vector(length=3)

conf_sup.c <- matrix(nrow=K/10, ncol=reps)


# Stage 1
for (n in 1:reps){

rand_ratio = 0.5
  
for(j in 1:S){
for (h in 1:T){
  y.c[h] = rbinom(1, round(nk.c[1,j]* Px.c[h]), Py.c[h,j])
}
Y.c[1,j]         = sum(y.c)
PostAlpha.c[1,j] = alpha.c + Y.c[1,j]
PostBeta.c[1,j]  = beta.c  + nk.c[1,j] -Y.c[1,j]
postmean.c[1,j]  = PostAlpha.c[1,j]/(PostAlpha.c[1,j]+PostBeta.c[1,j])
mle.c[1,j]       = sum(Y.c[1,j])/(nk.c[1,j])
}

for(i in 2:K){
  
nk.c[i, 1] = round(nk_.c * (1-rand_ratio.c))
nk.c[i, 2] = round(nk_.c * rand_ratio.c)

for(j in 1:S){
for (h in 1:T){
  y.c[h] = rbinom(1, round(nk.c[i,j]* Px.c[h]), Py.c[h,j])
}
Y.c[i,j]         = sum(y.c)
PostAlpha.c[i,j] = alpha.c + sum(Y.c[1:i,j], na.rm=T)
PostBeta.c[i,j]  = beta.c + sum(nk.c[1:i,j]-Y.c[1:i,j], na.rm=T)
postmean.c[i,j]  = PostAlpha.c[i,j]/(PostAlpha.c[i,j]+PostBeta.c[i,j])
mle.c[i,j]       = sum(Y.c[1:i,j])/sum(nk.c[1:i,j])


}

X1 <- rbeta(100, PostAlpha.c[i, 1], PostBeta.c[i, 1])
X2 <- rbeta(100, PostAlpha.c[i, 2], PostBeta.c[i, 2])
theta <- length(X1[X1>=X2])/100
thetas.c[n,i] <- theta
k = i/(K)
rand_ratio.c <- theta^k / (theta^k + (1-theta)^k)
rand_ratios.c[i] = rand_ratio

#making sure samples don't go to 0
if (rand_ratio.c>UPPER_LIMIT){
  rand_ratio.c = UPPER_LIMIT
}
if (rand_ratio.c<LOWER_LIMIT){
  rand_ratio.c = LOWER_LIMIT
}

if (i %% 10 == 0){
  pat_df <- data.frame(t=rep(1, i*nk_.c), o=rep(0, i*nk_.c))
  pat_df[1:sum(nk.c[1:i,1]), 't'] = 0
  pat_df[1:(sum(Y.c[1:i,1])), 'o'] = 1
  pat_df[(sum(nk.c[1:i,1])+1):(sum(nk.c[1:i,1])+sum(Y.c[1:i,2])+1), 'o'] = 1
  
  fit = glm(o~t,data=pat_df, family=binomial)
  conf_sup.c[i/10,n] = summary(fit)$coef['t','Pr(>|z|)']
}
}

nka.c[,,n] = nk.c
Ya.c[,,n] = Y.c
postmean_a.c[,,n] = postmean.c
mle_a.c[,,n] = mle.c

}



sum(Ya.c)/reps


Y_avg.c = matrix(nrow=K, ncol=S)
nk_avg.c = matrix(nrow=K, ncol=S)
Y_sd.c = vector(length=S)
nk_sd.c = vector(length=S)
for (j in 1:S){
  nk_avg.c[,j] = apply(nka.c[,j,], 1, mean)
  Y_avg.c[,j] = apply(Ya.c[,j,], 1, mean, na.rm=T)
  nk_sd.c[j] = sd(nka.c[,j,])
  Y_sd.c[j] = sd(Ya.c[,j,])
}
lines(nk_avg.c[seq(10, 100, 10),2]/nk_.c*100, lwd=4, type='b', col='brown')


nk_sd

sum(nk_avg.c[,1])
sum(nk_avg.c[,2])

sum(Y_avg.c[,1])
sum(Y_avg.c[,2])


#deaths at n=6450
sum(nk_avg.c[1:43,])
sum(Y_avg.c[1:43,])


conf_sup_avg = apply(conf_sup,1, mean)
plot(conf_sup_avg)
lines(conf_sup_avg)


power_avg = apply(conf_sup.c<0.05, 1, mean)
points(power_avg, lwd=4)
lines(power_avg, col = 'brown',lwd=4)
power_avg

plot(nk_avg[,1], ylim = c(0,60), col='red', xlab='Trial completion (%)', ylab='Patients allocated')
points(nk_avg[,2], col='blue')
legend("bottomleft",legend = c('no treatment', 'dexamethasone'),col = c("red", "blue"), lwd=4)

sum(Ya)/reps
#Y_avg = matrix(nrow=K, ncol=S)
#for (j in 1:S){
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

mle_avg = vector(length= 10) #vector(length=S)
mle_sd = vector(length=10)#vector(length=S)

treatment_effect = 0.257 - 0.229

for (i in 1:10){
  mle_avg[i] = mean(((mle_a.c[i*10,1,] - mle_a.c[i*10,2,])-treatment_effect)/treatment_effect)
  mle_sd[i] = sd(((mle_a.c[i*10,1,] - mle_a.c[i*10,2,])-treatment_effect)/treatment_effect)
}



upper = mle_avg + 1.96*mle_sd/sqrt(6450)
mid = mle_avg
lower = mle_avg - 1.96*mle_sd/sqrt(6450)

points(1:10-0.2, mid, ylim=c(-0.06,0.06), col= 'blue', lwd=6)
arrows(1:10-0.2, lower, 1:10-0.2, upper, code=3, angle=90, lty=8, col='blue', lwd=3)

legend('bottomleft', legend = c(expression('T'[f]), expression('RMC'[f])),col = c('blue','red'), 
       lty = c('dashed', 'dashed'), lwd=3)
abline(v= 4.3, lty='dashed')




# end of trial treatment effect bias


treatment_effect = 0.257 - 0.229
mle_t_avg = mean(((mle_a.c[43,1,] - mle_a.c[43,2,])-treatment_effect)/treatment_effect)
mle_t_sd = sd(((mle_a.c[43,1,] - mle_a.c[43,2,])-treatment_effect)/treatment_effect)


upper = mle_t_avg + 1.96*mle_t_sd/sqrt(6450)
mid = mle_t_avg
lower = mle_t_avg - 1.96*mle_t_sd/sqrt(6450)

points(13.5, mid, col = 'blue', lwd=6)
arrows(13.5, lower, 13.5, upper, code=3, angle=90, lty=8, col='blue', lwd=2)

#end of trial bias




mle_avg = vector(length=S)
mle_sd = vector(length=S)
for (j in 1:S){
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
