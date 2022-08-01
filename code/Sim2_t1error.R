
set.seed(20)

reps = 1000

K     = 100 #number of stages in trial
S     = 2
T     = 3


nk_.ce  = 80 #patients per stage
rand_ratio.ce  = 0.5
Nx.ce   = c(501, 1279, 324) + c(1034, 2604, 683)
Px.ce   = Nx.ce/sum(Nx.ce)
Py.ce   = matrix(c(0.140, 0.262, 0.414, 0.140, 0.262, 0.414), nrow=3)
#Py    = c(0.25,0.22) #probability of death on each treatment
#S     = length(Py[1,]) #Number of treatments
#T     = length(Px) #Number of care

alpha.ce       = 1 #uninformative beta
beta.ce        = 1 #uninformative beta
Y.ce           = matrix(nrow=K,ncol=S)
postmean.ce    = matrix(nrow=K,ncol=S)
mle.ce         = matrix(nrow=K,ncol=S)
PostAlpha.ce   = matrix(nrow=K,ncol=S)
PostBeta.ce    = matrix(nrow=K,ncol=S)

nk.ce          = matrix(round(nk_.ce*rand_ratio.ce), nrow=K,ncol=S)
rand_ratios.ce = vector(length=K)

UPPER_LIMIT <- 0.9
LOWER_LIMIT <- 0.1

nka.ce = array(dim = c(K, S, reps))
Ya.ce = array(dim = c(K, S, reps))
postmean_a.ce = array(dim = c(K, S, reps))
mle_a.ce = array(dim=c(K, S, reps))

thetas.ce = matrix(nrow=reps, ncol=K)

y.ce = vector(length=3)

conf_sup.ce <- matrix(nrow=K/10, ncol=reps)


# Stage 1
for (n in 1:reps){

rand_ratio.ce = 0.5
  
for(j in 1:S){
for (h in 1:T){
  y.ce[h] = rbinom(1, round(nk.ce[1,j]* Px.ce[h]), Py.ce[h,j])
}
Y.ce[1,j]         = sum(y.ce)
PostAlpha.ce[1,j] = alpha.ce + Y.ce[1,j]
PostBeta.ce[1,j]  = beta.ce  + nk.ce[1,j] -Y.ce[1,j]
postmean.ce[1,j]  = PostAlpha.ce[1,j]/(PostAlpha.ce[1,j]+PostBeta.ce[1,j])
mle.ce[1,j]       = sum(Y.ce[1,j])/(nk.ce[1,j])
}

for(i in 2:K){
  
nk.ce[i, 1] = round(nk_.ce * (1-rand_ratio.ce))
nk.ce[i, 2] = round(nk_.ce * rand_ratio.ce)

for(j in 1:S){
for (h in 1:T){
  y.ce[h] = rbinom(1, round(nk.ce[i,j]* Px.ce[h]), Py.ce[h,j])
}
Y.ce[i,j]         = sum(y.ce)
PostAlpha.ce[i,j] = alpha.ce + sum(Y.ce[1:i,j], na.rm=T)
PostBeta.ce[i,j]  = beta.ce + sum(nk.ce[1:i,j]-Y.ce[1:i,j], na.rm=T)
postmean.ce[i,j]  = PostAlpha.ce[i,j]/(PostAlpha.ce[i,j]+PostBeta.ce[i,j])
mle.ce[i,j]       = sum(Y.ce[1:i,j])/sum(nk.ce[1:i,j])


}

if (i > 28){

if (i %% 7 == 0){
X1 <- rbeta(100, PostAlpha.ce[i-28, 1], PostBeta.ce[i-28, 1])
X2 <- rbeta(100, PostAlpha.ce[i-28, 2], PostBeta.ce[i-28, 2])
theta <- length(X1[X1>=X2])/100
thetas.ce[n,i] <- theta
k = i/(K)
rand_ratio.ce <- theta^k / (theta^k + (1-theta)^k)
rand_ratios.ce[i] = rand_ratio.ce

#making sure samples don't go to 0
if (rand_ratio.ce>UPPER_LIMIT){
  rand_ratio.ce = UPPER_LIMIT
}
if (rand_ratio.ce<LOWER_LIMIT){
  rand_ratio.ce = LOWER_LIMIT
}
}
}
if (i %% 10 == 0){
  pat_df <- data.frame(t=rep(1, i*nk_.ce), o=rep(0, i*nk_.ce))
  pat_df[1:sum(nk.ce[1:i,1]), 't'] = 0
  pat_df[1:(sum(Y.ce[1:i,1])), 'o'] = 1
  pat_df[(sum(nk.ce[1:i,1])+1):(sum(nk.ce[1:i,1])+sum(Y.ce[1:i,2])+1), 'o'] = 1
  
  fit = glm(o~t,data=pat_df, family=binomial)
  conf_sup.ce[i/10,n] = summary(fit)$coef['t','Pr(>|z|)']
}
}

nka.ce[,,n] = nk.ce
Ya.ce[,,n] = Y.ce
postmean_a.ce[,,n] = postmean.ce
mle_a.ce[,,n] = mle.ce

}


sum(Ya.ce)/reps


Y_avg.ce = matrix(nrow=K, ncol=S)
nk_avg.ce = matrix(nrow=K, ncol=S)
Y_sd.ce = vector(length=S)
nk_sd.ce = vector(length=S)
for (j in 1:S){
  nk_avg.ce[,j] = apply(nka.ce[,j,], 1, mean)
  Y_avg.ce[,j] = apply(Ya.ce[,j,], 1, mean, na.rm=T)
  nk_sd.ce[j] = sd(nka.ce[,j,])
  Y_sd.ce[j] = sd(Ya.ce[,j,])
}



t1e_metric.ce = vector(length=10)

for (i in 1:10){
t1e_metric.ce[i] = mean(conf_sup.ce[i,]<0.05) #0.05 exactly!
}

t1e_metric.ce

plot(t1e_metric.ce*100, lwd=4, type='b', col='blue', ylim=c(0,15),
     axes=FALSE, xlab='Trial progress(%)', ylab='Family-wise error rate')
axis(1, seq(1,10), labels = seq(12.5, 125,12.5))
axis(2, seq(0, 15, 2.5), labels = seq(0, 15, 2.5))







conf_sup_avg = apply(conf_sup,1, mean)
plot(conf_sup_avg)
lines(conf_sup_avg)


power_avg = apply(conf_sup.ce<0.05, 1, mean)
power_avg = c(0, power_avg)
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
  mle_avg[i] = mean(((mle_a.ce[i*10,1,] - mle_a.ce[i*10,2,])-treatment_effect)/treatment_effect)
  mle_sd[i] = sd(((mle_a.ce[i*10,1,] - mle_a.ce[i*10,2,])-treatment_effect)/treatment_effect)
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
mle_t_avg = mean(((mle_a.ce[80,1,] - mle_a.ce[80,2,])-treatment_effect)/treatment_effect)
mle_t_sd = sd(((mle_a.ce[80,1,] - mle_a.ce[80,2,])-treatment_effect)/treatment_effect)


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
