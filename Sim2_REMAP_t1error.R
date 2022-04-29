
set.seed(20)

reps = 1000

K     = 100 #number of stages in trial
S     = 2
T     = 3

nk_.de   = 80 #patients per stage
Nx.de    = c(501, 1279, 324) + c(1034, 2604, 683)
Px.de    = Nx.de/sum(Nx.de)
Py.de   = matrix(c(0.140, 0.262, 0.414, 0.140, 0.262, 0.414), nrow=3)
#Py    = c(0.25,0.22) #probability of death on each treatment
#S     = length(Py[1,]) #Number of treatments
#T     = length(Px) #Number of care

rand_ratio.de <- rep(0.5, S)
pre_ratio.de <- vector(length=2)

alpha.de       = 1 #uninformative beta
beta.de        = 1 #uninformative beta
Y.de           = matrix(nrow=K,ncol=S)
postmean.de    = matrix(nrow=K,ncol=S)
mle.de         = matrix(nrow=K,ncol=S)
PostAlpha.de   = matrix(nrow=K,ncol=S)
PostBeta.de    = matrix(nrow=K,ncol=S)

nk.de          = matrix(round(nk_.de*rand_ratio.de), nrow=K,ncol=S)
rand_ratios.de = matrix(nrow=K, ncol=S)

UPPER_LIMIT <- 0.9
LOWER_LIMIT <- 0.1

nka.de = array(dim = c(K, S, reps))
Ya.de = array(dim = c(K, S, reps))
postmean_a.de = array(dim = c(K, S, reps))
mle_a.de = array(dim=c(K, S, reps))

X.de = matrix(nrow=100, ncol=S)

prob_optimal.de = vector(length=S)
prob_optimals.de = array(dim=c(reps, S, K))

y.de = vector(length=3)

conf_sup.de <- matrix(nrow=K/10, ncol=reps)

# Stage 1
for (n in 1:reps){

rand_ratio.de <- rep(0.5, S)
  
for(j in 1:S){
for (h in 1:T){
  y.de[h] = rbinom(1, round(nk.de[1,j]* Px.de[h]), Py.de[h,j])
}
Y.de[1,j]         = sum(y.de)
PostAlpha.de[1,j] = alpha.de + Y.de[1,j]
PostBeta.de[1,j]  = beta.de  + nk.de[1,j] -Y.de[1,j]
postmean.de[1,j]  = PostAlpha.de[1,j]/(PostAlpha.de[1,j]+PostBeta.de[1,j])
mle.de[1,j]       = sum(Y.de[1,j])/(nk.de[1,j])
}

for(i in 2:K){

for(j in 1:S){
  nk.de[i, j] = round(nk_.de *rand_ratio.de[j])
for (h in 1:T){
  y.de[h] = rbinom(1, round(nk.de[i,j]* Px.de[h]), Py.de[h,j])
}
Y.de[i,j]         = sum(y.de)
PostAlpha.de[i,j] = alpha.de + sum(Y.de[1:i,j], na.rm=T)
PostBeta.de[i,j]  = beta.de + sum(nk.de[1:i,j]-Y.de[1:i,j], na.rm=T)
postmean.de[i,j]  = PostAlpha.de[i,j]/(PostAlpha.de[i,j]+PostBeta.de[i,j])
mle.de[i,j]       = sum(Y.de[1:i,j])/sum(nk.de[1:i,j])
}
if (i > 28){
    
if (i %% 7 == 0){
for (j in 1:S){
  X.de[,j] = rbeta(100, PostAlpha.de[i-28, j], PostBeta.de[i-28, j])
}
optimum = apply(X.de, 1, which.min)
for (j in 1:S){
  prob_optimal.de[j] = mean(optimum == j)
  pre_ratio.de[j] = sqrt(prob_optimal.de[j]/(sum(nk.de[1:(i-28),j])+1))
}
prob_optimals.de[n,,i] = prob_optimal.de

rand_ratio.de = pre_ratio.de/sum(pre_ratio.de)
rand_ratios.de[i,] = rand_ratio.de

#making sure samples don't go to 0
for (j in 1:S){
if (rand_ratio.de[j]>1-(S-1)*LOWER_LIMIT){
  rand_ratio.de[j] = UPPER_LIMIT
}
if (rand_ratio.de[j]<LOWER_LIMIT){
  rand_ratio.de[j] = LOWER_LIMIT
}
}
}
}
if (i %% 10 == 0){
  pat_df <- data.frame(t=rep(1, i*nk_.de), o=rep(0, i*nk_.de))
  pat_df[1:sum(nk.de[1:i,1]), 't'] = 0
  pat_df[1:(sum(Y.de[1:i,1])), 'o'] = 1
  pat_df[(sum(nk.de[1:i,1])+1):(sum(nk.de[1:i,1])+sum(Y.de[1:i,2])+1), 'o'] = 1
  
  fit = glm(o~t,data=pat_df, family=binomial)
  conf_sup.de[i/10,n] = summary(fit)$coef['t','Pr(>|z|)']
}
}

nka.de[,,n] = nk.de
Ya.de[,,n] = Y.de
postmean_a.de[,,n] = postmean.de
mle_a.de[,,n] = mle.de
}

sum(Ya.de)/reps

mean(conf_sup.de[8,]<0.05) #0.05 exactly again!


Y_avg.de = matrix(nrow=K, ncol=S)
nk_avg.de = matrix(nrow=K, ncol=S)
Y_sd.de = vector(length=S)
nk_sd.de = vector(length=S)
for (j in 1:S){
  nk_avg.de[,j] = apply(nka.de[,j,], 1, mean)
  Y_avg.de[,j] = apply(Ya.de[,j,], 1, mean, na.rm=T)
  nk_sd.de[j] = sd(nka.de[,j,])
  Y_sd.de[j] = sd(Ya.de[,j,])
}



lines(nk_avg.de[c(1,1:10*10),2]/nk_.de*100, type='b',lwd=4, col='green')

legend('bottomright', legend= c(expression('T'[f]), expression('RMC'[f]), 'FeR',  'FuR'), 
       col=c('brown','green' ,'red' , 'blue'), lwd=4)

#legend('bottomright', legend= c('Adaptive randomization', 'Fixed randomization'), 
#       col=c('brown','red'), lwd=4)

abline(v=4.4, lty='dashed')


nk_sd

sum(nk_avg.de[1:80,1])
sum(nk_avg.de[1:80,2])

sum(Y_avg.de[,1])
sum(Y_avg.de[,2])

#deaths at n=6450
sum(nk_avg.de[1:80,])
sum(Y_avg.de[1:80,])


conf_sup_avg = apply(conf_sup,1, mean)
plot(conf_sup_avg)
lines(conf_sup_avg)


power_avg = apply(conf_sup.de<0.05, 1, mean)
power_avg = c(0, power_avg)
points(power_avg,lwd=4)
lines(power_avg, lwd=4,col = 'green')
abline(h=0.8, lty=4)
abline(h=0.9, lty=4)

legend('bottomright', legend= c('FeR', 'FuR', expression('RMC'[f]), expression('T'[f])),
       fill=c('blue', 'red', 'green', 'brown'), col= c('blue', 'red', 'green', 'brown'))

#abline(v=rec.prop, lty='dashed')




plot(nk_avg[,1], ylim = c(0,60), col='red', xlab='Trial completion (%)', ylab='Patients allocated')
points(nk_avg[,2], col='blue')
legend("bottomleft",legend = c('no treatment', 'dexamethasone'),col = c("red", "blue"), lwd=4)

sum(Ya)/reps

treatment_superior = vector(length=K)
for (i in 1:K){
  treatment_superior[i] = mean(postmean_a[i,1,] >= postmean_a[i,2,])
}

plot(treatment_superior)

thetas_avg = apply(prob_optimals[,2,], 2, mean)
plot(thetas_avg)
above_theta = prob_optimals[,2,] > 0.95
above_theta_avg = apply(above_theta, 2, mean)
plot(above_theta_avg)

#plot(postmean[,1],ylim=range(postmean),type="l",lwd=4,
#     ylab="Posterior mean",xlab="Interim Analysis")
#lines(postmean[,2],type="l",lwd=4,col="red")
#legend("bottomright",c("Treatment 1","Treatment 2"),col=c("black","red"),lwd=4,bty="n")

#lines(mle[,1],ylim=range(mle),type="l",lwd=4,lty=2)
#lines(mle[,2],type="l",lwd=4,col="red",lty=2)





#trial progress bias

treatment_effect = 0.257 - 0.229

for (i in 1:10){
  mle_avg[i] = mean(((mle_a.de[i*10,1,] - mle_a.de[i*10,2,])-treatment_effect)/treatment_effect)
  mle_sd[i] = sd(((mle_a.de[i*10,1,] - mle_a.de[i*10,2,])-treatment_effect)/treatment_effect)
}



upper = mle_avg + 1.96*mle_sd/sqrt(6450)
mid = mle_avg
lower = mle_avg - 1.96*mle_sd/sqrt(6450)

plot(1:10+0.2, mid, xlim = c(0.5, 10.5), ylim=c(min(lower)-0.01, 0.01), 
     xaxt = 'n', xlab = 'Trial Progress(%)', 
     ylab= 'Relative bias in treatment effect', col= 'red', lwd=6)
axis(1, seq(0,10,2), labels = seq(0, 125,25))
arrows(1:10+0.2, lower, 1:10+0.2, upper, code=3, angle=90, lty=8, col='red', lwd=3)
abline(h=0, lty='dashed')



# end of trial treatment effect bias

treatment_effect = 0.257 - 0.229

mle_t_avg = mean(((mle_a.de[80,1,] - mle_a.de[80,2,])-treatment_effect)/treatment_effect)
mle_t_sd = sd(((mle_a.de[80,1,] - mle_a.de[80,2,])-treatment_effect)/treatment_effect)



upper = mle_t_avg + 1.96*mle_t_sd/sqrt(6450)
mid = mle_t_avg
lower = mle_t_avg - 1.96*mle_t_sd/sqrt(6450)

points(14.5, mid, col = 'red', lwd=6)
arrows(14.5, lower, 14.5, upper, code=3, angle=90, lty=8, col='red', lwd=3)


axis(1, labels = c('Subgroup i', 'Subgroup ii','subgroup iii', 'Full cohort'), at = c(2,6,10,14))
legend('topleft', legend = c(expression('T'[s]), expression('RMC'[s])),
       col = c('blue','red'), lty = c('dashed','dashed'), lwd=2)
legend('topright', legend = c(expression('T'[f]), expression('RMC'[f])),
       col = c('blue','red'), lty = c('dashed','dashed'), lwd=2)

abline(v = 12)
abline(h=0, lty='dashed')

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
points(c(13.5, 14.5), bias, col=c_vector)
arrows(c(13.5,14.5), bias-2*mle_sd, c(13.5,14.5), bias+2*mle_sd, code=3, angle=90, lty=8, col=c_vector)
