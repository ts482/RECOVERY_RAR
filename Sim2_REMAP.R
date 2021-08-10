
set.seed(20)

reps = 1000

K     = 100 #number of stages in trial
S     = 2
T     = 3

nk_.d   = 150 #patients per stage
Nx.d    = c(501, 1279, 324) + c(1034, 2604, 683)
Px.d    = Nx.d/sum(Nx.d)
Py.d   = matrix(c(0.140, 0.262, 0.414, 0.178, 0.233, 0.293), nrow=3)
#Py    = c(0.25,0.22) #probability of death on each treatment
#S     = length(Py[1,]) #Number of treatments
#T     = length(Px) #Number of care

rand_ratio.d <- rep(0.5, S)
pre_ratio.d <- vector(length=2)

alpha.d       = 1 #uninformative beta
beta.d        = 1 #uninformative beta
Y.d           = matrix(nrow=K,ncol=S)
postmean.d    = matrix(nrow=K,ncol=S)
mle.d         = matrix(nrow=K,ncol=S)
PostAlpha.d   = matrix(nrow=K,ncol=S)
PostBeta.d    = matrix(nrow=K,ncol=S)

nk.d          = matrix(round(nk_.d*rand_ratio.d), nrow=K,ncol=S)
rand_ratios.d = matrix(nrow=K, ncol=S)

UPPER_LIMIT <- 0.9
LOWER_LIMIT <- 0.1

nka.d = array(dim = c(K, S, reps))
Ya.d = array(dim = c(K, S, reps))
postmean_a.d = array(dim = c(K, S, reps))
mle_a.d = array(dim=c(K, S, reps))

X.d = matrix(nrow=100, ncol=S)

prob_optimal.d = vector(length=S)
prob_optimals.d = array(dim=c(reps, S, K))

y.d = vector(length=3)

conf_sup.d <- matrix(nrow=K/10, ncol=reps)

# Stage 1
for (n in 1:reps){

rand_ratio.d <- rep(0.5, S)
  
for(j in 1:S){
for (h in 1:T){
  y.d[h] = rbinom(1, round(nk.d[1,j]* Px.d[h]), Py.d[h,j])
}
Y.d[1,j]         = sum(y.d)
PostAlpha.d[1,j] = alpha.d + Y.d[1,j]
PostBeta.d[1,j]  = beta.d  + nk.d[1,j] -Y.d[1,j]
postmean.d[1,j]  = PostAlpha.d[1,j]/(PostAlpha.d[1,j]+PostBeta.d[1,j])
mle.d[1,j]       = sum(Y.d[1,j])/(nk.d[1,j])
}

for(i in 2:K){

for(j in 1:S){
  nk.d[i, j] = round(nk_.d *rand_ratio.d[j])
for (h in 1:T){
  y.d[h] = rbinom(1, round(nk.d[i,j]* Px.d[h]), Py.d[h,j])
}
Y.d[i,j]         = sum(y.d)
PostAlpha.d[i,j] = alpha.d + sum(Y.d[1:i,j], na.rm=T)
PostBeta.d[i,j]  = beta.d + sum(nk.d[1:i,j]-Y.d[1:i,j], na.rm=T)
postmean.d[i,j]  = PostAlpha.d[i,j]/(PostAlpha.d[i,j]+PostBeta.d[i,j])
mle.d[i,j]       = sum(Y.d[1:i,j])/sum(nk.d[1:i,j])

X.d[,j] = rbeta(100, PostAlpha.d[i, j], PostBeta.d[i, j])
}

optimum = apply(X.d, 1, which.min)
for (j in 1:S){
  prob_optimal.d[j] = mean(optimum == j)
  pre_ratio.d[j] = sqrt(prob_optimal.d[j]/(sum(nk.d[1:i,j])+1))
}
prob_optimals.d[n,,i] = prob_optimal.d

rand_ratio.d = pre_ratio.d/sum(pre_ratio.d)
rand_ratios.d[i,] = rand_ratio.d

#making sure samples don't go to 0
for (j in 1:S){
if (rand_ratio.d[j]>1-(S-1)*LOWER_LIMIT){
  rand_ratio.d[j] = UPPER_LIMIT
}
if (rand_ratio.d[j]<LOWER_LIMIT){
  rand_ratio.d[j] = LOWER_LIMIT
}
}

if (i %% 10 == 0){
  pat_df <- data.frame(t=rep(1, i*nk_.d), o=rep(0, i*nk_.d))
  pat_df[1:sum(nk.d[1:i,1]), 't'] = 0
  pat_df[1:(sum(Y.d[1:i,1])), 'o'] = 1
  pat_df[(sum(nk.d[1:i,1])+1):(sum(nk.d[1:i,1])+sum(Y.d[1:i,2])+1), 'o'] = 1
  
  fit = glm(o~t,data=pat_df, family=binomial)
  conf_sup.d[i/10,n] = summary(fit)$coef['t','Pr(>|z|)']
}
}

nka.d[,,n] = nk.d
Ya.d[,,n] = Y.d
postmean_a.d[,,n] = postmean.d
mle_a.d[,,n] = mle.d
}

sum(Ya.d)/reps


Y_avg.d = matrix(nrow=K, ncol=S)
nk_avg.d = matrix(nrow=K, ncol=S)
Y_sd.d = vector(length=S)
nk_sd.d = vector(length=S)
for (j in 1:S){
  nk_avg.d[,j] = apply(nka.d[,j,], 1, mean)
  Y_avg.d[,j] = apply(Ya.d[,j,], 1, mean, na.rm=T)
  nk_sd.d[j] = sd(nka.d[,j,])
  Y_sd.d[j] = sd(Ya.d[,j,])
}
lines(nk_avg.d[seq(10, 100, 10),2]/nk_.d*100, type='b',lwd=4, col='green')

legend('bottomright', legend= c(expression('T'[f]), expression('RMC'[f]), 'FeR',  'FuR'), 
       fill=c('brown','green' ,'red' , 'blue'), col=c('brown','green' ,'red' , 'blue'), lwd=4, lty=1)

rec.prop = 6425/(nk_.d*K)*10
abline(v=rec.prop, lty='dashed')


nk_sd

sum(nk_avg.d[,1])
sum(nk_avg.d[,2])

sum(Y_avg.d[,1])
sum(Y_avg.d[,2])

#deaths at n=6450
sum(nk_avg.d[1:43,])
sum(Y_avg.d[1:43,])


conf_sup_avg = apply(conf_sup,1, mean)
plot(conf_sup_avg)
lines(conf_sup_avg)


power_avg = apply(conf_sup.d<0.05, 1, mean)
points(power_avg,lwd=4)
lines(power_avg, lwd=4,col = 'green')
abline(h=0.8, lty=4)
abline(h=0.9, lty=4)

legend('bottomright', legend= c('FeR', 'FuR', expression('RMC'[f]), expression('T'[f])),
       fill=c('blue', 'red', 'green', 'brown'), col= c('blue', 'red', 'green', 'brown'))

abline(v=rec.prop, lty='dashed')




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
  mle_avg[i] = mean(((mle_a.d[i*10,1,] - mle_a.d[i*10,2,])-treatment_effect)/treatment_effect)
  mle_sd[i] = sd(((mle_a.d[i*10,1,] - mle_a.d[i*10,2,])-treatment_effect)/treatment_effect)
}



upper = mle_avg + 1.96*mle_sd/sqrt(6450)
mid = mle_avg
lower = mle_avg - 1.96*mle_sd/sqrt(6450)

plot(1:10+0.2, mid, xlim = c(0.5, 10.5), ylim=c(min(lower)-0.01, max(upper)), 
     xaxt = 'n', xlab = 'Trial Progress(%)', 
     ylab= 'Relative bias in treatment effect', col= 'red', lwd=6)
axis(1, at = 1:10, labels = seq(10,100, 10))
arrows(1:10+0.2, lower, 1:10+0.2, upper, code=3, angle=90, lty=8, col='red', lwd=3)
abline(h=0, lty='dashed')



# end of trial treatment effect bias

treatment_effect = 0.257 - 0.229

mle_t_avg = mean(((mle_a.d[43,1,] - mle_a.d[43,2,])-treatment_effect)/treatment_effect)
mle_t_sd = sd(((mle_a.d[43,1,] - mle_a.d[43,2,])-treatment_effect)/treatment_effect)



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
