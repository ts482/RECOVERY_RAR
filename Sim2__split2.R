
set.seed(20)

reps = 1000

K     = 100 #number of stages in trial
S     = 2
T     = 3

nk_.g   = 150
Nx.g    = c(501, 1279, 324) + c(1034, 2604, 683)
Px.g    = Nx.g/sum(Nx.g)
nk_.g   = round(nk_.g*Px.g) #patients per stage
Py.g   = matrix(c(0.140, 0.262, 0.414, 0.178, 0.233, 0.293), nrow=3)
#Py    = c(0.25,0.22) #probability of death on each treatment
#S     = length(Py[1,]) #Number of treatments
#T     = length(Px) #Number of care

alpha.g       = 1 #uninformative beta
beta.g        = 1 #uninformative beta
Y.g           = array(dim=c(K,S,T))
postmean.g    = array(dim=c(K,S,T))
mle.g         = array(dim=c(K,S,T))
PostAlpha.g   = array(dim=c(K,S,T))
PostBeta.g    = array(dim=c(K,S,T))

nk.g          = array(dim=c(K,S,T))
for (h in 1:length(nk_.g)){
  nk.g[,,h] = round(nk_.g[h]/S)
}

rand_ratio.g = rep(0.5, T)
rand_ratios.g = matrix(nrow=K, ncol=T)

UPPER_LIMIT <- 0.9
LOWER_LIMIT <- 0.1

nka.g = array(dim = c(K, S, T,reps))
Ya.g = array(dim = c(K, S, T,reps))
postmean_a.g = array(dim = c(K, S, T,reps))
mle_a.g = array(dim=c(K, S,T, reps))

thetas.g = array(dim = c(reps,K,T))
#y = vector(length=3)

conf_sup.g <- array(dim=c(K/10, T,reps))

# Stage 1
for (n in 1:reps){

rand_ratio.g = rep(0.5, T)
  
for (h in 1:T){
for(j in 1:S){
y.g = rbinom(1, round(nk.g[1,j,h]), Py.g[h,j])
Y.g[1,j,h]         = y.g
PostAlpha.g[1,j,h] = alpha.g + Y.g[1,j,h]
PostBeta.g[1,j,h]  = beta.g  + nk.g[1,j,h] -Y.g[1,j,h]
postmean.g[1,j,h]  = PostAlpha.g[1,j,h]/(PostAlpha.g[1,j,h]+PostBeta.g[1,j,h])
mle.g[1,j,h]       = sum(Y.g[1,j,h])/(nk.g[1,j,h])
}
}

for(i in 2:K){

nk.g[i, 1,] = round(nk_.g * (1-rand_ratio.g))
nk.g[i, 2,] = round(nk_.g * rand_ratio.g)
for (h in 1:T){
for(j in 1:S){
  y.g = rbinom(1, round(nk.g[i,j,h]), Py.g[h,j])

Y.g[i,j,h]         = y.g
PostAlpha.g[i,j,h] = alpha.g + sum(Y.g[1:i,j,h], na.rm=T)
PostBeta.g[i,j,h]  = beta.g + sum(nk.g[1:i,j,h]-Y.g[1:i,j,h], na.rm=T)
postmean.g[i,j,h]  = PostAlpha.g[i,j,h]/(PostAlpha.g[i,j,h]+PostBeta.g[i,j,h])
mle.g[i,j,h]       = sum(Y.g[1:i,j,h])/sum(nk.g[1:i,j,h])


}

X1 <- rbeta(100, PostAlpha.g[i, 1,h], PostBeta.g[i, 1,h])
X2 <- rbeta(100, PostAlpha.g[i, 2,h], PostBeta.g[i, 2,h])
theta <- length(X1[X1>=X2])/100
thetas.g[n,i,h] <- theta
k = i/(K)
rand_ratio.g[h] <- theta^k / (theta^k + (1-theta)^k)
rand_ratios.g[i,h] = rand_ratio.g[h]

#making sure samples don't go to 0
if (rand_ratio.g[h]>UPPER_LIMIT){
  rand_ratio.g[h] = UPPER_LIMIT
}
if (rand_ratio.g[h]<LOWER_LIMIT){
  rand_ratio.g[h] = LOWER_LIMIT
}

if (i %% 10 == 0){
  pat_df <- data.frame(t=rep(1, i*nk_.g[h]), o=rep(0, i*nk_.g[h]))
  pat_df[1:sum(nk.g[1:i,1,h]), 't'] = 0
  pat_df[1:(sum(Y.g[1:i,1,h])), 'o'] = 1
  pat_df[(sum(nk.g[1:i,1,h])+1):(sum(nk.g[1:i,1,h])+sum(Y.g[1:i,2,h])), 'o'] = 1
  
  fit = glm(o~t,data=pat_df, family=binomial)
  conf_sup.g[i/10,h,n] = summary(fit)$coef['t','Pr(>|z|)']
}
}

}

nka.g[,,,n] = nk.g
Ya.g[,,,n] = Y.g
postmean_a.g[,,,n] = postmean.g
mle_a.g[,,,n] = mle.g
}

sum(Ya.g)/reps

nk_avg.g = array(dim=c(K,S,T))
Y_avg.g  = array(dim=c(K,S,T))

for (h in 1:T){
for (j in 1:S){
  nk_avg.g[,j,h] = apply(nka.g[,j,h,], 1, mean)
  Y_avg.g[,j,h]  = apply(Ya.g[,j,h,], 1, mean, na.rm=T)
}
}


for (h in 1:T){
print(sum(nk_avg.g[,1,h]))
print(sum(nk_avg.g[,2,h]))
}

sum(nk_avg.g)

for (h in 1:T){
  print(sum(Y_avg.g[,1,h]))
  print(sum(Y_avg.g[,2,h]))
}


#mortality at n=6450
sum(nk_avg.g[1:43,,])
sum(Y_avg.g[1:43,,])



for (h in 1:T){
  conf_sup_avg = apply(conf_sup[,h,],1, mean)
  plot(conf_sup_avg)
  lines(conf_sup_avg)
}


h=1

c_vector = c('red', 'blue', 'brown')

power_avg = apply(conf_sup.h[,h,]<0.05, 1, mean)
plot(power_avg, xlab='Trial Progress (%)', ylab = 'Power', ylim= c(0, 1), xaxt = 'n')

for (h in 1:T){
  power_avg = apply(conf_sup.g[,h,]<0.05, 1, mean)
  points(power_avg, lwd=2)
  lines(power_avg,lty='dashed', col=c_vector[h], lwd=3)
  
  power_avg = apply(conf_sup.h[,h,]<0.05, 1, mean)
  points(power_avg, lwd=2)
  lines(power_avg, col=c_vector[h], lwd=3)
  
  
  
}

#legend('bottomright', legend = c('REMAP-CAP', 'Tuning'), fill=c('red', 'blue'))

legend("bottomright",legend = c('subgroup iii', 'subgroup ii', 'subgroup i'),col = c_vector[c(3,2,1)], lwd=2)
legend('bottom', legend=c(expression('T'[s]), expression('RMC'[s])), lty = c('solid','dashed'))

axis(1, 1:10, 1:10*10)

abline(h=0.8, lty=4)
abline(h=0.9, lty=4)
abline(v=rec.prop, lty='dashed')






c_vector = c('red', 'blue', 'brown')
plot(nk_avg.g[1:10*10,2,1]/nk_.g[1]*100, ylim = c(0,100), col='red', xlab= 'Trial progress (%)', 
     ylab='Proportion of patients allocated dexamethasone(%)', type = 'b', xaxt='n')
axis(1, labels = 1:10*10, at = 1:10)

for (h in 1:T){
#plot(nk_avg[1:10*10,2,h], ylim = c(0,nk_[h]), col='red', xlab= 'trial progress (%)', 
#     ylab='number of patients')
points(nk_avg.g[1:10*10,2,h]/nk_.g[h]*100, col=c_vector[h], type = 'b', lwd=3)
#legend("bottom",legend = c('no treatment', 'dexamethasone'),col = c("red", "blue"), lwd=3)
#plot(nk_avg2[1:10*10,2,h], ylim = c(0,nk_[h]), col='red', xlab= 'trial progress (%)', 
#     ylab='number of patients')
points(nk_avg.h[1:10*10,2,h]/nk_.h[h]*100, col=c_vector[h], type = 'b', lty = 'dashed', lwd=3)
#legend("bottom",legend = c('no treatment', 'dexamethasone'),col = c("red", "blue"), lwd=4)
}

legend("right",legend = c('subgroup iii', 'subgroup ii','subgroup i'),col = c_vector[c(3,2,1)], lwd=4)
legend('bottomleft', legend = c(expression('T'[s]), expression('RMC'[s])), lwd =3, lty = c('solid','dashed'))
abline(v=rec.prop, lty='dashed')

sum(Ya)/reps
#Y_avg = matrix(nrow=K, ncol=S)
#for (j in 1:S){
#  Y_avg[,j] = apply(Ya[,j,], 1, mean)
#}

treatment_superior = matrix(nrow=K, ncol=T)
for (h in 1:T){
for (i in 1:K){
  treatment_superior[i,h] = mean(postmean_a[i,1,h,] >= postmean_a[i,2,h,])
}
}

for (h in 1:T){
plot(treatment_superior[,h])
}

for (h in 1:T){
thetas_avg = apply(thetas.g[,,h], 2, mean)
plot(thetas_avg)
}

for (h in 2:T){
above_theta = thetas.g > 0.95
above_theta_avg = apply(above_theta, 2, mean)
plot(above_theta_avg)
}
#plot(postmean[,1],ylim=range(postmean),type="l",lwd=4,
#     ylab="Posterior mean",xlab="Interim Analysis")
#lines(postmean[,2],type="l",lwd=4,col="red")
#legend("bottomright",c("Treatment 1","Treatment 2"),col=c("black","red"),lwd=4,bty="n")

#lines(mle[,1],ylim=range(mle),type="l",lwd=4,lty=2)
#lines(mle[,2],type="l",lwd=4,col="red",lty=2)

mle_avg = matrix(nrow=T, ncol=S)
mle_sd = matrix(nrow=T, ncol=S)

for (h in 1:T){
  for (j in 1:S){
    mle_avg[h,j] = mean(mle_a[100,j,h,])
    mle_sd[h,j] = sd(mle_a[100,j,h,])
  }
}  


bias = mle_avg - Py
c_vector =c('blue', 'red','yellow')
for (h in 1:T){
  points(c(9.5+(h-1)*18,10.5+(h-1)*18), bias[h,], col=c_vector)
  arrows(c(9.5+(h-1)*18,10.5+(h-1)*18), bias[h,]-2*mle_sd[h,], c(9.5+(h-1)*18,10.5+(h-1)*18), 
         bias[h,]+2*mle_sd[h,], code=3, angle=90, lty=8, col=c_vector)
}


#treatment effect bias at n=6450

mle_t_avg = vector(length=3)
mle_t_sd = vector(length=3)

treatment_effect = Py.g[,1] - Py.g[,2]

for (h in 1:T){
  mle_t_avg[h] = mean(((mle_a.g[43,1, h,] - mle_a.g[43,2,h,])- treatment_effect[h])/treatment_effect[h])
  mle_t_sd[h] = sd(((mle_a.g[43,1, h,] - mle_a.g[43,2,h,])- treatment_effect[h])/treatment_effect[h])
}


sample_size = round(6450*Px.g)


upper = mle_t_avg + 1.96*mle_t_sd/sqrt(sample_size)
mid = mle_t_avg
lower = mle_t_avg - 1.96*mle_t_sd/sqrt(sample_size)


for (h in 1:T){
  points(1.5+(h-1)*4, mid[h], col='blue', lwd=6)
  arrows(c(1.5+(h-1)*4), lower[h], c(1.5+(h-1)*4),
         upper[h], code=3, angle=90, lty=8, col='blue', lwd=2)
}

