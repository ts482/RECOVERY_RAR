
set.seed(20)

reps = 1000

K     = 100 #number of stages in trial
nk_.e   = 80 #patients per stage
Nx.e    = c(501, 1279, 324) + c(1034, 2604, 683)
Px.e    = Nx/sum(Nx)
nk_.e   = round(nk_*Px) #patients per stage
Py.e   = matrix(c(0.140, 0.262, 0.414, 0.178, 0.233, 0.293), nrow=3)
#Py    = c(0.25,0.22) #probability of death on each treatment
S     = length(Py[1,]) #Number of treatments
T     = length(Px) #Number of care


rand_ratio.e <- c(2/3, 1/3)
pre_ratio.e <- vector(length=2)

alpha.e       = 1 #uninformative beta
beta.e        = 1 #uninformative beta
Y.e           = array(dim=c(K,S,T))
postmean.e    = array(dim=c(K,S,T))
mle.e         = array(dim=c(K,S,T))
PostAlpha.e   = array(dim=c(K,S,T))
PostBeta.e    = array(dim=c(K,S,T))

nk.e          = matrix(nrow=K,ncol=S)
for (j in 1:S){
  nk.e[,j] = rand_ratio.e[j]*nk_.e
}

nk.e          = array(dim=c(K,S,T))
for (h in 1:length(nk_.e)){
  for (j in 1:S){
  nk[,j,h] = round(nk_.e[h] * rand_ratio.e[j])
  }
}
head(nk.e)

rand_ratios.e = matrix(nrow=K, ncol=S)

UPPER_LIMIT <- 0.9
LOWER_LIMIT <- 0.1

nka.e = array(dim = c(K, S, T,reps))
Ya.e = array(dim = c(K, S, T,reps))
postmean_a.e = array(dim = c(K, S, T,reps))
mle_a.e = array(dim=c(K, S,T, reps))

X.e = matrix(nrow=100, ncol=S)

y.e = vector(length=3)

prob_optimal.e = vector(length=S)
prob_optimals.e = array(dim=c(reps, S, K))


conf_sup3 <- array(dim=c(K/10, T,reps))
# Stage 1
for (n in 1:reps){
  for (h in 1:T){
    for(j in 1:S){
      y = rbinom(1, round(nk[1,j,h]), Py[h,j])
      Y[1,j,h]         = y
      PostAlpha[1,j,h] = alpha + Y[1,j,h]
      PostBeta[1,j,h]  = beta  + nk[1,j,h] -Y[1,j,h]
      postmean[1,j,h]  = PostAlpha[1,j,h]/(PostAlpha[1,j,h]+PostBeta[1,j,h])
      mle[1,j,h]       = sum(Y[1,j,h])/(nk[1,j,h])
    }
  }
  
for(i in 2:K){

#n_stand = rbinom(1, nk_, prob = c(rand_ratio[1], rand_ratio[2]))
#nk[i, ] = c(n_stand, nk_-n_stand)
  
for (h in 1:T){
  for(j in 1:S){
    y = rbinom(1, round(nk[i,j,h]), Py[h,j])
    
    Y[i,j,h]         = y
    PostAlpha[i,j,h] = alpha + sum(Y[1:i,j,h], na.rm=T)
    PostBeta[i,j,h]  = beta + sum(nk[1:i,j,h]-Y[1:i,j,h], na.rm=T)
    postmean[i,j,h]  = PostAlpha[i,j,h]/(PostAlpha[i,j,h]+PostBeta[i,j,h])
    mle[i,j,h]       = sum(Y[1:i,j,h])/sum(nk[1:i,j,h])
    
    
  }
  
#X[,j] = rbeta(100, PostAlpha[i, j], PostBeta[i, j])

#optimum = apply(X, 1, which.min)
#for (j in 1:S){
#  prob_optimal[j] = mean(optimum == j)
#}
#prob_optimals[n,,i] = prob_optimal

if (i %% 10 == 0){
  pat_df <- data.frame(t=rep(1, i*nk_[h]), o=rep(0, i*nk_[h]))
  pat_df[1:sum(nk[1:i,1,h]), 't'] = 0
  pat_df[1:(sum(Y[1:i,1,h])), 'o'] = 1
  pat_df[(sum(nk[1:i,1,h])+1):(sum(nk[1:i,1,h])+sum(Y[1:i,2,h])+1), 'o'] = 1
  
  fit = glm(o~t,data=pat_df, family=binomial)
  conf_sup3[i/10,h,n] = summary(fit)$coef['t','Pr(>|z|)']
  
}
}
}

  
nka[,,,n] = nk
Ya[,,,n] = Y
postmean_a[,,,n] = postmean
mle_a[,,,n] = mle
}

sum(Ya)/reps

Y_avg = matrix(nrow=K, ncol=S)
nk_avg = matrix(nrow=K, ncol=S)
Y_sd = vector(length=S)
nk_sd = vector(length=S)
for (j in 1:S){
  nk_avg[,j] = apply(nka[,j,], 1, mean)
  Y_avg[,j] = apply(Ya[,j,], 1, mean, na.rm=T)
  nk_sd[j] = sd(nka[,j,])
  Y_sd[j] = sd(Ya[,j,])
}
nk_sd

sum(nk_avg[,1])
sum(nk_avg[,2])

sum(Y_avg[,1])
sum(Y_avg[,2])



lines(nk_avg[seq(10, 100, 10),2], type='b', col='blue')
points(nk_avg[,2], col='blue')
legend("bottomleft",legend = c('no treatment', 'dexamethasone'),col = c("red", "blue"), lwd=4)

sum(Ya)/reps




conf_sup_avg = apply(conf_sup,1, mean)
plot(conf_sup_avg)
lines(conf_sup_avg)


power_avg = apply(conf_sup<0.05, 1, mean)
points(power_avg)
lines(power_avg, col='red')


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
  points(c(5.5+(h-1)*18,6.5+(h-1)*18), bias[h,], col=c_vector)
  arrows(c(5.5+(h-1)*18,6.5+(h-1)*18), bias[h,]-2*mle_sd[h,], c(5.5+(h-1)*18,6.5+(h-1)*18), 
         bias[h,]+2*mle_sd[h,], code=3, angle=90, lty=8, col=c_vector)
}
