rm(list=ls(all=TRUE))

set.seed(20)

reps = 100

K     = 100 #number of stages in trial
nk_   = 50 #patients per stage
Nx    = c(501, 1279, 324) + c(1034, 2604, 683)
Px    = Nx/sum(Nx)
Py   = matrix(c(0.140, 0.262, 0.414, 0.178, 0.233, 0.293), nrow=3)
#Py    = c(0.25,0.22) #probability of death on each treatment
S     = length(Py[1,]) #Number of treatments
T     = length(Px) #Number of care
Kx     = round(K*Px) #number of stages in trial
Ky = cumsum(Kx)

rand_ratio <- rep(0.5, S)
pre_ratio <- vector(length=2)

alpha       = 1 #uninformative beta
beta        = 1 #uninformative beta
Y           = matrix(nrow=K,ncol=S)
postmean    = matrix(nrow=K,ncol=S)
mle         = matrix(nrow=K,ncol=S)
PostAlpha   = matrix(nrow=K,ncol=S)
PostBeta    = matrix(nrow=K,ncol=S)

nk          = matrix(25, nrow=K,ncol=S)
rand_ratios = matrix(nrow=K, ncol=S)

UPPER_LIMIT <- 0.9
LOWER_LIMIT <- 0.1

nka = array(dim = c(K, S, reps))
Ya = array(dim = c(K, S, reps))
postmean_a = array(dim = c(K, S, reps))

X = matrix(nrow=100, ncol=S)

prob_optimal = vector(length=S)
prob_optimals = array(dim=c(reps, S, K))

start = 1
# Stage 1
for (h in 1:T){
for (n in start:Ky[h]){
for(j in 1:S){
y = rbinom(1, round(nk[1,j]* Px[h]), Py[h,j])
Y[1,j]         = y
PostAlpha[1,j] = alpha + Y[1,j]
PostBeta[1,j]  = beta  + nk[1,j] -Y[1,j]
postmean[1,j]  = PostAlpha[1,j]/(PostAlpha[1,j]+PostBeta[1,j])
mle[1,j]       = sum(Y[1,j])/(1*nk[1,j])
}

for(i in 2:K){

for(j in 1:S){
  nk[i, j] = round(nk_ *rand_ratio[j])
y = rbinom(1, round(nk[i,j]* Px[h]), Py[h,j])
Y[i,j]         = y
PostAlpha[i,j] = sum(Y[1:i,j], na.rm=T)
PostBeta[i,j]  = sum(nk[1:i,j]-Y[1:i,j], na.rm=T)
postmean[i,j]  = PostAlpha[i,j]/(PostAlpha[i,j]+PostBeta[i,j])
mle[i,j]       = sum(Y[1:i,j])/(i*nk[i,j])

X[,j] = rbeta(100, PostAlpha[i, j], PostBeta[i, j])
}

optimum = apply(X, 1, which.min)
for (j in 1:S){
  prob_optimal[j] = mean(optimum == j)
  pre_ratio[j] = sqrt(prob_optimal[j]/(sum(nk[1:i,j])+1))
}
prob_optimals[n,,i] = prob_optimal

rand_ratio = pre_ratio/sum(pre_ratio)
rand_ratios[i,] = rand_ratio

#making sure samples don't go to 0
for (j in 1:S){
if (rand_ratio[j]>1-(S-1)*LOWER_LIMIT){
  rand_ratio[j] = UPPER_LIMIT
}
if (rand_ratio[j]<LOWER_LIMIT){
  rand_ratio[j] = LOWER_LIMIT
}
}
}

nka[,,n] = nk
Ya[,,n] = Y
postmean_a[,,n] = postmean
}
start = Ky[h]+1
}

nk_avg = matrix(nrow=K, ncol=S)
for (j in 1:S){
  nk_avg[,j] = apply(nka[,j,85:100], 1, mean)
}

sum(nk_avg[,1])
sum(nk_avg[,2])

par(mfcol=c(3,1))
plot(nk_avg[,1], ylim = c(0,50), col='red')
points(nk_avg[,2], col='blue')
legend("bottomright",legend = c('no treatment', 'dexamethasone'),col = c("red", "blue"), lwd=4)

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
points(above_theta_avg, col='blue')

#plot(postmean[,1],ylim=range(postmean),type="l",lwd=4,
#     ylab="Posterior mean",xlab="Interim Analysis")
#lines(postmean[,2],type="l",lwd=4,col="red")
#legend("bottomright",c("Treatment 1","Treatment 2"),col=c("black","red"),lwd=4,bty="n")

#lines(mle[,1],ylim=range(mle),type="l",lwd=4,lty=2)
#lines(mle[,2],type="l",lwd=4,col="red",lty=2)
