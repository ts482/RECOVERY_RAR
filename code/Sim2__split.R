rm(list=ls(all=TRUE))

set.seed(20)

reps = 100

K = 100
nk_   = 50 #patients per stage
rand_ratio = 0.5
Nx    = c(501, 1279, 324) + c(1034, 2604, 683)
Px    = Nx/sum(Nx)
Py   = matrix(c(0.140, 0.262, 0.414, 0.178, 0.233, 0.293), nrow=3)
#Py    = c(0.25,0.22) #probability of death on each treatment
S     = length(Py[1,]) #Number of treatments
T     = length(Px) #Number of care
Kx     = round(K*Px) #number of stages in trial
Ky = cumsum(Kx)

alpha       = 1 #uninformative beta
beta        = 1 #uninformative beta
Y           = matrix(nrow=K, ncol=S)
postmean    = matrix(nrow=K, ncol=S)
mle         = matrix(nrow=K, ncol=S)
PostAlpha   = matrix(nrow=K, ncol=S)
PostBeta    = matrix(nrow=K, ncol=S)

nk          = matrix(25, nrow=K,ncol=S)
rand_ratios = vector(length=K)

UPPER_LIMIT <- 0.9
LOWER_LIMIT <- 0.1

nka = array(dim = c(K, S, reps))
Ya = array(dim = c(K, S, reps))
postmean_a = array(dim = c(K, S, reps))


thetas = matrix(nrow=reps, ncol=K)


#y = vector(length=3)

start = 1

# Stage 1
for (h in 1:T){
for (n in start:Ky[h]){
for(j in 1:S){
y = rbinom(1, nk[1,j], Py[h,j])
Y[1,j]         = y
PostAlpha[1,j] = alpha + Y[1,j]
PostBeta[1,j]  = beta  + nk[1,j] -Y[1,j]
postmean[1,j]  = PostAlpha[1,j]/(PostAlpha[1,j]+PostBeta[1,j])
mle[1,j]       = sum(Y[1,j])/(1*nk[1,j])

}
for(i in 2:K){

nk[i, 1] = round(nk_ * (1-rand_ratio))
nk[i, 2] = round(nk_ * rand_ratio)

for(j in 1:S){
y = rbinom(1, nk[i,j], Py[h,j])
Y[i,j]         = y
PostAlpha[i,j] = sum(Y[1:i,j], na.rm=T)
PostBeta[i,j]  = sum(nk[1:i,j]-Y[1:i,j], na.rm=T)
postmean[i,j]  = PostAlpha[i,j]/(PostAlpha[i,j]+PostBeta[i,j])
mle[i,j]       = sum(Y[1:i,j])/sum(nk[1:i,j])


}

X1 <- rbeta(100, PostAlpha[i, 1], PostBeta[i, 1])
X2 <- rbeta(100, PostAlpha[i, 2], PostBeta[i, 2])
theta <- length(X1[X1>=X2])/100
thetas[n,i] <- theta
k = i/(K)
rand_ratio <- theta^k / (theta^k + (1-theta)^k)
rand_ratios[i] = rand_ratio

#making sure samples don't go to 0
if (rand_ratio>UPPER_LIMIT){
  rand_ratio = UPPER_LIMIT
}
if (rand_ratio<LOWER_LIMIT){
  rand_ratio = LOWER_LIMIT
}
}
nka[,,n] = nk
Ya[,,n] = Y
postmean_a[,,n] = postmean
}

start = Ky[h]+1
}

Ky
nk_avg = matrix(nrow=K, ncol=S)
for (j in 1:S){
  nk_avg[,j] = apply(nka[,j,1:100], 1, mean)
}


sum(nk_avg[,1])
sum(nk_avg[,2])

par(mfcol=c(3,1))
plot(nk_avg[,1], ylim = c(0,50), col='red')
points(nk_avg[,2], col='blue')
legend("bottomright",legend = c('no treatment', 'dexamethasone'),col = c("red", "blue"), lwd=4)


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
points(above_theta_avg, col = 'red')

#plot(postmean[,1],ylim=range(postmean),type="l",lwd=4,
#     ylab="Posterior mean",xlab="Interim Analysis")
#lines(postmean[,2],type="l",lwd=4,col="red")
#legend("bottomright",c("Treatment 1","Treatment 2"),col=c("black","red"),lwd=4,bty="n")

#lines(mle[,1],ylim=range(mle),type="l",lwd=4,lty=2)
#lines(mle[,2],type="l",lwd=4,col="red",lty=2)
