par(mfrow=c(1,1))
barplot(c(1596, 1627, 1564, 1558, 1516, 1519), 
        names.arg=c('1:1 fixed', '2:1 fixed', 
                    'RMC cohort', 'Tuning cohort', 'RMC sub', 'tuning sub'), ylim=c(0,2000))
abline(h=c(1489,1671), lty='dashed')

mat = c(4844, 4804,3200, 2133)
matsd = c(11.1, 10, 4.0,3.8)
plot(1:length(mat), mat, ylab='number of patients given dexamethasone', xlab="randomisation algorithm",
     ylim=c(0,5000), pch=16, cex=1.5, xaxt='n')
axis(1, labels = c('RAR1', 'RAR2', '1:1 ratio', '2:1 ratio'), at = 1:length(mat))
arrows(1:length(mat), mat-matsd, 1:length(mat), mat+matsd, code=3, angle = 90, lwd=0.01)

#barplot(mat, horiz=T, xlab='number of patients', names.arg= c('Tuning', 'REMAPCAP','2:1 ratio', '1:1 ratio'),
#        legend.text = c('standard care', 'dexamethasone'))





#### barplots representing subgroup and full cohort algorithms
Nx.gr = c(501, 1279, 324) + c(1034, 2604, 683)
Px    = Nx.gr/sum(Nx.gr)
Px = c(Px, 1)

mat = matrix(data = c(937, 584, 1482, 2356, 335, 704, 2375, 4025), nrow=2, ncol=4)
barplot(mat, xlab='Randomisation algorithm and sub-group', ylab='number of patients', 
        names.arg= c(rep(expression('T'[s]), 3), expression('T'[f])),
        legend.text = c('standard care', 'dexamethasone'),  ylim=c(0,10000), width = c(0.07, 0.07, 0.07, 0.07)
        , space=c(3.5, 5, 5, 5), xlim = c(0,1.7), args.legend = list(x = 'topleft'))



#mat_sd = c(nk_sd.g[1,],nk_sd.c[1])
#arrows(c(0.275,0.7,1.125,1.55), mat[1,]-mat_sd, c(0.275,0.7,1.125,1.55), mat[1,]+mat_sd, code=3, angle=90, lty=8)


mat = matrix(data = c(923, 597, 1469, 2367, 307, 732, 2404, 3996), nrow=2, ncol=4)
barplot(mat, add=T,  names.arg= c(rep(expression('RMC'[s]), 3), expression('RMC'[f])),
          ylim=c(0,4500), width = c(0.07, 0.07, 0.07, 0.07),
         space=c(5, 5, 5,5), xlim = c(0,1.5))

#mat_sd = c(nk_sd.h[1,],nk_sd.d[1])
#arrows(0.1+c(0.275,0.7,1.125,1.55), mat[1,]-mat_sd, 0.1+c(0.275,0.7,1.125,1.55), 
#       mat[1,]+mat_sd, code=3, angle=90, lty=8)

for (p in 1:length(Px)){
  mat[,p] = Px[p]*c(2/3, 1/3)*6400
}

barplot(mat, add=T,  names.arg= rep('FuR', 4),
         ylim=c(0,4500), width = c(0.07, 0.07, 0.07, 0.07),
        space=c(2, 5, 5, 5), xlim = c(0,1.5))

for (p in 1:length(Px)){
  mat[,p] = Px[p]*c(1/2, 1/2)*6400
}

barplot(mat, add=T,  names.arg= rep('FeR', 4),
        ylim=c(0,4500), width = c(0.07, 0.07, 0.07, 0.07),
        space=c(0.5, 5, 5, 5), xlim = c(0,1.5))

text(x=0.2, y=3600, labels='subgroup i')
text(x=0.62, y=5500, labels='subgroup ii')
text(x=1.05, y=3000, labels='subgroup iii')
text(x= 1.5, y = 8000, labels = 'Full cohort')

plot(rep(0.52,100), xlab='trial progress', ylab='randomisation probability', col='red', 
     lwd=5, ylim=c(0,1))
points(rep(0.48, 100), col=c('blue'), lwd=5)
legend('bottomright', legend= c('dexamethasone', 'standard care'), fill=c('blue', 'red'))






###barplot for 4 arm trial
mat=  matrix(apply(nk_avg.l[1:80,], 2,sum), nrow=4,ncol=1)
barplot(mat, ylim=c(0,10000), xlim=c(0,8), 
        legend.text= c('standard care', 'dexamethasone', 'hydroxychloroquine','lopinavir'),
        names.arg = c('FeR'), col = c('red','blue','yellow','green'),
        xlab='Randomisation algorithm',ylab='Number of patients allocated to each arm')


mat=  matrix(apply(nk_avg.k[1:80,], 2,sum), nrow=4,ncol=1)
barplot(mat, col = c('red','blue','yellow','green'),
        add=T, space=c(1.5), names.arg = c('FuR'))

mat=  matrix(apply(nk_avg.i[1:80,], 2,sum), nrow=4,ncol=1)
barplot(mat, ylim=c(0,10000), xlim=c(0,6), col = c('red','blue','yellow','green'),
        add=T, space=c(3), names.arg = c(expression('T'[f])))

mat=  matrix(apply(nk_avg.m[1:80,], 2,sum), nrow=4,ncol=1)
barplot(mat, ylim=c(0,10000), xlim=c(0,6), col = c('red','blue','yellow','green'),
        add=T, space=c(4.5), names.arg = c(expression('RMC'[f])))















##### legacy code simulation study
#simulation study
set.seed(20)

#number randomised to each group during recovery

#order: no oxygen, oxygen, invasive mechanical ventilation
stand.npat <- c(1034, 2604, 683)
#stand.npat <- 4321
dex.npat <- c(501, 1279, 324)
#dex.npat <- 2104

#total number of patients in trial
pat_nums = c(stand.npat, dex.npat)
N_PATIENTS = sum(pat_nums)
n_pats.none <- sum(c(stand.npat[1], dex.npat[1]))
n_pats.ox <- sum(c(stand.npat[2], dex.npat[2]))
n_pats.vent <- sum(c(stand.npat[3], dex.npat[3]))

#ratios of randomisation
orig.ratios <- pat_nums[1:length(pat_nums)]/N_PATIENTS

#checking randomisation ratios
(stand.npat[1:3] - dex.npat[1:3])/stand.npat[1:3]
#all round up to 51/52% difference, therefore different ratios will not be accounted
treat.ratios <- c(sum(stand.npat), sum(dex.npat))/sum(c(stand.npat, dex.npat))

#ratio of care
care.ratios <- c(n_pats.none, n_pats.ox, n_pats.vent)/N_PATIENTS

#total probability of death in each arm
stand.pdeath <- c(0.140, 0.262, 0.414)
#stand.pdeath <- 0.257
dex.pdeath <- c(0.178, 0.233, 0.293)
#dex.pdeath <- 0.229

# odds ratio varies (0.91-1.55, 0.72-0.94, 0.51-0.81)
# potentially useful for best/worst case scenarios analysis

Py = matrix(nrow=3, ncol=2)
Py[,1] <- stand.pdeath
Py[,2] <- dex.pdeath
#Py = c(stand.pdeath, dex.pdeath)
S = dim(Py)[2]
T = dim(Py)[1]

## RECOVERY TRIAL SIMULATION

#initialising outcome vectors
rec.matrix <- matrix(nrow=N_PATIENTS, ncol = 3)


#looping through number of patients
for (i in 1:N_PATIENTS) {
  care = sample(1:T, size=1, prob = care.ratios) #1,2,3 for no treatment, oxygen, ventilation
  treatment = sample(1:S, size=1, prob= treat.ratios) #1 is standard, 2 is treatment
  outcome = sample(c(1, 0), size=1, prob = c(Py[care, treatment], 1-Py[care, treatment])) #1 is death, 0 is survive
  rec.matrix[i, 1:3] = c(care, treatment, outcome)
}

#outcome death metrics per treatment
mean(rec.matrix[rec.matrix[,2] == 1, 3])
mean(rec.matrix[rec.matrix[,2] == 2, 3])

#outcome death metrics per care per treatment
for (i in 1:T){
  for (j in 1:S){
    print(mean(rec.matrix[(rec.matrix[,1] == i & rec.matrix[,2] == j), 3]))
  }
}


## ADAPTIVE DESIGN SIMULATION

#intialising outcome vectors
ad.matrix <- matrix(nrow = N_PATIENTS, ncol = 3)

#number of patients required per group before 
#looking at adaptive randomisation
N_INIT = 30

#number of patients per randomisation ratio update
UPDATE_PER = 100

K = ceiling(N_PATIENTS / UPDATE_PER) #number of stages needed

alpha       = 1 #uninformative beta
beta        = 1 #uninformative beta
Y           = array(dim=c(K, S, T))
postmean    = array(dim=c(K, S, T))
mle         = array(dim=c(K, S, T))
PostAlpha   = array(dim=c(K, S, T))
PostBeta    = array(dim=c(K, S, T))
nk          = array(dim=c(K, S, T))
rand_ratios = matrix(nrow=K, ncol=T)

#tracking number of stages to update matrix
i = as.vector(rep(0,3))

#initial randomisation ratio
rand_ratio <- as.vector(rep(0.5,3))

#keeping track of last x value for filtering matrix
last_x = as.vector(rep(0,3))

#intialising vector to check sample size against threshold
threshold <- array(rep(F, S), dim=c(T, S))

#setting limits on randomisation ratios
UPPER_LIMIT = 0.9
LOWER_LIMIT = 0.1

X <- matrix(nrow=100, ncol=S)
#looping through same number of patients
for (x in 1:N_PATIENTS){
  
  care = sample(1:T, size=1, prob = care.ratios) #1,2,3 for no treatment, oxygen, ventilation
  for (h in 1:T){
    for (j in 1:S){
      threshold[h, j] <- sum((ad.matrix[,1] == h & ad.matrix[,2] == j), na.rm=T)>=N_INIT
    }
    #checking that both patient sets have passed initial threshold
    if (all(threshold[h,])){
      
      # updating randomisation every update_per times
      if (x%%UPDATE_PER==0){
        
        #record number of randomisation updates (but only once)
        i[h] <- i[h]+1
        
        #updating section of patients accounted for
        sub.matrix <- ad.matrix[last_x[h]:x,]
        
        #updating last x for next round
        last_x[h] <- x
        
        
        for (j in 1:S){
          #store number of patients randomised to this treatment
          nk[i[h], j, h] = sum(sub.matrix[,1] == h & sub.matrix[,2]==j, na.rm=T)
          
          #checking at least 1 patient randomised to this treatment
          if (nk[i[h],j,h]!=0){
          
            #update outcome in most recent stage
            Y[i[h],j, h] = sum(sub.matrix[(sub.matrix[,1] == h & sub.matrix[,2] == j), 3], na.rm=T)
            
          
            #step required in first randomisation
            if (i[h]<2){
              PostAlpha[1,j,h] = alpha + Y[1,j,h]
              PostBeta[1,j,h]  = beta  + nk[1, j,h]-Y[1,j,h]
              postmean[1,j,h]  = PostAlpha[1,j,h]/(PostAlpha[1,j,h]+PostBeta[1,j,h])
              mle[1,j,h]       = sum(Y[1,j,h])/(1*nk[1,j,h])
            } else {
              PostAlpha[i,j,h] = sum(Y[1:i[h],j,h], na.rm=T)
              PostBeta[i,j,h]  = sum(nk[1:i[h],j,h]-Y[1:i[h],j,h], na.rm=T)
              postmean[i,j,h]  = PostAlpha[i[h],j,h]/(PostAlpha[i[h],j,h]+PostBeta[i,j,h])
              mle[i,j,h]       = sum(Y[1:i[h],j,h])/(sum(nk[1:i[h],j,h], na.rm=T))
            }
          }
          X[,j] = rbeta(100, PostAlpha[i[h],j,h], PostBeta[i[h],j,h])
        }
        theta <- length(X[,1][X[,1]>=X[,2]])/100
        k = i[h]/(2*K)
        rand_ratio[h] <- theta^k / (theta^k + (1-theta)^k)
        rand_ratios[i,] = rand_ratio[h]
        
        #making sure samples don't go to 0
        if (rand_ratio[h]>UPPER_LIMIT){
          rand_ratio[h] = UPPER_LIMIT
        }
        if (rand_ratio[h]<LOWER_LIMIT){
          rand_ratio[h] = LOWER_LIMIT
        }
      }
    }
  }
  if (all(threshold[care,])){
    #rand_ratio <- 0.5    #use this for testing
    treatment = sample(c(1:S), size=1, prob=c(1-rand_ratio[care], rand_ratio[care]))
  } else {
    treatment = sample(c(1:S),size=1)
  }
  outcome = sample(c(1, 0), size=1, prob = c(Py[care, treatment], 1-Py[care, treatment])) #1 is death, 0 is survive
  ad.matrix[x, 1:3] = c(care, treatment, outcome)
}

#number of patients and death rte in each treatment
for (n in 1:S){
  print(sum(ad.matrix[,2]==n))
  print(mean(ad.matrix[ad.matrix[,2]==n,3]))
}

#number of patients and death ratein each care
for (n in 1:T){
  print(sum(ad.matrix[,1]==n))
  print(mean(ad.matrix[ad.matrix[,1]==n, 3]))
}

#sub-group analysis

for(j in 1:S){
  for (h in 1:T){
    print(sum(ad.matrix[,1]==h & ad.matrix[,2]==j))
    print(mean(ad.matrix[ad.matrix[,1]==h&ad.matrix[,2]==j, 3]))
  }
}

tail(rand_ratios)
#plotting
#plot(postmean[,1],type="l",ylim=range(postmean, na.rm=T),lwd=4,
#     ylab="Posterior mean",xlab="Interim Analysis")
#lines(postmean[,2],type="l",lwd=4,col="red")
#legend("bottomright",c("Treatment 1","Treatment 2"),col=c("black","red"),lwd=4,bty="n")
