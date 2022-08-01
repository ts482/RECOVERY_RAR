set.seed(20)

K = 100
Si = 4
reps = 1000

test_Py = c(0.5,0.5,0.5,0.5)

test.nk_ = 120

test.nk = matrix(nrow=K,ncol=Si)
test.Y = matrix(nrow=K, ncol=Si)

test.nka = array(dim = c(K, Si, reps))
test.Ya = array(dim = c(K, Si, reps))

test_t_metric = array(dim = c(Si-1, K/10, reps))


for (n in 1:reps){
  for (i in 1:K){
   for (j in 1:Si){
     test.nk[i,j] = test.nk_/Si
     test.Y[i,j] = rbinom(1,test.nk[i,j], test_Py[j])
   }
    if (i %% 10 == 0){
      pat_df <- data.frame(t=rep(1, sum(test.nk[1:i, c(1,2)])), 
                           o=rep(0, sum(test.nk[1:i, c(1,2)])))
      pat_df[1:sum(test.nk[1:i,1]), 't'] = 0
      pat_df[1:(sum(test.Y[1:i,1])), 'o'] = 1
      pat_df[(sum(test.nk[1:i,1])+1):(sum(test.nk[1:i,1])+sum(test.Y[1:i,2])+1), 'o'] = 1
      
      test_t_metric[1,i/10,n] = t.test(o~t, data = pat_df)$p.value
      
      
      #calculating power for hydroxy
      pat_df <- data.frame(t=rep(1, sum(test.nk[1:i, c(1,3)])), 
                           o=rep(0, sum(test.nk[1:i, c(1,3)])))
      pat_df[1:sum(test.nk[1:i,1]), 't'] = 0
      pat_df[1:(sum(test.Y[1:i,1])), 'o'] = 1
      pat_df[(sum(test.nk[1:i,1])+1):(sum(test.nk[1:i,1])+sum(test.Y[1:i,3])+1), 'o'] = 1
      
      test_t_metric[2,i/10,n] = t.test(o~t, data = pat_df)$p.value
      
      
      #calculating power for lopinavir
      pat_df <- data.frame(t=rep(1, sum(test.nk[1:i, c(1,4)])), 
                           o=rep(0, sum(test.nk[1:i, c(1,4)])))
      pat_df[1:sum(test.nk[1:i,1]), 't'] = 0
      pat_df[1:(sum(test.Y[1:i,1])), 'o'] = 1
      pat_df[(sum(test.nk[1:i,1])+1):(sum(test.nk[1:i,1])+sum(test.Y[1:i,4])+1), 'o'] = 1
      
      test_t_metric[3,i/10,n] = t.test(o~t, data = pat_df)$p.value
      
    }
  }
  test.nka[,,n] = test.nk
  test.Ya[,,n] = test.Y
}


test_t1e = matrix(nrow=10, ncol=3)
test_fwer = vector(length=10)
for (i in (1:10)){
  test_t1e[i,] = apply(test_t_metric[,i,]<0.05, 1, mean)
  test_fwer[i] = mean(apply(test_t_metric[,i,]<0.05/3, 2, any))
}

test_fwer
