
library(sampling); library(FNN); library(smotefamily)
library(MASS)

syspps=function(N_f,x){ 
  U=sample(N_f,N_f)
  xx=x[U]
  cx=cumsum(xx)
  r=runif(1)
  s=numeric()
  for(i in 1:N_f){
    if(cx[i]>=r){
      s=c(s,U[i])
      r=r+1
    }
  }
  return(s)
}

set.seed(14)

N = 3e4; M = 1e3
r1 = 0.7; r2 = 0.3
sig = matrix(c(  1, r1, r1, r1,
                r1,  1, r2, r2,
                r1, r2,  1, r2,
                r1, r2, r2,  1), 4, 4)

pop = as.data.frame(mvrnorm(n = N, mu = c(1,1,1,1), Sigma = sig))
colnames(pop) = c('y', 'x1', 'x2', 'x3')

scenario = expand.grid(P = c(0.01, 0.05, 0.1, 0.3, 0.5),
                       NP= c(0.1, 0.3, 0.5, 0.7))
Ind = c()

for (i in 1:nrow(scenario)) {
  if (scenario[i,1]<scenario[i,2]) {
    Ind = c(Ind,i)
  }
} # only look at NP>P

#i = m = 7

scen.fun <- function(i){
  
  n.P = scenario[i, 1]*N ; n.NP = scenario[i, 2]*N
  bias = mse = convergence = 0

for (m in 1:M) {
  
  # pps probability sample
  z = pop$x1
  c = (max(z)-50*min(z))/49
  s = sum(z+c)
  w = s/(z+c)/n.P
  z1 = (z+c)/s
  indx = syspps(N,z1*n.P)
  S = rep(F,N)
  S[indx] = T
  sample.P = pop[S, -1]
  W.P = w[S]
  
  # nonprobability sample 
  x = -pop$y*2
  f = function(theta) sum(exp(theta + x) / (1 + exp(theta + x))) - n.NP
  theta = uniroot(f, c(-100, 100))$root  
  includNP = exp(theta + x) / (1 + exp(theta + x))
  Sstar = as.logical(UPrandomsystematic(includNP))
  sample.NP = pop[Sstar,]
  W.NP = w[Sstar]
  
  # B set
  B = pop[Sstar+S == 1, -1]
  B$Z = Sstar[Sstar+S == 1]
  glmO = glm(Z ~ x1+x2+x3, data = B, family = "binomial")
  O = exp(predict(glmO, newdata = sample.NP))
  f18 = 1+(W.NP-1)/O
  
  #### classic SMOTE ####
  smote = SMOTE(B, B$Z) 
  OS_B = smote$data
  glmO = glm(Z ~ x1+x2+x3, data = OS_B, family = "binomial")
  O = exp(predict(glmO, newdata = sample.NP))
  f = mean(OS_B$Z)
  p = sum(1-B$Z)*f/(sum(B$Z)*(1-f))
  f18_OS = 1+(W.NP-1)*p/O
  
  
  #### srs smote ####
  B0 = S&!Sstar
  sample.P0 = pop[B0,-1]
  W.P0 = w[B0]
  w_reg = W.P0*sum(B$Z)/sum(W.P0) # regularized to have the same size with NP
  w_nn = trunc(w_reg) + UPrandomsystematic(w_reg - trunc(w_reg))-1 
  sub = which(w_nn>0)
  nn = knn.index(sample.P0, k = 5, algorithm = "kd_tree")
  
  syn_data = NULL
  
  for (i in sub) {
    nei = nn[i, sample.int(5, w_nn[i], replace = T)]
    P = matrix(unlist(sample.P0[i,   ]), w_nn[i], ncol(sample.P0), byrow = TRUE)
    Q = as.matrix(sample.P0[nei, ])
    new_data = P + runif(w_nn[i], min = 0, max = 1)*(Q - P)
    syn_data = rbind(new_data, syn_data)
  }
  
  syn_data = as.data.frame(syn_data)
  syn_data = rbind(syn_data, sample.P0)
  
  OS_B = rbind(syn_data, B[B$Z,-4]) 
  OS_B$Z = c(rep(0, nrow(syn_data)), rep(1, sum(B$Z)))
  
  glmO = glm(Z ~ x1+x2+x3, data = OS_B, family = "binomial")
  O = exp(predict(glmO, newdata = sample.NP))
  f = mean(OS_B$Z)
  cb0 = sum(B$Z)*(1-f)/f/sum(W.P0)
  f18_OS2 = 1+(W.NP-1)/(O*W.NP*cb0)
  
  #### under sampling ####
  max.class = names(which.max(table(B$Z)))
  B_maj = B[B$Z == max.class,]
  US_B = rbind(B_maj[sample(max(table(B$Z)), size = min(table(B$Z))),], B[B$Z != max.class,])
  glmO = glm(Z ~ x1+x2+x3, data = US_B, family = "binomial")
  O = exp(predict(glmO, newdata = sample.NP))
  f = mean(US_B$Z)
  p = sum(1-B$Z)*f/(sum(B$Z)*(1-f))
  f18_US = 1+(W.NP-1)*p/O
  
  #### mix ####
  US_B2 = rbind(B_maj[sample(max(table(B$Z)), size = 0.5*nrow(B)),], B[B$Z != max.class,])
  smote = SMOTE(US_B2, US_B2$Z) 
  OS_B = smote$data
  glmO = glm(Z ~ x1+x2+x3, data = OS_B, family = "binomial")
  O = exp(predict(glmO, newdata = sample.NP))
  f = mean(OS_B$Z)
  p = sum(1-B$Z)*f/(sum(B$Z)*(1-f))
  f18_mix = 1+(W.NP-1)*p/O
  
  
  
  est = c(mean(sample.NP$y), 
          sum(sample.NP$y*f18)/sum(f18), 
          sum(sample.NP$y*f18_OS)/sum(f18_OS),
          sum(sample.NP$y*f18_OS2)/sum(f18_OS2),
          sum(sample.NP$y*f18_US)/sum(f18_US),
          sum(sample.NP$y*f18_mix)/sum(f18_mix))
         
  bias =(est - mean(pop$y))/M + bias
  mse = (est - mean(pop$y))^2/M + mse
  
  convergence[m] = est[3]
 
}
  
  return(list(n.P/N, n.NP/N, mse, bias, convergence))
}

# parallel 

library(doRNG); library(doParallel)

registerDoParallel(cores = min(c(length(Ind), (detectCores()-1))))
registerDoRNG(13)

result <- foreach(i = Ind,
                  .packages = c('sampling', 'FNN', 'smotefamily', 'MASS'),
                  .combine = rbind) %dopar% scen.fun(i)

table = matrix(data = unlist(result[,c(1,2)]), nrow = length(Ind)) 
table = cbind(table, matrix(data = unlist(result[,3]), nrow = length(Ind), byrow = T)) 
table = cbind(table, matrix(data = unlist(result[,4]), nrow = length(Ind), byrow = T)) 

colnames(table) = c('nP', 'nNP', 
                    'mse_naive','mse_LSdW','mse_smote','mse_smote2','mse_us','mse_mix',
                    'bias_naive','bias_LSdW','bias_smote','bias_smote2','bias_us','bias_mix')

table

save.image("pps.RData")

# convergence
# plot(cumsum(result[,5]$result.1)/c(1:M))
# plot(cumsum(result[,5]$result.3)/c(1:M))
# plot(cumsum(result[,5]$result.6)/c(1:M))

# RB
RBtable = table
colnames(RBtable) = c('nP', 'nNP', 
                      '100*mse_np','100*mse_LSdW','100*mse_smote','100*mse_smote2','100*mse_us','100*mse_mix',
                      'RB_np','RB_LSdW','RB_smote','RB_smote2','RB_us','RB_mix')
RBtable[,3] = 100*table[,3]
RBtable[,4] = 100*table[,4]
RBtable[,5] = 100*table[,5]
RBtable[,6] = 100*table[,6]
RBtable[,7] = 100*table[,7]
RBtable[,8] = 100*table[,8]

RBtable[,9] = 100*table[,9]/mean(pop$y)
RBtable[,10] = 100*table[,10]/mean(pop$y)
RBtable[,11] = 100*table[,11]/mean(pop$y)
RBtable[,12] = 100*table[,12]/mean(pop$y)
RBtable[,13] = 100*table[,13]/mean(pop$y)
RBtable[,14] = 100*table[,14]/mean(pop$y)
round(RBtable,3)

capture.output(round(RBtable,3), file = 'pps.txt')

library("xtable")
RBtable = as.data.frame(RBtable)
Latex = xtable(RBtable, digits = c(0,c(rep(2,14))), caption = 'Estimates under pps')
Latex
 
