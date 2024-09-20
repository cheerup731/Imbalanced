
#library(sampling)
#library(themis)

library(here)

set.seed(14)
M = 300

### Load the population file
pop <- read.csv(file=here("synthdat.csv"),
                header=TRUE,
                sep=",")
pop$state <- as.factor(pop$state)
pop$indgrp <- as.factor(pop$indgrp)
pop = subset(pop,select=-c(id,stratum2))

N = nrow(pop)

scenario = expand.grid(P = c(0.001, 0.005),
                       NP= c(0.01, 0.03))
Ind = c()

for (i in 1:nrow(scenario)) {
  if (scenario[i,1]<scenario[i,2]) {
    Ind = c(Ind,i)
  }
} 


scen.fun <- function(i){
  
  n.P = scenario[i, 1]*N ; n.NP = scenario[i, 2]*N
  bias = mse = bias2 = mse2 = 0

for (m in 1:M) {
  
  # srs probability sample 
  indx = sample.int(N, size = n.P, replace = F)
  S = rep(F,N)
  S[indx] = T
  sample.P = pop[S, ]
  
  # nonprobability sample 
  includeNP = pop$pi.mar * n.NP/sum(pop$pi.mar)
  Sstar = as.logical(UPrandomsystematic(includeNP))
  sample.NP = pop[Sstar,]
  W.NP = rep(N/n.P, sum(Sstar))
  
  # B set
  B = pop[Sstar+S == 1, ]
  B$Z = Sstar[Sstar+S == 1]
  
  # estimate O
  glmO = glm(Z ~ size2  + indgrp, data = B, family=quasibinomial(link="logit"))
  O = exp(predict(glmO, newdata = sample.NP))
  
  f18 = 1+(W.NP-1)/O
  
  #### classic SMOTE ####
  B$Z <- as.factor(B$Z)
  OS_B = smotenc(B, var="Z") 
  OS_B$Z= as.numeric(OS_B$Z)-1
  glmO = glm(Z ~ size2  + indgrp, data = OS_B, family=quasibinomial(link="logit"))
  O = exp(predict(glmO, newdata = sample.NP))
  f = mean(OS_B$Z)
  B$Z <- as.numeric(B$Z)-1
  p = sum(1-B$Z)*f/(sum(B$Z)*(1-f))
  f18_OS = 1+(W.NP-1)*p/O
  
  #### under sampling ####
  max.class = names(which.max(table(B$Z)))
  B_maj = B[B$Z == max.class,]
  US_B = rbind(B_maj[sample(max(table(B$Z)), size = min(table(B$Z))),], B[B$Z != max.class,])
  glmO = glm(Z ~ size2  + indgrp, data = US_B, family=quasibinomial(link="logit"))
  O = exp(predict(glmO, newdata = sample.NP))
  f = mean(US_B$Z)
  p = sum(1-B$Z)*f/(sum(B$Z)*(1-f))
  f18_US = 1+(W.NP-1)*p/O
  
  #### mix ####
  US_B2 = rbind(B_maj[sample(max(table(B$Z)), size = 0.5*nrow(B)),], B[B$Z != max.class,])
  US_B2$Z = as.factor(US_B2$Z)
  OS_B = smotenc(US_B2, var="Z") 
  OS_B$Z = as.numeric(OS_B$Z)-1
  glmO = glm(Z ~ size2  + indgrp, data = OS_B, family=quasibinomial(link="logit"))
  O = exp(predict(glmO, newdata = sample.NP))
  f = mean(OS_B$Z)
  p = sum(1-B$Z)*f/(sum(B$Z)*(1-f))
  f18_mix = 1+(W.NP-1)*p/O
  
  
  #### estimates ####
  
  est = c(mean(sample.NP$ovt), 
          sum(sample.NP$ovt*f18)/sum(f18), 
          sum(sample.NP$ovt*f18_OS)/sum(f18_OS),
          sum(sample.NP$ovt*f18_US)/sum(f18_US),
          sum(sample.NP$ovt*f18_mix)/sum(f18_mix))
  
  bias =(est - mean(pop$ovt))/M + bias
  mse = (est - mean(pop$ovt))^2/M + mse
  
  est2 = c(mean(sample.NP$earnings),
           sum(sample.NP$earnings*f18)/sum(f18), 
           sum(sample.NP$earnings*f18_OS)/sum(f18_OS),
           sum(sample.NP$earnings*f18_US)/sum(f18_US),
           sum(sample.NP$earnings*f18_mix)/sum(f18_mix))
  
  bias2 =(est2 - mean(pop$earnings))/M + bias2
  mse2 = (est2 - mean(pop$earnings))^2/M + mse2
  
}
  return(list(n.P/N, n.NP/N, mse, bias, mse2, bias2))
}
  



# parallel

library(doRNG); library(doParallel)

registerDoParallel(cores = (detectCores()-1))
registerDoRNG(13)

result <- foreach(i = Ind,
                  .packages = c('sampling', 'themis'),
                  .combine = rbind) %dopar% scen.fun(i)

table = matrix(data = unlist(result[,c(1,2)]), nrow = length(Ind))
table = cbind(table, matrix(data = unlist(result[,3]), nrow = length(Ind), byrow = T))
table = cbind(table, matrix(data = unlist(result[,4]), nrow = length(Ind), byrow = T))
table = cbind(table, matrix(data = unlist(result[,5]), nrow = length(Ind), byrow = T))
table = cbind(table, matrix(data = unlist(result[,6]), nrow = length(Ind), byrow = T))

colnames(table) = c('nP', 'nNP',
                    'mse_naive','mse_LSdW','mse_smote','mse_us','mse_mix',
                    'bias_naive','bias_LSdW','bias_smote','bias_us','bias_mix',
                    'mse_naive','mse_LSdW','mse_smote','mse_us','mse_mix',
                    'bias_naive','bias_LSdW','bias_smote','bias_us','bias_mix')

#table

save.image("app_srs.RData")


# RB
RBtable = table[,c(1:22)]
colnames(RBtable) = c('nP', 'nNP',
                    'mse_np','mse_LSdW','mse_smote','mse_us','mse_mix',
                    'RB_np','RB_LSdW','RB_smote','RB_us', 'RB_mix',
                    'mse_np','mse_LSdW','mse_smote','mse_us','mse_mix',
                    'RB_np','RB_LSdW','RB_smote','RB_us', 'RB_mix')
RBtable[,3] = 1e-2*table[,3]
RBtable[,4] = 1e-2*table[,4]
RBtable[,5] = 1e-2*table[,5]
RBtable[,6] = 1e-2*table[,6]
RBtable[,7] = 1e-2*table[,7]

RBtable[,8] = 100*table[,8]/mean(pop$ovt)
RBtable[,9] = 100*table[,9]/mean(pop$ovt)
RBtable[,10] = 100*table[,10]/mean(pop$ovt)
RBtable[,11] = 100*table[,11]/mean(pop$ovt)
RBtable[,12] = 100*table[,12]/mean(pop$ovt)

RBtable[,13] = 1e-5*table[,13]
RBtable[,14] = 1e-5*table[,14]
RBtable[,15] = 1e-5*table[,15]
RBtable[,16] =1e-5*table[,16]
RBtable[,17] = 1e-5*table[,17]

RBtable[,18] = 100*table[,18]/mean(pop$earnings)
RBtable[,19] = 100*table[,19]/mean(pop$earnings)
RBtable[,20] = 100*table[,20]/mean(pop$earnings)
RBtable[,21] = 100*table[,21]/mean(pop$earnings)
RBtable[,22] = 100*table[,22]/mean(pop$earnings)

round(RBtable,2)

capture.output(round(RBtable,3), file = 'app_srs.txt')


library("xtable")
RBtable = as.data.frame(RBtable)
RBtable[,c(1,2)] = round(RBtable[,c(1,2)]*100, 1)
Latex = xtable(RBtable[, c(1:7)], digits = c(0,1,0, c(rep(3,5))), caption = 'Ovt mse estimates under srs')
Latex

Latex = xtable(RBtable[, c(1,2, 13:17)], digits = c(0,1,0, c(rep(3,5))), caption = 'earnings mse estimates under srs')
Latex

Latex = xtable(RBtable[, c(1,2, 8:12)], digits = c(0,1,0, c(rep(3,5))), caption = 'Ovt RB estimates under srs')
Latex

Latex = xtable(RBtable[, c(1,2, 18:22)], digits = c(0,1,0, c(rep(3,5))), caption = 'earnings RB estimates under srs')
Latex