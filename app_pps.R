
library(sampling);library(FNN);library(themis)
  
library(here)

set.seed(14)

M = 300

### Load the population file
pop <- read.csv(file=here("synthdat.csv"),
                header=TRUE,
                sep=",")

### Calculate the number of observations in each stratum2 category
stratum2 <- as.data.frame(table(pop$stratum2, dnn=list("stratum2")),responseName="Nh")

### Merge the stratum2 counts with the population file
pop <- merge(x=pop,
             y=stratum2,
             by="stratum2",
             all.x=TRUE)

### Load the allocation file
alloc <- read.csv(file=here("ref alloc.csv"),
                  header=TRUE,
                  sep=",")

### Add on the stratum allocation to population file
pop <- merge(x=pop,
             y=alloc,
             by="stratum2",
             all.x=TRUE)

pop$deswgt <- pop$Nh/pop$sampalloc

### Place sample allocation into a vector
alloc <- alloc[order(alloc$stratum2),]
sampalloc <- alloc$sampalloc

pop$state <- as.factor(pop$state)
pop$indgrp <- as.factor(pop$indgrp)

N = nrow(pop)

scenario = c(0.03, 0.05, 0.07, 0.10, 0.30)


scen.fun <- function(i){
  
  n.NP = scenario[i]*N
  bias = mse = bias2 = mse2 = 0

for (m in 1:M) {
  
  # Draw a Stratified SRS probability sample 
  indx = sampling::strata(data=pop,
                          stratanames="stratum2",
                          size=sampalloc,
                          method="srswor")
  
  S = rep(F,N)
  S[indx$ID_unit] = T
  sample.P = pop[S, ]
  n.P <- nrow(sample.P) # P sample size
  
  # nonprobability sample 
  includeNP = pop$pi.mar * n.NP/sum(pop$pi.mar)
  Sstar = as.logical(UPrandomsystematic(includeNP))
  sample.NP = pop[Sstar, ]
  
  ### Set delta.i = 1 if observation in sample.P is also in sample.NP
  sample.P$delta.i <- sample.P$id %in% sample.NP$id
  
  W.NP = pop[Sstar,"deswgt"]  # Design weights for P sample, for the NP units
  
  # B set includes units in P or NP sample, not both
  # Separate out the CE and sampled sector units
  B = pop[Sstar+S == 1, -1]
  B$Z = Sstar[Sstar+S == 1]
  
  B.samp <- B[B$deswgt > 1,]
  B.ce <- B[B$deswgt == 1,]
  
  # estimate O
  glmO = glm(Z ~ size2  + indgrp + state , data = B.samp, family=quasibinomial(link="logit"))
  
  ### Separate out NP sample into the CE and sampled sectors
  sample.NP.samp <- sample.NP[sample.NP$deswgt > 1,]
  sample.NP.ce <- sample.NP[sample.NP$deswgt == 1,]
  W.NP.samp <- sample.NP.samp$deswgt
  
  ### Apply model to the NP sample that have weights > 1
  O = exp(predict(glmO, newdata = sample.NP.samp))
  f18 = 1+(W.NP.samp-1)/O
  
  
  #### classic SMOTE - on the B.samp portion only ####
  B.samp$Z <- as.factor(B.samp$Z)
  B.samp <- subset(B.samp,select=-c(id,Nh,sampalloc,pi.mar,ovt))
  OS_B = smotenc(B.samp,var="Z")
  OS_B$Z <- as.numeric(OS_B$Z)-1

  B.samp$Z <- as.numeric(B.samp$Z)-1
  glmO = glm(Z ~ size2  + indgrp + state, data = OS_B, family=quasibinomial(link="logit"))
  O = exp(predict(glmO, newdata = sample.NP.samp))
  f = mean(OS_B$Z)
  p = sum(1-B.samp$Z)*f/(sum(B.samp$Z)*(1-f))
  f18_OS = 1+(W.NP.samp-1)*p/O
 
  
  
  #### srs smote ####
  B0 <- B.samp$Z == 0
  sample.P0 = B.samp[B0, ]
  W.P0 = B.samp$deswgt[B0]
  w_reg = W.P0*sum(B.samp$Z)/sum(W.P0) # regularized to have the same size with NP
  w_nn = trunc(w_reg) + UPrandomsystematic(w_reg - trunc(w_reg))-1 
  sub = which(w_nn>0)
  nn = knn.index(sample.P0[,-c(1:3)], k = 5, algorithm = "kd_tree")
  
  syn_data = NULL
  
  for (i in sub) {
    nei = nn[i, sample.int(5, w_nn[i], replace = T)]
    P = matrix(unlist(sample.P0[i, -c(1:3)]), w_nn[i], ncol(sample.P0)-3, byrow = TRUE)
    Q = as.matrix(sample.P0[nei, -c(1:3)])
    sq = P + runif(w_nn[i], min = 0, max = 1)*(Q - P)
    Q_cat = sample.P0[nei, c(1:3)]
    new_data = cbind(Q_cat, sq)
    syn_data = rbind(new_data, syn_data)
  }
  
  syn_data = rbind(syn_data, sample.P0)
  OS_B = rbind(syn_data, B.samp[B.samp$Z == 1,]) 
  glmO = glm(Z ~ size2  + indgrp + state, data = OS_B, family=quasibinomial(link="logit"))
  O = exp(predict(glmO, newdata = sample.NP.samp))
  f = mean(OS_B$Z)
  cb0 = sum(B.samp$Z)*(1-f)/f/sum(W.P0)
  f18_OS2 = 1+(W.NP.samp-1)/(O*W.NP.samp*cb0)
  
  
  
  #### under sampling ####
  max.class = names(which.max(table(B.samp$Z)))
  B_maj = B.samp[B.samp$Z == max.class,]
  US_B = rbind(B_maj[sample(max(table(B.samp$Z)), size = min(table(B.samp$Z))),], B.samp[B.samp$Z != max.class,])
  glmO = glm(Z ~ size2  + indgrp + state, data = US_B, family=quasibinomial(link="logit"))
  O = exp(predict(glmO, newdata = sample.NP.samp))
  f = mean(US_B$Z)
  p = sum(1-B.samp$Z)*f/(sum(B.samp$Z)*(1-f))
  f18_US = 1+(W.NP.samp-1)*p/O
  
  
  
  #### mix ####
  US_B2 = rbind(B_maj[sample(max(table(B.samp$Z)), size = 0.5*nrow(B.samp)),], B.samp[B.samp$Z != max.class,])
  US_B2$Z = as.factor(US_B2$Z)
  OS_B = smotenc(US_B2,var="Z")
  OS_B$Z = as.numeric(OS_B$Z)-1
  glmO = glm(Z ~ size2  + indgrp + state, data = OS_B, family=quasibinomial(link="logit"))
  O = exp(predict(glmO, newdata = sample.NP.samp))
  f = mean(OS_B$Z)
  p = sum(1-B.samp$Z)*f/(sum(B.samp$Z)*(1-f))
  f18_mix = 1+(W.NP.samp-1)*p/O
  
  
  
  ###############################
  
  ### Calculate estimates for CE sector
  
  sample.P.ce <- sample.P[sample.P$deswgt == 1,]
  cemodel <- glm(delta.i ~ size2  + indgrp + state, data = sample.P.ce, family=quasibinomial(link="logit"))

  ### Apply model to the P sample CE sector

  sample.P.ce$prop <- predict(cemodel, newdata=sample.P.ce, type="response")

  ### Limit the CE units in the P sample to those in the NP sample
  sample.P.ce.bd <- sample.P.ce[sample.P.ce$delta.i == 1,]
  ce.wgt <- 1/sample.P.ce.bd$prop
  


  
  #### estimates ####
  ce.pro = nrow(sample.P.ce)/N
  
  est = c( mean(sample.NP$ovt), 
           ((1-ce.pro)*sum(sample.NP.samp$ovt*f18)/sum(f18)         + ce.pro*sum(sample.P.ce.bd$ovt*ce.wgt)/sum(ce.wgt)),
           ((1-ce.pro)*sum(sample.NP.samp$ovt*f18_OS)/sum(f18_OS)   + ce.pro*sum(sample.P.ce.bd$ovt*ce.wgt)/sum(ce.wgt)),
           ((1-ce.pro)*sum(sample.NP.samp$ovt*f18_OS2)/sum(f18_OS2) + ce.pro*sum(sample.P.ce.bd$ovt*ce.wgt)/sum(ce.wgt)),
           ((1-ce.pro)*sum(sample.NP.samp$ovt*f18_US)/sum(f18_US)   + ce.pro*sum(sample.P.ce.bd$ovt*ce.wgt)/sum(ce.wgt)),
           ((1-ce.pro)*sum(sample.NP.samp$ovt*f18_mix)/sum(f18_mix) + ce.pro*sum(sample.P.ce.bd$ovt*ce.wgt)/sum(ce.wgt)))
  
  bias =(est - mean(pop$ovt))/M + bias
  mse = (est - mean(pop$ovt))^2/M + mse
  
  
  est2 = c(mean(sample.NP$earnings),
           ((1-ce.pro)*sum(sample.NP.samp$earnings*f18)/sum(f18)         + ce.pro*sum(sample.P.ce.bd$earnings*ce.wgt)/sum(ce.wgt)),
           ((1-ce.pro)*sum(sample.NP.samp$earnings*f18_OS)/sum(f18_OS)   + ce.pro*sum(sample.P.ce.bd$earnings*ce.wgt)/sum(ce.wgt)),
           ((1-ce.pro)*sum(sample.NP.samp$earnings*f18_OS2)/sum(f18_OS2) + ce.pro*sum(sample.P.ce.bd$earnings*ce.wgt)/sum(ce.wgt)),
           ((1-ce.pro)*sum(sample.NP.samp$earnings*f18_US)/sum(f18_US)   + ce.pro*sum(sample.P.ce.bd$earnings*ce.wgt)/sum(ce.wgt)),
           ((1-ce.pro)*sum(sample.NP.samp$earnings*f18_mix)/sum(f18_mix) + ce.pro*sum(sample.P.ce.bd$earnings*ce.wgt)/sum(ce.wgt)))
  
  bias2 =(est2 - mean(pop$earnings))/M + bias2
  mse2 = (est2 - mean(pop$earnings))^2/M + mse2
 
}
  
  return(list(n.P/N, n.NP/N, mse, bias, mse2, bias2))
  
}



# parallel

library(doRNG); library(doParallel)

registerDoParallel(cores = min(c(length(scenario), (detectCores()-1))))
registerDoRNG(13)

result <- foreach(i = c(1:length(scenario)),
                  .packages = c('sampling', 'themis', 'FNN'),
                  .combine = rbind) %dopar% scen.fun(i)

table = matrix(data = unlist(result[,c(1,2)]), nrow = length(scenario))
table = cbind(table, matrix(data = unlist(result[,3]), nrow = length(scenario), byrow = T))
table = cbind(table, matrix(data = unlist(result[,4]), nrow = length(scenario), byrow = T))
table = cbind(table, matrix(data = unlist(result[,5]), nrow = length(scenario), byrow = T))
table = cbind(table, matrix(data = unlist(result[,6]), nrow = length(scenario), byrow = T))

colnames(table) = c('nP', 'nNP',
                    'mse_naive','mse_LSdW','mse_smote','mse_smote2','mse_us','mse_mix',
                    'bias_naive','bias_LSdW','bias_smote','bias_smote2','bias_us','bias_mix',
                    'mse_naive','mse_LSdW','mse_smote','mse_smote2','mse_us','mse_mix',
                    'bias_naive','bias_LSdW','bias_smote','bias_smote2','bias_us','bias_mix')

save.image("app_pps.RData")

#table
# RB
RBtable = table[,c(1:26)]
colnames(RBtable) = c('nP', 'nNP',
                      'mse_np','mse_LSdW','mse_smote','mse_smote2','mse_us','mse_mix',
                      'RB_np','RB_LSdW','RB_smote','RB_smote2','RB_us', 'RB_mix',
                      'mse_np','mse_LSdW','mse_smote','mse_smote2','mse_us','mse_mix',
                      'RB_np','RB_LSdW','RB_smote','RB_smote2','RB_us', 'RB_mix')

RBtable[,3] = 1e-2*table[,3]
RBtable[,4] = 1e-2*table[,4]
RBtable[,5] = 1e-2*table[,5]
RBtable[,6] = 1e-2*table[,6]
RBtable[,7] = 1e-2*table[,7]
RBtable[,8] = 1e-2*table[,8]

RBtable[,9] = 100*table[,9]/mean(pop$ovt)
RBtable[,10] = 100*table[,10]/mean(pop$ovt)
RBtable[,11] = 100*table[,11]/mean(pop$ovt)
RBtable[,12] = 100*table[,12]/mean(pop$ovt)
RBtable[,13] = 100*table[,13]/mean(pop$ovt)
RBtable[,14] = 100*table[,14]/mean(pop$ovt)

RBtable[,15] = 1e-5*table[,15]
RBtable[,16] = 1e-5*table[,16]
RBtable[,17] = 1e-5*table[,17]
RBtable[,18] = 1e-5*table[,18]
RBtable[,19] = 1e-5*table[,19]
RBtable[,20] = 1e-5*table[,20]

RBtable[,21] = 100*table[,21]/mean(pop$earnings)
RBtable[,22] = 100*table[,22]/mean(pop$earnings)
RBtable[,23] = 100*table[,23]/mean(pop$earnings)
RBtable[,24] = 100*table[,24]/mean(pop$earnings)
RBtable[,25] = 100*table[,25]/mean(pop$earnings)
RBtable[,26] = 100*table[,26]/mean(pop$earnings)

round(RBtable,2)

capture.output(round(RBtable,3), file = 'app_pps.txt')


library("xtable")
RBtable = as.data.frame(RBtable)
RBtable[,c(1,2)] = round(RBtable[,c(1,2)]*100)
Latex = xtable(RBtable[, c(2:8)], digits = c(0,0, c(rep(3,6))), caption = 'Ovt mse estimates under pps')
Latex

Latex = xtable(RBtable[, c(2, 15:20)], digits = c(0,0, c(rep(3,6))), caption = 'earnings mse estimates under pps')
Latex

Latex = xtable(RBtable[, c(2,9:14)], digits = c(0,0, c(rep(3,6))), caption = 'Ovt RB estimates under pps')
Latex

Latex = xtable(RBtable[, c(2, 21:26)], digits = c(0,0, c(rep(3,6))), caption = 'earnings RB estimates under pps')
Latex

