library(exchanger)
library(comparator)
library(clevr)
library(progress)
library(Matrix)
library(testit)
library("BRL")
source("~/Documents/FlexRL/recordlinkage.r")
Rcpp:::sourceCpp("~/Documents/FlexRL/functions.cpp")

nbrSimu = 500

ExchangerNSamples = 2000
ExchangerNBurnin = 1000

MatchMindNIter = 300
MatchMindNBurnin = 150 
MatchMindNMStepIter = 25
MatchMindNcutBurninStEM = 150

didnotworkExchanger = 0

NdataA = 300
NdataB = 400
Nlinks = 200

Delta = matrix(0, nrow=NdataA, ncol=NdataB)
for (l in 1:Nlinks)
{
  Delta[l,l]=1
}

XA = data.frame()
XB = data.frame()
pSameH.varA = list(c(), c(), c(), c(), c())
pSameH.varB = list(c(), c(), c(), c(), c())
pivs_stable = c(TRUE, TRUE, TRUE, TRUE, FALSE)
pivs = c("V1", "V2", "V3", "V4", "V5")

Nval       = c(    15,     16,     17,     18,     25)
ptypos     = c(  0.03,   0.03,   0.03,   0.03,   0.03)
pmissing1  = c(  0.01,   0.01,   0.01,   0.01,   0.01)
pmissing2  = c(  0.01,   0.01,   0.01,   0.01,   0.01)

NRecords = c(NdataA, NdataB)

results_f1score = data.frame( matrix(0, nrow=nbrSimu, ncol=5) )
colnames(results_f1score) = c("ExchangerSelf", "Exchanger", "BRL", "MatchMind", "MatchMindstable")
results_recall = data.frame( matrix(0, nrow=nbrSimu, ncol=5) )
colnames(results_recall) = c("ExchangerSelf", "Exchanger", "BRL", "MatchMind", "MatchMindstable")
results_precision = data.frame( matrix(0, nrow=nbrSimu, ncol=5) )
colnames(results_precision) = c("ExchangerSelf", "Exchanger", "BRL", "MatchMind", "MatchMindstable")
results_FN = data.frame( matrix(0, nrow=nbrSimu, ncol=5) )
colnames(results_FN) = c("ExchangerSelf", "Exchanger", "BRL", "MatchMind", "MatchMindstable")
results_FP = data.frame( matrix(0, nrow=nbrSimu, ncol=5) )
colnames(results_FP) = c("ExchangerSelf", "Exchanger", "BRL", "MatchMind", "MatchMindstable")
results_TP = data.frame( matrix(0, nrow=nbrSimu, ncol=5) )
colnames(results_TP) = c("ExchangerSelf", "Exchanger", "BRL", "MatchMind", "MatchMindstable")
results_distance = data.frame( matrix(0, nrow=nbrSimu, ncol=5) )
colnames(results_distance) = c("ExchangerSelf", "Exchanger", "BRL", "MatchMind", "MatchMindstable")
results_exchangerdedup = data.frame( matrix(0, nrow=nbrSimu, ncol=1) )
colnames(results_exchangerdedup) = c("Exchanger")
results_gamma = data.frame( matrix(0, nrow=nbrSimu, ncol=MatchMindNIter-MatchMindNcutBurninStEM) )
results_omega_param = data.frame( matrix(0, nrow=nbrSimu, ncol=MatchMindNIter-MatchMindNcutBurninStEM) )
results_omega_probaEstimate = data.frame( matrix(0, nrow=nbrSimu, ncol=NdataA) )
results_omega_timesEstimate = data.frame( matrix(0, nrow=nbrSimu, ncol=NdataA) )
results_omega_probaTrue = data.frame( matrix(0, nrow=nbrSimu, ncol=Nlinks) )
results_omega_timesTrue = data.frame( matrix(0, nrow=nbrSimu, ncol=Nlinks) )
resultsstable_gamma = data.frame( matrix(0, nrow=nbrSimu, ncol=MatchMindNIter-MatchMindNcutBurninStEM) )
# overall data generated
simu_NA_A = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simu_NA_A) = pivs
simu_NA_B = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simu_NA_B) = pivs
simu_uniquevalues = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simu_uniquevalues) = pivs
simu_truelinks_agree = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simu_truelinks_agree) = pivs
simu_truelinks_unstable_change = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simu_truelinks_unstable_change) = pivs
# NA for exchanger results
simuExchanger_NAA_TP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuExchanger_NAA_TP) = pivs
simuExchanger_NAA_FP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuExchanger_NAA_FP) = pivs
simuExchanger_NAA_FN = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuExchanger_NAA_FN) = pivs
simuExchanger_NAA_linked = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuExchanger_NAA_linked) = pivs
simuExchanger_NAB_TP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuExchanger_NAB_TP) = pivs
simuExchanger_NAB_FP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuExchanger_NAB_FP) = pivs
simuExchanger_NAB_FN = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuExchanger_NAB_FN) = pivs
simuExchanger_NAB_linked = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuExchanger_NAB_linked) = pivs
# agreements for exchanger results
simuExchanger_agree_TP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuExchanger_agree_TP) = pivs
simuExchanger_agree_FP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuExchanger_agree_FP) = pivs
simuExchanger_agree_FN = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuExchanger_agree_FN) = pivs
simuExchanger_agree_linked = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuExchanger_agree_linked) = pivs
# changes for exchanger results
simuExchanger_unstable_change_TP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuExchanger_unstable_change_TP) = pivs
simuExchanger_unstable_change_FP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuExchanger_unstable_change_FP) = pivs
simuExchanger_unstable_change_FN = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuExchanger_unstable_change_FN) = pivs
simuExchanger_unstable_change_linked = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuExchanger_unstable_change_linked) = pivs
# NA for BRL results
simuBRL_NAA_TP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuBRL_NAA_TP) = pivs
simuBRL_NAA_FP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuBRL_NAA_FP) = pivs
simuBRL_NAA_FN = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuBRL_NAA_FN) = pivs
simuBRL_NAA_linked = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuBRL_NAA_linked) = pivs
simuBRL_NAB_TP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuBRL_NAB_TP) = pivs
simuBRL_NAB_FP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuBRL_NAB_FP) = pivs
simuBRL_NAB_FN = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuBRL_NAB_FN) = pivs
simuBRL_NAB_linked = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuBRL_NAB_linked) = pivs
# agreements for BRL results
simuBRL_agree_TP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuBRL_agree_TP) = pivs
simuBRL_agree_FP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuBRL_agree_FP) = pivs
simuBRL_agree_FN = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuBRL_agree_FN) = pivs
simuBRL_agree_linked = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuBRL_agree_linked) = pivs
# changes for BRL results
simuBRL_unstable_change_TP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuBRL_unstable_change_TP) = pivs
simuBRL_unstable_change_FP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuBRL_unstable_change_FP) = pivs
simuBRL_unstable_change_FN = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuBRL_unstable_change_FN) = pivs
simuBRL_unstable_change_linked = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuBRL_unstable_change_linked) = pivs
# NA for our results
simuMatchMind_NAA_TP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMind_NAA_TP) = pivs
simuMatchMind_NAA_FP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMind_NAA_FP) = pivs
simuMatchMind_NAA_FN = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMind_NAA_FN) = pivs
simuMatchMind_NAA_linked = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMind_NAA_linked) = pivs
simuMatchMind_NAB_TP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMind_NAB_TP) = pivs
simuMatchMind_NAB_FP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMind_NAB_FP) = pivs
simuMatchMind_NAB_FN = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMind_NAB_FN) = pivs
simuMatchMind_NAB_linked = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMind_NAB_linked) = pivs
# agreements for our results
simuMatchMind_agree_TP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMind_agree_TP) = pivs
simuMatchMind_agree_FP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMind_agree_FP) = pivs
simuMatchMind_agree_FN = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMind_agree_FN) = pivs
simuMatchMind_agree_linked = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMind_agree_linked) = pivs
# changes for our results
simuMatchMind_unstable_change_TP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMind_unstable_change_TP) = pivs
simuMatchMind_unstable_change_FP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMind_unstable_change_FP) = pivs
simuMatchMind_unstable_change_FN = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMind_unstable_change_FN) = pivs
simuMatchMind_unstable_change_linked = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMind_unstable_change_linked) = pivs
# NA for our results STABLE
simuMatchMindstable_NAA_TP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMindstable_NAA_TP) = pivs
simuMatchMindstable_NAA_FP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMindstable_NAA_FP) = pivs
simuMatchMindstable_NAA_FN = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMindstable_NAA_FN) = pivs
simuMatchMindstable_NAA_linked = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMindstable_NAA_linked) = pivs
simuMatchMindstable_NAB_TP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMindstable_NAB_TP) = pivs
simuMatchMindstable_NAB_FP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMindstable_NAB_FP) = pivs
simuMatchMindstable_NAB_FN = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMindstable_NAB_FN) = pivs
simuMatchMindstable_NAB_linked = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMindstable_NAB_linked) = pivs
# agreements for our results
simuMatchMindstable_agree_TP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMindstable_agree_TP) = pivs
simuMatchMindstable_agree_FP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMindstable_agree_FP) = pivs
simuMatchMindstable_agree_FN = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMindstable_agree_FN) = pivs
simuMatchMindstable_agree_linked = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMindstable_agree_linked) = pivs
# changes for our results
simuMatchMindstable_unstable_change_TP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMindstable_unstable_change_TP) = pivs
simuMatchMindstable_unstable_change_FP = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMindstable_unstable_change_FP) = pivs
simuMatchMindstable_unstable_change_FN = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMindstable_unstable_change_FN) = pivs
simuMatchMindstable_unstable_change_linked = data.frame( matrix(0, nrow=nbrSimu, ncol=length(pivs)) )
colnames(simuMatchMindstable_unstable_change_linked) = pivs

for(iterSimu in 1:nbrSimu){
  
  pivs_stable = c(TRUE, TRUE, TRUE, TRUE, FALSE)
  
  ### DATA CREATION
  
  for(i in 1:2)
  {	
    dataSet=c()
    for(u in 1:length(Nval))
    {
      # Simulate different probabilities for the true values (if desired)
      # TODOTEST
      xp =  exp(0.1 *(0:(Nval[u]-1)))
      # xp = exp(0 *(0:(Nval[u]-1)))
      probx = xp/sum(xp)
      
      dataSet = cbind(dataSet, sample(1:Nval[u], NRecords[i], replace=TRUE, prob=probx))
    }
    dataSet = as.data.frame(dataSet)
    names(dataSet) = c(paste("V", 1:(ncol(dataSet)), sep=""))
    assign(paste("dataSet", i, sep=""),dataSet  )
  }
  
  # Add overlapping units
  # match = rep(FALSE,nrow(dataSet1))
  # match[1:Nlinks] = TRUE
  dataSet2[1:Nlinks,] = dataSet1[1:Nlinks,]
  
  # Add typos
  for(x in 1:ncol(dataSet1))
  { 
    biased = as.logical(rbinom(nrow(dataSet1),1,ptypos[x]))
    if(sum(biased)>0)
      dataSet1[,x][biased] = sapply(dataSet1[,x][biased], function(i) sample((1:Nval[x])[-c(i)], 1))
  }  
  
  for(x in 1:ncol(dataSet2))
  { 
    biased = as.logical(rbinom(nrow(dataSet2),1,ptypos[x]))
    if(sum(biased)>0)
      dataSet2[,x][biased] = sapply(dataSet2[,x][biased], function(i) sample((1:Nval[x])[-c(i)], 1))
  } 
  
  # Add missings
  for(x in 1:ncol(dataSet1))
  { 
    biased = as.logical(rbinom(nrow(dataSet1),1,pmissing1[x]))
    if(sum(biased)>0)
      dataSet1[,x][biased] = NA
  }  
  
  for(x in 1:ncol(dataSet2))
  { 
    biased = as.logical(rbinom(nrow(dataSet2),1,pmissing2[x]))
    if(sum(biased)>0)
      dataSet2[,x][biased] = NA
  }  
  
  # Add registration time for the survival model
  
  regTimeA = runif(nrow(dataSet1), 0, 1)
  regTimeB = runif(nrow(dataSet2), 1, 2)
  time_difference = regTimeB[1:Nlinks] - regTimeA[1:Nlinks]
  
  # TODOTEST
  regTimeA = runif(nrow(dataSet1), 0, 3)
  regTimeB = runif(nrow(dataSet2), 3, 6)
  time_difference = regTimeB[1:Nlinks] - regTimeA[1:Nlinks]
  # Create a moving process for a potential unstable variable (proba of same values decreases with time)
  proba_same_H = exp( - 0.28 * (time_difference) ) # 0.28
  # plot( sort(proba_same_H), ylim=c(0,1), xlab = "true records pairs", ylab = sprintf("true model for proba of no change in %s", pivs[5]) )
  # plot( time_difference, proba_same_H, ylim=c(0,1), xlab = "time difference", ylab = sprintf("true model for proba of no change in %s", pivs[5]) )
  
  dataSet1$change = FALSE
  dataSet2$change = FALSE
  
  x = which(!pivs_stable)
  for(i in 1:Nlinks)
  {
    is_not_moving = rbinom(1, 1, proba_same_H[i])
    if(!is_not_moving)
    {
      # the value MUST be different
      dataSet2[i,x] = sample((1:Nval[x])[-c(dataSet1[i,x])], 1)
      dataSet2[i,"change"] = TRUE
    }
  }  
  A = dataSet1
  B = dataSet2
  
  # Recode the pivs
  levels_pivs = lapply(pivs, function(x) levels(factor(as.character(c(A[,x], B[,x])))))
  
  for(i in 1:length(pivs))
  {
    A[,pivs[i]] = as.numeric(factor(as.character(A[,pivs[i]]), levels=levels_pivs[[i]]))
    B[,pivs[i]] = as.numeric(factor(as.character(B[,pivs[i]]), levels=levels_pivs[[i]]))
  }
  
  nvalues = sapply(levels_pivs, length)
  
  A$localID = 1:nrow(A)
  B$localID = 1:nrow(B)
  
  A$source = "A"
  B$source = "B"
  
  # PREPARE DATA FOR EXCHANGER
  RLdata = rbind(A[,pivs], B[,pivs])
  rownames(RLdata) = 1:nrow(RLdata)
  
  RLdata = rbind(A, B)
  rownames(RLdata) = 1:nrow(RLdata)
  true_pairs = cbind( rownames(RLdata[(RLdata$source=="A")&(RLdata$localID %in% 1:Nlinks),]),
                      rownames(RLdata[(RLdata$source=="B")&(RLdata$localID %in% 1:Nlinks),]) )
  
  # check that data are encoded (without 0 now)
  # for(i in 1:length(pivs)){
  #   print( sort(unique(RLdata[,pivs[i]])) )
  # }
  
  # STORY TELLING ABOUT THE DATA
  simu_NA_A[iterSimu,] = ( colSums(is.na(A)) / nrow(A) )[pivs]
  simu_NA_B[iterSimu,] = ( colSums(is.na(B)) / nrow(B) )[pivs]
  simu_uniquevalues[iterSimu,] = nvalues
  
  simu_truelinks_agree_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simu_truelinks_agree_tmp) = pivs
  simu_truelinks_unstable_change_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simu_truelinks_unstable_change_tmp) = pivs
  for( i in 1:Nlinks ){
    entityA = A[i,]
    entityB = B[i,]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        simu_truelinks_agree_tmp = rbind(simu_truelinks_agree_tmp, entityA[,pivs] == entityB[,pivs])
        simu_truelinks_unstable_change_tmp = rbind(simu_truelinks_unstable_change_tmp, entityB[,"change"]*(!pivs_stable))
      }
    }
  }
  simu_truelinks_agree_tmp = colSums( simu_truelinks_agree_tmp, na.rm = TRUE ) / nrow( simu_truelinks_agree_tmp )
  simu_truelinks_unstable_change_tmp = colSums( simu_truelinks_unstable_change_tmp, na.rm = TRUE ) / nrow( simu_truelinks_unstable_change_tmp )
  simu_truelinks_agree[iterSimu,] = simu_truelinks_agree_tmp
  simu_truelinks_unstable_change[iterSimu,] = simu_truelinks_unstable_change_tmp
  
  ### LAUNCH EXCHANGER (MARCHANT, STEORTS, ...)
  unif_prior <- BetaRV(1, 1)
  attr_params <- list(
    V1 = CategoricalAttribute(distort_prob_prior = unif_prior),
    V2 = CategoricalAttribute(distort_prob_prior = unif_prior),
    V3 = CategoricalAttribute(distort_prob_prior = unif_prior),
    V4 = CategoricalAttribute(distort_prob_prior = unif_prior),
    V5 = CategoricalAttribute(distort_prob_prior = unif_prior)
  )
  clust_prior <- PitmanYorRP(alpha = GammaRV(1, 1), d = BetaRV(1, 1))
  
  model <- exchanger(RLdata, attr_params, clust_prior)
  
  result <- run_inference(model, n_samples=ExchangerNSamples, thin_interval=10, burnin_interval=ExchangerNBurnin)
  pred_clust <- smp_clusters(result)
  n_records <- nrow(RLdata)
  pred_pairs <- clusters_to_pairs(pred_clust)
  if ( nrow(pred_pairs)<=1 ){
    print("RE RUN EXCHANGER DID NOT WORK")
    didnotworkExchanger = didnotworkExchanger + 1
    next
  }
  
  measures <- eval_report_pairs(true_pairs, pred_pairs, num_pairs=n_records*(n_records-1)/2)
  
  DeltaExchanger = matrix(0, nrow=nrow(A), ncol=nrow(B))
  dedup = 0
  for (l in 1:nrow(pred_pairs))
  {
    idx_row_RLdata_A = as.integer(pred_pairs[l,1])
    idx_row_RLdata_B = as.integer(pred_pairs[l,2])
    idxA = RLdata[rownames(RLdata) == idx_row_RLdata_A, ]$localID
    idxB = RLdata[rownames(RLdata) == idx_row_RLdata_B, ]$localID
    
    if( (idxA <= nrow(DeltaExchanger)) & (idxB <= ncol(DeltaExchanger)) ){
      DeltaExchanger[idxA,idxB] = 1
    }else{
      dedup = dedup + 1
    }
  }
  # image(x=1:nrow(A), y=1:nrow(B), z=as.matrix(DeltaExchanger), xlab="obs. in A", ylab="obs. in B")
  # title(main = c("Estimated linkage matrix fitting the data", "Exchanger method", sprintf("%s linked record pairs", sum(DeltaExchanger))), font.main = 4)
  
  # STORY TELLING ABOUT THE DATA LINKED!!!
  
  linkedpairs = which(DeltaExchanger>0.5, arr.ind=TRUE)
  linkedpairsA = linkedpairs[,1]
  linkedpairsB = linkedpairs[,2]
  simuExchanger_NAA_linked[iterSimu,] = ( colSums(is.na(A[linkedpairsA,])) / length(linkedpairsA) )[pivs]
  simuExchanger_NAB_linked[iterSimu,] = ( colSums(is.na(B[linkedpairsB,])) / length(linkedpairsB) )[pivs]
  
  simuExchanger_agree_linked_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuExchanger_agree_linked_tmp) = pivs
  simuExchanger_unstable_change_linked_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuExchanger_unstable_change_linked_tmp) = pivs
  for( i in 1:nrow(linkedpairs) ){
    entityA = A[linkedpairsA[i],]
    entityB = B[linkedpairsB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        simuExchanger_agree_linked_tmp = rbind(simuExchanger_agree_linked_tmp, entityA[,pivs] == entityB[,pivs])
        simuExchanger_unstable_change_linked_tmp = rbind(simuExchanger_unstable_change_linked_tmp, entityB[,"change"]*(!pivs_stable))
      }
    }
  }
  simuExchanger_agree_linked_tmp = colSums( simuExchanger_agree_linked_tmp, na.rm = TRUE ) / nrow( simuExchanger_agree_linked_tmp )
  simuExchanger_unstable_change_linked_tmp = colSums( simuExchanger_unstable_change_linked_tmp, na.rm = TRUE ) / nrow( simuExchanger_unstable_change_linked_tmp )
  simuExchanger_agree_linked[iterSimu,] = simuExchanger_agree_linked_tmp
  simuExchanger_unstable_change_linked[iterSimu,] = simuExchanger_unstable_change_linked_tmp
  
  linkedTP = which((DeltaExchanger>0.5) & (Delta==1), arr.ind=TRUE)
  linkedTPA = linkedTP[,1]
  linkedTPB = linkedTP[,2]
  simuExchanger_NAA_TP[iterSimu,] = ( colSums(is.na(A[linkedTPA,])) / length(linkedTPA) )[pivs]
  simuExchanger_NAB_TP[iterSimu,] = ( colSums(is.na(B[linkedTPB,])) / length(linkedTPB) )[pivs]
  
  simuExchanger_agree_TP_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuExchanger_agree_TP_tmp) = pivs
  simuExchanger_unstable_change_TP_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuExchanger_unstable_change_TP_tmp) = pivs
  for( i in 1:nrow(linkedTP) ){
    entityA = A[linkedTPA[i],]
    entityB = B[linkedTPB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        simuExchanger_agree_TP_tmp = rbind(simuExchanger_agree_TP_tmp, entityA[,pivs] == entityB[,pivs])
        simuExchanger_unstable_change_TP_tmp = rbind(simuExchanger_unstable_change_TP_tmp, entityB[,"change"]*(!pivs_stable))
      }
    }
  }
  simuExchanger_agree_TP_tmp = colSums( simuExchanger_agree_TP_tmp, na.rm = TRUE ) / nrow( simuExchanger_agree_TP_tmp )
  simuExchanger_unstable_change_TP_tmp = colSums( simuExchanger_unstable_change_TP_tmp, na.rm = TRUE ) / nrow( simuExchanger_unstable_change_TP_tmp )
  simuExchanger_agree_TP[iterSimu,] = simuExchanger_agree_TP_tmp
  simuExchanger_unstable_change_TP[iterSimu,] = simuExchanger_unstable_change_TP_tmp
  
  linkedFP = which((DeltaExchanger>0.5) & (Delta==0), arr.ind=TRUE)
  linkedFPA = linkedFP[,1]
  linkedFPB = linkedFP[,2]
  simuExchanger_NAA_FP[iterSimu,] = ( colSums(is.na(A[linkedFPA,])) / length(linkedFPA) )[pivs]
  simuExchanger_NAB_FP[iterSimu,] = ( colSums(is.na(B[linkedFPB,])) / length(linkedFPB) )[pivs]
  
  simuExchanger_agree_FP_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuExchanger_agree_FP_tmp) = pivs
  simuExchanger_unstable_change_FP_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuExchanger_unstable_change_FP_tmp) = pivs
  for( i in 1:nrow(linkedFP) ){
    entityA = A[linkedFPA[i],]
    entityB = B[linkedFPB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        simuExchanger_agree_FP_tmp = rbind(simuExchanger_agree_FP_tmp, entityA[,pivs] == entityB[,pivs])
        simuExchanger_unstable_change_FP_tmp = rbind(simuExchanger_unstable_change_FP_tmp, entityB[,"change"]*(!pivs_stable))
      }
    }
  }
  simuExchanger_agree_FP_tmp = colSums( simuExchanger_agree_FP_tmp, na.rm = TRUE ) / nrow( simuExchanger_agree_FP_tmp )
  simuExchanger_unstable_change_FP_tmp = colSums( simuExchanger_unstable_change_FP_tmp, na.rm = TRUE ) / nrow( simuExchanger_unstable_change_FP_tmp )
  simuExchanger_agree_FP[iterSimu,] = simuExchanger_agree_FP_tmp
  simuExchanger_unstable_change_FP[iterSimu,] = simuExchanger_unstable_change_FP_tmp
  
  linkedFN = which((DeltaExchanger<0.5) & (Delta==1), arr.ind=TRUE)
  linkedFNA = linkedFN[,1]
  linkedFNB = linkedFN[,2]
  simuExchanger_NAA_FN[iterSimu,] = ( colSums(is.na(A[linkedFNA,])) / length(linkedFNA) )[pivs]
  simuExchanger_NAB_FN[iterSimu,] = ( colSums(is.na(B[linkedFNB,])) / length(linkedFNB) )[pivs]
  
  simuExchanger_agree_FN_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuExchanger_agree_FN_tmp) = pivs
  simuExchanger_unstable_change_FN_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuExchanger_unstable_change_FN_tmp) = pivs
  for( i in 1:nrow(linkedFN) ){
    entityA = A[linkedFNA[i],]
    entityB = B[linkedFNB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        simuExchanger_agree_FN_tmp = rbind(simuExchanger_agree_FN_tmp, entityA[,pivs] == entityB[,pivs])
        simuExchanger_unstable_change_FN_tmp = rbind(simuExchanger_unstable_change_FN_tmp, entityB[,"change"]*(!pivs_stable))
      }
    }
  }
  simuExchanger_agree_FN_tmp = colSums( simuExchanger_agree_FN_tmp, na.rm = TRUE ) / nrow( simuExchanger_agree_FN_tmp )
  simuExchanger_unstable_change_FN_tmp = colSums( simuExchanger_unstable_change_FN_tmp, na.rm = TRUE ) / nrow( simuExchanger_unstable_change_FN_tmp )
  simuExchanger_agree_FN[iterSimu,] = simuExchanger_agree_FN_tmp
  simuExchanger_unstable_change_FN[iterSimu,] = simuExchanger_unstable_change_FN_tmp
  
  truepositive = sum( (DeltaExchanger>0.5) & (Delta==1) )
  falsepositive = sum( (DeltaExchanger>0.5) & (Delta==0) ) + dedup
  falsenegative = sum( (DeltaExchanger<0.5) & (Delta==1) )
  precision = truepositive / (truepositive + falsepositive)
  recall = truepositive / (truepositive + falsenegative)
  f1score = 2 * (precision * recall) / (precision + recall)
  
  results_f1score[iterSimu, "Exchanger"] = f1score
  results_recall[iterSimu, "Exchanger"] = recall
  results_precision[iterSimu, "Exchanger"] = precision
  results_FN[iterSimu, "Exchanger"] = falsenegative
  results_FP[iterSimu, "Exchanger"] = falsepositive
  results_TP[iterSimu, "Exchanger"] = truepositive
  results_distance[iterSimu, "Exchanger"] = sqrt(sum((DeltaExchanger - Delta)**2))
  results_exchangerdedup[iterSimu, "Exchanger"] = dedup
  
  results_f1score[iterSimu, "ExchangerSelf"] = measures$f1score
  results_recall[iterSimu, "ExchangerSelf"] = measures$recall
  results_precision[iterSimu, "ExchangerSelf"] = measures$precision
  results_FN[iterSimu, "ExchangerSelf"] = falsenegative
  results_FP[iterSimu, "ExchangerSelf"] = falsepositive
  results_TP[iterSimu, "ExchangerSelf"] = truepositive
  results_distance[iterSimu, "ExchangerSelf"] = sqrt(sum((DeltaExchanger - Delta)**2))
  
  # posterior_DeltaExchanger = posterior_DeltaExchanger + DeltaExchanger
  
  ### LAUNCH BRL (SADINLE)
  
  DeltaBRL = matrix(0, nrow=nrow(A), ncol=nrow(B))
  Zhat <- BRL(B, A, flds=pivs, types=c("bi","bi","bi","bi","bi"), nIter=1000) # breaks = c(0, .25, .5)
  n1 <- nrow(B)
  idxA = which( Zhat <= n1 ) # index in A
  idxB = Zhat[ Zhat <= n1 ] # index in B
  for (l in 1:length(idxA))
  {
    DeltaBRL[idxA[l], idxB[l]] = 1
  }
  # image(x=1:nrow(A), y=1:nrow(B), z=as.matrix(DeltaBRL), xlab="obs. in A", ylab="obs. in B")
  # title(main = c("Estimated linkage matrix fitting the data", "BRL method (Sadinle)", sprintf("%s linked record pairs", sum(DeltaBRL))), font.main = 4)
  
  # STORY TELLING ABOUT THE DATA LINKED!!!
  
  linkedpairs = which(DeltaBRL>0.5, arr.ind=TRUE)
  linkedpairsA = linkedpairs[,1]
  linkedpairsB = linkedpairs[,2]
  simuBRL_NAA_linked[iterSimu,] = ( colSums(is.na(A[linkedpairsA,])) / length(linkedpairsA) )[pivs]
  simuBRL_NAB_linked[iterSimu,] = ( colSums(is.na(B[linkedpairsB,])) / length(linkedpairsB) )[pivs]
  
  simuBRL_agree_linked_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuBRL_agree_linked_tmp) = pivs
  simuBRL_unstable_change_linked_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuBRL_unstable_change_linked_tmp) = pivs
  for( i in 1:nrow(linkedpairs) ){
    entityA = A[linkedpairsA[i],]
    entityB = B[linkedpairsB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        simuBRL_agree_linked_tmp = rbind(simuBRL_agree_linked_tmp, entityA[,pivs] == entityB[,pivs])
        simuBRL_unstable_change_linked_tmp = rbind(simuBRL_unstable_change_linked_tmp, entityB[,"change"]*(!pivs_stable))
      }
    }
  }
  simuBRL_agree_linked_tmp = colSums( simuBRL_agree_linked_tmp, na.rm = TRUE ) / nrow( simuBRL_agree_linked_tmp )
  simuBRL_unstable_change_linked_tmp = colSums( simuBRL_unstable_change_linked_tmp, na.rm = TRUE ) / nrow( simuBRL_unstable_change_linked_tmp )
  simuBRL_agree_linked[iterSimu,] = simuBRL_agree_linked_tmp
  simuBRL_unstable_change_linked[iterSimu,] = simuBRL_unstable_change_linked_tmp
  
  linkedTP = which((DeltaBRL>0.5) & (Delta==1), arr.ind=TRUE)
  linkedTPA = linkedTP[,1]
  linkedTPB = linkedTP[,2]
  simuBRL_NAA_TP[iterSimu,] = ( colSums(is.na(A[linkedTPA,])) / length(linkedTPA) )[pivs]
  simuBRL_NAB_TP[iterSimu,] = ( colSums(is.na(B[linkedTPB,])) / length(linkedTPB) )[pivs]
  
  simuBRL_agree_TP_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuBRL_agree_TP_tmp) = pivs
  simuBRL_unstable_change_TP_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuBRL_unstable_change_TP_tmp) = pivs
  for( i in 1:nrow(linkedTP) ){
    entityA = A[linkedTPA[i],]
    entityB = B[linkedTPB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        simuBRL_agree_TP_tmp = rbind(simuBRL_agree_TP_tmp, entityA[,pivs] == entityB[,pivs])
        simuBRL_unstable_change_TP_tmp = rbind(simuBRL_unstable_change_TP_tmp, entityB[,"change"]*(!pivs_stable))
      }
    }
  }
  simuBRL_agree_TP_tmp = colSums( simuBRL_agree_TP_tmp, na.rm = TRUE ) / nrow( simuBRL_agree_TP_tmp )
  simuBRL_unstable_change_TP_tmp = colSums( simuBRL_unstable_change_TP_tmp, na.rm = TRUE ) / nrow( simuBRL_unstable_change_TP_tmp )
  simuBRL_agree_TP[iterSimu,] = simuBRL_agree_TP_tmp
  simuBRL_unstable_change_TP[iterSimu,] = simuBRL_unstable_change_TP_tmp
  
  linkedFP = which((DeltaBRL>0.5) & (Delta==0), arr.ind=TRUE)
  linkedFPA = linkedFP[,1]
  linkedFPB = linkedFP[,2]
  simuBRL_NAA_FP[iterSimu,] = ( colSums(is.na(A[linkedFPA,])) / length(linkedFPA) )[pivs]
  simuBRL_NAB_FP[iterSimu,] = ( colSums(is.na(B[linkedFPB,])) / length(linkedFPB) )[pivs]
  
  simuBRL_agree_FP_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuBRL_agree_FP_tmp) = pivs
  simuBRL_unstable_change_FP_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuBRL_unstable_change_FP_tmp) = pivs
  for( i in 1:nrow(linkedFP) ){
    entityA = A[linkedFPA[i],]
    entityB = B[linkedFPB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        simuBRL_agree_FP_tmp = rbind(simuBRL_agree_FP_tmp, entityA[,pivs] == entityB[,pivs])
        simuBRL_unstable_change_FP_tmp = rbind(simuBRL_unstable_change_FP_tmp, entityB[,"change"]*(!pivs_stable))
      }
    }
  }
  simuBRL_agree_FP_tmp = colSums( simuBRL_agree_FP_tmp, na.rm = TRUE ) / nrow( simuBRL_agree_FP_tmp )
  simuBRL_unstable_change_FP_tmp = colSums( simuBRL_unstable_change_FP_tmp, na.rm = TRUE ) / nrow( simuBRL_unstable_change_FP_tmp )
  simuBRL_agree_FP[iterSimu,] = simuBRL_agree_FP_tmp
  simuBRL_unstable_change_FP[iterSimu,] = simuBRL_unstable_change_FP_tmp
  
  linkedFN = which((DeltaBRL<0.5) & (Delta==1), arr.ind=TRUE)
  linkedFNA = linkedFN[,1]
  linkedFNB = linkedFN[,2]
  simuBRL_NAA_FN[iterSimu,] = ( colSums(is.na(A[linkedFNA,])) / length(linkedFNA) )[pivs]
  simuBRL_NAB_FN[iterSimu,] = ( colSums(is.na(B[linkedFNB,])) / length(linkedFNB) )[pivs]
  
  simuBRL_agree_FN_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuBRL_agree_FN_tmp) = pivs
  simuBRL_unstable_change_FN_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuBRL_unstable_change_FN_tmp) = pivs
  for( i in 1:nrow(linkedFN) ){
    entityA = A[linkedFNA[i],]
    entityB = B[linkedFNB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        simuBRL_agree_FN_tmp = rbind(simuBRL_agree_FN_tmp, entityA[,pivs] == entityB[,pivs])
        simuBRL_unstable_change_FN_tmp = rbind(simuBRL_unstable_change_FN_tmp, entityB[,"change"]*(!pivs_stable))
      }
    }
  }
  simuBRL_agree_FN_tmp = colSums( simuBRL_agree_FN_tmp, na.rm = TRUE ) / nrow( simuBRL_agree_FN_tmp )
  simuBRL_unstable_change_FN_tmp = colSums( simuBRL_unstable_change_FN_tmp, na.rm = TRUE ) / nrow( simuBRL_unstable_change_FN_tmp )
  simuBRL_agree_FN[iterSimu,] = simuBRL_agree_FN_tmp
  simuBRL_unstable_change_FN[iterSimu,] = simuBRL_unstable_change_FN_tmp
  
  truepositive = sum( (DeltaBRL>0.5) & (Delta==1) )
  falsepositive = sum( (DeltaBRL>0.5) & (Delta==0) )
  falsenegative = sum( (DeltaBRL<0.5) & (Delta==1) )
  precision = truepositive / (truepositive + falsepositive)
  recall = truepositive / (truepositive + falsenegative)
  f1score = 2 * (precision * recall) / (precision + recall)
  
  results_f1score[iterSimu, "BRL"] = f1score
  results_recall[iterSimu, "BRL"] = recall
  results_precision[iterSimu, "BRL"] = precision
  results_FN[iterSimu, "BRL"] = falsenegative
  results_FP[iterSimu, "BRL"] = falsepositive
  results_TP[iterSimu, "BRL"] = truepositive
  results_distance[iterSimu, "BRL"] = sqrt(sum((DeltaBRL - Delta)**2))
  
  # posterior_DeltaBRL = posterior_DeltaBRL + DeltaBRL
  
  ### LAUNCH MATCHMIND (OURS)
  
  encodedA = A
  encodedB = B
  
  # Value for missing data
  encodedA[,pivs][ is.na(encodedA[,pivs]) ] = 0
  encodedB[,pivs][ is.na(encodedB[,pivs]) ] = 0
  
  dataSimu = list( pivsA        = encodedA[, pivs], 
                   pivsB        = encodedB[, pivs], 
                   nvalues      = nvalues,
                   pivs_stable  = pivs_stable, 
                   regTimeA     = regTimeA, 
                   regTimeB     = regTimeB,
                   XA           = XA,
                   XB           = XB,
                   pSameH.varA  = pSameH.varA,
                   pSameH.varB  = pSameH.varB,
                   include_time = TRUE)
  
  fit = stEM(  data          = dataSimu, 
               nIter         = MatchMindNIter, 
               nBurnin       = MatchMindNBurnin, 
               MStepIter     = MatchMindNMStepIter,
               cutBurninStEM = MatchMindNcutBurninStEM,
               trace      = 1 )
  
  linkedpairs = which(fit$Delta>0.5, arr.ind=TRUE)
  linkedpairsA = linkedpairs[,1]
  linkedpairsB = linkedpairs[,2]
  simuMatchMind_NAA_linked[iterSimu,] = ( colSums(is.na(A[linkedpairsA,])) / length(linkedpairsA) )[pivs]
  simuMatchMind_NAB_linked[iterSimu,] = ( colSums(is.na(B[linkedpairsB,])) / length(linkedpairsB) )[pivs]
  
  simuMatchMind_agree_linked_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuMatchMind_agree_linked_tmp) = pivs
  simuMatchMind_unstable_change_linked_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuMatchMind_unstable_change_linked_tmp) = pivs
  for( i in 1:nrow(linkedpairs) ){
    entityA = A[linkedpairsA[i],]
    entityB = B[linkedpairsB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        simuMatchMind_agree_linked_tmp = rbind(simuMatchMind_agree_linked_tmp, entityA[,pivs] == entityB[,pivs])
        simuMatchMind_unstable_change_linked_tmp = rbind(simuMatchMind_unstable_change_linked_tmp, entityB[,"change"]*(!pivs_stable))
      }
    }
  }
  simuMatchMind_agree_linked_tmp = colSums( simuMatchMind_agree_linked_tmp, na.rm = TRUE ) / nrow( simuMatchMind_agree_linked_tmp )
  simuMatchMind_unstable_change_linked_tmp = colSums( simuMatchMind_unstable_change_linked_tmp, na.rm = TRUE ) / nrow( simuMatchMind_unstable_change_linked_tmp )
  simuMatchMind_agree_linked[iterSimu,] = simuMatchMind_agree_linked_tmp
  simuMatchMind_unstable_change_linked[iterSimu,] = simuMatchMind_unstable_change_linked_tmp
  
  linkedTP = which((fit$Delta>0.5) & (Delta==1), arr.ind=TRUE)
  linkedTPA = linkedTP[,1]
  linkedTPB = linkedTP[,2]
  simuMatchMind_NAA_TP[iterSimu,] = ( colSums(is.na(A[linkedTPA,])) / length(linkedTPA) )[pivs]
  simuMatchMind_NAB_TP[iterSimu,] = ( colSums(is.na(B[linkedTPB,])) / length(linkedTPB) )[pivs]
  
  simuMatchMind_agree_TP_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuMatchMind_agree_TP_tmp) = pivs
  simuMatchMind_unstable_change_TP_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuMatchMind_unstable_change_TP_tmp) = pivs
  for( i in 1:nrow(linkedTP) ){
    entityA = A[linkedTPA[i],]
    entityB = B[linkedTPB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        simuMatchMind_agree_TP_tmp = rbind(simuMatchMind_agree_TP_tmp, entityA[,pivs] == entityB[,pivs])
        simuMatchMind_unstable_change_TP_tmp = rbind(simuMatchMind_unstable_change_TP_tmp, entityB[,"change"]*(!pivs_stable))
      }
    }
  }
  simuMatchMind_agree_TP_tmp = colSums( simuMatchMind_agree_TP_tmp, na.rm = TRUE ) / nrow( simuMatchMind_agree_TP_tmp )
  simuMatchMind_unstable_change_TP_tmp = colSums( simuMatchMind_unstable_change_TP_tmp, na.rm = TRUE ) / nrow( simuMatchMind_unstable_change_TP_tmp )
  simuMatchMind_agree_TP[iterSimu,] = simuMatchMind_agree_TP_tmp
  simuMatchMind_unstable_change_TP[iterSimu,] = simuMatchMind_unstable_change_TP_tmp
  
  linkedFP = which((fit$Delta>0.5) & (Delta==0), arr.ind=TRUE)
  linkedFPA = linkedFP[,1]
  linkedFPB = linkedFP[,2]
  simuMatchMind_NAA_FP[iterSimu,] = ( colSums(is.na(A[linkedFPA,])) / length(linkedFPA) )[pivs]
  simuMatchMind_NAB_FP[iterSimu,] = ( colSums(is.na(B[linkedFPB,])) / length(linkedFPB) )[pivs]
  
  simuMatchMind_agree_FP_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuMatchMind_agree_FP_tmp) = pivs
  simuMatchMind_unstable_change_FP_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuMatchMind_unstable_change_FP_tmp) = pivs
  for( i in 1:nrow(linkedFP) ){
    entityA = A[linkedFPA[i],]
    entityB = B[linkedFPB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        simuMatchMind_agree_FP_tmp = rbind(simuMatchMind_agree_FP_tmp, entityA[,pivs] == entityB[,pivs])
        simuMatchMind_unstable_change_FP_tmp = rbind(simuMatchMind_unstable_change_FP_tmp, entityB[,"change"]*(!pivs_stable))
      }
    }
  }
  simuMatchMind_agree_FP_tmp = colSums( simuMatchMind_agree_FP_tmp, na.rm = TRUE ) / nrow( simuMatchMind_agree_FP_tmp )
  simuMatchMind_unstable_change_FP_tmp = colSums( simuMatchMind_unstable_change_FP_tmp, na.rm = TRUE ) / nrow( simuMatchMind_unstable_change_FP_tmp )
  simuMatchMind_agree_FP[iterSimu,] = simuMatchMind_agree_FP_tmp
  simuMatchMind_unstable_change_FP[iterSimu,] = simuMatchMind_unstable_change_FP_tmp
  
  linkedFN = which((fit$Delta<0.5) & (Delta==1), arr.ind=TRUE)
  linkedFNA = linkedFN[,1]
  linkedFNB = linkedFN[,2]
  simuMatchMind_NAA_FN[iterSimu,] = ( colSums(is.na(A[linkedFNA,])) / length(linkedFNA) )[pivs]
  simuMatchMind_NAB_FN[iterSimu,] = ( colSums(is.na(B[linkedFNB,])) / length(linkedFNB) )[pivs]
  
  simuMatchMind_agree_FN_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuMatchMind_agree_FN_tmp) = pivs
  simuMatchMind_unstable_change_FN_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuMatchMind_unstable_change_FN_tmp) = pivs
  for( i in 1:nrow(linkedFN) ){
    entityA = A[linkedFNA[i],]
    entityB = B[linkedFNB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        simuMatchMind_agree_FN_tmp = rbind(simuMatchMind_agree_FN_tmp, entityA[,pivs] == entityB[,pivs])
        simuMatchMind_unstable_change_FN_tmp = rbind(simuMatchMind_unstable_change_FN_tmp, entityB[,"change"]*(!pivs_stable))
      }
    }
  }
  simuMatchMind_agree_FN_tmp = colSums( simuMatchMind_agree_FN_tmp, na.rm = TRUE ) / nrow( simuMatchMind_agree_FN_tmp )
  simuMatchMind_unstable_change_FN_tmp = colSums( simuMatchMind_unstable_change_FN_tmp, na.rm = TRUE ) / nrow( simuMatchMind_unstable_change_FN_tmp )
  simuMatchMind_agree_FN[iterSimu,] = simuMatchMind_agree_FN_tmp
  simuMatchMind_unstable_change_FN[iterSimu,] = simuMatchMind_unstable_change_FN_tmp
  
  truepositive = sum( (fit$Delta>0.5) & (Delta==1) )
  falsepositive = sum( (fit$Delta>0.5) & (Delta==0) )
  falsenegative = sum( (fit$Delta<0.5) & (Delta==1) )
  precision = truepositive / (truepositive + falsepositive)
  recall = truepositive / (truepositive + falsenegative)
  f1score = 2 * (precision * recall) / (precision + recall)
  
  results_f1score[iterSimu, "MatchMind"] = f1score
  results_recall[iterSimu, "MatchMind"] = recall
  results_precision[iterSimu, "MatchMind"] = precision
  results_FN[iterSimu, "MatchMind"] = falsenegative
  results_FP[iterSimu, "MatchMind"] = falsepositive
  results_TP[iterSimu, "MatchMind"] = truepositive
  results_distance[iterSimu, "MatchMind"] = sqrt(sum((fit$Delta - Delta)**2))
  
  fit_gamma = (fit$gamma)[(MatchMindNcutBurninStEM+1):MatchMindNIter , drop=FALSE]
  results_gamma[iterSimu,] = fit_gamma # 1 vector of length iter post burnin
  
  fit_omega = lapply(fit$omega, function(x) x[(MatchMindNcutBurninStEM+1):MatchMindNIter, , drop=FALSE])
  unstable_piv = which(!pivs_stable)
  omega_avg = lapply(fit_omega, function(x) apply(x, 2, mean))
  links = which(fit$Delta>0.5, arr.ind=TRUE)
  times = abs(dataSimu$regTimeB[links[,2]] - dataSimu$regTimeA[links[,1]])
  Xomega = cbind( times, dataSimu$XA[links[,1], dataSimu$pSameH.varA[[unstable_piv]], drop=FALSE], dataSimu$XB[links[,2], dataSimu$pSameH.varB[[unstable_piv]], drop=FALSE] )
  
  
  results_omega_param[iterSimu,] = fit_omega[[unstable_piv]] # 1 vector of length: iter post burnin
  omega_probaEstimate = exp( - as.matrix(Xomega) %*% omega_avg[[unstable_piv]] )
  results_omega_probaEstimate[iterSimu,] = append( omega_probaEstimate, rep(0, NdataA - length(omega_probaEstimate)) ) # 1 vector of length linked records + filled with 0
  results_omega_timesEstimate[iterSimu,] = append( times, rep(0, NdataA - length(omega_probaEstimate)) ) # 1 vector of length linked records + filled with 0
  results_omega_probaTrue[iterSimu,] = proba_same_H # 1 vector of length Nlinks
  results_omega_timesTrue[iterSimu,] = time_difference # 1 vector of length Nlinks
  # plot( sort(proba_same_H), ylim=c(0,1), xlab = "true records pairs", ylab = sprintf("true model for proba of no change in %s", pivs[5]) )
  # plot( time_difference, proba_same_H, ylim=c(0,1), xlab = "time difference", ylab = sprintf("true model for proba of no change in %s", pivs[5]) )
  
  ###
  pivs_stable  = c(TRUE,TRUE,TRUE,TRUE,TRUE)
  dataSimu = list( pivsA        = encodedA[, pivs], 
                   pivsB        = encodedB[, pivs], 
                   nvalues      = nvalues,
                   pivs_stable  = pivs_stable, 
                   regTimeA     = regTimeA, 
                   regTimeB     = regTimeB,
                   XA           = XA,
                   XB           = XB,
                   pSameH.varA  = pSameH.varA,
                   pSameH.varB  = pSameH.varB,
                   include_time = TRUE)
  
  fit = stEM(  data          = dataSimu, 
               nIter         = MatchMindNIter, 
               nBurnin       = MatchMindNBurnin, 
               MStepIter     = MatchMindNMStepIter,
               cutBurninStEM = MatchMindNcutBurninStEM,
               trace      = 1 )
  
  linkedpairs = which(fit$Delta>0.5, arr.ind=TRUE)
  linkedpairsA = linkedpairs[,1]
  linkedpairsB = linkedpairs[,2]
  simuMatchMindstable_NAA_linked[iterSimu,] = ( colSums(is.na(A[linkedpairsA,])) / length(linkedpairsA) )[pivs]
  simuMatchMindstable_NAB_linked[iterSimu,] = ( colSums(is.na(B[linkedpairsB,])) / length(linkedpairsB) )[pivs]
  
  simuMatchMind_agree_linked_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuMatchMind_agree_linked_tmp) = pivs
  simuMatchMind_unstable_change_linked_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuMatchMind_unstable_change_linked_tmp) = pivs
  for( i in 1:nrow(linkedpairs) ){
    entityA = A[linkedpairsA[i],]
    entityB = B[linkedpairsB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        simuMatchMind_agree_linked_tmp = rbind(simuMatchMind_agree_linked_tmp, entityA[,pivs] == entityB[,pivs])
        simuMatchMind_unstable_change_linked_tmp = rbind(simuMatchMind_unstable_change_linked_tmp, entityB[,"change"]*(!pivs_stable))
      }
    }
  }
  simuMatchMind_agree_linked_tmp = colSums( simuMatchMind_agree_linked_tmp, na.rm = TRUE ) / nrow( simuMatchMind_agree_linked_tmp )
  simuMatchMind_unstable_change_linked_tmp = colSums( simuMatchMind_unstable_change_linked_tmp, na.rm = TRUE ) / nrow( simuMatchMind_unstable_change_linked_tmp )
  simuMatchMindstable_agree_linked[iterSimu,] = simuMatchMind_agree_linked_tmp
  simuMatchMindstable_unstable_change_linked[iterSimu,] = simuMatchMind_unstable_change_linked_tmp
  
  linkedTP = which((fit$Delta>0.5) & (Delta==1), arr.ind=TRUE)
  linkedTPA = linkedTP[,1]
  linkedTPB = linkedTP[,2]
  simuMatchMindstable_NAA_TP[iterSimu,] = ( colSums(is.na(A[linkedTPA,])) / length(linkedTPA) )[pivs]
  simuMatchMindstable_NAB_TP[iterSimu,] = ( colSums(is.na(B[linkedTPB,])) / length(linkedTPB) )[pivs]
  
  simuMatchMind_agree_TP_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuMatchMind_agree_TP_tmp) = pivs
  simuMatchMind_unstable_change_TP_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuMatchMind_unstable_change_TP_tmp) = pivs
  for( i in 1:nrow(linkedTP) ){
    entityA = A[linkedTPA[i],]
    entityB = B[linkedTPB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        simuMatchMind_agree_TP_tmp = rbind(simuMatchMind_agree_TP_tmp, entityA[,pivs] == entityB[,pivs])
        simuMatchMind_unstable_change_TP_tmp = rbind(simuMatchMind_unstable_change_TP_tmp, entityB[,"change"]*(!pivs_stable))
      }
    }
  }
  simuMatchMind_agree_TP_tmp = colSums( simuMatchMind_agree_TP_tmp, na.rm = TRUE ) / nrow( simuMatchMind_agree_TP_tmp )
  simuMatchMind_unstable_change_TP_tmp = colSums( simuMatchMind_unstable_change_TP_tmp, na.rm = TRUE ) / nrow( simuMatchMind_unstable_change_TP_tmp )
  simuMatchMindstable_agree_TP[iterSimu,] = simuMatchMind_agree_TP_tmp
  simuMatchMindstable_unstable_change_TP[iterSimu,] = simuMatchMind_unstable_change_TP_tmp
  
  linkedFP = which((fit$Delta>0.5) & (Delta==0), arr.ind=TRUE)
  linkedFPA = linkedFP[,1]
  linkedFPB = linkedFP[,2]
  simuMatchMindstable_NAA_FP[iterSimu,] = ( colSums(is.na(A[linkedFPA,])) / length(linkedFPA) )[pivs]
  simuMatchMindstable_NAB_FP[iterSimu,] = ( colSums(is.na(B[linkedFPB,])) / length(linkedFPB) )[pivs]
  
  simuMatchMind_agree_FP_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuMatchMind_agree_FP_tmp) = pivs
  simuMatchMind_unstable_change_FP_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuMatchMind_unstable_change_FP_tmp) = pivs
  for( i in 1:nrow(linkedFP) ){
    entityA = A[linkedFPA[i],]
    entityB = B[linkedFPB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        simuMatchMind_agree_FP_tmp = rbind(simuMatchMind_agree_FP_tmp, entityA[,pivs] == entityB[,pivs])
        simuMatchMind_unstable_change_FP_tmp = rbind(simuMatchMind_unstable_change_FP_tmp, entityB[,"change"]*(!pivs_stable))
      }
    }
  }
  simuMatchMind_agree_FP_tmp = colSums( simuMatchMind_agree_FP_tmp, na.rm = TRUE ) / nrow( simuMatchMind_agree_FP_tmp )
  simuMatchMind_unstable_change_FP_tmp = colSums( simuMatchMind_unstable_change_FP_tmp, na.rm = TRUE ) / nrow( simuMatchMind_unstable_change_FP_tmp )
  simuMatchMindstable_agree_FP[iterSimu,] = simuMatchMind_agree_FP_tmp
  simuMatchMindstable_unstable_change_FP[iterSimu,] = simuMatchMind_unstable_change_FP_tmp
  
  linkedFN = which((fit$Delta<0.5) & (Delta==1), arr.ind=TRUE)
  linkedFNA = linkedFN[,1]
  linkedFNB = linkedFN[,2]
  simuMatchMindstable_NAA_FN[iterSimu,] = ( colSums(is.na(A[linkedFNA,])) / length(linkedFNA) )[pivs]
  simuMatchMindstable_NAB_FN[iterSimu,] = ( colSums(is.na(B[linkedFNB,])) / length(linkedFNB) )[pivs]
  
  simuMatchMind_agree_FN_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuMatchMind_agree_FN_tmp) = pivs
  simuMatchMind_unstable_change_FN_tmp = data.frame( matrix(0, nrow=0, ncol=length(pivs)) )
  colnames(simuMatchMind_unstable_change_FN_tmp) = pivs
  for( i in 1:nrow(linkedFN) ){
    entityA = A[linkedFNA[i],]
    entityB = B[linkedFNB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        simuMatchMind_agree_FN_tmp = rbind(simuMatchMind_agree_FN_tmp, entityA[,pivs] == entityB[,pivs])
        simuMatchMind_unstable_change_FN_tmp = rbind(simuMatchMind_unstable_change_FN_tmp, entityB[,"change"]*(!pivs_stable))
      }
    }
  }
  simuMatchMind_agree_FN_tmp = colSums( simuMatchMind_agree_FN_tmp, na.rm = TRUE ) / nrow( simuMatchMind_agree_FN_tmp )
  simuMatchMind_unstable_change_FN_tmp = colSums( simuMatchMind_unstable_change_FN_tmp, na.rm = TRUE ) / nrow( simuMatchMind_unstable_change_FN_tmp )
  simuMatchMindstable_agree_FN[iterSimu,] = simuMatchMind_agree_FN_tmp
  simuMatchMindstable_unstable_change_FN[iterSimu,] = simuMatchMind_unstable_change_FN_tmp
  
  truepositive = sum( (fit$Delta>0.5) & (Delta==1) )
  falsepositive = sum( (fit$Delta>0.5) & (Delta==0) )
  falsenegative = sum( (fit$Delta<0.5) & (Delta==1) )
  precision = truepositive / (truepositive + falsepositive)
  recall = truepositive / (truepositive + falsenegative)
  f1score = 2 * (precision * recall) / (precision + recall)
  
  results_f1score[iterSimu, "MatchMindstable"] = f1score
  results_recall[iterSimu, "MatchMindstable"] = recall
  results_precision[iterSimu, "MatchMindstable"] = precision
  results_FN[iterSimu, "MatchMindstable"] = falsenegative
  results_FP[iterSimu, "MatchMindstable"] = falsepositive
  results_TP[iterSimu, "MatchMindstable"] = truepositive
  results_distance[iterSimu, "MatchMindstable"] = sqrt(sum((fit$Delta - Delta)**2))
  
  fit_gamma = (fit$gamma)[(MatchMindNcutBurninStEM+1):MatchMindNIter , drop=FALSE]
  resultsstable_gamma[iterSimu,] = fit_gamma # 1 vector of length iter post burnin
  
  # no unstability here
  
  write.csv(simu_NA_A, "datasetSimu_NA_A.csv")
  write.csv(simu_NA_B, "datasetSimu_NA_B.csv")
  write.csv(simu_uniquevalues, "datasetSimu_uniquevalues.csv")
  write.csv(simu_truelinks_agree, "datasetSimu_truelinks_agree.csv")
  write.csv(simu_truelinks_unstable_change, "datasetSimu_truelinks_unstable_change.csv")
  
  write.csv(simuExchanger_NAA_TP, "datasetSimuExchanger_NAA_TP.csv")
  write.csv(simuExchanger_NAA_FP, "datasetSimuExchanger_NAA_FP.csv")
  write.csv(simuExchanger_NAA_FN, "datasetSimuExchanger_NAA_FN.csv")
  write.csv(simuExchanger_NAA_linked, "datasetSimuExchanger_NAA_linked.csv")
  write.csv(simuExchanger_NAB_TP, "datasetSimuExchanger_NAB_TP.csv")
  write.csv(simuExchanger_NAB_FP, "datasetSimuExchanger_NAB_FP.csv")
  write.csv(simuExchanger_NAB_FN, "datasetSimuExchanger_NAB_FN.csv")
  write.csv(simuExchanger_NAB_linked, "datasetSimuExchanger_NAB_linked.csv")
  write.csv(simuExchanger_agree_TP, "datasetSimuExchanger_agree_TP.csv")
  write.csv(simuExchanger_agree_FP, "datasetSimuExchanger_agree_FP.csv")
  write.csv(simuExchanger_agree_FN, "datasetSimuExchanger_agree_FN.csv")
  write.csv(simuExchanger_agree_linked, "datasetSimuExchanger_agree_linked.csv")
  write.csv(simuExchanger_unstable_change_TP, "datasetSimuExchanger_unstable_change_TP.csv")
  write.csv(simuExchanger_unstable_change_FP, "datasetSimuExchanger_unstable_change_FP.csv")
  write.csv(simuExchanger_unstable_change_FN, "datasetSimuExchanger_unstable_change_FN.csv")
  write.csv(simuExchanger_unstable_change_linked, "datasetSimuExchanger_unstable_change_linked.csv")
  
  write.csv(simuBRL_NAA_TP, "datasetSimuBRL_NAA_TP.csv")
  write.csv(simuBRL_NAA_FP, "datasetSimuBRL_NAA_FP.csv")
  write.csv(simuBRL_NAA_FN, "datasetSimuBRL_NAA_FN.csv")
  write.csv(simuBRL_NAA_linked, "datasetSimuBRL_NAA_linked.csv")
  write.csv(simuBRL_NAB_TP, "datasetSimuBRL_NAB_TP.csv")
  write.csv(simuBRL_NAB_FP, "datasetSimuBRL_NAB_FP.csv")
  write.csv(simuBRL_NAB_FN, "datasetSimuBRL_NAB_FN.csv")
  write.csv(simuBRL_NAB_linked, "datasetSimuBRL_NAB_linked.csv")
  write.csv(simuBRL_agree_TP, "datasetSimuBRL_agree_TP.csv")
  write.csv(simuBRL_agree_FP, "datasetSimuBRL_agree_FP.csv")
  write.csv(simuBRL_agree_FN, "datasetSimuBRL_agree_FN.csv")
  write.csv(simuBRL_agree_linked, "datasetSimuBRL_agree_linked.csv")
  write.csv(simuBRL_unstable_change_TP, "datasetSimuBRL_unstable_change_TP.csv")
  write.csv(simuBRL_unstable_change_FP, "datasetSimuBRL_unstable_change_FP.csv")
  write.csv(simuBRL_unstable_change_FN, "datasetSimuBRL_unstable_change_FN.csv")
  write.csv(simuBRL_unstable_change_linked, "datasetSimuBRL_unstable_change_linked.csv")
  
  write.csv(simuMatchMind_NAA_TP, "datasetSimuMatchMind_NAA_TP.csv")
  write.csv(simuMatchMind_NAA_FP, "datasetSimuMatchMind_NAA_FP.csv")
  write.csv(simuMatchMind_NAA_FN, "datasetSimuMatchMind_NAA_FN.csv")
  write.csv(simuMatchMind_NAA_linked, "datasetSimuMatchMind_NAA_linked.csv")
  write.csv(simuMatchMind_NAB_TP, "datasetSimuMatchMind_NAB_TP.csv")
  write.csv(simuMatchMind_NAB_FP, "datasetSimuMatchMind_NAB_FP.csv")
  write.csv(simuMatchMind_NAB_FN, "datasetSimuMatchMind_NAB_FN.csv")
  write.csv(simuMatchMind_NAB_linked, "datasetSimuMatchMind_NAB_linked.csv")
  write.csv(simuMatchMind_agree_TP, "datasetSimuMatchMind_agree_TP.csv")
  write.csv(simuMatchMind_agree_FP, "datasetSimuMatchMind_agree_FP.csv")
  write.csv(simuMatchMind_agree_FN, "datasetSimuMatchMind_agree_FN.csv")
  write.csv(simuMatchMind_agree_linked, "datasetSimuMatchMind_agree_linked.csv")
  write.csv(simuMatchMind_unstable_change_TP, "datasetSimuMatchMind_unstable_change_TP.csv")
  write.csv(simuMatchMind_unstable_change_FP, "datasetSimuMatchMind_unstable_change_FP.csv")
  write.csv(simuMatchMind_unstable_change_FN, "datasetSimuMatchMind_unstable_change_FN.csv")
  write.csv(simuMatchMind_unstable_change_linked, "datasetSimuMatchMind_unstable_change_linked.csv")
  
  write.csv(simuMatchMindstable_NAA_TP, "datasetSimuMatchMindstable_NAA_TP.csv")
  write.csv(simuMatchMindstable_NAA_FP, "datasetSimuMatchMindstable_NAA_FP.csv")
  write.csv(simuMatchMindstable_NAA_FN, "datasetSimuMatchMindstable_NAA_FN.csv")
  write.csv(simuMatchMindstable_NAA_linked, "datasetSimuMatchMindstable_NAA_linked.csv")
  write.csv(simuMatchMindstable_NAB_TP, "datasetSimuMatchMindstable_NAB_TP.csv")
  write.csv(simuMatchMindstable_NAB_FP, "datasetSimuMatchMindstable_NAB_FP.csv")
  write.csv(simuMatchMindstable_NAB_FN, "datasetSimuMatchMindstable_NAB_FN.csv")
  write.csv(simuMatchMindstable_NAB_linked, "datasetSimuMatchMindstable_NAB_linked.csv")
  write.csv(simuMatchMindstable_agree_TP, "datasetSimuMatchMindstable_agree_TP.csv")
  write.csv(simuMatchMindstable_agree_FP, "datasetSimuMatchMindstable_agree_FP.csv")
  write.csv(simuMatchMindstable_agree_FN, "datasetSimuMatchMindstable_agree_FN.csv")
  write.csv(simuMatchMindstable_agree_linked, "datasetSimuMatchMindstable_agree_linked.csv")
  write.csv(simuMatchMindstable_unstable_change_TP, "datasetSimuMatchMindstable_unstable_change_TP.csv")
  write.csv(simuMatchMindstable_unstable_change_FP, "datasetSimuMatchMindstable_unstable_change_FP.csv")
  write.csv(simuMatchMindstable_unstable_change_FN, "datasetSimuMatchMindstable_unstable_change_FN.csv")
  write.csv(simuMatchMindstable_unstable_change_linked, "datasetSimuMatchMindstable_unstable_change_linked.csv")
  
  write.csv(results_f1score, "datasetSimu_results_f1score.csv")
  write.csv(results_recall, "datasetSimu_results_recall.csv")
  write.csv(results_precision, "datasetSimu_results_precision.csv")
  write.csv(results_FN, "datasetSimu_results_FN.csv")
  write.csv(results_FP, "datasetSimu_results_FP.csv")
  write.csv(results_TP, "datasetSimu_results_TP.csv")
  write.csv(results_distance, "datasetSimu_results_distance.csv")
  write.csv(results_exchangerdedup, "datasetSimu_results_exchangerdedup.csv")
  
  write.csv(results_gamma, "datasetSimu_results_gamma.csv")
  write.csv(results_omega_param, "datasetSimu_results_omega_param.csv")
  write.csv(results_omega_probaEstimate, "datasetSimu_results_omega_probaEstimate.csv")
  write.csv(results_omega_timesEstimate, "datasetSimu_results_omega_timesEstimate.csv")
  write.csv(results_omega_probaTrue, "datasetSimu_results_omega_probaTrue.csv")
  write.csv(results_omega_timesTrue, "datasetSimu_results_omega_timesTrue.csv")
  
  write.csv(resultsstable_gamma, "datasetSimu_resultsstable_gamma.csv")
  
}

pdf("datasetSimu_eta_foronce.pdf")
par(mfrow=c(2,3))
# eta
eta = fit$eta
for(k in 1:length(eta))
{
  vec = c("eta 1")
  plot(eta[[k]][,1], ylim=c(0,1), type="l", col=1, xlab = "MC-EM Iterations", ylab = sprintf("PIV %s", pivs[k]))
  for(j in 2:ncol(eta[[k]]))
  {
    lines(eta[[k]][,j], ylim=c(0,1), col=j)
    vec <- append(vec, sprintf("eta %s", j))
  }
  legend("topright", legend=vec, col=seq_len(ncol(eta[[k]])), lty=1:4, cex=0.8)
}
dev.off()

pdf("datasetSimu_phi_foronce.pdf")
par(mfrow=c(2,3))
# phi
phi = fit$phi
for(k in 1:length(phi))
{
  plot(phi[[k]][,1], ylim=c(0,1), type="l", col=1, xlab = "MC-EM Iterations", ylab = sprintf("PIV %s", pivs[k]), lty=1)
  lines(phi[[k]][,2], col=2, lty=2)
  lines(phi[[k]][,3], col=3, lty=3)
  legend("left", legend=c("agree", "missings A", "missings B"), col=c(1,2,3), lty=1:3, cex=0.8)
}
dev.off()