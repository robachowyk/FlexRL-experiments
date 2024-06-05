################################################################################
#                 RUN SIMULATIONS FOR THE RECORD LINKAGE TASK                  #
#         AND COMPARE FLEXRL WITH BRL, EXCHANGER AND A NAIVE APPROACH          #
################################################################################

################################################################################
#                                 IMPORTS                                      #
################################################################################
library(exchanger)
library(comparator)
library(clevr)
library(progress)
library(Matrix)
library(testit)
library("BRL")
source("../recordlinkage.r")
Rcpp:::sourceCpp("../functions.cpp")

################################################################################
#                          FIX AMOUNT OF ITERATIONS                            #
################################################################################
Nsimu = 2

ExchangerNSamples = 20000
ExchangerNBurnin = 10000

FlexRLUnstableStEMIter     = 100
FlexRLUnstableStEMBurnin   = 75
FlexRLUnstableGibbsIter    = 100
FlexRLUnstableGibbsBurnin  = 50

FlexRLStableStEMIter     = 100
FlexRLStableStEMBurnin   = 75
FlexRLStableGibbsIter    = 100
FlexRLStableGibbsBurnin  = 50

# BRL by default is 1000 samples and 10% = 100 burn-in

################################################################################
#                                SIMULATE DATA                                 #
################################################################################
NdataA = 800
NdataB = 1000
Nlinks = 500

Delta = matrix(0, nrow=NdataA, ncol=NdataB)
for (l in 1:Nlinks)
{
  Delta[l,l]=1
}

Nval       = c(      6,       7,       8,       9,      15)
Pmistakes  = c(   0.02,    0.02,    0.02,    0.02,    0.00)
PmissingA  = c(  0.007,   0.007,   0.007,   0.007,   0.007)
PmissingB  = c(  0.007,   0.007,   0.007,   0.007,   0.007)

NRecords = c(NdataA, NdataB)

PIVs_config = list( V1 = list(stable = TRUE),
                    V2 = list(stable = TRUE),
                    V3 = list(stable = TRUE),
                    V4 = list(stable = TRUE),
                    V5 = list(stable = FALSE, conditionalHazard = FALSE, pSameH.cov.A = c(), pSameH.cov.B = c()) )

PIVs = names(PIVs_config)

PIVs_stable = sapply(PIVs_config, function(x) x$stable)

newDirectory = sprintf("Simulations %s", Sys.time())
dir.create(newDirectory)

################################################################################
#                       CREATE TABLES TO SAVE RESULTS                          #
################################################################################
methods = c("Exchanger", "BRL", "FlexRL_instability4V5", "FlexRL_noinstability", "Naive")

Nmethods                                                = length(methods)
results_f1score                                         = data.frame( matrix(0, nrow=Nsimu, ncol=Nmethods) )
colnames(results_f1score)                               = methods
results_recall                                          = data.frame( matrix(0, nrow=Nsimu, ncol=Nmethods) )
colnames(results_recall)                                = methods
results_precision                                       = data.frame( matrix(0, nrow=Nsimu, ncol=Nmethods) )
colnames(results_precision)                             = methods
results_FN                                              = data.frame( matrix(0, nrow=Nsimu, ncol=Nmethods) )
colnames(results_FN)                                    = methods
results_FP                                              = data.frame( matrix(0, nrow=Nsimu, ncol=Nmethods) )
colnames(results_FP)                                    = methods
results_TP                                              = data.frame( matrix(0, nrow=Nsimu, ncol=Nmethods) )
colnames(results_TP)                                    = methods
results_MatrixDistance                                  = data.frame( matrix(0, nrow=Nsimu, ncol=Nmethods) )
colnames(results_MatrixDistance)                        = methods

# overall data generated
simu_NA_A                                               = data.frame( matrix(0, nrow=Nsimu, ncol=length(PIVs)) )
colnames(simu_NA_A)                                     = PIVs
simu_NA_B                                               = data.frame( matrix(0, nrow=Nsimu, ncol=length(PIVs)) )
colnames(simu_NA_B)                                     = PIVs
simu_truelinks_NA_A                                     = data.frame( matrix(0, nrow=Nsimu, ncol=length(PIVs)) )
colnames(simu_truelinks_NA_A)                           = PIVs
simu_truelinks_NA_B                                     = data.frame( matrix(0, nrow=Nsimu, ncol=length(PIVs)) )
colnames(simu_truelinks_NA_B)                           = PIVs
simu_uniquevalues                                       = data.frame( matrix(0, nrow=Nsimu, ncol=length(PIVs)) )
colnames(simu_uniquevalues)                             = PIVs
simu_truelinks_agree                                    = data.frame( matrix(0, nrow=Nsimu, ncol=length(PIVs)) )
colnames(simu_truelinks_agree)                          = PIVs
simu_truelinks_unstable_change                          = data.frame( matrix(0, nrow=Nsimu, ncol=length(PIVs)) )
colnames(simu_truelinks_unstable_change)                = PIVs

results_unstable_gamma                                  = data.frame( matrix(0, nrow=Nsimu, ncol=FlexRLUnstableStEMIter) )# -FlexRLUnstableStEMBurnin) )
results_unstable_phi_agree_V1                           = data.frame( matrix(0, nrow=Nsimu, ncol=FlexRLUnstableStEMIter) )# -FlexRLUnstableStEMBurnin) )
results_unstable_phi_agree_V2                           = data.frame( matrix(0, nrow=Nsimu, ncol=FlexRLUnstableStEMIter) )# -FlexRLUnstableStEMBurnin) )
results_unstable_phi_agree_V3                           = data.frame( matrix(0, nrow=Nsimu, ncol=FlexRLUnstableStEMIter) )# -FlexRLUnstableStEMBurnin) )
results_unstable_phi_agree_V4                           = data.frame( matrix(0, nrow=Nsimu, ncol=FlexRLUnstableStEMIter) )# -FlexRLUnstableStEMBurnin) )
results_unstable_phi_agree_V5                           = data.frame( matrix(0, nrow=Nsimu, ncol=FlexRLUnstableStEMIter) )# -FlexRLUnstableStEMBurnin) )
results_unstable_alpha_param                            = data.frame( matrix(0, nrow=Nsimu, ncol=FlexRLUnstableStEMIter) )# -FlexRLUnstableStEMBurnin) )
results_unstable_alpha_probaEstimate                    = data.frame( matrix(0, nrow=Nsimu, ncol=NdataA) )
results_unstable_alpha_timesEstimate                    = data.frame( matrix(0, nrow=Nsimu, ncol=NdataA) )
results_unstable_alpha_probaTrue                        = data.frame( matrix(0, nrow=Nsimu, ncol=Nlinks) )
results_unstable_alpha_timesTrue                        = data.frame( matrix(0, nrow=Nsimu, ncol=Nlinks) )

results_stable_gamma                                    = data.frame( matrix(0, nrow=Nsimu, ncol=FlexRLStableStEMIter) )# -FlexRLStableStEMBurnin) )
results_stable_phi_agree_V1                             = data.frame( matrix(0, nrow=Nsimu, ncol=FlexRLStableStEMIter) )# -FlexRLStableStEMBurnin) )
results_stable_phi_agree_V2                             = data.frame( matrix(0, nrow=Nsimu, ncol=FlexRLStableStEMIter) )# -FlexRLStableStEMBurnin) )
results_stable_phi_agree_V3                             = data.frame( matrix(0, nrow=Nsimu, ncol=FlexRLStableStEMIter) )# -FlexRLStableStEMBurnin) )
results_stable_phi_agree_V4                             = data.frame( matrix(0, nrow=Nsimu, ncol=FlexRLStableStEMIter) )# -FlexRLStableStEMBurnin) )
results_stable_phi_agree_V5                             = data.frame( matrix(0, nrow=Nsimu, ncol=FlexRLStableStEMIter) )# -FlexRLStableStEMBurnin) )

for(i in 1:Nmethods){
  method = methods[i]
  
  # NA
  simu_NAA_TP                                           = data.frame( matrix(0, nrow=Nsimu, ncol=length(PIVs)) )
  colnames(simu_NAA_TP)                                 = PIVs
  simu_NAA_FP                                           = data.frame( matrix(0, nrow=Nsimu, ncol=length(PIVs)) )
  colnames(simu_NAA_FP)                                 = PIVs
  simu_NAA_FN                                           = data.frame( matrix(0, nrow=Nsimu, ncol=length(PIVs)) )
  colnames(simu_NAA_FN)                                 = PIVs
  simu_NAA_linked                                       = data.frame( matrix(0, nrow=Nsimu, ncol=length(PIVs)) )
  colnames(simu_NAA_linked)                             = PIVs
  simu_NAB_TP                                           = data.frame( matrix(0, nrow=Nsimu, ncol=length(PIVs)) )
  colnames(simu_NAB_TP)                                 = PIVs
  simu_NAB_FP                                           = data.frame( matrix(0, nrow=Nsimu, ncol=length(PIVs)) )
  colnames(simu_NAB_FP)                                 = PIVs
  simu_NAB_FN                                           = data.frame( matrix(0, nrow=Nsimu, ncol=length(PIVs)) )
  colnames(simu_NAB_FN)                                 = PIVs
  simu_NAB_linked                                       = data.frame( matrix(0, nrow=Nsimu, ncol=length(PIVs)) )
  colnames(simu_NAB_linked)                             = PIVs
  if(method == "Exchanger"){
    assign( paste(method, "simu_NARL_TP",                  sep="_"), simu_NAA_TP )
    assign( paste(method, "simu_NARL_FP",                  sep="_"), simu_NAA_FP )
    assign( paste(method, "simu_NARL_FN",                  sep="_"), simu_NAA_FN )
    assign( paste(method, "simu_NARL_linked",              sep="_"), simu_NAA_linked )
  }else{
    assign( paste(method, "simu_NAA_TP",                  sep="_"), simu_NAA_TP )
    assign( paste(method, "simu_NAA_FP",                  sep="_"), simu_NAA_FP )
    assign( paste(method, "simu_NAA_FN",                  sep="_"), simu_NAA_FN )
    assign( paste(method, "simu_NAA_linked",              sep="_"), simu_NAA_linked )
    assign( paste(method, "simu_NAB_TP",                  sep="_"), simu_NAB_TP )
    assign( paste(method, "simu_NAB_FP",                  sep="_"), simu_NAB_FP )
    assign( paste(method, "simu_NAB_FN",                  sep="_"), simu_NAB_FN )
    assign( paste(method, "simu_NAB_linked",              sep="_"), simu_NAB_linked )
  }
  
  # AGREEMENTS
  simu_agree_TP                                         = data.frame( matrix(0, nrow=Nsimu, ncol=length(PIVs)) )
  colnames(simu_agree_TP)                               = PIVs
  simu_agree_FP                                         = data.frame( matrix(0, nrow=Nsimu, ncol=length(PIVs)) )
  colnames(simu_agree_FP)                               = PIVs
  simu_agree_FN                                         = data.frame( matrix(0, nrow=Nsimu, ncol=length(PIVs)) )
  colnames(simu_agree_FN)                               = PIVs
  simu_agree_linked                                     = data.frame( matrix(0, nrow=Nsimu, ncol=length(PIVs)) )
  colnames(simu_agree_linked)                           = PIVs
  assign( paste(method, "simu_agree_TP",                sep="_"), simu_agree_TP )
  assign( paste(method, "simu_agree_FP",                sep="_"), simu_agree_FP )
  assign( paste(method, "simu_agree_FN",                sep="_"), simu_agree_FN )
  assign( paste(method, "simu_agree_linked",            sep="_"), simu_agree_linked )
  
  # CHANGES
  simu_unstable_change_TP                               = data.frame( matrix(0, nrow=Nsimu, ncol=length(PIVs)) )
  colnames(simu_unstable_change_TP)                     = PIVs
  simu_unstable_change_FP                               = data.frame( matrix(0, nrow=Nsimu, ncol=length(PIVs)) )
  colnames(simu_unstable_change_FP)                     = PIVs
  simu_unstable_change_FN                               = data.frame( matrix(0, nrow=Nsimu, ncol=length(PIVs)) )
  colnames(simu_unstable_change_FN)                     = PIVs
  simu_unstable_change_linked                           = data.frame( matrix(0, nrow=Nsimu, ncol=length(PIVs)) )
  colnames(simu_unstable_change_linked)                 = PIVs
  assign( paste(method, "simu_unstable_change_TP",      sep="_"), simu_unstable_change_TP )
  assign( paste(method, "simu_unstable_change_FP",      sep="_"), simu_unstable_change_FP )
  assign( paste(method, "simu_unstable_change_FN",      sep="_"), simu_unstable_change_FN )
  assign( paste(method, "simu_unstable_change_linked",  sep="_"), simu_unstable_change_linked )
}

################################################################################
#                              LAUNCH SIMULATIONS                              #
################################################################################
for(iter in 1:Nsimu){
  
  PIVs_config = list( V1 = list(stable = TRUE),
                      V2 = list(stable = TRUE),
                      V3 = list(stable = TRUE),
                      V4 = list(stable = TRUE),
                      V5 = list(stable = FALSE, conditionalHazard = FALSE, pSameH.cov.A = c(), pSameH.cov.B = c()) )
  PIVs = names(PIVs_config)
  PIVs_stable = sapply(PIVs_config, function(x) x$stable)
  Data              = DataCreation(PIVs_config, PIVs, PIVs_stable, Nval, NRecords, Nlinks, Pmistakes, PmissingA, PmissingB, UnstableNoMistake=TRUE)
  A                 = Data$A
  B                 = Data$B
  Nvalues           = Data$Nvalues
  TimeDifference    = Data$TimeDifference
  proba_same_H      = Data$proba_same_H
  
  # PREPARE DATA FOR EXCHANGER
  RLdata            = rbind(A[,PIVs], B[,PIVs])
  rownames(RLdata)  = 1:nrow(RLdata)
  RLdata            = rbind(A, B)
  rownames(RLdata)  = 1:nrow(RLdata)
  true_pairs        = cbind( rownames(RLdata[(RLdata$source=="A")&(RLdata$localID %in% 1:Nlinks),]),
                             rownames(RLdata[(RLdata$source=="B")&(RLdata$localID %in% 1:Nlinks),]) )
  
  # STORY TELLING ABOUT THE DATA
  simu_NA_A[iter,] = ( colSums(is.na(A)) / nrow(A) )[PIVs]
  simu_NA_B[iter,] = ( colSums(is.na(B)) / nrow(B) )[PIVs]
  simu_truelinks_NA_A[iter,] = ( colSums(is.na(A[1:Nlinks,])) / Nlinks )[PIVs]
  simu_truelinks_NA_B[iter,] = ( colSums(is.na(B[1:Nlinks,])) / Nlinks )[PIVs]
  simu_uniquevalues[iter,] = Nvalues
  simu_truelinks_agree_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(simu_truelinks_agree_tmp) = PIVs
  simu_truelinks_unstable_change_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(simu_truelinks_unstable_change_tmp) = PIVs
  for( i in 1:Nlinks ){
    entityA = A[i,]
    entityB = B[i,]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        simu_truelinks_agree_tmp = rbind(simu_truelinks_agree_tmp, entityA[,PIVs] == entityB[,PIVs])
        simu_truelinks_unstable_change_tmp = rbind(simu_truelinks_unstable_change_tmp, entityB[,"change"])
      }
    }
  }
  simu_truelinks_agree_tmp = colSums( simu_truelinks_agree_tmp, na.rm = TRUE ) / nrow( simu_truelinks_agree_tmp )
  simu_truelinks_unstable_change_tmp = colSums( simu_truelinks_unstable_change_tmp, na.rm = TRUE ) / nrow( simu_truelinks_unstable_change_tmp )
  simu_truelinks_agree[iter,] = simu_truelinks_agree_tmp
  simu_truelinks_unstable_change[iter,] = simu_truelinks_unstable_change_tmp
  
  ### LAUNCH EXCHANGER
  distort_prior <- BetaRV(1, 4)
  distort_prior_0 <- BetaRV(1, 100)
  attr_params <- list(
                        V1 = CategoricalAttribute(distort_prob_prior = distort_prior),
                        V2 = CategoricalAttribute(distort_prob_prior = distort_prior),
                        V3 = CategoricalAttribute(distort_prob_prior = distort_prior),
                        V4 = CategoricalAttribute(distort_prob_prior = distort_prior),
                        V5 = CategoricalAttribute(distort_prob_prior = distort_prior_0) )
  clust_prior <- PitmanYorRP(alpha = GammaRV(1, 0.01), d = BetaRV(1, 1))
  model <- exchanger(RLdata, attr_params, clust_prior)
  result <- run_inference(model, n_samples=ExchangerNSamples, thin_interval=10, burnin_interval=ExchangerNBurnin)
  pred_clust <- smp_clusters(result)
  n_records <- nrow(RLdata)
  pred_pairs <- clusters_to_pairs(pred_clust)
  comb_pairs <- rbind(true_pairs, pred_pairs)
  true_pairs <- comb_pairs[seq_len(nrow(true_pairs)),]
  pred_pairs <- comb_pairs[nrow(true_pairs) + seq_len(nrow(pred_pairs)),]
  df_pred_pairs <- as.data.frame(canonicalize_pairs(pred_pairs, ordered = FALSE))
  df_true_pairs <- as.data.frame(canonicalize_pairs(true_pairs, ordered = FALSE))
  df_true_pairs$match = rep(TRUE, times=nrow(df_true_pairs))
  df_pred_pairs$pred_match = rep(TRUE, times=nrow(df_pred_pairs))
  merged_pairs = merge(df_true_pairs, df_pred_pairs, by=c("V1", "V2"), all=TRUE)
  merged_pairs[is.na(merged_pairs)] = FALSE
  prediction = factor(merged_pairs$pred_match, levels = c(TRUE, FALSE))
  truth = factor(merged_pairs$match, levels = c(TRUE, FALSE))
  CT = table(prediction, truth, dnn = c("Prediction", "Truth"))
  ### SAVE RESULTS
  # STORY TELLING
  # LINKED
  Exchanger_simu_NARL_linked[iter,] = ( ( colSums(is.na(RLdata[df_pred_pairs[,1],])) + colSums(is.na(RLdata[df_pred_pairs[,2],])) ) / nrow(df_pred_pairs) )[PIVs]
  Exchanger_simu_agree_linked_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(Exchanger_simu_agree_linked_tmp) = PIVs
  Exchanger_simu_unstable_change_linked_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(Exchanger_simu_unstable_change_linked_tmp) = PIVs
  for( i in 1:nrow(df_pred_pairs) ){
    entity1 = RLdata[df_pred_pairs[i,1],]
    entity2 = RLdata[df_pred_pairs[i,2],]
    if(nrow(entity1)>0){
      if(nrow(entity2)>0){
        Exchanger_simu_agree_linked_tmp = rbind(Exchanger_simu_agree_linked_tmp, entity1[,PIVs] == entity2[,PIVs])
        Exchanger_simu_unstable_change_linked_tmp = rbind(Exchanger_simu_unstable_change_linked_tmp, entity1[,"change"]+entity2[,"change"])
      }
    }
  }
  Exchanger_simu_agree_linked_tmp = colSums( Exchanger_simu_agree_linked_tmp, na.rm = TRUE ) / nrow( Exchanger_simu_agree_linked_tmp )
  Exchanger_simu_unstable_change_linked_tmp = colSums( Exchanger_simu_unstable_change_linked_tmp, na.rm = TRUE ) / nrow( Exchanger_simu_unstable_change_linked_tmp )
  Exchanger_simu_agree_linked[iter,] = Exchanger_simu_agree_linked_tmp
  Exchanger_simu_unstable_change_linked[iter,] = Exchanger_simu_unstable_change_linked_tmp
  # TP
  linkedTP = merged_pairs[merged_pairs$match==TRUE & merged_pairs$pred_match==TRUE, ]
  Exchanger_simu_NARL_TP[iter,] = ( ( colSums(is.na(RLdata[linkedTP[,1],])) + colSums(is.na(RLdata[linkedTP[,2],])) ) / nrow(linkedTP) )[PIVs]
  Exchanger_simu_agree_TP_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(Exchanger_simu_agree_TP_tmp) = PIVs
  Exchanger_simu_unstable_change_TP_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(Exchanger_simu_unstable_change_TP_tmp) = PIVs
  for( i in 1:nrow(linkedTP) ){
    entity1 = RLdata[ linkedTP[i,1], ]
    entity2 = RLdata[ linkedTP[i,2], ]
    if(nrow(entity1)>0){
      if(nrow(entity2)>0){
        Exchanger_simu_agree_TP_tmp = rbind(Exchanger_simu_agree_TP_tmp, entity1[,PIVs] == entity2[,PIVs])
        Exchanger_simu_unstable_change_TP_tmp = rbind(Exchanger_simu_unstable_change_TP_tmp, entity1[,"change"]+entity2[,"change"])
      }
    }
  }
  Exchanger_simu_agree_TP_tmp = colSums( Exchanger_simu_agree_TP_tmp, na.rm = TRUE ) / nrow( Exchanger_simu_agree_TP_tmp )
  Exchanger_simu_unstable_change_TP_tmp = colSums( Exchanger_simu_unstable_change_TP_tmp, na.rm = TRUE ) / nrow( Exchanger_simu_unstable_change_TP_tmp )
  Exchanger_simu_agree_TP[iter,] = Exchanger_simu_agree_TP_tmp
  Exchanger_simu_unstable_change_TP[iter,] = Exchanger_simu_unstable_change_TP_tmp
  # FP
  linkedFP = merged_pairs[merged_pairs$match==FALSE & merged_pairs$pred_match==TRUE, ]
  Exchanger_simu_NARL_FP[iter,] = ( ( colSums(is.na(RLdata[linkedFP[,1],])) + colSums(is.na(RLdata[linkedFP[,2],])) ) / nrow(linkedFP) )[PIVs]
  Exchanger_simu_agree_FP_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(Exchanger_simu_agree_FP_tmp) = PIVs
  Exchanger_simu_unstable_change_FP_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(Exchanger_simu_unstable_change_FP_tmp) = PIVs
  for( i in 1:nrow(linkedFP) ){
    entity1 = RLdata[ linkedFP[i,1], ]
    entity2 = RLdata[ linkedFP[i,2], ]
    if(nrow(entity1)>0){
      if(nrow(entity2)>0){
        Exchanger_simu_agree_FP_tmp = rbind(Exchanger_simu_agree_FP_tmp, entity1[,PIVs] == entity2[,PIVs])
        Exchanger_simu_unstable_change_FP_tmp = rbind(Exchanger_simu_unstable_change_FP_tmp, entity1[,"change"]+entity2[,"change"])
      }
    }
  }
  Exchanger_simu_agree_FP_tmp = colSums( Exchanger_simu_agree_FP_tmp, na.rm = TRUE ) / nrow( Exchanger_simu_agree_FP_tmp )
  Exchanger_simu_unstable_change_FP_tmp = colSums( Exchanger_simu_unstable_change_FP_tmp, na.rm = TRUE ) / nrow( Exchanger_simu_unstable_change_FP_tmp )
  Exchanger_simu_agree_FP[iter,] = Exchanger_simu_agree_FP_tmp
  Exchanger_simu_unstable_change_FP[iter,] = Exchanger_simu_unstable_change_FP_tmp
  # FN
  linkedFN = merged_pairs[merged_pairs$match==TRUE & merged_pairs$pred_match==FALSE, ]
  Exchanger_simu_NARL_FN[iter,] = ( ( colSums(is.na(RLdata[linkedFN[,1],])) + colSums(is.na(RLdata[linkedFN[,2],])) ) / nrow(linkedFN) )[PIVs]
  Exchanger_simu_agree_FN_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(Exchanger_simu_agree_FN_tmp) = PIVs
  Exchanger_simu_unstable_change_FN_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(Exchanger_simu_unstable_change_FN_tmp) = PIVs
  for( i in 1:nrow(linkedFN) ){
    entity1 = RLdata[ linkedFN[i,1], ]
    entity2 = RLdata[ linkedFN[i,2], ]
    if(nrow(entity1)>0){
      if(nrow(entity2)>0){
        Exchanger_simu_agree_FN_tmp = rbind(Exchanger_simu_agree_FN_tmp, entity1[,PIVs] == entity2[,PIVs])
        Exchanger_simu_unstable_change_FN_tmp = rbind(Exchanger_simu_unstable_change_FN_tmp, entity1[,"change"]+entity2[,"change"])
      }
    }
  }
  Exchanger_simu_agree_FN_tmp = colSums( Exchanger_simu_agree_FN_tmp, na.rm = TRUE ) / nrow( Exchanger_simu_agree_FN_tmp )
  Exchanger_simu_unstable_change_FN_tmp = colSums( Exchanger_simu_unstable_change_FN_tmp, na.rm = TRUE ) / nrow( Exchanger_simu_unstable_change_FN_tmp )
  Exchanger_simu_agree_FN[iter,] = Exchanger_simu_agree_FN_tmp
  Exchanger_simu_unstable_change_FN[iter,] = Exchanger_simu_unstable_change_FN_tmp
  # METRICS
  tp <- CT["TRUE", "TRUE"]
  fp <- CT["TRUE", "FALSE"]
  fn <- CT["FALSE", "TRUE"]
  precision = tp / (tp + fp)
  recall = tp / (tp + fn)
  f1score = 2 * (precision * recall) / (precision + recall)
  results_f1score[iter, "Exchanger"] = f1score
  results_precision[iter, "Exchanger"] = precision
  results_recall[iter, "Exchanger"] = recall
  results_FN[iter, "Exchanger"] = fn
  results_FP[iter, "Exchanger"] = fp
  results_TP[iter, "Exchanger"] = tp
  
  ### LAUNCH BRL
  Zhat <- BRL(B, A, flds=PIVs, types=c("bi","bi","bi","bi","bi"), nIter=1000) # breaks = c(0, .25, .5)
  ### SAVE RESULTS
  DeltaBRL = matrix(0, nrow=nrow(A), ncol=nrow(B))
  n1 <- nrow(B)
  idxA = which( Zhat <= n1 ) # index in A
  idxB = Zhat[ Zhat <= n1 ] # index in B
  for (l in 1:length(idxA))
  {
    DeltaBRL[idxA[l], idxB[l]] = 1
  }
  # STORY TELLING
  # LINKED
  linkedpairs = which(DeltaBRL>0.5, arr.ind=TRUE)
  linkedpairsA = linkedpairs[,1]
  linkedpairsB = linkedpairs[,2]
  BRL_simu_NAA_linked[iter,] = ( colSums(is.na(A[linkedpairsA,])) / length(linkedpairsA) )[PIVs]
  BRL_simu_NAB_linked[iter,] = ( colSums(is.na(B[linkedpairsB,])) / length(linkedpairsB) )[PIVs]
  BRL_simu_agree_linked_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(BRL_simu_agree_linked_tmp) = PIVs
  BRL_simu_unstable_change_linked_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(BRL_simu_unstable_change_linked_tmp) = PIVs
  for( i in 1:nrow(linkedpairs) ){
    entityA = A[linkedpairsA[i],]
    entityB = B[linkedpairsB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        BRL_simu_agree_linked_tmp = rbind(BRL_simu_agree_linked_tmp, entityA[,PIVs] == entityB[,PIVs])
        BRL_simu_unstable_change_linked_tmp = rbind(BRL_simu_unstable_change_linked_tmp, entityB[,"change"])
      }
    }
  }
  BRL_simu_agree_linked_tmp = colSums( BRL_simu_agree_linked_tmp, na.rm = TRUE ) / nrow( BRL_simu_agree_linked_tmp )
  BRL_simu_unstable_change_linked_tmp = colSums( BRL_simu_unstable_change_linked_tmp, na.rm = TRUE ) / nrow( BRL_simu_unstable_change_linked_tmp )
  BRL_simu_agree_linked[iter,] = BRL_simu_agree_linked_tmp
  BRL_simu_unstable_change_linked[iter,] = BRL_simu_unstable_change_linked_tmp
  # TP
  linkedTP = which((DeltaBRL>0.5) & (Delta==1), arr.ind=TRUE)
  linkedTPA = linkedTP[,1]
  linkedTPB = linkedTP[,2]
  BRL_simu_NAA_TP[iter,] = ( colSums(is.na(A[linkedTPA,])) / length(linkedTPA) )[PIVs]
  BRL_simu_NAB_TP[iter,] = ( colSums(is.na(B[linkedTPB,])) / length(linkedTPB) )[PIVs]
  BRL_simu_agree_TP_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(BRL_simu_agree_TP_tmp) = PIVs
  BRL_simu_unstable_change_TP_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(BRL_simu_unstable_change_TP_tmp) = PIVs
  for( i in 1:nrow(linkedTP) ){
    entityA = A[linkedTPA[i],]
    entityB = B[linkedTPB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        BRL_simu_agree_TP_tmp = rbind(BRL_simu_agree_TP_tmp, entityA[,PIVs] == entityB[,PIVs])
        BRL_simu_unstable_change_TP_tmp = rbind(BRL_simu_unstable_change_TP_tmp, entityB[,"change"])
      }
    }
  }
  BRL_simu_agree_TP_tmp = colSums( BRL_simu_agree_TP_tmp, na.rm = TRUE ) / nrow( BRL_simu_agree_TP_tmp )
  BRL_simu_unstable_change_TP_tmp = colSums( BRL_simu_unstable_change_TP_tmp, na.rm = TRUE ) / nrow( BRL_simu_unstable_change_TP_tmp )
  BRL_simu_agree_TP[iter,] = BRL_simu_agree_TP_tmp
  BRL_simu_unstable_change_TP[iter,] = BRL_simu_unstable_change_TP_tmp
  # FP
  linkedFP = which((DeltaBRL>0.5) & (Delta==0), arr.ind=TRUE)
  linkedFPA = linkedFP[,1]
  linkedFPB = linkedFP[,2]
  BRL_simu_NAA_FP[iter,] = ( colSums(is.na(A[linkedFPA,])) / length(linkedFPA) )[PIVs]
  BRL_simu_NAB_FP[iter,] = ( colSums(is.na(B[linkedFPB,])) / length(linkedFPB) )[PIVs]
  BRL_simu_agree_FP_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(BRL_simu_agree_FP_tmp) = PIVs
  BRL_simu_unstable_change_FP_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(BRL_simu_unstable_change_FP_tmp) = PIVs
  for( i in 1:nrow(linkedFP) ){
    entityA = A[linkedFPA[i],]
    entityB = B[linkedFPB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        BRL_simu_agree_FP_tmp = rbind(BRL_simu_agree_FP_tmp, entityA[,PIVs] == entityB[,PIVs])
        BRL_simu_unstable_change_FP_tmp = rbind(BRL_simu_unstable_change_FP_tmp, entityB[,"change"])
      }
    }
  }
  BRL_simu_agree_FP_tmp = colSums( BRL_simu_agree_FP_tmp, na.rm = TRUE ) / nrow( BRL_simu_agree_FP_tmp )
  BRL_simu_unstable_change_FP_tmp = colSums( BRL_simu_unstable_change_FP_tmp, na.rm = TRUE ) / nrow( BRL_simu_unstable_change_FP_tmp )
  BRL_simu_agree_FP[iter,] = BRL_simu_agree_FP_tmp
  BRL_simu_unstable_change_FP[iter,] = BRL_simu_unstable_change_FP_tmp
  # FN
  linkedFN = which((DeltaBRL<0.5) & (Delta==1), arr.ind=TRUE)
  linkedFNA = linkedFN[,1]
  linkedFNB = linkedFN[,2]
  BRL_simu_NAA_FN[iter,] = ( colSums(is.na(A[linkedFNA,])) / length(linkedFNA) )[PIVs]
  BRL_simu_NAB_FN[iter,] = ( colSums(is.na(B[linkedFNB,])) / length(linkedFNB) )[PIVs]
  BRL_simu_agree_FN_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(BRL_simu_agree_FN_tmp) = PIVs
  BRL_simu_unstable_change_FN_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(BRL_simu_unstable_change_FN_tmp) = PIVs
  for( i in 1:nrow(linkedFN) ){
    entityA = A[linkedFNA[i],]
    entityB = B[linkedFNB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        BRL_simu_agree_FN_tmp = rbind(BRL_simu_agree_FN_tmp, entityA[,PIVs] == entityB[,PIVs])
        BRL_simu_unstable_change_FN_tmp = rbind(BRL_simu_unstable_change_FN_tmp, entityB[,"change"])
      }
    }
  }
  BRL_simu_agree_FN_tmp = colSums( BRL_simu_agree_FN_tmp, na.rm = TRUE ) / nrow( BRL_simu_agree_FN_tmp )
  BRL_simu_unstable_change_FN_tmp = colSums( BRL_simu_unstable_change_FN_tmp, na.rm = TRUE ) / nrow( BRL_simu_unstable_change_FN_tmp )
  BRL_simu_agree_FN[iter,] = BRL_simu_agree_FN_tmp
  BRL_simu_unstable_change_FN[iter,] = BRL_simu_unstable_change_FN_tmp
  # METRICS
  truepositive = sum( (DeltaBRL>0.5) & (Delta==1) )
  falsepositive = sum( (DeltaBRL>0.5) & (Delta==0) )
  falsenegative = sum( (DeltaBRL<0.5) & (Delta==1) )
  precision = truepositive / (truepositive + falsepositive)
  recall = truepositive / (truepositive + falsenegative)
  f1score = 2 * (precision * recall) / (precision + recall)
  results_f1score[iter, "BRL"] = f1score
  results_recall[iter, "BRL"] = recall
  results_precision[iter, "BRL"] = precision
  results_FN[iter, "BRL"] = falsenegative
  results_FP[iter, "BRL"] = falsepositive
  results_TP[iter, "BRL"] = truepositive
  results_MatrixDistance[iter, "BRL"] = sqrt(sum((DeltaBRL - Delta)**2))
  
  ### LAUNCH FlexRL with Instability for V5
  encodedA = A
  encodedB = B
  # Value for missing data
  encodedA[,PIVs][ is.na(encodedA[,PIVs]) ] = 0
  encodedB[,PIVs][ is.na(encodedB[,PIVs]) ] = 0
  # Data
  PIVs_config = list( V1 = list(stable = TRUE),
                      V2 = list(stable = TRUE),
                      V3 = list(stable = TRUE),
                      V4 = list(stable = TRUE),
                      V5 = list(stable = FALSE, conditionalHazard = FALSE, pSameH.cov.A = c(), pSameH.cov.B = c()) )
  PIVs = names(PIVs_config)
  PIVs_stable = sapply(PIVs_config, function(x) x$stable)
  
  ### CHECK THE NAIVE
  DeltaNaive = matrix(0, nrow=nrow(A), ncol=nrow(B))
  # A not missing avec all B
  isNotMissingA = apply(encodedA[,PIVs]!=0, 1, all)
  A_PIVs_notMissing = encodedA[isNotMissingA,PIVs]
  A_PIVs_notMissing_ID = encodedA[isNotMissingA,"localID"]
  UA = sspaste2(as.matrix(A_PIVs_notMissing))
  UB = sspaste2(as.matrix(encodedB[,PIVs]))
  valuesU = unique(c(UA,UB))
  UA = as.numeric(factor(UA,levels=valuesU))
  UB = as.numeric(factor(UB,levels=valuesU)) 
  tmpA = F2(UA, length(valuesU))
  tmpB = F2(UB, length(valuesU))
  select = F33(tmpA, tmpB, length(tmpA))
  if(nrow(select)>1){
    for (l in 1:nrow(select))
    {
      idxA = as.integer(select[l,1])
      idxA = A_PIVs_notMissing_ID[idxA]
      idxB = as.integer(select[l,2])
      if( (idxA <= nrow(DeltaNaive)) & (idxB <= ncol(DeltaNaive)) ){
        DeltaNaive[idxA,idxB] = 1
      }
    }
  }
  # A missing avec all B
  for(k in 1:length(PIVs)){
    isMissingA_k = encodedA[,k]==0
    A_PIVs_k_Missing = encodedA[isMissingA_k,PIVs]
    A_PIVs_k_Missing_ID = encodedA[isMissingA_k,"localID"]
    UA = sspaste2(as.matrix(A_PIVs_k_Missing[,-k]))
    UB = sspaste2(as.matrix(encodedB[,PIVs][,-k]))
    valuesU = unique(c(UA,UB))
    UA = as.numeric(factor(UA,levels=valuesU))
    UB = as.numeric(factor(UB,levels=valuesU)) 
    tmpA = F2(UA, length(valuesU))
    tmpB = F2(UB, length(valuesU))
    select = F33(tmpA, tmpB, length(tmpA))
    if(nrow(select)>1){
      for (l in 1:nrow(select))
      {
        idxA = as.integer(select[l,1])
        idxA = A_PIVs_k_Missing_ID[idxA]
        idxB = as.integer(select[l,2])
        if( (idxA <= nrow(DeltaNaive)) & (idxB <= ncol(DeltaNaive)) ){
          DeltaNaive[idxA,idxB] = 1
        }
      }
    }
  }
  # B not missing avec all A
  isNotMissingB = apply(encodedB[,PIVs]!=0, 1, all)
  B_PIVs_notMissing = encodedB[isNotMissingB,PIVs]
  B_PIVs_notMissing_ID = encodedB[isNotMissingB,"localID"]
  UB = sspaste2(as.matrix(B_PIVs_notMissing))
  UA = sspaste2(as.matrix(encodedA[,PIVs]))
  valuesU = unique(c(UA,UB))
  UA = as.numeric(factor(UA,levels=valuesU))
  UB = as.numeric(factor(UB,levels=valuesU)) 
  tmpA = F2(UA, length(valuesU))
  tmpB = F2(UB, length(valuesU))
  select = F33(tmpA, tmpB, length(tmpA))
  if(nrow(select)>1){
    for (l in 1:nrow(select))
    {
      idxA = as.integer(select[l,1])
      idxB = as.integer(select[l,2])
      idxB = B_PIVs_notMissing_ID[idxB]
      if( (idxA <= nrow(DeltaNaive)) & (idxB <= ncol(DeltaNaive)) ){
        DeltaNaive[idxA,idxB] = 1
      }
    }
  }
  # B missing avec all A
  for(k in 1:length(PIVs)){
    isMissingB_k = encodedB[,k]==0
    B_PIVs_k_Missing = encodedB[isMissingB_k,PIVs]
    B_PIVs_k_Missing_ID = encodedB[isMissingB_k,"localID"]
    UB = sspaste2(as.matrix(B_PIVs_k_Missing[,-k]))
    UA = sspaste2(as.matrix(encodedA[,PIVs][,-k]))
    valuesU = unique(c(UA,UB))
    UA = as.numeric(factor(UA,levels=valuesU))
    UB = as.numeric(factor(UB,levels=valuesU)) 
    tmpA = F2(UA, length(valuesU))
    tmpB = F2(UB, length(valuesU))
    select = F33(tmpA, tmpB, length(tmpA))
    if(nrow(select)>1){
      for (l in 1:nrow(select))
      {
        idxA = as.integer(select[l,1])
        idxB = as.integer(select[l,2])
        idxB = B_PIVs_k_Missing_ID[idxB]
        if( (idxA <= nrow(DeltaNaive)) & (idxB <= ncol(DeltaNaive)) ){
          DeltaNaive[idxA,idxB] = 1
        }
      }
    }
  }
  # STORY TELLING
  # LINKED
  linkedpairs = which(DeltaNaive>0.5, arr.ind=TRUE)
  linkedpairsA = linkedpairs[,1]
  linkedpairsB = linkedpairs[,2]
  Naive_simu_NAA_linked[iter,] = ( colSums(is.na(A[linkedpairsA,])) / length(linkedpairsA) )[PIVs]
  Naive_simu_NAB_linked[iter,] = ( colSums(is.na(B[linkedpairsB,])) / length(linkedpairsB) )[PIVs]
  Naive_simu_agree_linked_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(Naive_simu_agree_linked_tmp) = PIVs
  Naive_simu_unstable_change_linked_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(Naive_simu_unstable_change_linked_tmp) = PIVs
  for( i in 1:nrow(linkedpairs) ){
    entityA = A[linkedpairsA[i],]
    entityB = B[linkedpairsB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        Naive_simu_agree_linked_tmp = rbind(Naive_simu_agree_linked_tmp, entityA[,PIVs] == entityB[,PIVs])
        Naive_simu_unstable_change_linked_tmp = rbind(Naive_simu_unstable_change_linked_tmp, entityB[,"change"])
      }
    }
  }
  Naive_simu_agree_linked_tmp = colSums( Naive_simu_agree_linked_tmp, na.rm = TRUE ) / nrow( Naive_simu_agree_linked_tmp )
  Naive_simu_unstable_change_linked_tmp = colSums( Naive_simu_unstable_change_linked_tmp, na.rm = TRUE ) / nrow( Naive_simu_unstable_change_linked_tmp )
  Naive_simu_agree_linked[iter,] = Naive_simu_agree_linked_tmp
  Naive_simu_unstable_change_linked[iter,] = Naive_simu_unstable_change_linked_tmp
  # TP
  linkedTP = which((DeltaNaive>0.5) & (Delta==1), arr.ind=TRUE)
  linkedTPA = linkedTP[,1]
  linkedTPB = linkedTP[,2]
  Naive_simu_NAA_TP[iter,] = ( colSums(is.na(A[linkedTPA,])) / length(linkedTPA) )[PIVs]
  Naive_simu_NAB_TP[iter,] = ( colSums(is.na(B[linkedTPB,])) / length(linkedTPB) )[PIVs]
  Naive_simu_agree_TP_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(Naive_simu_agree_TP_tmp) = PIVs
  Naive_simu_unstable_change_TP_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(Naive_simu_unstable_change_TP_tmp) = PIVs
  for( i in 1:nrow(linkedTP) ){
    entityA = A[linkedTPA[i],]
    entityB = B[linkedTPB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        Naive_simu_agree_TP_tmp = rbind(Naive_simu_agree_TP_tmp, entityA[,PIVs] == entityB[,PIVs])
        Naive_simu_unstable_change_TP_tmp = rbind(Naive_simu_unstable_change_TP_tmp, entityB[,"change"])
      }
    }
  }
  Naive_simu_agree_TP_tmp = colSums( Naive_simu_agree_TP_tmp, na.rm = TRUE ) / nrow( Naive_simu_agree_TP_tmp )
  Naive_simu_unstable_change_TP_tmp = colSums( Naive_simu_unstable_change_TP_tmp, na.rm = TRUE ) / nrow( Naive_simu_unstable_change_TP_tmp )
  Naive_simu_agree_TP[iter,] = Naive_simu_agree_TP_tmp
  Naive_simu_unstable_change_TP[iter,] = Naive_simu_unstable_change_TP_tmp
  # FP
  linkedFP = which((DeltaNaive>0.5) & (Delta==0), arr.ind=TRUE)
  linkedFPA = linkedFP[,1]
  linkedFPB = linkedFP[,2]
  Naive_simu_NAA_FP[iter,] = ( colSums(is.na(A[linkedFPA,])) / length(linkedFPA) )[PIVs]
  Naive_simu_NAB_FP[iter,] = ( colSums(is.na(B[linkedFPB,])) / length(linkedFPB) )[PIVs]
  Naive_simu_agree_FP_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(Naive_simu_agree_FP_tmp) = PIVs
  Naive_simu_unstable_change_FP_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(Naive_simu_unstable_change_FP_tmp) = PIVs
  for( i in 1:nrow(linkedFP) ){
    entityA = A[linkedFPA[i],]
    entityB = B[linkedFPB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        Naive_simu_agree_FP_tmp = rbind(Naive_simu_agree_FP_tmp, entityA[,PIVs] == entityB[,PIVs])
        Naive_simu_unstable_change_FP_tmp = rbind(Naive_simu_unstable_change_FP_tmp, entityB[,"change"])
      }
    }
  }
  Naive_simu_agree_FP_tmp = colSums( Naive_simu_agree_FP_tmp, na.rm = TRUE ) / nrow( Naive_simu_agree_FP_tmp )
  Naive_simu_unstable_change_FP_tmp = colSums( Naive_simu_unstable_change_FP_tmp, na.rm = TRUE ) / nrow( Naive_simu_unstable_change_FP_tmp )
  Naive_simu_agree_FP[iter,] = Naive_simu_agree_FP_tmp
  Naive_simu_unstable_change_FP[iter,] = Naive_simu_unstable_change_FP_tmp
  # FN
  linkedFN = which((DeltaNaive<0.5) & (Delta==1), arr.ind=TRUE)
  linkedFNA = linkedFN[,1]
  linkedFNB = linkedFN[,2]
  Naive_simu_NAA_FN[iter,] = ( colSums(is.na(A[linkedFNA,])) / length(linkedFNA) )[PIVs]
  Naive_simu_NAB_FN[iter,] = ( colSums(is.na(B[linkedFNB,])) / length(linkedFNB) )[PIVs]
  Naive_simu_agree_FN_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(Naive_simu_agree_FN_tmp) = PIVs
  Naive_simu_unstable_change_FN_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(Naive_simu_unstable_change_FN_tmp) = PIVs
  for( i in 1:nrow(linkedFN) ){
    entityA = A[linkedFNA[i],]
    entityB = B[linkedFNB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        Naive_simu_agree_FN_tmp = rbind(Naive_simu_agree_FN_tmp, entityA[,PIVs] == entityB[,PIVs])
        Naive_simu_unstable_change_FN_tmp = rbind(Naive_simu_unstable_change_FN_tmp, entityB[,"change"])
      }
    }
  }
  Naive_simu_agree_FN_tmp = colSums( Naive_simu_agree_FN_tmp, na.rm = TRUE ) / nrow( Naive_simu_agree_FN_tmp )
  Naive_simu_unstable_change_FN_tmp = colSums( Naive_simu_unstable_change_FN_tmp, na.rm = TRUE ) / nrow( Naive_simu_unstable_change_FN_tmp )
  Naive_simu_agree_FN[iter,] = Naive_simu_agree_FN_tmp
  Naive_simu_unstable_change_FN[iter,] = Naive_simu_unstable_change_FN_tmp
  # METRICS
  truepositive = sum( (DeltaNaive>0.5) & (Delta==1) )
  falsepositive = sum( (DeltaNaive>0.5) & (Delta==0) )
  falsenegative = sum( (DeltaNaive<0.5) & (Delta==1) )
  precision = truepositive / (truepositive + falsepositive)
  recall = truepositive / (truepositive + falsenegative)
  f1score = 2 * (precision * recall) / (precision + recall)
  results_f1score[iter, "Naive"] = f1score
  results_recall[iter, "Naive"] = recall
  results_precision[iter, "Naive"] = precision
  results_FN[iter, "Naive"] = falsenegative
  results_FP[iter, "Naive"] = falsepositive
  results_TP[iter, "Naive"] = truepositive
  results_MatrixDistance[iter, "Naive"] = sqrt(sum((DeltaNaive - Delta)**2))
  
  ### LAUNCH FlexRL with Instability for V5
  dataSimu = list( A             = encodedA,
                   B             = encodedB, 
                   PIVs_config   = PIVs_config,
                   PIVs          = PIVs,
                   PIVs_stable   = PIVs_stable,
                   Nvalues       = Nvalues,
                   newDirectory  = newDirectory )
  
  fit = stEM(  data               = dataSimu,
               StEMIter           = FlexRLUnstableStEMIter,
               StEMBurnin         = FlexRLUnstableStEMBurnin,
               GibbsIter          = FlexRLUnstableGibbsIter,
               GibbsBurnin        = FlexRLUnstableGibbsBurnin,
               UnstableNoMistake  = TRUE  )
  
  ### SAVE RESULTS
  # fit$Delta
  # STORY TELLING
  # LINKED
  linkedpairs = which(fit$Delta>0.5, arr.ind=TRUE)
  linkedpairsA = linkedpairs[,1]
  linkedpairsB = linkedpairs[,2]
  FlexRL_instability4V5_simu_NAA_linked[iter,] = ( colSums(is.na(A[linkedpairsA,])) / length(linkedpairsA) )[PIVs]
  FlexRL_instability4V5_simu_NAB_linked[iter,] = ( colSums(is.na(B[linkedpairsB,])) / length(linkedpairsB) )[PIVs]
  FlexRL_instability4V5_simu_agree_linked_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(FlexRL_instability4V5_simu_agree_linked_tmp) = PIVs
  FlexRL_instability4V5_simu_unstable_change_linked_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(FlexRL_instability4V5_simu_unstable_change_linked_tmp) = PIVs
  for( i in 1:nrow(linkedpairs) ){
    entityA = A[linkedpairsA[i],]
    entityB = B[linkedpairsB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        FlexRL_instability4V5_simu_agree_linked_tmp = rbind(FlexRL_instability4V5_simu_agree_linked_tmp, entityA[,PIVs] == entityB[,PIVs])
        FlexRL_instability4V5_simu_unstable_change_linked_tmp = rbind(FlexRL_instability4V5_simu_unstable_change_linked_tmp, entityB[,"change"])
      }
    }
  }
  FlexRL_instability4V5_simu_agree_linked_tmp = colSums( FlexRL_instability4V5_simu_agree_linked_tmp, na.rm = TRUE ) / nrow( FlexRL_instability4V5_simu_agree_linked_tmp )
  FlexRL_instability4V5_simu_unstable_change_linked_tmp = colSums( FlexRL_instability4V5_simu_unstable_change_linked_tmp, na.rm = TRUE ) / nrow( FlexRL_instability4V5_simu_unstable_change_linked_tmp )
  FlexRL_instability4V5_simu_agree_linked[iter,] = FlexRL_instability4V5_simu_agree_linked_tmp
  FlexRL_instability4V5_simu_unstable_change_linked[iter,] = FlexRL_instability4V5_simu_unstable_change_linked_tmp
  # TP
  linkedTP = which((fit$Delta>0.5) & (Delta==1), arr.ind=TRUE)
  linkedTPA = linkedTP[,1]
  linkedTPB = linkedTP[,2]
  FlexRL_instability4V5_simu_NAA_TP[iter,] = ( colSums(is.na(A[linkedTPA,])) / length(linkedTPA) )[PIVs]
  FlexRL_instability4V5_simu_NAB_TP[iter,] = ( colSums(is.na(B[linkedTPB,])) / length(linkedTPB) )[PIVs]
  FlexRL_instability4V5_simu_agree_TP_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(FlexRL_instability4V5_simu_agree_TP_tmp) = PIVs
  FlexRL_instability4V5_simu_unstable_change_TP_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(FlexRL_instability4V5_simu_unstable_change_TP_tmp) = PIVs
  for( i in 1:nrow(linkedTP) ){
    entityA = A[linkedTPA[i],]
    entityB = B[linkedTPB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        FlexRL_instability4V5_simu_agree_TP_tmp = rbind(FlexRL_instability4V5_simu_agree_TP_tmp, entityA[,PIVs] == entityB[,PIVs])
        FlexRL_instability4V5_simu_unstable_change_TP_tmp = rbind(FlexRL_instability4V5_simu_unstable_change_TP_tmp, entityB[,"change"])
      }
    }
  }
  FlexRL_instability4V5_simu_agree_TP_tmp = colSums( FlexRL_instability4V5_simu_agree_TP_tmp, na.rm = TRUE ) / nrow( FlexRL_instability4V5_simu_agree_TP_tmp )
  FlexRL_instability4V5_simu_unstable_change_TP_tmp = colSums( FlexRL_instability4V5_simu_unstable_change_TP_tmp, na.rm = TRUE ) / nrow( FlexRL_instability4V5_simu_unstable_change_TP_tmp )
  FlexRL_instability4V5_simu_agree_TP[iter,] = FlexRL_instability4V5_simu_agree_TP_tmp
  FlexRL_instability4V5_simu_unstable_change_TP[iter,] = FlexRL_instability4V5_simu_unstable_change_TP_tmp
  # FP
  linkedFP = which((fit$Delta>0.5) & (Delta==0), arr.ind=TRUE)
  linkedFPA = linkedFP[,1]
  linkedFPB = linkedFP[,2]
  FlexRL_instability4V5_simu_NAA_FP[iter,] = ( colSums(is.na(A[linkedFPA,])) / length(linkedFPA) )[PIVs]
  FlexRL_instability4V5_simu_NAB_FP[iter,] = ( colSums(is.na(B[linkedFPB,])) / length(linkedFPB) )[PIVs]
  FlexRL_instability4V5_simu_agree_FP_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(FlexRL_instability4V5_simu_agree_FP_tmp) = PIVs
  FlexRL_instability4V5_simu_unstable_change_FP_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(FlexRL_instability4V5_simu_unstable_change_FP_tmp) = PIVs
  for( i in 1:nrow(linkedFP) ){
    entityA = A[linkedFPA[i],]
    entityB = B[linkedFPB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        FlexRL_instability4V5_simu_agree_FP_tmp = rbind(FlexRL_instability4V5_simu_agree_FP_tmp, entityA[,PIVs] == entityB[,PIVs])
        FlexRL_instability4V5_simu_unstable_change_FP_tmp = rbind(FlexRL_instability4V5_simu_unstable_change_FP_tmp, entityB[,"change"])
      }
    }
  }
  FlexRL_instability4V5_simu_agree_FP_tmp = colSums( FlexRL_instability4V5_simu_agree_FP_tmp, na.rm = TRUE ) / nrow( FlexRL_instability4V5_simu_agree_FP_tmp )
  FlexRL_instability4V5_simu_unstable_change_FP_tmp = colSums( FlexRL_instability4V5_simu_unstable_change_FP_tmp, na.rm = TRUE ) / nrow( FlexRL_instability4V5_simu_unstable_change_FP_tmp )
  FlexRL_instability4V5_simu_agree_FP[iter,] = FlexRL_instability4V5_simu_agree_FP_tmp
  FlexRL_instability4V5_simu_unstable_change_FP[iter,] = FlexRL_instability4V5_simu_unstable_change_FP_tmp
  # FN
  linkedFN = which((fit$Delta<0.5) & (Delta==1), arr.ind=TRUE)
  linkedFNA = linkedFN[,1]
  linkedFNB = linkedFN[,2]
  FlexRL_instability4V5_simu_NAA_FN[iter,] = ( colSums(is.na(A[linkedFNA,])) / length(linkedFNA) )[PIVs]
  FlexRL_instability4V5_simu_NAB_FN[iter,] = ( colSums(is.na(B[linkedFNB,])) / length(linkedFNB) )[PIVs]
  FlexRL_instability4V5_simu_agree_FN_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(FlexRL_instability4V5_simu_agree_FN_tmp) = PIVs
  FlexRL_instability4V5_simu_unstable_change_FN_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(FlexRL_instability4V5_simu_unstable_change_FN_tmp) = PIVs
  for( i in 1:nrow(linkedFN) ){
    entityA = A[linkedFNA[i],]
    entityB = B[linkedFNB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        FlexRL_instability4V5_simu_agree_FN_tmp = rbind(FlexRL_instability4V5_simu_agree_FN_tmp, entityA[,PIVs] == entityB[,PIVs])
        FlexRL_instability4V5_simu_unstable_change_FN_tmp = rbind(FlexRL_instability4V5_simu_unstable_change_FN_tmp, entityB[,"change"])
      }
    }
  }
  FlexRL_instability4V5_simu_agree_FN_tmp = colSums( FlexRL_instability4V5_simu_agree_FN_tmp, na.rm = TRUE ) / nrow( FlexRL_instability4V5_simu_agree_FN_tmp )
  FlexRL_instability4V5_simu_unstable_change_FN_tmp = colSums( FlexRL_instability4V5_simu_unstable_change_FN_tmp, na.rm = TRUE ) / nrow( FlexRL_instability4V5_simu_unstable_change_FN_tmp )
  FlexRL_instability4V5_simu_agree_FN[iter,] = FlexRL_instability4V5_simu_agree_FN_tmp
  FlexRL_instability4V5_simu_unstable_change_FN[iter,] = FlexRL_instability4V5_simu_unstable_change_FN_tmp
  # METRICS
  truepositive = sum( (fit$Delta>0.5) & (Delta==1) )
  falsepositive = sum( (fit$Delta>0.5) & (Delta==0) )
  falsenegative = sum( (fit$Delta<0.5) & (Delta==1) )
  precision = truepositive / (truepositive + falsepositive)
  recall = truepositive / (truepositive + falsenegative)
  f1score = 2 * (precision * recall) / (precision + recall)
  results_f1score[iter, "FlexRL_instability4V5"] = f1score
  results_recall[iter, "FlexRL_instability4V5"] = recall
  results_precision[iter, "FlexRL_instability4V5"] = precision
  results_FN[iter, "FlexRL_instability4V5"] = falsenegative
  results_FP[iter, "FlexRL_instability4V5"] = falsepositive
  results_TP[iter, "FlexRL_instability4V5"] = truepositive
  results_MatrixDistance[iter, "FlexRL_instability4V5"] = sqrt(sum((fit$Delta - Delta)**2))
  # Gamma: proba to form a link
  fit_gamma = fit$gamma # (fit$gamma)[(FlexRLStableStEMBurnin+1):FlexRLStableStEMIter, drop=FALSE]
  results_unstable_gamma[iter,] = fit_gamma # 1 vector of length iter post burnin
  # Phi: agreements between registered and truth for links
  fit_phi = fit$phi # lapply(fit$phi, function(x) x[(FlexRLStableStEMBurnin+1):FlexRLStableStEMIter, , drop=FALSE])
  results_unstable_phi_agree_V1[iter,] = fit_phi[[1]][,1] # 1 vector of length iter post burnin
  results_unstable_phi_agree_V2[iter,] = fit_phi[[2]][,1]
  results_unstable_phi_agree_V3[iter,] = fit_phi[[3]][,1]
  results_unstable_phi_agree_V4[iter,] = fit_phi[[4]][,1]
  results_unstable_phi_agree_V5[iter,] = fit_phi[[5]][,1]
  # Alpha: parameter for the survival model for instability
  fit_alpha = fit$alpha # lapply(fit$alpha, function(x) x[(FlexRLStableStEMBurnin+1):FlexRLStableStEMIter, , drop=FALSE])
  unstablePIV = which(!PIVs_stable)
  alpha_avg = lapply(fit_alpha, function(x) apply(x, 2, mean))
  links = which(fit$Delta>0.5, arr.ind=TRUE)
  times = abs(B[links[,2], "date"] - A[links[,1], "date"])
  Xalpha = cbind( times,
                  A[links[,1], PIVs_config[[unstablePIV]]$pSameH.cov.A, drop=FALSE],
                  B[links[,2], PIVs_config[[unstablePIV]]$pSameH.cov.B, drop=FALSE] )
  results_unstable_alpha_param[iter,] = fit_alpha[[unstablePIV]] # 1 vector of length: iter post burnin
  alpha_probaEstimate = exp( - as.matrix(Xalpha) %*% alpha_avg[[unstablePIV]] )
  results_unstable_alpha_probaEstimate[iter,] = append( alpha_probaEstimate, rep(0, NdataA - length(alpha_probaEstimate)) ) # 1 vector of length linked records + filled with 0
  results_unstable_alpha_timesEstimate[iter,] = append( times, rep(0, NdataA - length(alpha_probaEstimate)) ) # 1 vector of length linked records + filled with 0
  results_unstable_alpha_probaTrue[iter,] = proba_same_H # 1 vector of length Nlinks
  results_unstable_alpha_timesTrue[iter,] = TimeDifference # 1 vector of length Nlinks
  
  ### LAUNCH FlexRL without Instability
  PIVs_config = list( V1 = list(stable = TRUE),
                      V2 = list(stable = TRUE),
                      V3 = list(stable = TRUE),
                      V4 = list(stable = TRUE),
                      V5 = list(stable = TRUE) )
  PIVs = names(PIVs_config)
  PIVs_stable = sapply(PIVs_config, function(x) x$stable)
  
  dataSimu = list( A             = encodedA,
                   B             = encodedB, 
                   PIVs_config   = PIVs_config,
                   PIVs          = PIVs,
                   PIVs_stable   = PIVs_stable,
                   Nvalues       = Nvalues,
                   newDirectory  = newDirectory )
  
  fit = stEM(  data               = dataSimu,
               StEMIter           = FlexRLStableStEMIter,
               StEMBurnin         = FlexRLStableStEMBurnin,
               GibbsIter          = FlexRLStableGibbsIter,
               GibbsBurnin        = FlexRLStableGibbsBurnin   )
  
  ### SAVE RESULTS
  # STORY TELLING
  # LINKED
  linkedpairs = which(fit$Delta>0.5, arr.ind=TRUE)
  linkedpairsA = linkedpairs[,1]
  linkedpairsB = linkedpairs[,2]
  FlexRL_noinstability_simu_NAA_linked[iter,] = ( colSums(is.na(A[linkedpairsA,])) / length(linkedpairsA) )[PIVs]
  FlexRL_noinstability_simu_NAB_linked[iter,] = ( colSums(is.na(B[linkedpairsB,])) / length(linkedpairsB) )[PIVs]
  FlexRL_noinstability_simu_agree_linked_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(FlexRL_noinstability_simu_agree_linked_tmp) = PIVs
  FlexRL_noinstability_simu_unstable_change_linked_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(FlexRL_noinstability_simu_unstable_change_linked_tmp) = PIVs
  for( i in 1:nrow(linkedpairs) ){
    entityA = A[linkedpairsA[i],]
    entityB = B[linkedpairsB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        FlexRL_noinstability_simu_agree_linked_tmp = rbind(FlexRL_noinstability_simu_agree_linked_tmp, entityA[,PIVs] == entityB[,PIVs])
        FlexRL_noinstability_simu_unstable_change_linked_tmp = rbind(FlexRL_noinstability_simu_unstable_change_linked_tmp, entityB[,"change"])
      }
    }
  }
  FlexRL_noinstability_simu_agree_linked_tmp = colSums( FlexRL_noinstability_simu_agree_linked_tmp, na.rm = TRUE ) / nrow( FlexRL_noinstability_simu_agree_linked_tmp )
  FlexRL_noinstability_simu_unstable_change_linked_tmp = colSums( FlexRL_noinstability_simu_unstable_change_linked_tmp, na.rm = TRUE ) / nrow( FlexRL_noinstability_simu_unstable_change_linked_tmp )
  FlexRL_noinstability_simu_agree_linked[iter,] = FlexRL_noinstability_simu_agree_linked_tmp
  FlexRL_noinstability_simu_unstable_change_linked[iter,] = FlexRL_noinstability_simu_unstable_change_linked_tmp
  # TP
  linkedTP = which((fit$Delta>0.5) & (Delta==1), arr.ind=TRUE)
  linkedTPA = linkedTP[,1]
  linkedTPB = linkedTP[,2]
  FlexRL_noinstability_simu_NAA_TP[iter,] = ( colSums(is.na(A[linkedTPA,])) / length(linkedTPA) )[PIVs]
  FlexRL_noinstability_simu_NAB_TP[iter,] = ( colSums(is.na(B[linkedTPB,])) / length(linkedTPB) )[PIVs]
  FlexRL_noinstability_simu_agree_TP_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(FlexRL_noinstability_simu_agree_TP_tmp) = PIVs
  FlexRL_noinstability_simu_unstable_change_TP_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(FlexRL_noinstability_simu_unstable_change_TP_tmp) = PIVs
  for( i in 1:nrow(linkedTP) ){
    entityA = A[linkedTPA[i],]
    entityB = B[linkedTPB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        FlexRL_noinstability_simu_agree_TP_tmp = rbind(FlexRL_noinstability_simu_agree_TP_tmp, entityA[,PIVs] == entityB[,PIVs])
        FlexRL_noinstability_simu_unstable_change_TP_tmp = rbind(FlexRL_noinstability_simu_unstable_change_TP_tmp, entityB[,"change"])
      }
    }
  }
  FlexRL_noinstability_simu_agree_TP_tmp = colSums( FlexRL_noinstability_simu_agree_TP_tmp, na.rm = TRUE ) / nrow( FlexRL_noinstability_simu_agree_TP_tmp )
  FlexRL_noinstability_simu_unstable_change_TP_tmp = colSums( FlexRL_noinstability_simu_unstable_change_TP_tmp, na.rm = TRUE ) / nrow( FlexRL_noinstability_simu_unstable_change_TP_tmp )
  FlexRL_noinstability_simu_agree_TP[iter,] = FlexRL_noinstability_simu_agree_TP_tmp
  FlexRL_noinstability_simu_unstable_change_TP[iter,] = FlexRL_noinstability_simu_unstable_change_TP_tmp
  # FP
  linkedFP = which((fit$Delta>0.5) & (Delta==0), arr.ind=TRUE)
  linkedFPA = linkedFP[,1]
  linkedFPB = linkedFP[,2]
  FlexRL_noinstability_simu_NAA_FP[iter,] = ( colSums(is.na(A[linkedFPA,])) / length(linkedFPA) )[PIVs]
  FlexRL_noinstability_simu_NAB_FP[iter,] = ( colSums(is.na(B[linkedFPB,])) / length(linkedFPB) )[PIVs]
  FlexRL_noinstability_simu_agree_FP_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(FlexRL_noinstability_simu_agree_FP_tmp) = PIVs
  FlexRL_noinstability_simu_unstable_change_FP_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(FlexRL_noinstability_simu_unstable_change_FP_tmp) = PIVs
  for( i in 1:nrow(linkedFP) ){
    entityA = A[linkedFPA[i],]
    entityB = B[linkedFPB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        FlexRL_noinstability_simu_agree_FP_tmp = rbind(FlexRL_noinstability_simu_agree_FP_tmp, entityA[,PIVs] == entityB[,PIVs])
        FlexRL_noinstability_simu_unstable_change_FP_tmp = rbind(FlexRL_noinstability_simu_unstable_change_FP_tmp, entityB[,"change"])
      }
    }
  }
  FlexRL_noinstability_simu_agree_FP_tmp = colSums( FlexRL_noinstability_simu_agree_FP_tmp, na.rm = TRUE ) / nrow( FlexRL_noinstability_simu_agree_FP_tmp )
  FlexRL_noinstability_simu_unstable_change_FP_tmp = colSums( FlexRL_noinstability_simu_unstable_change_FP_tmp, na.rm = TRUE ) / nrow( FlexRL_noinstability_simu_unstable_change_FP_tmp )
  FlexRL_noinstability_simu_agree_FP[iter,] = FlexRL_noinstability_simu_agree_FP_tmp
  FlexRL_noinstability_simu_unstable_change_FP[iter,] = FlexRL_noinstability_simu_unstable_change_FP_tmp
  # FN
  linkedFN = which((fit$Delta<0.5) & (Delta==1), arr.ind=TRUE)
  linkedFNA = linkedFN[,1]
  linkedFNB = linkedFN[,2]
  FlexRL_noinstability_simu_NAA_FN[iter,] = ( colSums(is.na(A[linkedFNA,])) / length(linkedFNA) )[PIVs]
  FlexRL_noinstability_simu_NAB_FN[iter,] = ( colSums(is.na(B[linkedFNB,])) / length(linkedFNB) )[PIVs]
  FlexRL_noinstability_simu_agree_FN_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(FlexRL_noinstability_simu_agree_FN_tmp) = PIVs
  FlexRL_noinstability_simu_unstable_change_FN_tmp = data.frame( matrix(0, nrow=0, ncol=length(PIVs)) )
  colnames(FlexRL_noinstability_simu_unstable_change_FN_tmp) = PIVs
  for( i in 1:nrow(linkedFN) ){
    entityA = A[linkedFNA[i],]
    entityB = B[linkedFNB[i],]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        FlexRL_noinstability_simu_agree_FN_tmp = rbind(FlexRL_noinstability_simu_agree_FN_tmp, entityA[,PIVs] == entityB[,PIVs])
        FlexRL_noinstability_simu_unstable_change_FN_tmp = rbind(FlexRL_noinstability_simu_unstable_change_FN_tmp, entityB[,"change"])
      }
    }
  }
  FlexRL_noinstability_simu_agree_FN_tmp = colSums( FlexRL_noinstability_simu_agree_FN_tmp, na.rm = TRUE ) / nrow( FlexRL_noinstability_simu_agree_FN_tmp )
  FlexRL_noinstability_simu_unstable_change_FN_tmp = colSums( FlexRL_noinstability_simu_unstable_change_FN_tmp, na.rm = TRUE ) / nrow( FlexRL_noinstability_simu_unstable_change_FN_tmp )
  FlexRL_noinstability_simu_agree_FN[iter,] = FlexRL_noinstability_simu_agree_FN_tmp
  FlexRL_noinstability_simu_unstable_change_FN[iter,] = FlexRL_noinstability_simu_unstable_change_FN_tmp
  # METRICS
  truepositive = sum( (fit$Delta>0.5) & (Delta==1) )
  falsepositive = sum( (fit$Delta>0.5) & (Delta==0) )
  falsenegative = sum( (fit$Delta<0.5) & (Delta==1) )
  precision = truepositive / (truepositive + falsepositive)
  recall = truepositive / (truepositive + falsenegative)
  f1score = 2 * (precision * recall) / (precision + recall)
  results_f1score[iter, "FlexRL_noinstability"] = f1score
  results_recall[iter, "FlexRL_noinstability"] = recall
  results_precision[iter, "FlexRL_noinstability"] = precision
  results_FN[iter, "FlexRL_noinstability"] = falsenegative
  results_FP[iter, "FlexRL_noinstability"] = falsepositive
  results_TP[iter, "FlexRL_noinstability"] = truepositive
  results_MatrixDistance[iter, "FlexRL_noinstability"] = sqrt(sum((fit$Delta - Delta)**2))
  # Gamma: proba to form a link
  fit_gamma = fit$gamma # (fit$gamma)[(FlexRLStableStEMBurnin+1):FlexRLStableStEMIter, drop=FALSE]
  results_stable_gamma[iter,] = fit_gamma # 1 vector of length iter post burnin
  # Phi: agreements between registered and truth for links
  fit_phi = fit$phi # lapply(fit$phi, function(x) x[(FlexRLStableStEMBurnin+1):FlexRLStableStEMIter, , drop=FALSE])
  results_stable_phi_agree_V1[iter,] = fit_phi[[1]][,1] # 1 vector of length iter post burnin
  results_stable_phi_agree_V2[iter,] = fit_phi[[2]][,1]
  results_stable_phi_agree_V3[iter,] = fit_phi[[3]][,1]
  results_stable_phi_agree_V4[iter,] = fit_phi[[4]][,1]
  results_stable_phi_agree_V5[iter,] = fit_phi[[5]][,1]
  
  write.csv(simu_NA_A, file.path(newDirectory, "datasetSimu_NA_A.csv"))
  write.csv(simu_NA_B, file.path(newDirectory, "datasetSimu_NA_B.csv"))
  write.csv(simu_truelinks_NA_A, file.path(newDirectory, "datasetSimu_truelinks_NA_A.csv"))
  write.csv(simu_truelinks_NA_B, file.path(newDirectory, "datasetSimu_truelinks_NA_B.csv"))
  write.csv(simu_uniquevalues, file.path(newDirectory, "datasetSimu_uniquevalues.csv"))
  write.csv(simu_truelinks_agree, file.path(newDirectory, "datasetSimu_truelinks_agree.csv"))
  write.csv(simu_truelinks_unstable_change, file.path(newDirectory, "datasetSimu_truelinks_unstable_change.csv"))
  
  write.csv(Exchanger_simu_NARL_TP, file.path(newDirectory, "datasetSimuExchanger_NARL_TP.csv"))
  write.csv(Exchanger_simu_NARL_FP, file.path(newDirectory, "datasetSimuExchanger_NARL_FP.csv"))
  write.csv(Exchanger_simu_NARL_FN, file.path(newDirectory, "datasetSimuExchanger_NARL_FN.csv"))
  write.csv(Exchanger_simu_NARL_linked, file.path(newDirectory, "datasetSimuExchanger_NARL_linked.csv"))
  write.csv(Exchanger_simu_agree_TP, file.path(newDirectory, "datasetSimuExchanger_agree_TP.csv"))
  write.csv(Exchanger_simu_agree_FP, file.path(newDirectory, "datasetSimuExchanger_agree_FP.csv"))
  write.csv(Exchanger_simu_agree_FN, file.path(newDirectory, "datasetSimuExchanger_agree_FN.csv"))
  write.csv(Exchanger_simu_agree_linked, file.path(newDirectory, "datasetSimuExchanger_agree_linked.csv"))
  write.csv(Exchanger_simu_unstable_change_TP, file.path(newDirectory, "datasetSimuExchanger_unstable_change_TP.csv"))
  write.csv(Exchanger_simu_unstable_change_FP, file.path(newDirectory, "datasetSimuExchanger_unstable_change_FP.csv"))
  write.csv(Exchanger_simu_unstable_change_FN, file.path(newDirectory, "datasetSimuExchanger_unstable_change_FN.csv"))
  write.csv(Exchanger_simu_unstable_change_linked, file.path(newDirectory, "datasetSimuExchanger_unstable_change_linked.csv"))
  
  write.csv(BRL_simu_NAA_TP, file.path(newDirectory, "datasetSimuBRL_NAA_TP.csv"))
  write.csv(BRL_simu_NAA_FP, file.path(newDirectory, "datasetSimuBRL_NAA_FP.csv"))
  write.csv(BRL_simu_NAA_FN, file.path(newDirectory, "datasetSimuBRL_NAA_FN.csv"))
  write.csv(BRL_simu_NAA_linked, file.path(newDirectory, "datasetSimuBRL_NAA_linked.csv"))
  write.csv(BRL_simu_NAB_TP, file.path(newDirectory, "datasetSimuBRL_NAB_TP.csv"))
  write.csv(BRL_simu_NAB_FP, file.path(newDirectory, "datasetSimuBRL_NAB_FP.csv"))
  write.csv(BRL_simu_NAB_FN, file.path(newDirectory, "datasetSimuBRL_NAB_FN.csv"))
  write.csv(BRL_simu_NAB_linked, file.path(newDirectory, "datasetSimuBRL_NAB_linked.csv"))
  write.csv(BRL_simu_agree_TP, file.path(newDirectory, "datasetSimuBRL_agree_TP.csv"))
  write.csv(BRL_simu_agree_FP, file.path(newDirectory, "datasetSimuBRL_agree_FP.csv"))
  write.csv(BRL_simu_agree_FN, file.path(newDirectory, "datasetSimuBRL_agree_FN.csv"))
  write.csv(BRL_simu_agree_linked, file.path(newDirectory, "datasetSimuBRL_agree_linked.csv"))
  write.csv(BRL_simu_unstable_change_TP, file.path(newDirectory, "datasetSimuBRL_unstable_change_TP.csv"))
  write.csv(BRL_simu_unstable_change_FP, file.path(newDirectory, "datasetSimuBRL_unstable_change_FP.csv"))
  write.csv(BRL_simu_unstable_change_FN, file.path(newDirectory, "datasetSimuBRL_unstable_change_FN.csv"))
  write.csv(BRL_simu_unstable_change_linked, file.path(newDirectory, "datasetSimuBRL_unstable_change_linked.csv"))
  
  write.csv(Naive_simu_NAA_TP, file.path(newDirectory, "datasetSimuNaive_NAA_TP.csv"))
  write.csv(Naive_simu_NAA_FP, file.path(newDirectory, "datasetSimuNaive_NAA_FP.csv"))
  write.csv(Naive_simu_NAA_FN, file.path(newDirectory, "datasetSimuNaive_NAA_FN.csv"))
  write.csv(Naive_simu_NAA_linked, file.path(newDirectory, "datasetSimuNaive_NAA_linked.csv"))
  write.csv(Naive_simu_NAB_TP, file.path(newDirectory, "datasetSimuNaive_NAB_TP.csv"))
  write.csv(Naive_simu_NAB_FP, file.path(newDirectory, "datasetSimuNaive_NAB_FP.csv"))
  write.csv(Naive_simu_NAB_FN, file.path(newDirectory, "datasetSimuNaive_NAB_FN.csv"))
  write.csv(Naive_simu_NAB_linked, file.path(newDirectory, "datasetSimuNaive_NAB_linked.csv"))
  write.csv(Naive_simu_agree_TP, file.path(newDirectory, "datasetSimuNaive_agree_TP.csv"))
  write.csv(Naive_simu_agree_FP, file.path(newDirectory, "datasetSimuNaive_agree_FP.csv"))
  write.csv(Naive_simu_agree_FN, file.path(newDirectory, "datasetSimuNaive_agree_FN.csv"))
  write.csv(Naive_simu_agree_linked, file.path(newDirectory, "datasetSimuNaive_agree_linked.csv"))
  write.csv(Naive_simu_unstable_change_TP, file.path(newDirectory, "datasetSimuNaive_unstable_change_TP.csv"))
  write.csv(Naive_simu_unstable_change_FP, file.path(newDirectory, "datasetSimuNaive_unstable_change_FP.csv"))
  write.csv(Naive_simu_unstable_change_FN, file.path(newDirectory, "datasetSimuNaive_unstable_change_FN.csv"))
  write.csv(Naive_simu_unstable_change_linked, file.path(newDirectory, "datasetSimuNaive_unstable_change_linked.csv"))
  
  write.csv(FlexRL_instability4V5_simu_NAA_TP, file.path(newDirectory, "datasetSimuFlexRL_instability4V5_NAA_TP.csv"))
  write.csv(FlexRL_instability4V5_simu_NAA_FP, file.path(newDirectory, "datasetSimuFlexRL_instability4V5_NAA_FP.csv"))
  write.csv(FlexRL_instability4V5_simu_NAA_FN, file.path(newDirectory, "datasetSimuFlexRL_instability4V5_NAA_FN.csv"))
  write.csv(FlexRL_instability4V5_simu_NAA_linked, file.path(newDirectory, "datasetSimuFlexRL_instability4V5_NAA_linked.csv"))
  write.csv(FlexRL_instability4V5_simu_NAB_TP, file.path(newDirectory, "datasetSimuFlexRL_instability4V5_NAB_TP.csv"))
  write.csv(FlexRL_instability4V5_simu_NAB_FP, file.path(newDirectory, "datasetSimuFlexRL_instability4V5_NAB_FP.csv"))
  write.csv(FlexRL_instability4V5_simu_NAB_FN, file.path(newDirectory, "datasetSimuFlexRL_instability4V5_NAB_FN.csv"))
  write.csv(FlexRL_instability4V5_simu_NAB_linked, file.path(newDirectory, "datasetSimuFlexRL_instability4V5_NAB_linked.csv"))
  write.csv(FlexRL_instability4V5_simu_agree_TP, file.path(newDirectory, "datasetSimuFlexRL_instability4V5_agree_TP.csv"))
  write.csv(FlexRL_instability4V5_simu_agree_FP, file.path(newDirectory, "datasetSimuFlexRL_instability4V5_agree_FP.csv"))
  write.csv(FlexRL_instability4V5_simu_agree_FN, file.path(newDirectory, "datasetSimuFlexRL_instability4V5_agree_FN.csv"))
  write.csv(FlexRL_instability4V5_simu_agree_linked, file.path(newDirectory, "datasetSimuFlexRL_instability4V5_agree_linked.csv"))
  write.csv(FlexRL_instability4V5_simu_unstable_change_TP, file.path(newDirectory, "datasetSimuFlexRL_instability4V5_unstable_change_TP.csv"))
  write.csv(FlexRL_instability4V5_simu_unstable_change_FP, file.path(newDirectory, "datasetSimuFlexRL_instability4V5_unstable_change_FP.csv"))
  write.csv(FlexRL_instability4V5_simu_unstable_change_FN, file.path(newDirectory, "datasetSimuFlexRL_instability4V5_unstable_change_FN.csv"))
  write.csv(FlexRL_instability4V5_simu_unstable_change_linked, file.path(newDirectory, "datasetSimuFlexRL_instability4V5_unstable_change_linked.csv"))
  
  write.csv(FlexRL_noinstability_simu_NAA_TP, file.path(newDirectory, "datasetSimuFlexRL_noinstability_NAA_TP.csv"))
  write.csv(FlexRL_noinstability_simu_NAA_FP, file.path(newDirectory, "datasetSimuFlexRL_noinstability_NAA_FP.csv"))
  write.csv(FlexRL_noinstability_simu_NAA_FN, file.path(newDirectory, "datasetSimuFlexRL_noinstability_NAA_FN.csv"))
  write.csv(FlexRL_noinstability_simu_NAA_linked, file.path(newDirectory, "datasetSimuFlexRL_noinstability_NAA_linked.csv"))
  write.csv(FlexRL_noinstability_simu_NAB_TP, file.path(newDirectory, "datasetSimuFlexRL_noinstability_NAB_TP.csv"))
  write.csv(FlexRL_noinstability_simu_NAB_FP, file.path(newDirectory, "datasetSimuFlexRL_noinstability_NAB_FP.csv"))
  write.csv(FlexRL_noinstability_simu_NAB_FN, file.path(newDirectory, "datasetSimuFlexRL_noinstability_NAB_FN.csv"))
  write.csv(FlexRL_noinstability_simu_NAB_linked, file.path(newDirectory, "datasetSimuFlexRL_noinstability_NAB_linked.csv"))
  write.csv(FlexRL_noinstability_simu_agree_TP, file.path(newDirectory, "datasetSimuFlexRL_noinstability_agree_TP.csv"))
  write.csv(FlexRL_noinstability_simu_agree_FP, file.path(newDirectory, "datasetSimuFlexRL_noinstability_agree_FP.csv"))
  write.csv(FlexRL_noinstability_simu_agree_FN, file.path(newDirectory, "datasetSimuFlexRL_noinstability_agree_FN.csv"))
  write.csv(FlexRL_noinstability_simu_agree_linked, file.path(newDirectory, "datasetSimuFlexRL_noinstability_agree_linked.csv"))
  write.csv(FlexRL_noinstability_simu_unstable_change_TP, file.path(newDirectory, "datasetSimuFlexRL_noinstability_unstable_change_TP.csv"))
  write.csv(FlexRL_noinstability_simu_unstable_change_FP, file.path(newDirectory, "datasetSimuFlexRL_noinstability_unstable_change_FP.csv"))
  write.csv(FlexRL_noinstability_simu_unstable_change_FN, file.path(newDirectory, "datasetSimuFlexRL_noinstability_unstable_change_FN.csv"))
  write.csv(FlexRL_noinstability_simu_unstable_change_linked, file.path(newDirectory, "datasetSimuFlexRL_noinstability_unstable_change_linked.csv"))
  
  write.csv(results_f1score, file.path(newDirectory, "datasetSimu_results_f1score.csv"))
  write.csv(results_recall, file.path(newDirectory, "datasetSimu_results_recall.csv"))
  write.csv(results_precision, file.path(newDirectory, "datasetSimu_results_precision.csv"))
  write.csv(results_FN, file.path(newDirectory, "datasetSimu_results_FN.csv"))
  write.csv(results_FP, file.path(newDirectory, "datasetSimu_results_FP.csv"))
  write.csv(results_TP, file.path(newDirectory, "datasetSimu_results_TP.csv"))
  write.csv(results_MatrixDistance, file.path(newDirectory, "datasetSimu_results_MatrixDistance.csv"))
  
  write.csv(results_unstable_gamma, file.path(newDirectory, "datasetSimu_results_unstable_gamma.csv"))
  write.csv(results_unstable_phi_agree_V1, file.path(newDirectory, "datasetSimu_results_unstable_phi_agree_V1.csv"))
  write.csv(results_unstable_phi_agree_V2, file.path(newDirectory, "datasetSimu_results_unstable_phi_agree_V2.csv"))
  write.csv(results_unstable_phi_agree_V3, file.path(newDirectory, "datasetSimu_results_unstable_phi_agree_V3.csv"))
  write.csv(results_unstable_phi_agree_V4, file.path(newDirectory, "datasetSimu_results_unstable_phi_agree_V4.csv"))
  write.csv(results_unstable_phi_agree_V5, file.path(newDirectory, "datasetSimu_results_unstable_phi_agree_V5.csv"))
  write.csv(results_unstable_alpha_param, file.path(newDirectory, "datasetSimu_results_unstable_alpha_param.csv"))
  write.csv(results_unstable_alpha_probaEstimate, file.path(newDirectory, "datasetSimu_results_unstable_alpha_probaEstimate.csv"))
  write.csv(results_unstable_alpha_timesEstimate, file.path(newDirectory, "datasetSimu_results_unstable_alpha_timesEstimate.csv"))
  write.csv(results_unstable_alpha_probaTrue, file.path(newDirectory, "datasetSimu_results_unstable_alpha_probaTrue.csv"))
  write.csv(results_unstable_alpha_timesTrue, file.path(newDirectory, "datasetSimu_results_unstable_alpha_timesTrue.csv"))
  
  write.csv(results_stable_gamma, file.path(newDirectory, "datasetSimu_results_stable_gamma.csv"))
  write.csv(results_stable_phi_agree_V1, file.path(newDirectory, "datasetSimu_results_stable_phi_agree_V1.csv"))
  write.csv(results_stable_phi_agree_V2, file.path(newDirectory, "datasetSimu_results_stable_phi_agree_V2.csv"))
  write.csv(results_stable_phi_agree_V3, file.path(newDirectory, "datasetSimu_results_stable_phi_agree_V3.csv"))
  write.csv(results_stable_phi_agree_V4, file.path(newDirectory, "datasetSimu_results_stable_phi_agree_V4.csv"))
  write.csv(results_stable_phi_agree_V5, file.path(newDirectory, "datasetSimu_results_stable_phi_agree_V5.csv"))
  
}
