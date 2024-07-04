################################################################################
#               RUN SHIW APPLICATION FOR THE RECORD LINKAGE TASK               #
#           AND COMPARE THE METHODS ON REGIONAL SUBSETS OF THE DATA            #
################################################################################

################################################################################
#                                 IMPORTS                                      #
################################################################################
DF = read.csv("SHIWData/SHIWdata.csv")
DF = DF[,c("SESSO","PAR","ANASCI","STACIV","IREG","ANNO","ID","STUDIO","NASCREG")]
# sesso = sex
# par = household position
# anasci = year of birth
# staciv = marital status
# ireg = regional code
# studio = educational qualification
# nasreg = region of birth

library(exchanger)
library(comparator)
library(clevr)
library(progress)
library(Matrix)
library(testit)
library(BRL)
source("../recordlinkage.r")
Rcpp:::sourceCpp("../functions.cpp")

################################################################################
#                       CREATE A TABLE TO SAVE RESULTS                         #
################################################################################

newDirectory = sprintf("SHIW Subsets Regions %s", Sys.time())
dir.create(newDirectory)

results = data.frame(matrix(0, nrow=6, ncol=4))
colnames(results) = c("FlexRL", "Naive","Exchanger","BRL")
rownames(results) = c("F1Score", "FDR", "Sensitivity", "FN", "FP", "TP")

################################################################################
#                               PREPARE THE DATA                               #
################################################################################

PIVs_config = list( SESSO      = list(stable = TRUE),
                    ANASCI     = list(stable = TRUE),
                    STACIV     = list(stable = TRUE),
                    STUDIO     = list(stable = TRUE),
                    NASCREG    = list(stable = TRUE) )

PIVs = names(PIVs_config)

PIVs_stable = sapply(PIVs_config, function(x) x$stable)

A = DF[DF$ANNO==2020,]
B = DF[DF$ANNO==2016,]

################################################################################
#                       RUN THE METHODS ON EACH REGION                         #
################################################################################

unique_reg = sort(unique(c(A$IREG,B$IREG)))

for(i in 1:length(unique_reg)){
  
  A = DF[DF$ANNO==2020,]
  B = DF[DF$ANNO==2016,]
  
  A = A[ !is.na(A$IREG) & (A$IREG==region), ]
  B = B[ !is.na(B$IREG) & (B$IREG==region), ]
  
  if( length(duplicated(A$ID))>0 ){
    A = A[ !A$ID %in% A$ID[duplicated(A$ID)] , ]
  }
  if( length(duplicated(B$ID))>0 ){
    B = B[ !B$ID %in% B$ID[duplicated(B$ID)] , ]
  }
  
  ### FILTER THE DATA ON THE INTERSECTING SUPPORT OF THE PIVS
  for(i in 1:length(PIVs)){
    intersect_support_piv = intersect( unique(A[,PIVs[i]]), unique(B[,PIVs[i]]) )
    A = A[A[,PIVs[i]] %in% intersect_support_piv,]
    B = B[B[,PIVs[i]] %in% intersect_support_piv,]
  }
  
  ### A IS THE SMALLEST FILE
  nbrRecordsA = nrow(A)
  nbrRecordsB = nrow(B)
  nbrRecordsA < nbrRecordsB
  if(nbrRecordsA >= nbrRecordsB){
    tmpA = A
    tmpB = B
    B = tmpA
    A = tmpB
    nbrRecordsA = nrow(A)
    nbrRecordsB = nrow(B)
  }
  
  linkedID = intersect(A$ID, B$ID)
  Nlinks = length( linkedID )
  
  ### CREATE THE TRUE DELTA WITH TRUE IDENTIFIERS
  DeltaTrue = matrix(0, nrow=nbrRecordsA, ncol=nbrRecordsB)
  for (i in 1:Nlinks)
  {
    id = linkedID[i]
    idA = which(A$ID == id)
    idB = which(B$ID == id)
    DeltaTrue[idA,idB] = 1
  }

  ### MISSING IN THE WHOLE DATASETS (ENCODED WITH NA)
  NA_A = ( colSums(is.na(A)) / nrow(A) )[PIVs]
  NA_B = ( colSums(is.na(B)) / nrow(B) )[PIVs]
  
  ### MISSING IN THE TRUELINKS (ENCODED WITH NA)
  NA_A_true = ( colSums(is.na(A[A$ID %in% linkedID,])) / Nlinks )[PIVs]
  NA_B_true = ( colSums(is.na(B[B$ID %in% linkedID,])) / Nlinks )[PIVs]
  
  ### UNIQUE VALUES IN THE WHOLE DATASETS (ENCODED WITH NA)
  unique_values = sapply( rbind(A[,PIVs],B[,PIVs]), function(x) length(unique(x[!(is.na(x))])) )
  
  ### AGREEMENTS IN THE TRUE LINKS
  recapAgreementsTrueLinks = data.frame(matrix(0, nrow=0, ncol=length(PIVs)))
  colnames(recapAgreementsTrueLinks) = PIVs
  for( i in 1:Nlinks ){
    matricule_id = linkedID[i]
    entityA = A[A$ID == matricule_id,]
    entityB = B[B$ID == matricule_id,]
    if(nrow(entityA)>0){
      if(nrow(entityB)>0){
        recapAgreementsTrueLinks = rbind(recapAgreementsTrueLinks, entityA[,PIVs] == entityB[,PIVs])
      }
    }
  }
  recapAgreementsTrueLinks = colSums( recapAgreementsTrueLinks, na.rm = TRUE ) / nrow( recapAgreementsTrueLinks )
  
  ### STORY TELLING
  df = data.frame( rbind(NA_A, NA_B, NA_A_true, NA_B_true, recapAgreementsTrueLinks, unique_values, nbrRecordsA, nbrRecordsB, Nlinks) )
  colnames(df) = PIVs
  rownames(df) = c("NaN in A", "NaN in B", "NaN in A true", "NaN in B true", "agreements btw A and B true links", "unique values", "size A", "size B", "Nlinks")
  write.csv(df, file.path(newDirectory, sprintf("datasetSHIW_recapstory_%s.csv", region)))
                      
  ### PROCESSING FOR TRUE LINKS
  A$localID = 1:nrow(A)
  B$localID = 1:nrow(B)
  
  A$source = "A"
  B$source = "B"
  
  ### GROUP THE DATA FOR EXCHANGER
  RLdata = rbind(A, B)
  rownames(RLdata) = 1:nrow(RLdata)
  true_pairs = cbind( rownames(RLdata[(RLdata$source=="A")&(RLdata$ID %in% linkedID),]),
                      rownames(RLdata[(RLdata$source=="B")&(RLdata$ID %in% linkedID),]) )
  
  ### EXCHANGER
  distort_prior <- BetaRV(1, 4)
  attr_params <- list(
    SESSO   = CategoricalAttribute(distort_prior,
                                   distort_dist_prior = DirichletProcess(GammaRV(2, 1e-4)),
                                   entity_dist_prior = DirichletRV(1.0)),
    ANASCI  = CategoricalAttribute(distort_prob_prior = distort_prior,
                                   distort_dist_prior = DirichletProcess(GammaRV(2, 1e-4)),
                                   entity_dist_prior = DirichletRV(1.0)),
    STACIV  = CategoricalAttribute(distort_prior,
                                   distort_dist_prior = DirichletProcess(GammaRV(2, 1e-4)),
                                   entity_dist_prior = DirichletRV(1.0)),
    STUDIO  = CategoricalAttribute(distort_prior,
                                   distort_dist_prior = DirichletProcess(GammaRV(2, 1e-4)),
                                   entity_dist_prior = DirichletRV(1.0)),
    NASCREG = CategoricalAttribute(distort_prior,
                                   distort_dist_prior = DirichletProcess(GammaRV(2, 1e-4)),
                                   entity_dist_prior = DirichletRV(1.0))
  )
  clust_prior <- PitmanYorRP(alpha = GammaRV(1, .01), d = BetaRV(1, 1))
  model <- exchanger(RLdata, attr_params, clust_prior)
  result <- exchanger::run_inference(model, n_samples=20, thin_interval=10, burnin_interval=10)
  pred_clust <- smp_clusters(result)
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
  V1_value = merged_pairs[(merged_pairs$pred_match)&(!merged_pairs$match),"V1"]
  V2_value = merged_pairs[(merged_pairs$pred_match)&(!merged_pairs$match),"V2"]
  deduplication_count = sum( RLdata[V1_value,"source"] == RLdata[V2_value,"source"] )
  prediction = factor(merged_pairs$pred_match, levels = c(TRUE, FALSE))
  truth = factor(merged_pairs$match, levels = c(TRUE, FALSE))
  CT = table(prediction, truth, dnn = c("Prediction", "Truth"))
  ### PERFORMANCE METRICS
  tp <- CT["TRUE", "TRUE"]
  fp <- CT["TRUE", "FALSE"]
  fp = fp - deduplication_count
  fn <- CT["FALSE", "TRUE"]
  precision = tp / (tp + fp)
  recall = tp / (tp + fn)
  f1score = 2 * (precision * recall) / (precision + recall)
  results["F1Score", "Exchanger"] = f1score
  results["FDR", "Exchanger"] = 1 - precision
  results["Sensitivity", "Exchanger"] = recall
  results["FN", "Exchanger"] = fn
  results["FP", "Exchanger"] = fp
  results["TP", "Exchanger"] = tp
  results["dedup", "Exchanger"] = deduplication_count
  
  ### BRL
  Zhat <- BRL(B, A, flds=PIVs, types=c("bi","bi","bi","bi","bi"), nIter=1000) # breaks = c(0, .25, .5)
  n1 <- nrow(B)
  idxA = which( Zhat <= n1 ) # index in A
  idxB = Zhat[ Zhat <= n1 ] # index in B
  DeltaBRL = matrix(0, nrow=nrow(A), ncol=nrow(B))
  for (l in 1:length(idxA))
  {
    DeltaBRL[idxA[l], idxB[l]] = 1
  }
  linkedpairs = which(DeltaBRL>0.5, arr.ind=TRUE)
  linkedpairsA = linkedpairs[,1]
  linkedpairsB = linkedpairs[,2]
  linkedTP = which((DeltaBRL>0.5) & (DeltaTrue==1), arr.ind=TRUE)
  linkedTPA = linkedTP[,1]
  linkedTPB = linkedTP[,2]
  linkedFP = which((DeltaBRL>0.5) & (DeltaTrue==0), arr.ind=TRUE)
  linkedFPA = linkedFP[,1]
  linkedFPB = linkedFP[,2]
  linkedFN = which((DeltaBRL<0.5) & (DeltaTrue==1), arr.ind=TRUE)
  linkedFNA = linkedFN[,1]
  linkedFNB = linkedFN[,2]
  ### PERFORMANCE METRICS
  linkedpairs = data.frame(which(DeltaBRL>0.5, arr.ind=TRUE))
  linkedpairs = do.call(paste, c(linkedpairs[,c("row","col")], list(sep="_")))
  true_pairs = data.frame(which(Delta==1, arr.ind=TRUE))
  true_pairs = do.call(paste, c(true_pairs[,c("row","col")], list(sep="_")))
  truepositive = length( intersect(linkedpairs,true_pairs) )
  falsepositive = length( setdiff(linkedpairs, true_pairs) )
  falsenegative = length( setdiff(true_pairs, linkedpairs)  )
  precision = truepositive / (truepositive + falsepositive)
  recall = truepositive / (truepositive + falsenegative)
  f1score = 2 * (precision * recall) / (precision + recall)
  results["F1Score","BRL"] = f1score
  results["FDR","BRL"] = 1 - precision
  results["Sensitivity","BRL"] = recall
  results["FN","BRL"] = falsenegative
  results["FP","BRL"] = falsepositive
  results["TP","BRL"] = truepositive
  
  encodedA = A
  encodedB = B
  levels_PIVs = lapply(PIVs, function(x) levels(factor(as.character(c(encodedA[,x], encodedB[,x])))))
  for(i in 1:length(PIVs))
  {
    encodedA[,PIVs[i]] = as.numeric(factor(as.character(encodedA[,PIVs[i]]), levels=levels_PIVs[[i]]))
    encodedB[,PIVs[i]] = as.numeric(factor(as.character(encodedB[,PIVs[i]]), levels=levels_PIVs[[i]]))
  }
  nvalues = sapply(levels_PIVs, length)
  encodedA[,PIVs][ is.na(encodedA[,PIVs]) ] = 0
  encodedB[,PIVs][ is.na(encodedB[,PIVs]) ] = 0
  
  ### LAUNCH SIMPLISTIC APPROACH
  DeltaNaive = matrix(0, nrow=nrow(A), ncol=nrow(B))
  # A NOT MISSING THAT MATCH WITH B
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
  if(nrow(select>0)){
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
  # A MISSING THAT MATCH WITH B
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
    if(nrow(select>0)){
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
  # B NOT MISSING THAT MATCH WITH A
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
  if(nrow(select>0)){
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
  # B MISSING THAT MATCH WITH A
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
    if(nrow(select>0)){
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
  linkedpairs = which(DeltaNaive>0.5, arr.ind=TRUE)
  linkedpairsA = linkedpairs[,1]
  linkedpairsB = linkedpairs[,2]
  linkedTP = which((DeltaNaive>0.5) & (DeltaTrue==1), arr.ind=TRUE)
  linkedTPA = linkedTP[,1]
  linkedTPB = linkedTP[,2]
  linkedFP = which((DeltaNaive>0.5) & (DeltaTrue==0), arr.ind=TRUE)
  linkedFPA = linkedFP[,1]
  linkedFPB = linkedFP[,2]
  linkedFN = which((DeltaNaive<0.5) & (DeltaTrue==1), arr.ind=TRUE)
  linkedFNA = linkedFN[,1]
  linkedFNB = linkedFN[,2]
  ### PERFORMANCE METRICS
  truepositive = sum( (DeltaNaive>0.5) & (DeltaTrue==1) )
  falsepositive = sum( (DeltaNaive>0.5) & (DeltaTrue==0) )
  falsenegative = sum( (DeltaNaive<0.5) & (DeltaTrue==1) )
  precision = truepositive / (truepositive + falsepositive)
  recall = truepositive / (truepositive + falsenegative)
  f1score = 2 * (precision * recall) / (precision + recall)
  results["F1Score","Naive"] = f1score
  results["FDR","Naive"] = 1 - precision
  results["Sensitivity","Naive"] = recall
  results["FN","Naive"] = falsenegative
  results["FP","Naive"] = falsepositive
  results["TP","Naive"] = truepositive
  
  ### LAUNCH FLEXRL (all PIVs stable)
  dataSimu = list( A             = encodedA,
                   B             = encodedB,
                   PIVs_config   = PIVs_config,
                   PIVs          = PIVs,
                   PIVs_stable   = PIVs_stable,
                   Nvalues       = nvalues,
                   newDirectory  = newDirectory )

  fit = stEM(  data               = dataSimu,
               StEMIter           = 2,
               StEMBurnin         = 1,
               GibbsIter          = 2,
               GibbsBurnin        = 1  )

  linkedpairs = which(fit$Delta>0.5, arr.ind=TRUE)
  linkedpairsA = linkedpairs[,1]
  linkedpairsB = linkedpairs[,2]
  linkedTP = which((fit$Delta>0.5) & (DeltaTrue==1), arr.ind=TRUE)
  linkedTPA = linkedTP[,1]
  linkedTPB = linkedTP[,2]
  linkedFP = which((fit$Delta>0.5) & (DeltaTrue==0), arr.ind=TRUE)
  linkedFPA = linkedFP[,1]
  linkedFPB = linkedFP[,2]
  linkedFN = which(((fit$Delta<0.5) & (DeltaTrue==1)), arr.ind=TRUE)
  linkedFNA = linkedFN[,1]
  linkedFNB = linkedFN[,2]
  ### PERFORMANCE METRICS
  truepositive = sum( (fit$Delta>0.5) & (DeltaTrue==1) )
  falsepositive = sum( (fit$Delta>0.5) & (DeltaTrue==0) )
  falsenegative = sum( (fit$Delta<0.5) & (DeltaTrue==1) )
  precision = truepositive / (truepositive + falsepositive)
  recall = truepositive / (truepositive + falsenegative)
  f1score = 2 * (precision * recall) / (precision + recall)
  results["F1Score","FlexRL"] = f1score
  results["FDR","FlexRL"] = 1 - precision
  results["Sensitivity","FlexRL"] = recall
  results["FN","FlexRL"] = falsenegative
  results["FP","FlexRL"] = falsepositive
  results["TP","FlexRL"] = truepositive
  
  write.csv(results, file.path(newDirectory, sprintf("datasetSHIW_results_%s.csv", region)))

}
