################################################################################
#              RUN NLTCS APPLICATION FOR THE RECORD LINKAGE TASK               #
#          AND COMPARE FLEXRL WITH A NAIVE APPROACH ON THE FULL DATA           #
################################################################################

################################################################################
#                                 IMPORTS                                      #
################################################################################
library(progress)
library(Matrix)
library(MASS)
library(testit)
source("../recordlinkage.r")
Rcpp:::sourceCpp("../functions.cpp")

################################################################################
#                       CREATE A TABLE TO SAVE RESULTS                         #
################################################################################

newDirectory = sprintf("NLTCS Full 19821994 %s", Sys.time())
dir.create(newDirectory)

results = data.frame(matrix(0, nrow=6, ncol=2))
colnames(results) = c("FlexRL", "Naive")
rownames(results) = c("F1Score", "FDR", "Sensitivity", "FN", "FP", "TP")

################################################################################
#                   IMPORT ENCODED DATA (WE ARE NOT ALLOWED                    #
#           TO PUBLISH THE ORIGINAL DATA, THOUGH IT DOES NOT AFFECT            #
#     THE RL PROCESSES SINCE THE ORIGINAL DATA ARE CATEGORICAL VARIABLES)      #
################################################################################

B = read.table("NLTCSData/encodedNLTCS1982.txt")
A = read.table("NLTCSData/encodedNLTCS1994.txt") 

################################################################################
#  STATE AND REGIONAL CODE MAY BE MERGED AS ONE PIV TO AVOID INCONSISTENCIES   #
################################################################################

A$statereg = do.call(paste, c(A[,c("state","reg")], list(sep="_")))
B$statereg = do.call(paste, c(B[,c("state","reg")], list(sep="_")))
levels_statereg = levels(factor(as.character(c(A$statereg, B$statereg))))
A$statereg = as.numeric(factor(as.character(A$statereg), levels=levels_statereg))
B$statereg = as.numeric(factor(as.character(B$statereg), levels=levels_statereg))

################################################################################
#                               PREPARE THE DATA                               #
################################################################################

PIVs_config = list( sex      = list(stable = TRUE),
                    dob_yy   = list(stable = TRUE),
                    dob_mm   = list(stable = TRUE),
                    state    = list(stable = TRUE),
                    reg      = list(stable = TRUE) )

PIVs = names(PIVs_config)

PIVs_stable = sapply(PIVs_config, function(x) x$stable)

### FILTER THE DATA ON THE INTERSECTING SUPPORT OF THE PIVS
for(i in 1:length(PIVs)){
  intersect_support_piv = intersect( unique(A[,PIVs[i]]), unique(B[,PIVs[i]]) )
  A = A[A[,PIVs[i]] %in% c(NA,intersect_support_piv),]
  B = B[B[,PIVs[i]] %in% c(NA,intersect_support_piv),]
}

A = A[ !A$seq %in% A$seq[duplicated(A$seq)] , ]
B = B[ !B$seq %in% B$seq[duplicated(B$seq)] , ]

linkedID = intersect(A$seq, B$seq)
Nlinks = length( linkedID )

### A IS THE SMALLEST FILE
nbrRecordsA = nrow(A)
nbrRecordsB = nrow(B)
nbrRecordsA < nbrRecordsB

### PROCESSING FOR TRUE LINKS
A$localID = 1:nrow(A)
B$localID = 1:nrow(B)

A$source = "A"
B$source = "B"

### CREATE THE TRUE DELTA WITH TRUE IDENTIFIERS
Delta = matrix(0, nrow=nbrRecordsA, ncol=nbrRecordsB)
for (i in 1:Nlinks)
{
  id = linkedID[i]
  idA = which(A$seq == id)
  idB = which(B$seq == id)
  Delta[idA,idB] = 1
}

### MISSING IN THE WHOLE DATASETS (ENCODED WITH NA)
NA_A = ( colSums(is.na(A)) / nrow(A) )[PIVs]
NA_B = ( colSums(is.na(B)) / nrow(B) )[PIVs]

### MISSING IN THE TRUELINKS (ENCODED WITH NA)
NA_A_true = ( colSums(is.na(A[A$seq %in% linkedID,])) / Nlinks )[PIVs]
NA_B_true = ( colSums(is.na(B[B$seq %in% linkedID,])) / Nlinks )[PIVs]

### UNIQUE VALUES IN THE WHOLE DATASETS (ENCODED WITH NA)
unique_values = sapply( rbind(A[,PIVs],B[,PIVs]), function(x) length(unique(x[!(is.na(x))])) )

### AGREEMENTS IN THE TRUE LINKS
recapAgreementsTrueLinks = data.frame(matrix(0, nrow=0, ncol=length(PIVs)))
colnames(recapAgreementsTrueLinks) = PIVs
for( i in 1:Nlinks ){
  matricule_id = linkedID[i]
  entityA = A[A$seq == matricule_id,]
  entityB = B[B$seq == matricule_id,]
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
write.csv(df, file.path(newDirectory, "datasetNLTCS_recapstory.csv"))

### ENCODE THE DATA
encodedA = A
encodedB = B
levels_PIVs = lapply(PIVs, function(x) levels(factor(as.character(c(encodedA[,x], encodedB[,x])))))
for(i in 1:length(PIVs))
{
  encodedA[,PIVs[i]] = as.numeric(factor(as.character(encodedA[,PIVs[i]]), levels=levels_PIVs[[i]]))
  encodedB[,PIVs[i]] = as.numeric(factor(as.character(encodedB[,PIVs[i]]), levels=levels_PIVs[[i]]))
}
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
### PERFORMANCE METRICS
linkedpairs = data.frame(which(DeltaNaive>0.5, arr.ind=TRUE))
linkedpairs = do.call(paste, c(linkedpairs[,c("row","col")], list(sep="_")))
true_pairs = data.frame(which(Delta==1, arr.ind=TRUE))
true_pairs = do.call(paste, c(true_pairs[,c("row","col")], list(sep="_")))
truepositive = length( intersect(linkedpairs,true_pairs) )
falsepositive = length( setdiff(linkedpairs, true_pairs) )
falsenegative = length( setdiff(true_pairs, linkedpairs)  )
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
PIVs_config = list( sex      = list(stable = TRUE),
                    dob_yy   = list(stable = TRUE),
                    dob_mm   = list(stable = TRUE),
                    statereg = list(stable = TRUE) )

PIVs = names(PIVs_config)

PIVs_stable = sapply(PIVs_config, function(x) x$stable)

for(i in 1:length(PIVs)){
  intersect_support_piv = intersect( unique(A[,PIVs[i]]), unique(B[,PIVs[i]]) )
  A = A[A[,PIVs[i]] %in% c(NA,intersect_support_piv),]
  B = B[B[,PIVs[i]] %in% c(NA,intersect_support_piv),]
}
### ENCODE THE DATA
encodedA = A
encodedB = B
levels_PIVs = lapply(PIVs, function(x) levels(factor(as.character(c(encodedA[,x], encodedB[,x])))))
for(i in 1:length(PIVs))
{
  encodedA[,PIVs[i]] = as.numeric(factor(as.character(encodedA[,PIVs[i]]), levels=levels_PIVs[[i]]))
  encodedB[,PIVs[i]] = as.numeric(factor(as.character(encodedB[,PIVs[i]]), levels=levels_PIVs[[i]]))
}
encodedA[,PIVs][ is.na(encodedA[,PIVs]) ] = 0
encodedB[,PIVs][ is.na(encodedB[,PIVs]) ] = 0

nvalues = sapply(levels_PIVs, length)

dataSimu = list( A             = encodedA,
                 B             = encodedB, 
                 PIVs_config   = PIVs_config,
                 PIVs          = PIVs,
                 PIVs_stable   = PIVs_stable,
                 Nvalues       = nvalues,
                 newDirectory  = newDirectory )

fit = stEM(  data               = dataSimu,
             StEMIter           = 3, 
             StEMBurnin         = 1, 
             GibbsIter          = 2, 
             GibbsBurnin        = 1  )

### PERFORMANCE METRICS
linkedpairs = data.frame(which(fit$Delta>0.5, arr.ind=TRUE))
linkedpairs = do.call(paste, c(linkedpairs[,c("row","col")], list(sep="_")))
true_pairs = data.frame(which(Delta==1, arr.ind=TRUE))
true_pairs = do.call(paste, c(true_pairs[,c("row","col")], list(sep="_")))
truepositive = length( intersect(linkedpairs, true_pairs) )
falsepositive = length( setdiff(linkedpairs, true_pairs) )
falsenegative = length( setdiff(true_pairs, linkedpairs)  )
precision = truepositive / (truepositive + falsepositive)
recall = truepositive / (truepositive + falsenegative)
f1score = 2 * (precision * recall) / (precision + recall)
results["F1Score","FlexRL"] = f1score
results["FDR","FlexRL"] = 1 - precision
results["Sensitivity","FlexRL"] = recall
results["FN","FlexRL"] = falsenegative
results["FP","FlexRL"] = falsepositive
results["TP","FlexRL"] = truepositive

write.csv(results, file.path(newDirectory, "datasetNLTCS_results.csv"))

fit_phi = fit$phi

write.csv(fit$gamma, file.path(newDirectory, "datasetNLTCS_results_gamma.csv"))
write.csv(fit$alpha, file.path(newDirectory, "datasetSHIW_results_alpha.csv"))
write.csv(fit_phi[[1]][,1], file.path(newDirectory, "datasetNLTCS_results_phi_agree_V1sex.csv"))
write.csv(fit_phi[[2]][,1], file.path(newDirectory, "datasetNLTCS_results_phi_agree_V2dob_yy.csv"))
write.csv(fit_phi[[3]][,1], file.path(newDirectory, "datasetNLTCS_results_phi_agree_V3dob_mm.csv"))
write.csv(fit_phi[[4]][,1], file.path(newDirectory, "datasetNLTCS_results_phi_agree_V4statereg.csv"))

### SAVE THE VECTOR OF LINKAGE PROBABILITIES (TOO BIG TO SAVE ALL THE ALMOST 0 PROBA)
DeltaVector = as.vector(fit$Delta)
NewDeltaVector = DeltaVector[DeltaVector>0.0001]
nonrpz = length( DeltaVector[DeltaVector<=0.0001] )
rpzonly = length( NewDeltaVector )
write.csv(NewDeltaVector, file.path(newDirectory, sprintf("DeltaVector_WithOnly%sProba_Missing%sAlmostZero.csv", rpzonly, nonrpz)))
