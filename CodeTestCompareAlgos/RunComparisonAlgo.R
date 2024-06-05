################################################################################
#       COMPARE OLD MICHEL'S CODE, KAYA'S VERSION AND KAYA'S CPP VERSION       #
#                 ON A REGIONAL SUBSET OF THE NLTCS DATA SETS                  #
################################################################################


################################################################################
#                                 IMPORTS                                      #
################################################################################
library(numDeriv)  
library(progress)
library(Matrix)
library(MASS) 
library(plyr)
library(testit)

################################################################################
#                          FIX AMOUNT OF ITERATIONS                            #
################################################################################
StEMBurninIterations     = 0
StEMIterations           = 20
GibbsIterations          = 50
GibbsBurninIterations    = 25

################################################################################
#                          USE NLTCS REGIONAL SUBSET                           #
################################################################################
B = read.table("NLTCS1982.txt")
A = read.table("NLTCS1994.txt")

PIVs_config = list( sex      = list(stable = TRUE),
                    dob_yy   = list(stable = TRUE),
                    dob_mm   = list(stable = TRUE),
                    state    = list(stable = TRUE) )

PIVs = names(PIVs_config)

PIVs_stable = sapply(PIVs_config, function(x) x$stable)

## FILTER THE DATA ON THE INTERSECTING SUPPORT OF THE PIVS
for(i in 1:length(PIVs)){
  intersect_support_piv = intersect( unique(A[,PIVs[i]]), unique(B[,PIVs[i]]) )
  A = A[ A[,PIVs[i]] %in% c(NA, intersect_support_piv), ]
  B = B[ B[,PIVs[i]] %in% c(NA, intersect_support_piv), ]
}

A = A[ !A$seq %in% A$seq[duplicated(A$seq)] , ]
B = B[ !B$seq %in% B$seq[duplicated(B$seq)] , ]

A = A[ !is.na(A$reg) & A$reg==31 , ]
B = B[ !is.na(B$reg) & B$reg==31 , ]

nbrRecordsA = nrow(A)
nbrRecordsB = nrow(B)
sprintf("Is file A smaller than file B? %s", nbrRecordsA < nbrRecordsB)

linkedID = intersect(A$seq, B$seq)
Nlinks = length( linkedID )
sprintf("Proportion of links (as a fraction of smallest file): %s", Nlinks / nbrRecordsA)

Delta = matrix(0, nrow=nbrRecordsA, ncol=nbrRecordsB)
for (i in 1:Nlinks)
{
  id = linkedID[i]
  idA = which(A$seq == id)
  idB = which(B$seq == id)
  Delta[idA,idB] = 1
}

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

results = data.frame(matrix(0, nrow=6, ncol=3))
colnames(results) = c("FlexRL michel okt", "FlexRL kaya v1", "FlexRL kaya cpp")
rownames(results) = c("F1Score", "FDR", "Sensitivity", "FN", "FP", "TP")

################################################################################
#                   OLD MICHEL'S CODE VERSION OKTOBER 2023                     #
################################################################################
source("recordlinkage_new.r", local = RLokt <- new.env())
Rcpp:::sourceCpp("functions_new.cpp")

set.seed(123)

pSameH.varA = list(c(), c(), c(), c())
pSameH.varB = list(c(), c(), c(), c())
XA = data.frame()
XB = data.frame()

data = list(pivsA        = encodedA[,PIVs], 
            pivsB        = encodedB[,PIVs], 
            nvalues      = nvalues, 
            levels_pivs  = levels_PIVs, 
            pivs_stable  = PIVs_stable, 
            regTimeA     = rep(1,nrow(A)), 
            regTimeB     = rep(1,nrow(B)),
            XA           = XA,
            XB           = XB,
            pSameH.varA  = pSameH.varA,
            pSameH.varB  = pSameH.varB)

fit = RLokt$stEM(data, nIter=StEMIterations, nBurnin=GibbsBurninIterations, MStepIter=GibbsIterations-GibbsBurninIterations, trace=1)

linked_pairs = data.frame(which(fit$Delta>0.5, arr.ind=TRUE))
linked_pairs = do.call(paste, c(linked_pairs[,c("row","col")], list(sep="_")))
true_pairs = data.frame(which(Delta==1, arr.ind=TRUE))
true_pairs = do.call(paste, c(true_pairs[,c("row","col")], list(sep="_")))
truepositive = length( intersect(linked_pairs,true_pairs) )
falsepositive = length( setdiff(linked_pairs, true_pairs) )
falsenegative = length( setdiff(true_pairs, linked_pairs)  )
precision = truepositive / (truepositive + falsepositive)
sensitivity = truepositive / (truepositive + falsenegative)
f1score = 2 * (precision * recall) / (precision + recall)
results["F1Score","FlexRL michel okt"] = f1score
results["FDR","FlexRL michel okt"] = 1 - precision
results["Sensitivity","FlexRL michel okt"] = sensitivity
results["FN","FlexRL michel okt"] = falsenegative
results["FP","FlexRL michel okt"] = falsepositive
results["TP","FlexRL michel okt"] = truepositive
print(results)
print(fit$gamma[(StEMIterations-10):StEMIterations])

################################################################################
#                               KAYA'S VERSION 1                               #
################################################################################
source("recordlinkagekv1.r", local = RLkv1 <- new.env())
Rcpp:::sourceCpp("functions.cpp")

set.seed(123)

newDirectory = sprintf("RL kaya v1 results directory %s", Sys.time())
dir.create(newDirectory)

dataSimu = list( A             = encodedA,
                 B             = encodedB,
                 PIVs_config   = PIVs_config,
                 PIVs          = PIVs,
                 PIVs_stable   = PIVs_stable,
                 Nvalues       = nvalues,
                 newDirectory  = newDirectory )

fitkv1 = RLkv1$stEM(  data                   = dataSimu,
                          StEMIter           = StEMIterations,
                          StEMBurnin         = StEMBurninIterations,
                          GibbsIter          = GibbsIterations,
                          GibbsBurnin        = GibbsBurninIterations  )

linked_pairs = data.frame(which(fitkv1$Delta>0.5, arr.ind=TRUE))
linked_pairs = do.call(paste, c(linked_pairs[,c("row","col")], list(sep="_")))
true_pairs = data.frame(which(Delta==1, arr.ind=TRUE))
true_pairs = do.call(paste, c(true_pairs[,c("row","col")], list(sep="_")))
truepositive = length( intersect(linked_pairs,true_pairs) )
falsepositive = length( setdiff(linked_pairs, true_pairs) )
falsenegative = length( setdiff(true_pairs, linked_pairs)  )
precision = truepositive / (truepositive + falsepositive)
sensitivity = truepositive / (truepositive + falsenegative)
f1score = 2 * (precision * recall) / (precision + recall)
results["F1Score","FlexRL kaya v1"] = f1score
results["FDR","FlexRL kaya v1"] = 1 - precision
results["Sensitivity","FlexRL kaya v1"] = sensitivity
results["FN","FlexRL kaya v1"] = falsenegative
results["FP","FlexRL kaya v1"] = falsepositive
results["TP","FlexRL kaya v1"] = truepositive
print(results)
print(fitkv1$gamma[(StEMIterations-10):StEMIterations])

################################################################################
#                    KAYA'S VERSION 2 WITH MORE CPP CODE                       #
################################################################################
source("recordlinkageCppv2.r", local = RLkv2cpp <- new.env())
Rcpp:::sourceCpp("functions.cpp")

set.seed(123)

newDirectory = sprintf("RL kaya v2 cpp results directory %s", Sys.time())
dir.create(newDirectory)

dataSimu = list( A             = encodedA,
                 B             = encodedB,
                 PIVs_config   = PIVs_config,
                 PIVs          = PIVs,
                 PIVs_stable   = PIVs_stable,
                 Nvalues       = nvalues,
                 newDirectory  = newDirectory )

fitkv2cpp = RLkv2cpp$stEM(  data      = dataSimu,
                   StEMIter           = StEMIterations,
                   StEMBurnin         = StEMBurninIterations,
                   GibbsIter          = GibbsIterations,
                   GibbsBurnin        = GibbsBurninIterations  )

linked_pairs = data.frame(which(fitkv2cpp$Delta>0.5, arr.ind=TRUE))
linked_pairs = do.call(paste, c(linked_pairs[,c("row","col")], list(sep="_")))
true_pairs = data.frame(which(Delta==1, arr.ind=TRUE))
true_pairs = do.call(paste, c(true_pairs[,c("row","col")], list(sep="_")))
truepositive = length( intersect(linked_pairs,true_pairs) )
falsepositive = length( setdiff(linked_pairs, true_pairs) )
falsenegative = length( setdiff(true_pairs, linked_pairs)  )
precision = truepositive / (truepositive + falsepositive)
sensitivity = truepositive / (truepositive + falsenegative)
f1score = 2 * (precision * recall) / (precision + recall)
results["F1Score","FlexRL kaya cpp"] = f1score
results["FDR","FlexRL kaya cpp"] = 1 - precision
results["Sensitivity","FlexRL kaya cpp"] = sensitivity
results["FN","FlexRL kaya cpp"] = falsenegative
results["FP","FlexRL kaya cpp"] = falsepositive
results["TP","FlexRL kaya cpp"] = truepositive
print(results)
print(fitkv2cpp$gamma[(StEMIterations-10):StEMIterations])