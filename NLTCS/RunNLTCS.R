### Run StEM on full NLTCS data 

library(progress)
library(Matrix)
source("../recordlinkage.r")
Rcpp:::sourceCpp("../functions.cpp")

results = data.frame( matrix(0, nrow=7, ncol=1) )
rownames(results) = c("F1Score", "Recall", "Precision", "FN", "FP", "TP", "MatrixDistance")
colnames(results) = c("Ours")

B = read.table("completeA1982.txt")
A = read.table("completeB1994.txt") # A has to be the smallest dataset: 1994

# B$helperdo1 = rep(NA,nrow(B))
# B$helperdo2 = rep(NA,nrow(B))
# A$helperpresent = rep(NA,nrow(A))

linkedID = intersect(A$seq, B$seq)
Nlinks = length( linkedID )

# can re arrange the data:
# A_linked = A[A$seq %in% linkedID,]
# A_linked <- A_linked[order(A_linked$seq),]
# B_linked = B[B$seq %in% linkedID,]
# B_linked <- B_linked[order(B_linked$seq),]
# A_nonlinked = A[ ! (A$seq %in% linkedID) , ]
# B_nonlinked = B[ ! (B$seq %in% linkedID) , ]
# A = rbind(A_linked, A_nonlinked)
# B = rbind(B_linked, B_nonlinked)
Delta = matrix(0, nrow=nrow(A), ncol=nrow(B))
for (i in 1:Nlinks)
{
  id = linkedID[i]
  idA = which(A$seq == id)
  idB = which(B$seq == id)
  Delta[idA,idB] = 1
}

# try to find some unstability model

# head(A_linked)
# head(B_linked)
# change_occured = (A_linked$reg != B_linked$reg)
# A_linked$HelperTypeChange = paste(A_linked$helpertype, B_linked$helpertype)
# A_linked$HelperTypeChange = as.numeric(factor(A_linked$HelperTypeChange, levels=levels(factor(A_linked$HelperTypeChange))))
# A_linked[change_occured,"HelperTypeChange"]
# B_linked[change_occured,"helperpresent"] - A_linked[change_occured,"helperpresent"]
# hist(abs(B_linked[change_occured,"helpertype"] - A_linked[change_occured,"helpertype"]))
# hist(abs(B_linked[!change_occured,"helpertype"] - A_linked[!change_occured,"helpertype"]))
# B_linked[change_occured,"helperpresent"]
# B_linked[!change_occured,"helperpresent"]

pivs = c("sex", "dob_yy", "dob_mm", "dob_dd", "reg", "state")
pivs_stable = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)

# STORY TELLING ABOUT THE DATA
NA_A = ( colSums((A)==0) / nrow(A) )[pivs]
NA_B = ( colSums((B)==0) / nrow(B) )[pivs]
recap = data.frame(matrix(0, nrow=0, ncol=length(pivs)))
colnames(recap) = pivs
for( i in 1:Nlinks ){
  matricule_id = linkedID[i]
  entityA = A[A$seq == matricule_id,]
  entityB = B[B$seq == matricule_id,]
  if(nrow(entityA)>0){
    if(nrow(entityB)>0){
      recap = rbind(recap, entityA[,pivs] == entityB[,pivs])
    }
  }
}
df = data.frame( rbind(NA_A, NA_B, colSums( recap, na.rm = TRUE ) / nrow( recap ), sapply( rbind(A[,pivs],B[,pivs]), function(x) length(unique(x[!((x)==0)])) )) )
colnames(df) = pivs
rownames(df) = c("NaN in A", "NaN in B", "agreements btw A and B", "unique values")
write.csv(df, "datasetNLCTS_recapstory.csv")
nvalues = df["unique values",]

XA = data.frame()
XB = data.frame()
regTimeA = rep(1982, nrow(A))
regTimeB = rep(1994, nrow(B))
pSameH.varA = list(c(), c(), c(), c(), c(), c())
pSameH.varB = list(c(), c(), c(), c(), c(), c())

dataSimu = list( pivsA        = A[, pivs], 
                 pivsB        = B[, pivs], 
                 nvalues      = as.integer(nvalues),
                 pivs_stable  = pivs_stable, 
                 regTimeA     = regTimeA, 
                 regTimeB     = regTimeB,
                 XA           = XA,
                 XB           = XB,
                 pSameH.varA  = pSameH.varA,
                 pSameH.varB  = pSameH.varB,
                 include_time = FALSE)

# fit = stEM( data       = dataSimu, 
#             nIter      = 300, 
#             nBurnin    = 150, 
#             MStepIter  = 25,
#             cutBurninStEM = 150,
#             trace      = 1 )

fit = stEM( data       = dataSimu, 
            nIter      = 2, 
            nBurnin    = 1, 
            MStepIter  = 1,
            cutBurninStEM = 1,
            trace      = 1 )

# save story about linked pairs, NA, agreements
linkedpairs = which(fit$Delta>0.5, arr.ind=TRUE)
linkedpairsA = linkedpairs[,1]
linkedpairsB = linkedpairs[,2]
NA_A = ( colSums((A[linkedpairsA,])==0) / length(linkedpairsA) )[pivs]
NA_B = ( colSums((B[linkedpairsB,])==0) / length(linkedpairsB) )[pivs]
recap = data.frame(matrix(0, nrow=0, ncol=length(pivs)))
colnames(recap) = pivs
for( i in 1:nrow(linkedpairs) ){
  entityA = A[linkedpairsA[i],]
  entityB = B[linkedpairsB[i],]
  if(nrow(entityA)>0){
    if(nrow(entityB)>0){
      recap = rbind(recap, entityA[,pivs] == entityB[,pivs])
    }
  }
}
df = data.frame( rbind(NA_A, NA_B, colSums( recap, na.rm = TRUE ) / nrow( recap )) )
colnames(df) = pivs
rownames(df) = c("NaN in A", "NaN in B", "agreements btw A and B")
write.csv(df, "datasetNLTCS_recaplinkedpairs_Ours.csv")

# save story about TP links
linkedTP = which((fit$Delta>0.5) & (Delta==1), arr.ind=TRUE)
linkedTPA = linkedTP[,1]
linkedTPB = linkedTP[,2]
NA_A = ( colSums((A[linkedTPA,])==0) / length(linkedTPA) )[pivs]
NA_B = ( colSums((B[linkedTPB,])==0) / length(linkedTPB) )[pivs]
recap = data.frame(matrix(0, nrow=0, ncol=length(pivs)))
colnames(recap) = pivs
for( i in 1:nrow(linkedTP) ){
  entityA = A[linkedTPA[i],]
  entityB = B[linkedTPB[i],]
  if(nrow(entityA)>0){
    if(nrow(entityB)>0){
      recap = rbind(recap, entityA[,pivs] == entityB[,pivs])
    }
  }
}
df = data.frame( rbind(NA_A, NA_B, colSums( recap, na.rm = TRUE ) / nrow( recap )) )
colnames(df) = pivs
rownames(df) = c("NaN in A", "NaN in B", "agreements btw A and B")
write.csv(df, "datasetNLTCS_recaplinkedTP_Ours.csv")

# save story about FP links
linkedFP = which((fit$Delta>0.5) & (Delta==0), arr.ind=TRUE)
linkedFPA = linkedFP[,1]
linkedFPB = linkedFP[,2]
NA_A = ( colSums((A[linkedFPA,])==0) / length(linkedFPA) )[pivs]
NA_B = ( colSums((B[linkedFPB,])==0) / length(linkedFPB) )[pivs]
recap = data.frame(matrix(0, nrow=0, ncol=length(pivs)))
colnames(recap) = pivs
for( i in 1:nrow(linkedFP) ){
  entityA = A[linkedFPA[i],]
  entityB = B[linkedFPB[i],]
  if(nrow(entityA)>0){
    if(nrow(entityB)>0){
      recap = rbind(recap, entityA[,pivs] == entityB[,pivs])
    }
  }
}
df = data.frame( rbind(NA_A, NA_B, colSums( recap, na.rm = TRUE ) / nrow( recap )) )
colnames(df) = pivs
rownames(df) = c("NaN in A", "NaN in B", "agreements btw A and B")
write.csv(df, "datasetNLTCS_recaplinkedFP_Ours.csv")

# save story about FN links
linkedFN = which((fit$Delta<0.5) & (Delta==1), arr.ind=TRUE)
linkedFNA = linkedFN[,1]
linkedFNB = linkedFN[,2]
NA_A = ( colSums((A[linkedFNA,])==0) / length(linkedFNA) )[pivs]
NA_B = ( colSums((B[linkedFNB,])==0) / length(linkedFNB) )[pivs]
recap = data.frame(matrix(0, nrow=0, ncol=length(pivs)))
colnames(recap) = pivs
for( i in 1:nrow(linkedFN) ){
  entityA = A[linkedFNA[i],]
  entityB = B[linkedFNB[i],]
  if(nrow(entityA)>0){
    if(nrow(entityB)>0){
      recap = rbind(recap, entityA[,pivs] == entityB[,pivs])
    }
  }
}
df = data.frame( rbind(NA_A, NA_B, colSums( recap, na.rm = TRUE ) / nrow( recap )) )
colnames(df) = pivs
rownames(df) = c("NaN in A", "NaN in B", "agreements btw A and B")
write.csv(df, "datasetNLTCS_recaplinkedFN_Ours.csv")

pdf("datasetNLTCS_gamma.pdf")
par(mfrow=c(1,1))
# gamma
gamma = fit$gamma
plot(gamma, type="l", ylim=c(0,1), xlab = "MC-EM Iterations", ylab = "gamma")
abline(h = Nlinks/nrow(A), col="red")
dev.off()

pdf("datasetNLTCS_eta.pdf")
par(mar=c(5, 4, 4, 8), mfrow=c(2,3), xpd=TRUE)
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
  # legend("topright", inset=c(-0.2, 0), legend=c("df1","df2"), pch=c(1,3), title="Data")
  legend("topright", inset=c(-0.8, 0), legend=vec, col=seq_len(ncol(eta[[k]])), lty=1:4, cex=0.8)
} 
dev.off()

if(any(!dataSimu$pivs_stable)){
  pdf("datasetNLTCS_omega.pdf")
  par(mfrow=c(2,3))
  # omega
  omega = fit$omega
  for(k in 1:length(omega))
  {
    if(!dataSimu$pivs_stable[k]){
      if(dataSimu$include_time){
        vec = c("time difference")
        plot(omega[[k]][,1], ylim=c(min(omega[[k]])-0.5, max(omega[[k]])+0.5), type="l", col=1, xlab = "MC-EM Iterations", ylab = sprintf("Survival model coef. for unstable PIV %s", pivs[k]))
        if(ncol(omega[[k]])>=2){
          vec <- append(vec, pSameH.varA[[k]])
          vec <- append(vec, pSameH.varB[[k]])
          for(c in 2:ncol(omega[[k]]))
          {
            lines(log(omega[[k]][,c]), col=c)
          }
        }
        legend("topright", legend=vec, col=c(1:ncol(omega[[k]])), lty=1:4, cex=0.8)
      }else{
        vec = c(pSameH.varA[[k]][1])
        plot(omega[[k]][,1], ylim=c(min(omega[[k]])-0.5, max(omega[[k]])+0.5), type="l", col=1, xlab = "MC-EM Iterations", ylab = sprintf("Survival model coef. for unstable PIV %s", pivs[k]))
        if(ncol(omega[[k]])>=2){
          vec <- append(vec, pSameH.varA[[k]])
          vec <- append(vec, pSameH.varB[[k]])
          for(c in 2:ncol(omega[[k]]))
          {
            lines((omega[[k]][,c]), col=c)
          }
        }
        legend("topright", legend=vec, col=c(1:ncol(omega[[k]])), lty=1:4, cex=0.8)
      }
    }
  }
  omega_avg = lapply(omega, function(x) apply(x[(nrow(x)-25):nrow(x),], 2, mean))
  for(k in 1:length(pivs)){
    if(!pivs_stable[[k]]){
      omega_unstable = omega_avg[[k]]
      links = which(fit$Delta>0.5, arr.ind=TRUE)
      times = abs(dataSimu$regTimeB[links[,2]] - dataSimu$regTimeA[links[,1]])
      if(dataSimu$include_time){
        Xomega = cbind( times, dataSimu$XA[links[,1], dataSimu$pSameH.varA[[k]], drop=FALSE], dataSimu$XB[links[,2], dataSimu$pSameH.varB[[k]], drop=FALSE] )
      }else{
        Xomega = cbind( dataSimu$XA[links[,1], dataSimu$pSameH.varA[[k]], drop=FALSE], dataSimu$XB[links[,2], dataSimu$pSameH.varB[[k]], drop=FALSE] )
      }
      plot( sort( exp( - as.matrix(Xomega) %*% omega_unstable ) ), ylim=c(0,1), xlab = "linked records", ylab = sprintf("estimated proba of no change in %s", pivs[k]) )
      # plot( sort(proba_same_H), ylim=c(0,1), xlab = "true records pairs", ylab = sprintf("true model for proba of no change in %s", pivs[k]) )
    }
  }
  dev.off()
}

pdf("datasetNLTCS_phi.pdf")
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

truepositive = sum( (fit$Delta>0.5) & (Delta==1) )
falsepositive = sum( (fit$Delta>0.5) & (Delta==0) )
falsenegative = sum( (fit$Delta<0.5) & (Delta==1) )
precision = truepositive / (truepositive + falsepositive)
recall = truepositive / (truepositive + falsenegative)
f1score = 2 * (precision * recall) / (precision + recall)

results["F1Score", "Ours"] = f1score
results["Recall", "Ours"] = recall
results["Precision", "Ours"] = precision
results["FN", "Ours"] = falsenegative
results["FP", "Ours"] = falsepositive
results["TP", "Ours"] = truepositive
results["MatrixDistance", "Ours"] = sqrt( sum( (fit$Delta - Delta)**2 ) )

write.csv(results, "datasetNLTCS_results_Ours.csv")