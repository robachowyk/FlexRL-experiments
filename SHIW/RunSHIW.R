### Run StEM on full SHIW data 

library(progress)
library(Matrix)
source("../recordlinkage.r")
Rcpp:::sourceCpp("../functions.cpp")

results = data.frame( matrix(0, nrow=7, ncol=1) )
rownames(results) = c("F1Score", "Recall", "Precision", "FN", "FP", "TP", "MatrixDistance")
colnames(results) = c("Ours")

DF = read.csv("comp.fixed.csv")

DF = DF[,c("NQUEST","ANNO","NORD","NORDC","SESSO","IPROV","PAR","ANASCI","STACIV","ETA","IREG","AREA5","AREA3","PERC")]
B = DF[DF$ANNO==2008,]
A = DF[DF$ANNO==2010,]

pivs = c("SESSO", "IREG", "PAR", "ANASCI", "STACIV")
pivs_stable = c(TRUE, TRUE, TRUE, TRUE, TRUE)

B$ID = do.call(paste, c(B[,c("NQUEST", "NORDC")], list(sep="_")))
A$ID = do.call(paste, c(A[,c("NQUEST", "NORDC")], list(sep="_")))

linkedID = intersect(B$ID, A$ID)
Nlinks = length( linkedID )

Delta = matrix(0, nrow=nrow(A), ncol=nrow(B))
for (i in 1:Nlinks)
{
  id = linkedID[i]
  idA = which(A$ID == id)
  idB = which(B$ID == id)
  Delta[idA,idB] = 1
}

# STORY TELLING ABOUT THE DATA
NA_A = ( colSums(is.na(A)) / nrow(A) )[pivs]
NA_B = ( colSums(is.na(B)) / nrow(B) )[pivs]
recap = data.frame(matrix(0, nrow=0, ncol=length(pivs)))
colnames(recap) = pivs
for( i in 1:Nlinks ){
  id = linkedID[i]
  entityA = A[A$ID == id,]
  entityB = B[B$ID == id,]
  if(nrow(entityA)>0){
    if(nrow(entityB)>0){
      recap = rbind(recap, entityA[,pivs] == entityB[,pivs])
    }
  }
}
df = data.frame( rbind(NA_A, NA_B, colSums( recap, na.rm = TRUE ) / nrow( recap ), sapply( rbind(A[,pivs],B[,pivs]), function(x) length(unique(x[!(is.na(x))])) )) )
colnames(df) = pivs
rownames(df) = c("NaN in A", "NaN in B", "agreements btw A and B", "unique values")
write.csv(df, "datasetSHIW_recapstory.csv")

encodedA = A
encodedB = B

levels_pivs = lapply(pivs, function(x) levels(factor(as.character(c(encodedA[,x], encodedB[,x])))))

for(i in 1:length(pivs))
{
  encodedA[,pivs[i]] = as.numeric(factor(as.character(encodedA[,pivs[i]]), levels=levels_pivs[[i]]))
  encodedB[,pivs[i]] = as.numeric(factor(as.character(encodedB[,pivs[i]]), levels=levels_pivs[[i]]))
}
nvalues = sapply(levels_pivs, length)

encodedA[,pivs][ is.na(encodedA[,pivs]) ] = 0
encodedB[,pivs][ is.na(encodedB[,pivs]) ] = 0

XA = data.frame() 
XB = data.frame()
regTimeA = rep(2010, nrow(encodedA))
regTimeB = rep(2008, nrow(encodedB))
pSameH.varA = list(c(), c(), c(), c(), c()) 
pSameH.varB = list(c(), c(), c(), c(), c())

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
NA_A = ( colSums(is.na(A[linkedpairsA,])) / length(linkedpairsA) )[pivs]
NA_B = ( colSums(is.na(B[linkedpairsB,])) / length(linkedpairsB) )[pivs]
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
write.csv(df, "datasetSHIW_recaplinkedpairs_Ours.csv")

# save story about TP links
linkedTP = which((fit$Delta>0.5) & (Delta==1), arr.ind=TRUE)
linkedTPA = linkedTP[,1]
linkedTPB = linkedTP[,2]
NA_A = ( colSums(is.na(A[linkedTPA,])) / length(linkedTPA) )[pivs]
NA_B = ( colSums(is.na(B[linkedTPB,])) / length(linkedTPB) )[pivs]
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
write.csv(df, "datasetSHIW_recaplinkedTP_Ours.csv")

# save story about FP links
linkedFP = which((fit$Delta>0.5) & (Delta==0), arr.ind=TRUE)
linkedFPA = linkedFP[,1]
linkedFPB = linkedFP[,2]
NA_A = ( colSums(is.na(A[linkedFPA,])) / length(linkedFPA) )[pivs]
NA_B = ( colSums(is.na(B[linkedFPB,])) / length(linkedFPB) )[pivs]
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
write.csv(df, "datasetSHIW_recaplinkedFP_Ours.csv")

# save story about FN links
linkedFN = which((fit$Delta<0.5) & (Delta==1), arr.ind=TRUE)
linkedFNA = linkedFN[,1]
linkedFNB = linkedFN[,2]
NA_A = ( colSums(is.na(A[linkedFNA,])) / length(linkedFNA) )[pivs]
NA_B = ( colSums(is.na(B[linkedFNB,])) / length(linkedFNB) )[pivs]
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
write.csv(df, "datasetSHIW_recaplinkedFN_Ours.csv")

pdf("datasetSHIW_gamma.pdf")
par(mfrow=c(1,1))
# gamma
gamma = fit$gamma
plot(gamma, type="l", ylim=c(0,1), xlab = "MC-EM Iterations", ylab = "gamma")
abline(h = Nlinks/nrow(A), col="red")
dev.off()

pdf("datasetSHIW_eta.pdf")
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
  pdf("datasetSHIW_omega.pdf")
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

pdf("datasetSHIW_phi.pdf")
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

write.csv(results, "datasetSHIW_results_Ours.csv")