DataCreation = function(PIVs_config, PIVs, PIVs_stable, Nval, NRecords, Nlinks, Pmistakes, PmissingA, PmissingB, UnstableNoMistake=FALSE){
  
  for(i in 1:2)
  {	
    dataSet=c()
    for(u in 1:length(Nval))
    {
      # Simulate different probabilities for the true values (xp = exp(0 *(0:(Nval[u]-1))) for uniform distribution of PIVs)
      xp            = exp(0.25 *(0:(Nval[u]-1)))
      probx         = xp/sum(xp)
      dataSet       = cbind(dataSet, sample(1:Nval[u], NRecords[i], replace=TRUE, prob=probx))
    }
    dataSet         = as.data.frame(dataSet)
    names(dataSet)  = PIVs
    assign( paste("dataSet", i, sep=""), dataSet )
  }
  
  # Add overlapping units
  dataSet2[1:Nlinks,] = dataSet1[1:Nlinks,]
  
  # Add typos
  for(x in 1:ncol(dataSet1))
  { 
    biased = as.logical(rbinom(nrow(dataSet1),1,Pmistakes[x]))
    if(sum(biased)>0)
      dataSet1[,x][biased] = sapply(dataSet1[,x][biased], function(i) sample((1:Nval[x])[-c(i)], 1))
  }  
  
  for(x in 1:ncol(dataSet2))
  { 
    biased = as.logical(rbinom(nrow(dataSet2),1,Pmistakes[x]))
    if(sum(biased)>0)
      dataSet2[,x][biased] = sapply(dataSet2[,x][biased], function(i) sample((1:Nval[x])[-c(i)], 1))
  } 
  
  # Add missings
  for(x in 1:ncol(dataSet1))
  { 
    biased = as.logical(rbinom(nrow(dataSet1),1,PmissingA[x]))
    if(sum(biased)>0)
      dataSet1[,x][biased] = NA
  }  
  
  for(x in 1:ncol(dataSet2))
  { 
    biased = as.logical(rbinom(nrow(dataSet2),1,PmissingB[x]))
    if(sum(biased)>0)
      dataSet2[,x][biased] = NA
  }
  
  # Add registration time
  dataSet1$date                     = runif(nrow(dataSet1), 0, 3)
  dataSet2$date                     = runif(nrow(dataSet2), 3, 6)
  if(!UnstableNoMistake){
    # There is instability and mistakes so: to avoid identifiability problems:
    # there should be null registration time difference in links 
    nullTimeDiff                      = as.integer(Nlinks/2)
    dataSet1[1:nullTimeDiff, "date"]  = runif(nullTimeDiff, 0.00, 0.01)
    dataSet2[1:nullTimeDiff, "date"]  = runif(nullTimeDiff, 0.00, 0.01)
  }
  
  TimeDifference                    = abs( dataSet2[1:Nlinks, "date"] - dataSet1[1:Nlinks, "date"] )
  proba_same_H                      = exp( - 0.28 * (TimeDifference) )
  
  dataSet1$change                   = FALSE
  dataSet2$change                   = FALSE
  
  unstablePIV = which(!PIVs_stable)
  for(i in 1:Nlinks)
  {
    is_not_moving           = rbinom(1, 1, proba_same_H[i])
    if(!is_not_moving)
    {
      # the value MUST be different
      dataSet2[i,x]         = sample((1:Nval[x])[-c(dataSet1[i,x])], 1)
      dataSet2[i,"change"]  = TRUE
    }
  }
  # plot( sort(proba_same_H, decreasing=TRUE), ylim=c(0,1), xlab = "true records pairs", ylab = sprintf("proba of no change in %s", PIVs[unstablePIV]) )
  # plot( TimeDifference, proba_same_H, ylim=c(0,1), xlab = "time difference", ylab = sprintf("proba of no change in %s", PIVs[unstablePIV]) )
  
  A = dataSet1
  B = dataSet2
  
  # Recode the PIVs
  levels_PIVs = lapply(PIVs, function(x) levels(factor(as.character(c(A[,x], B[,x])))))
  
  for(i in 1:length(PIVs))
  {
    A[,PIVs[i]] = as.numeric(factor(as.character(A[,PIVs[i]]), levels=levels_PIVs[[i]]))
    B[,PIVs[i]] = as.numeric(factor(as.character(B[,PIVs[i]]), levels=levels_PIVs[[i]]))
  }
  
  Nvalues       = sapply(levels_PIVs, length)
  
  A$localID     = 1:nrow(A)
  B$localID     = 1:nrow(B)
  
  A$source      = "A"
  B$source      = "B"
  
  list(A=A, B=B, Nvalues=Nvalues, TimeDifference=TimeDifference, proba_same_H=proba_same_H)
}

createOmegadata <- function(ncoef, stable){
  if (!stable){ return ( data.frame(matrix(nrow = 0, ncol = ncoef)) ) } }

logPossibleConfig = function(Brecords,sumD)
{
  return = 0
  if(sumD>0)
    return = sum( log(Brecords:(Brecords-sumD+1))  )
  return 
}

loglik = function(LLL, LLA, LLB, D, links, sumRowD, sumColD, gamma)
{
  # Sanity check for using logPossibleConfig(.) below
  if(ncol(D) - sum(D) + 1 <= 0){
    print("The number of records in B has to be >= to the number of linked records.")
    print("Number of records in B:")
    print(ncol(D))
    print("Number of linked records:")
    print(sum(D))
    print("The problem may come from the fact that file A is bigger than file B, which should not happen.")
    print("Number of records in A:")
    print(nrow(D))
  }
  logPossD = sum(log(gamma) * sumRowD + log(1-gamma) * (1-sumRowD)) - logPossibleConfig(ncol(D),sum(D))
  logPossD + sum(LLA[sumRowD==0]) + sum(LLB[sumColD==0]) + sum(LLL[links])
}

simulateH_old = function(data, D, linksR, omegaData, sumRowD, sumColD, eta, phi)
{
  truepivsA = array(NA, dim(data$A[,data$PIVs]))
  truepivsB = array(NA, dim(data$B[,data$PIVs]))
  
  nonlinkedA = sumRowD==0
  nonlinkedB = sumColD==0
  
  # print(linksR)
  
  for(k in 1:length(data$Nvalues))
  {
    phi_k_A = phi[[k]][c(1,2)]
    phi_k_B = phi[[k]][c(1,3)]
    truepivsA[nonlinkedA,k] = sampleNL(G=data$A[nonlinkedA,data$PIVs][,k], eta=eta[[k]], phi=phi_k_A)
    truepivsB[nonlinkedB,k] = sampleNL(G=data$B[nonlinkedB,data$PIVs][,k], eta=eta[[k]], phi=phi_k_B)
    # print(head(truepivsA[nonlinkedA,k]))
    
    if(nrow(linksR)>0)
    {
      # Choice sets and comparisons (this is just some bookkeeping)
      choice_set = expand.grid(1:data$Nvalues[k], 1:data$Nvalues[k])
      # print(head(choice_set),15)
      choice_equal = as.numeric(choice_set[,2] == choice_set[,1])
      # print(head(choice_equal),15)
      # Prepare P(HA)
      eta_choice = eta[[k]][choice_set[,1]]
      
      # If the variable is stable: times are 0 such that pSame = exp(-0) = 1
      if(data$PIVs_stable[k])
      {
        survivalpSameH = rep(1, nrow(linksR))
      }else{
        survivalpSameH = exp(-omegaData[[k]])
      }

      out = sampleL(   GA=data$A[linksR[,1],data$PIVs][,k], 
                       GB=data$B[linksR[,2],data$PIVs][,k], 
                       survivalpSameH=survivalpSameH,
                       choice_set=as.matrix(choice_set), 
                       choice_equal=choice_equal, 
                       nval = data$Nvalues[k], 
                       phikA=phi_k_A, 
                       phikB=phi_k_B,
                       eta=eta_choice)
      
      truepivsA[linksR[,1],k] = choice_set[out,1]
      truepivsB[linksR[,2],k] = choice_set[out,2]
    }
  }
  list(truepivsA=truepivsA, truepivsB=truepivsB)
}

simulateH = function(data, D, links, omegaData, sumRowD, sumColD, eta, phi)
{
  nonlinkedA = sumRowD==0
  nonlinkedB = sumColD==0
  
  # 2 alterantives to sample H in cpp / similar / sampleH is faster
  # if (nrow(links) > 0){
  #   times = abs(data$regTimeB[links[,2]+1] - data$regTimeA[links[,1]+1])
  # }else{
  #   times = 0
  # }
  # truePIVs = sampleH_old( nA=dim(data$pivsA), nB=dim(data$pivsB), times=times, Xomegaomega=omegaData, pivs_stable=data$pivs_stable, pivsA=data$pivsA, pivsB=data$pivsB, links=links+1, nvalues=data$nvalues, D=D, nonlinkedA=nonlinkedA, nonlinkedB=nonlinkedB, eta=eta, omega=omega, phi=phi )
  truePIVs = sampleH(nA=dim(data$A[,data$PIVs]), nB=dim(data$B[,data$PIVs]), links=links, omegaData=omegaData, pivs_stable=data$PIVs_stable, pivsA=data$A[,data$PIVs], pivsB=data$B[,data$PIVs], nvalues=data$Nvalues, D=D, nonlinkedA=nonlinkedA, nonlinkedB=nonlinkedB, eta=eta, phi=phi)

  list(truepivsA=truePIVs$truepivsA, truepivsB=truePIVs$truepivsB)
}

simulateD = function(data, D, linksR, sumRowD, sumColD, truepivsA, truepivsB, gamma, eta, expalpha, phi)
{
  #Determine which observation pairs can be matches (or not) based on the true values of the PIVS
  UA = sspaste2(as.matrix(truepivsA[,data$PIVs_stable]))#do.call(paste, c(as.data.frame(truepivsA), list(sep="_")))
  UB = sspaste2(as.matrix(truepivsB[,data$PIVs_stable]))#do.call(paste, c(as.data.frame(truepivsB), list(sep="_")))
  # print(head(UA))
  tmpUA = UA
  tmpUB = UB
  valuesU = unique(c(UA,UB))
  # print(sort(valuesU))
  UA = as.numeric(factor(UA,levels=valuesU))
  UB = as.numeric(factor(UB,levels=valuesU)) 
  tmpA = F2(UA, length(valuesU))
  tmpB = F2(UB, length(valuesU))
  select = F33(tmpA, tmpB, length(tmpA))
  #print(select)
  # print(select)
  pLink = rep(gamma, nrow(data$A[,data$PIVs]))
  # What should we add/substract from the loglikelihood if an observation is a not linked
  LLA = rep(0,nrow(data$A[,data$PIVs]))
  for(k in 1:length(data$Nvalues))
  {
    logpTrue = log(eta[[k]])[truepivsA[,k]]
    pMissingA = phi[[k]][2]
    pTypoA = (1-pMissingA) * (1-phi[[k]][1]) / (data$Nvalues[k]-1)
    pAgreeA = (1-pMissingA) * phi[[k]][1]
    # Contribution to the likelihood
    contr = rep(pAgreeA, nrow(data$A[,data$PIVs]))
    contr[data$A[,data$PIVs][,k] != truepivsA[,k]] = pTypoA 
    contr[data$A[,data$PIVs][,k] == 0] = pMissingA
    LLA = LLA + logpTrue + log(contr)
  } 
  # What should we add/substract from the loglikelihood if an observation is a not linked
  LLB = rep(0,nrow(data$B[,data$PIVs])) 
  for(k in 1:length(data$Nvalues))
  {
    logpTrue = log(eta[[k]])[truepivsB[,k]]
    pMissingB = phi[[k]][3]
    pTypoB = (1-pMissingB) * (1-phi[[k]][1]) / (data$Nvalues[k]-1)
    pAgreeB = (1-pMissingB) * phi[[k]][1]
    # Contribution to the likelihood
    contr = rep(pAgreeB, nrow(data$B[,data$PIVs]))
    contr[data$B[,data$PIVs][,k] != truepivsB[,k]] = pTypoB
    contr[data$B[,data$PIVs][,k] == 0] = pMissingB
    LLB = LLB + logpTrue + log(contr)
  } 
  # What do we add/substract from the loglikelihood if a pair is linked
  LLL = Matrix(0, nrow=nrow(data$A[,data$PIVs]), ncol=nrow(data$B[,data$PIVs]))
  for(k in 1:length(data$Nvalues))
  {
    HA = truepivsA[ select[,1],k ]
    HB = truepivsB[ select[,2],k ]
    logpTrue = log(eta[[k]])[HA]
    pMissingA = phi[[k]][2]
    pTypoA = (1-pMissingA) * (1-phi[[k]][1]) / (data$Nvalues[k]-1)
    pAgreeA = (1-pMissingA) * phi[[k]][1]
    pMissingB = phi[[k]][3]
    pTypoB = (1-pMissingB) * (1-phi[[k]][1]) / (data$Nvalues[k]-1)
    pAgreeB = (1-pMissingB) * phi[[k]][1]
    # Contribution to the likelihood of linked observation from A
    helpA = rep(pAgreeA, length(HA))
    helpA[data$A[,data$PIVs][select[,1],k] != HA] = pTypoA
    helpA[data$A[,data$PIVs][select[,1],k] == 0] = pMissingA
    # Contribution to the likelihood of linked observation from B
    helpB = rep(pAgreeB, length(HB))
    helpB[data$B[,data$PIVs][select[,2],k] != HB] = pTypoB
    helpB[data$B[,data$PIVs][select[,2],k] == 0] = pMissingB
    LLL[select] = LLL[select] + logpTrue + log(helpA) + log(helpB)
    # Add unstable part if unstable
    if(!data$PIVs_stable[k])
    {
      # MODIF NEW
      #times = rep(1, nrow(select))
      times = abs(data$B[select[,2], "date"] - data$A[select[,1], "date"])
      times = cbind( times )
      cov = cbind( data$A[select[,1], data$PIVs_config[[k]]$pSameH.cov.A, drop=FALSE],
                   data$B[select[,2], data$PIVs_config[[k]]$pSameH.cov.B, drop=FALSE] ) # ??? TODO
      cumulativeHazard_k = as.matrix(times) %*% expalpha[[k]]
      pSameH = exp(- cumulativeHazard_k)
      helpH = pSameH^(HA==HB) * ((1-pSameH)/(data$Nvalues[k]-1))^(HA!=HB)
      LLL[select] = LLL[select] + log(helpH)
    }
  }  
  LLL[select][is.na(LLL[select])] = Inf
  # Complete data likelihood
  LL0 = loglik(LLL=LLL, LLA=LLA, LLB=LLB, D=D, links=linksR, sumRowD=sumRowD, sumColD=sumColD, gamma=gamma)
  # Single run through D
  Dsample = sampleD(S=as.matrix(select),
                    LLA=LLA, 
                    LLB=LLB, 
                    LLL=LLL,
                    gamma=pLink, 
                    loglik=LL0, 
                    D=D, 
                    nlinkrec=as.integer(sum(D)), 
                    sumRowD=sumRowD>0, 
                    sumColD=sumColD>0)
  linksR = Dsample$links+1
  # Sanity check: does it give the same likelihood?
  if (round(Dsample$loglik, digits = 3) != round(loglik( LLL=LLL, LLA=LLA, LLB=LLB, D=Dsample$D, links=linksR, sumRowD=Dsample$sumRowD, sumColD=Dsample$sumColD, gamma=pLink), digits = 3))
  {
    print( "Sanity check failed" )
    print( "Log likelihood associated with new Delta:" )
    print( Dsample$loglik )
    print( "Log likelihood computed on the new Delta:" )
    print( loglik( LLL=LLL, LLA=LLA, LLB=LLB, D=Dsample$D, links=linksR, sumRowD=Dsample$sumRowD, sumColD=Dsample$sumColD, gamma=pLink) )
  }
  Dsample
}

loglikOmega = function(logexpalpha, dataOmega)
{
  S = exp( - as.matrix(dataOmega[,!names(dataOmega) %in% c("equal")]) %*% exp(logexpalpha) )
  -sum(log(S) * dataOmega$equal + log(1-S) * !dataOmega$equal)
} 

#
# dataSimu = list( A             = A,
#                  B             = B, 
#                  PIVs_config   = PIVs_config,
#                  PIVs          = PIVs,
#                  PIVs_stable   = PIVs_stable,
#                  Nvalues       = Nvalues
# )
# if false in stable: and no cov specified, simple model with time ==> implies include time
# if false in stable: and cov specified, we add the covs to time! impleis time in any case
# StEMIter     = 
# StEMBurnin   = 
# GibbsIter    = 
# GibbsBurnin  = 
# a revoir: recompter les burnin differemment
# V5 = list(stable = FALSE, conditionalHazard = FALSE, pSameH.cov.A = c(), pSameH.cov.B = c()))

stEM = function(data, StEMIter, StEMBurnin, GibbsIter, GibbsBurnin, trace=1, UnstableNoMistake=FALSE)
{
  
  instability = any(!data$PIVs_stable)
  if(instability){
    unstablePIVs = which(!data$PIVs_stable)
    instability_conditionalHazard = unlist(sapply(data$PIVs_config, function(x) x$conditionalHazard))
    conditionalHazard = any(instability_conditionalHazard)
    if(conditionalHazard){
      for(k in 1:length(unstablePIVs)){
        unstablePIV = unstablePIVs[k]
        if( length(data$PIVs_config[[unstablePIV]]$pSameH.cov.A)>0 )
          assert( "Some variables from A to include in the survival model for unstable PIVs do not exist.", data$PIVs_config[[unstablePIV]]$pSameH.cov.A %in% colnames(A) )
        if( length(data$PIVs_config[[unstablePIV]]$pSameH.cov.B)>0 )
          assert( "Some variables from B to include in the survival model for unstable PIVs do not exist.", data$PIVs_config[[unstablePIV]]$pSameH.cov.B %in% colnames(B) )
      }
    }
  }
  # now assume conditional hazard = FALSE (none of the usntable has conditional hazard function, only baseline), to work on!
  
  nGibbsIter = GibbsIter - GibbsBurnin
  assert( "Number of iterations for StEM should be positive.", StEMIter - StEMBurnin > 0 )
  assert( "Number of iterations for StEM Gibbs sampler should be positive.", nGibbsIter > 0 )
  
  # linePID = system("ps -C R")
  # lineV = system("ps v")
  # write(linePID, file="memorytrackSimuRL.txt", append=TRUE)
  # write(lineV, file="memorytrackSimuRL.txt", append=TRUE)
  
  # Parameters for PIVs
  
  # Parameter for the probability for a pair to be linked
  gamma = 0.05
  
  # Parameters for the distributions of the true values
  eta = lapply(data$Nvalues, function(x) rep(1/x,x))
  
  # Parameters for the survival model describing pSameH for unstable PIVs (over time and potentially more covariates)
  # Number of coefficients: covariates from A + covariates from B + time
  ncoef = lapply(seq_along(data$PIVs_stable), function(idx) if(data$PIVs_stable[idx]){ 0 }else{ ncol(data$A[, data$PIVs_config[[idx]]$pSameH.cov.A, drop=FALSE]) + ncol(data$B[, data$PIVs_config[[idx]]$pSameH.cov.B, drop=FALSE]) + 1 } )
  expalpha = lapply(seq_along(data$PIVs_stable), function(idx) if(data$PIVs_stable[idx]){ c(-Inf) }else{ rep(0.05, ncoef[idx]) })
  
  # Parameters for the registration errors (agreement, missing in A, missing in B)
  phi = lapply(data$Nvalues, function(x)  c(0.99,0.01,0.01))
  
  NmissingA = lapply(seq_along(data$Nvalues), function(k) sum(data$A[,data$PIVs][,k]==0))
  NmissingB = lapply(seq_along(data$Nvalues), function(k) sum(data$B[,data$PIVs][,k]==0))
  
  gamma.iter = array(NA, c(StEMIter, length(gamma)))
  eta.iter = lapply(data$Nvalues, function(x) array(NA, c(StEMIter, x)))
  phi.iter = lapply(data$Nvalues, function(x) array(NA, c(StEMIter, 3)))
  alpha.iter = lapply(ncoef, function(x) array(NA, c(StEMIter, x)))
  
  gamma.iter[1,] = gamma
  
  for(k in 1:length(data$Nvalues))
  {
    eta.iter[[k]][1,] = eta[[k]]
    alpha.iter[[k]][1,] = expalpha[[k]]
    phi.iter[[k]][1,] = phi[[k]]
  }
  
  # print(gamma)
  # print(eta)
  # print(phi)
  
  time.iter=c()
  if(trace==1)
    pb = progress_bar$new(format = "Running StEM algorithm [:bar] :percent in :elapsed",       total = StEMIter, clear = FALSE, width= 60)
  
  if(trace==2)
    cat("-----------------------\n")

  
  # MC-EM iteration
  for(iter in 2:StEMIter)
  {
    # print("ITER STEM")
    tijdM = Sys.time()
    
    # linePID = system("ps -C R")
    # lineV = system("ps v")
    # write(linePID, file="memorytrackSimuRL.txt", append=TRUE)
    # write(lineV, file="memorytrackSimuRL.txt", append=TRUE)
    
    D = Matrix(0, nrow=nrow(data$A), ncol=nrow(data$B), sparse=TRUE)
    linksR = which(D==1, arr.ind=TRUE)
    linksCpp = linksR
    sumRowD = rowSums(D)
    sumColD = colSums(D)
    cumulativeHazard = vector(mode='list', length=length(data$PIVs))
    nlinkrec = 0
    
    # if burnin value is 0, the algorithm will explore the necessary burnin for the number of linked records to stagnate
    countBurnin = 15
    Burnin_total = 0
    
    # Burn-in period Gibbs sampler
    if(GibbsBurnin == 0){
      pb_burnin_explore = progress_bar$new(format = "(:spin)     Burn-in period: :whatburnin     Number of linked records: :whatlinkrec", total = 300, clear = FALSE, width= 100)
      while(countBurnin != 0)
      {
        Burnin_total = Burnin_total + 1
        # 2 alterantives to simulate H in R / similar / simulateH is faster
        newTruePivs = simulateH_old(data=data, D=D, linksR=linksR, omegaData=omegaData, sumRowD=sumRowD, sumColD=sumColD, eta=eta, phi=phi)
        #newTruePivs = simulateH(data=data, D=D, links=linksCpp, omegaData=cumulativeHazard, sumRowD=sumRowD, sumColD=sumColD, eta=eta, phi=phi)
        truepivsA = newTruePivs$truepivsA
        truepivsB = newTruePivs$truepivsB
        Dsample = simulateD(data=data, D=D, linksR=linksR, sumRowD=sumRowD, sumColD=sumColD, truepivsA=truepivsA, truepivsB=truepivsB, gamma=gamma, eta=eta, expalpha=expalpha, phi=phi)
        
        D = Dsample$D
        linksCpp = Dsample$links
        linksR = linksCpp + 1
        sumRowD = Dsample$sumRowD
        sumColD = Dsample$sumColD
        loglikelihood = Dsample$loglik
        new_nlinkred = Dsample$nlinkrec

        if( abs(new_nlinkred - nlinkrec) < 10 ){ 
          countBurnin = countBurnin - 1
        }else{
          countBurnin = 15
          nlinkrec = new_nlinkred
        }
        
        if(instability){
          if(nrow(linksR)>0){
            # MODIF NEW
            #times = rep(1, nrow(links))
            times = abs(data$B[linksR[,2], "date"] - data$A[linksR[,1], "date"])
            times = cbind( times )
            cumulativeHazard = vector(mode='list', length=length(data$PIVs))
            if(conditionalHazard){
              for(k in 1:length(data$PIVs)){
                if(!data$PIVs_stable[k]){
                  cov = cbind( data$A[linksR[,1], data$PIVs_config[[k]]$pSameH.cov.A, drop=FALSE],
                               data$B[linksR[,2], data$PIVs_config[[k]]$pSameH.cov.B, drop=FALSE] )
                  cumulativeHazard[[k]] = as.matrix(times) %*% expalpha[[k]]
                  # what if we need to incorporate cov????? TODO
                }
              }
            }else{
              for(k in 1:length(data$PIVs)){
                if(!data$PIVs_stable[k]){
                  cumulativeHazard[[k]] = as.matrix(times) %*% expalpha[[k]]
                }
              }
            }
          }
        }
        pb_burnin_explore$tick(tokens = list(whatburnin = Burnin_total, whatlinkrec = nlinkrec))
        if(Burnin_total>=300){
          countBurnin = 0
          cat("The burn-in period exceeded 300 iterations, this is the default maximum burn-in in that case.\nIf you want to increase the burn-in period fix it with the parameter 'GibbsBurnin' of the stEM function.")
          pb_burnin_explore$terminate()
        }
      }
    }else{
      for(j in 1:GibbsBurnin)
      {
        # print("ITER GIBSS")
        # 2 alterantives to simulate H in R / similar / simulateH is faster
        # print("call simuH")
        newTruePivs = simulateH_old(data=data, D=D, linksR=linksR, omegaData=omegaData, sumRowD=sumRowD, sumColD=sumColD, eta=eta, phi=phi)
        #newTruePivs = simulateH(data=data, D=D, links=linksCpp, omegaData=cumulativeHazard, sumRowD=sumRowD, sumColD=sumColD, eta=eta, phi=phi)
        truepivsA = newTruePivs$truepivsA
        truepivsB = newTruePivs$truepivsB
        # print("call simuD")
        Dsample = simulateD(data=data, D=D, linksR=linksR, sumRowD=sumRowD, sumColD=sumColD, truepivsA=truepivsA, truepivsB=truepivsB, gamma=gamma, eta=eta, expalpha=expalpha, phi=phi)
        
        D = Dsample$D
        linksCpp = Dsample$links
        linksR = linksCpp + 1
        sumRowD = Dsample$sumRowD
        sumColD = Dsample$sumColD
        loglikelihood = Dsample$loglik
        nlinkred = Dsample$nlinkrec
        # print(linksR)
        
        if(instability){
          if(nrow(linksR)>0){
            # MODIF NEW
            #times = rep(1, nrow(links))
            times = abs(data$B[linksR[,2], "date"] - data$A[linksR[,1], "date"])
            times = cbind( times )
            cumulativeHazard = vector(mode='list', length=length(data$PIVs))
            if(conditionalHazard){
              for(k in 1:length(data$PIVs)){
                if(!data$PIVs_stable[k]){
                  cov = cbind( data$A[linksR[,1], data$PIVs_config[[k]]$pSameH.cov.A, drop=FALSE],
                               data$B[linksR[,2], data$PIVs_config[[k]]$pSameH.cov.B, drop=FALSE] )
                  cumulativeHazard[[k]] = as.matrix(times) %*% expalpha[[k]]
                  # what if we need to incorporate cov????? TODO
                }
              }
            }else{
              for(k in 1:length(data$PIVs)){
                if(!data$PIVs_stable[k]){
                  cumulativeHazard[[k]] = as.matrix(times) %*% expalpha[[k]]
                }
              }
            }
          }
        }
      }
    }

    # Administration for the M-step
    Vgamma = c()
    Veta = lapply(data$Nvalues, function(x) c())
    Vomega = mapply(createOmegadata, ncoef, data$PIVs_stable, SIMPLIFY=FALSE)
    Vphi = lapply(data$Nvalues, function(x) c())
    
    for(j in 1:nGibbsIter)
    {
      # 2 alterantives to simulate H in R / similar / simulateH is faster
      newTruePivs = simulateH_old(data=data, D=D, linksR=linksR, omegaData=omegaData, sumRowD=sumRowD, sumColD=sumColD, eta=eta, phi=phi)
      #newTruePivs = simulateH(data=data, D=D, links=linksCpp, omegaData=cumulativeHazard, sumRowD=sumRowD, sumColD=sumColD, eta=eta, phi=phi)
      truepivsA = newTruePivs$truepivsA
      truepivsB = newTruePivs$truepivsB
      # if(iter==29){
      #   print(eta)
      #   print(expalpha)
      #   print(phi)
      #   print(truepivsA)
      #   print(truepivsB)
      # }
      Dsample = simulateD(data=data, D=D, linksR=linksR, sumRowD=sumRowD, sumColD=sumColD, truepivsA=truepivsA, truepivsB=truepivsB, gamma=gamma, eta=eta, expalpha=expalpha, phi=phi)
      
      D = Dsample$D
      linksCpp = Dsample$links
      linksR = linksCpp + 1
      # if(iter==29){
      #   print(linksR)
      # }
      sumRowD = Dsample$sumRowD
      sumColD = Dsample$sumColD
      loglikelihood = Dsample$loglik
      nlinkrec = Dsample$nlinkrec
      
      if(instability){
        if(nrow(linksR)>0){
          # MODIF NEW
          #times = rep(1, nrow(links))
          times = abs(data$B[linksR[,2], "date"] - data$A[linksR[,1], "date"])
          times = cbind( times )
          cumulativeHazard = vector(mode='list', length=length(data$PIVs))
          if(conditionalHazard){
            for(k in 1:length(data$PIVs)){
              if(!data$PIVs_stable[k]){
                equal = truepivsA[linksR[,1],k] == truepivsB[linksR[,2],k]
                cov = cbind( data$A[linksR[,1], data$PIVs_config[[k]]$pSameH.cov.A, drop=FALSE],
                             data$B[linksR[,2], data$PIVs_config[[k]]$pSameH.cov.B, drop=FALSE] )
                Vtmp_k = cbind( equal,
                                times )
                Vomega[[k]] = rbind(Vomega[[k]], Vtmp_k)
                cumulativeHazard[[k]] = as.matrix(times) %*% expalpha[[k]]
                # what if we need to incorporate cov????? TODO
              }
            }
          }else{
            for(k in 1:length(data$PIVs)){
              if(!data$PIVs_stable[k]){
                equal = truepivsA[linksR[,1],k] == truepivsB[linksR[,2],k]
                Vtmp_k = cbind( equal,
                                times )
                Vomega[[k]] = rbind(Vomega[[k]], Vtmp_k)
                cumulativeHazard[[k]] = as.matrix(times) %*% expalpha[[k]]
              }
            }
          }
        }
      }
      
      # Update gamma
      Vgamma[j] = sum(D)
      if(sum(D)==0){
        print( "loglikelihood:" )
        print(loglikelihood)
        print( "gamma so far:" )
        print(gamma)
        print( "In the last iteration no link has been made." )
        print( "This is probably due to a difference in the support of some PIV between file A and file B." )
      }
      
      # Update eta
      for(k in 1:length(data$Nvalues))
      {
        facpivsA = factor(truepivsA[,k], levels=1:data$Nvalues[k])
        facpivsB = factor(truepivsB[,k], levels=1:data$Nvalues[k])
        Veta[[k]] = rbind(Veta[[k]],
                          table(facpivsA[sumRowD==0]) + table(facpivsB[sumColD==0]) + table(facpivsA[sumRowD==1]))
      }

      # Update phi: count agreements / missings
      for(k in 1:length(data$Nvalues)){
        Vphi[[k]] = rbind( Vphi[[k]], 
                           c( sum(truepivsA[,k]==data$A[,data$PIVs][,k]) + sum(truepivsB[,k]==data$B[,data$PIVs][,k]), 
                              NmissingA[[k]],
                              NmissingB[[k]]) )
      }
    }
    
    # Calculate new parameters gamma/eta/phi/omega 
    
    # New gamma
    gamma = sum(Vgamma) / (nGibbsIter*nrow(truepivsA))
    # if(iter==29){
    #   print(gamma)
    # }
    
    # New eta
    for(k in 1:length(data$Nvalues))
      eta[[k]] = colSums(Veta[[k]])/sum(Veta[[k]])
    
    # New omega
    if(nrow(linksR)>0){
      for(k in 1:length(data$Nvalues))
      {
        if(!data$PIVs_stable[k])
        {
          omegaInit = rep(0.05, ncoef[k])
          expalpha[[k]] = exp(nlminb(omegaInit, loglikOmega, control=list(trace=FALSE), dataOmega=Vomega[[k]])$par)
          # MODIF NEW
          # IREG 10%, PAR 10%, state 23%, reg 23%
          # if(data$PIVs[k]=="PAR"){
          #   expalpha[[k]] = 0.11
          # }
          # if(data$PIVs[k]=="IREG"){
          #   expalpha[[k]] = 0.05
          # }
          # if(data$PIVs[k]=="state"){
          #   expalpha[[k]] = 0.22
          # }
          # if(data$PIVs[k]=="reg"){
          #   expalpha[[k]] = 0.11
          # }
          # if(data$PIVs[k]=="IREG"){
          #   omegaInit = rep(0.05, ncoef[k])
          #   expalpha[[k]] = exp(nlminb(omegaInit, loglikOmega, control=list(trace=FALSE), dataOmega=Vomega[[k]])$par)
          # }
        }
      } 
    }
    
    # New phi
    for(k in 1:length(data$Nvalues))
    {
      Ntotal = nGibbsIter * (nrow(data$A) + nrow(data$B))
      phi[[k]][1] = sum(Vphi[[k]][,1]) / (Ntotal - sum(Vphi[[k]][,2]) - sum(Vphi[[k]][,3])) # max( 0.9, sum(Vphi[[k]][,1]) / (Ntotal - sum(Vphi[[k]][,2]) - sum(Vphi[[k]][,3])) )
      phi[[k]][2] = sum(Vphi[[k]][,2]) / (nGibbsIter*nrow(data$A))
      phi[[k]][3] = sum(Vphi[[k]][,3]) / (nGibbsIter*nrow(data$A))
    } 
    if(UnstableNoMistake){
      # There is instability but no mistake so: no identifiability problem
      # there should be null registration time difference in links 
      if(instability){
        for(i in 1:length(unstablePIVs)){
          phi[[ unstablePIVs[i] ]][1] = 1 # proba agreement set to 1 ie proba mistake at 0 for unstable PIV
        }
      }
    }
    
    # Administration
    
    gamma.iter[iter,] = gamma

    for(k in 1:length(data$Nvalues))
    {
      eta.iter[[k]][iter,] = eta[[k]]
      alpha.iter[[k]][iter,] = expalpha[[k]]
      phi.iter[[k]][iter,] = phi[[k]]
    }
    
    if(trace==1)
      pb$tick()
    
    if(trace==2)
    {
      cat("-------------------------------------------------------------------------\n")
      cat("Iter:", iter, "\n")
      
      time.iter = time.iter + as.numeric(difftime(Sys.time(), tijdM, units="min"))
      cat("Approx time remaining:\t", (StEMIter-iter) * time.iter/iter , "mins\n")
      
      cat("gamma:\t", gamma, "\n")
      cat("phi:")
      for(i in 1:length(data$Nvalues))
        cat("\t\t", phi[[i]], "\n")
      
      cat("eta:")
      for(i in 1:length(data$Nvalues))
        cat("\t\t",  eta[[i]], "\n") 
      
      cat("omega:")
      for(i in 1:length(data$Nvalues))
        cat("\t\t",  expalpha[[i]], "\n") 
    }
    
    # save current environment, global variables and local variables
    save.image(file=file.path(data$newDirectory, 'myEnvironment.RData'))
    save(iter, gamma, eta, expalpha, phi, gamma.iter, eta.iter, phi.iter, alpha.iter, file=file.path(data$newDirectory, 'myEnvironmentLocal.RData'))
  }

  if(trace==1)
    pb$terminate()
  
  Delta = matrix(0, nrow=nrow(data$A), ncol=nrow(data$B))
  
  gamma_avg = apply(gamma.iter, 2, function(x) mean(x[StEMBurnin:StEMIter , drop=FALSE]))
  eta_avg = lapply(eta.iter, function(x) apply(x[StEMBurnin:StEMIter, , drop=FALSE], 2, mean))
  expalpha_avg = lapply(alpha.iter, function(x) apply(x[StEMBurnin:StEMIter, , drop=FALSE], 2, mean))
  phi_avg = lapply(phi.iter, function(x) apply(x[StEMBurnin:StEMIter, , drop=FALSE], 2, mean))
  
  pbfinal = progress_bar$new(format = "Drawing Delta          [:bar] :percent in :elapsed",       total = 1000, clear = FALSE, width= 60)
  
  for(m in 1:1000)
  {
    
    # linePID = system("ps -C R")
    # lineV = system("ps v")
    # write(linePID, file="memorytrackSimuRL.txt", append=TRUE)
    # write(lineV, file="memorytrackSimuRL.txt", append=TRUE)
    
    # 2 alterantives to simulate H in R / similar / simulateH is faster
    newTruePivs = simulateH_old(data=data, D=D, linksR=linksR, omegaData=omegaData, sumRowD=sumRowD, sumColD=sumColD, eta=eta_avg, phi=phi_avg)
    #newTruePivs = simulateH(data=data, D=D, links=linksCpp, omegaData=cumulativeHazard, sumRowD=sumRowD, sumColD=sumColD, eta=eta_avg, phi=phi_avg)
    truepivsA = newTruePivs$truepivsA
    truepivsB = newTruePivs$truepivsB
    Dsample = simulateD(data=data, D=D, linksR=linksR, sumRowD=sumRowD, sumColD=sumColD, truepivsA=truepivsA, truepivsB=truepivsB, gamma=gamma_avg, eta=eta_avg, expalpha=expalpha_avg, phi=phi_avg)

    D = Dsample$D
    linksCpp = Dsample$links
    linksR = linksCpp+1
    sumRowD = Dsample$sumRowD
    sumColD = Dsample$sumColD
    loglikelihood = Dsample$loglik
    nlinkrec = Dsample$nlinkrec
    
    if(instability){
      if(nrow(linksR)>0){
        # MODIF NEW
        #times = rep(1, nrow(links))
        times = abs(data$B[linksR[,2], "date"] - data$A[linksR[,1], "date"])
        times = cbind( times )
        cumulativeHazard = vector(mode='list', length=length(data$PIVs))
        if(conditionalHazard){
          for(k in 1:length(data$PIVs)){
            if(!data$PIVs_stable[k]){
              cov = cbind( data$A[linksR[,1], data$PIVs_config[[k]]$pSameH.cov.A, drop=FALSE],
                           data$B[linksR[,2], data$PIVs_config[[k]]$pSameH.cov.B, drop=FALSE] )
              cumulativeHazard[[k]] = as.matrix(times) %*% expalpha_avg[[k]]
              # what if we need to incorporate cov????? TODO
            }
          }
        }else{
          for(k in 1:length(data$PIVs)){
            if(!data$PIVs_stable[k]){
              cumulativeHazard[[k]] = as.matrix(times) %*% expalpha_avg[[k]]
            }
          }
        }
      }
    }

    Delta = Delta + D
    pbfinal$tick()
  }
  Delta = Delta / 1000
  #print(truepivsA)
  #print(truepivsB)
  
  pbfinal$terminate()
  
  # browseURL('https://www.youtube.com/watch?v=NTa6Xbzfq1U')
  
  list(Delta=Delta, gamma=gamma.iter, eta=eta.iter, alpha=alpha.iter, phi=phi.iter)
}