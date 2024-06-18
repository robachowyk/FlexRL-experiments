##### FOR THE NLTCS DATA

library(MASS)

write.matrix(data,"data_sparseMat.txt",sep="\t")  

# set directory
deltahat = c(read.csv("DeltaVector_only41110_miss195212378.csv"))
# > length(test[test<=0.0001])
# [1] 195220714
# > length(test[test>0.0001])
# [1] 32774
histogram = hist(deltahat$x)
histogram$counts
histogram$counts[1] = histogram$counts[1] + 195212378
histogram$counts
histogram$counts = log10(histogram$counts)
plot(histogram)

# estimated FDR
fdrhat = c() # 1 - true predictions / (total predicitons)
for(threshold in seq(0.5, 1, length=100)){
  fdrhat = c( fdrhat, 1 - sum( deltahat$x[deltahat$x>threshold] ) / sum( deltahat$x>threshold ) )
}
plot(seq(0.5, 1, length=100), fdrhat, ylab="hat FDR", xlab="threshold", type="l")

# true FDR = 1 - precision
# displayed for the different thresholds
lambda = c(0.5,0.6,0.7,0.8,0.9)
trueFDR = c(0.10, 0.078, 0.069, 0.066, 0.04)

################################################################################
#           GRAPH FDR SENSITIVITY FLEX RL COMPLETE NLTCS DATA SETS             #
################################################################################

# FDR HAT for FlexRL_state
deltahat = c(read.csv("DeltaVector_only50896_miss195202592.csv"))
histogram = hist(deltahat$x)
histogram$counts
histogram$counts[1] = histogram$counts[1] + 195202592
histogram$counts
histogram$counts = log10(histogram$counts)
plot(histogram)
fdrhat = c() # 1 - true predictions / (total predicitons)
for(threshold in seq(0.5, 1, length=100)){
  fdrhat = c( fdrhat, 1 - sum( deltahat$x[deltahat$x>threshold] ) / sum( deltahat$x>threshold ) )
}
fdrhat
which(fdrhat<0.10)
seq(0.5, 1, length=100)[73]
# fdrhat < 10%: 0.0977 --> threshold lambda = 0.86

# FDR HAT for FlexRL_state_reg
deltahat = c(read.csv("DeltaVector_only32774_miss195220714.csv"))
histogram = hist(deltahat$x)
histogram$counts
histogram$counts[1] = histogram$counts[1] + 195220714
histogram$counts
histogram$counts = log10(histogram$counts)
plot(histogram)
fdrhat = c() # 1 - true predictions / (total predicitons)
for(threshold in seq(0.5, 1, length=100)){
  fdrhat = c( fdrhat, 1 - sum( deltahat$x[deltahat$x>threshold] ) / sum( deltahat$x>threshold ) )
}
fdrhat
which(fdrhat<0.10)
seq(0.5, 1, length=100)[28]
# fdrhat < 10%: 0.0996 --> threshold lambda = 0.64

# FDR HAT for FlexRL_statereg_glued
deltahat = c(read.csv("DeltaVector_only98336_miss188832887.csv"))
histogram = hist(deltahat$x)
histogram$counts
histogram$counts[1] = histogram$counts[1] + 188832887
histogram$counts
histogram$counts = log10(histogram$counts)
plot(histogram)
fdrhat = c() # 1 - true predictions / (total predicitons)
for(threshold in seq(0.5, 1, length=100)){
  fdrhat = c( fdrhat, 1 - sum( deltahat$x[deltahat$x>threshold] ) / sum( deltahat$x>threshold ) )
}
fdrhat
which(fdrhat<0.10)
seq(0.5, 1, length=100)[70]
# fdrhat < 10%: 0.099 --> threshold lambda = 0.85

tikz('testfdrplot.tex', width = 3.25, height = 3.25)
# New Graph, true sensitivity and FDR for complete NLTCS
lambda_labels_FlexRL_statereg_glued = c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9)
true_sensitivity_FlexRL_statereg_glued = c(0.21647406611126, 0.19067454985219, 0.162456328943832, 0.134909970438054, 0.113141628594464, 0.0889545821015856, 0.0686643375436711, 0.0427304488040849, 0.0210964794410105)
true_fdr_FlexRL_statereg_glued = c(0.109452736318408, 0.0990476190476191, 0.0923423423423423, 0.0897552130553038, 0.0827886710239651, 0.0805555555555556, 0.0606617647058824, 0.0701754385964912, 0.0710059171597633)
TP_FlexRL_FlexRL_statereg_glued = c(1611, 1419, 1209, 1004, 842, 662, 511, 318, 157) # , 28
plot(true_sensitivity_FlexRL_statereg_glued, true_fdr_FlexRL_statereg_glued, type="l", col=2, xlim=c(0,1), ylim=c(0,1), xlab="Sensitivity", ylab="FDR", lty=2)

lambda_labels_FlexRL_statereg = c(0.5, 0.6, 0.7, 0.8, 0.9)
true_sensitivity_FlexRL_statereg = c(0.4269574, 0.3775617, 0.3384130, 0.3030741, 0.2673410)
true_fdr_FlexRL_statereg = c(0.4869771, 0.4342520, 0.3887043, 0.3481209, 0.3064076)
lines(true_sensitivity_FlexRL_statereg, true_fdr_FlexRL_statereg, type="l", col=3, lty=3)

lambda_labels_FlexRL_stateonly = c(0.5, 0.6, 0.7, 0.8, 0.9)
true_sensitivity_FlexRL_stateonly = c(0.1845770, 0.1382028, 0.09774041, 0.05977404, 0.01523910)
true_fdr_FlexRL_stateonly = c(0.1016624, 0.08759757, 0.07347447, 0.05601660, 0.05691057)
lines(true_sensitivity_FlexRL_stateonly, true_fdr_FlexRL_stateonly, type="l", col=4, lty=1)

points(0.904864283794679, 0.641865659735149, pch = 12, col=2) # simplistic_statereg_glued
points(0.884655806621125, 0.641865659735149, pch = 0, col=3) # simplistic_state_reg
points(0.904755648975302, 0.684964091304149, pch = 15, col=4) # simplistic_state

points(0.21647406611126, 0.109452736318408, pch = 10, col=2) # FlexRL_statereg_glued
points(0.426957435627956, 0.486977111286504, pch = 1, col=3) # FlexRL_state_reg
points(0.184576983709932, 0.101662404092072, pch = 16, col=4) # FlexRL_state

points(0.04, 0.07, pch = 8, col=2) # FlexRLcontrolled_statereg_glued
points(0.36, 0.41, pch = 8, col=3) # FlexRLcontrolled_state_reg
points(0.038, 0.056445, pch = 8, col=4) # FlexRLcontrolled_state
dev.off()

legend("bottomright", legend = c("simplisitc", "FlexRL default 0.5", "FlexRL hatFDR controlled < 10%","state and reg glued", "state and reg", "state only"),
       pch = c(15, 16, 8, 10, 1,16), 
       lty = c(NA, NA, NA, 2, 3, 1),
       col = c(1,1,1,2,3,4), 
       text.col = c(1,1,1,2,3,4))
dev.off()


##### FOR THE SHIW DATA

# set directory
deltahat = c(read.csv("DeltaVector_only42127_miss245267938.csv"))

# estimated FDR
fdrhat = c() # 1 - true predictions / (total predicitons)
for(threshold in seq(0.5, 1, length=100)){
  fdrhat = c( fdrhat, 1 - sum( deltahat$x[deltahat$x>threshold] ) / sum( deltahat$x>threshold ) )
}
plot(seq(0.5, 1, length=100), fdrhat, ylab="hat FDR", xlab="threshold", type="l")
fdrhat
which(fdrhat<0.10)
seq(0.5, 1, length=100)[23]
# fdrhat < 10%: 0.0977 --> threshold lambda = 0.61