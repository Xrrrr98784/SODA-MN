standard.data.preparation=function(data.path,code.path,file.name,format){
#read standard

library(CluMSID)
path <- paste(data.path,file.name,sep="")
setwd(path)
ms1.standard <- read.csv('standard.csv',header = TRUE, sep = ",", quote = "\"",
                         dec = ".", fill = TRUE, comment.char = "")


ms2.standard <- grep(pattern = format, dir(path), value = TRUE)
standard.spec <- list()

for (i in 1:length(ms2.standard)) {
  standard.spec[[i]] <- extractMS2spectra(MSfile=ms2.standard[i], min_peaks = 2, recalibrate_precursor = FALSE,RTlims = NULL)
}

#Find the most matching standard for each single standard file
setwd(code.path)
source('remove.node.duplicates.R')
mz.tol <- 10
rt.tol <- 60
standard.MS1.MS2 <- data.frame()
for (i in 1:nrow(ms1.standard)){
  for (j in 1:length(standard.spec[[i]])){
    mz.ms1 <- ms1.standard[i,1]
    rt.ms1 <- ms1.standard[i,2]
    mz.ms2 <- standard.spec[[i]][[j]]@precursor
    rt.ms2 <- standard.spec[[i]][[j]]@rt
    spec_t <- max(standard.spec[[i]][[j]]@spectrum[,2])
    score.mz <- abs(mz.ms1-mz.ms2)*10e5/mz.ms1
    score.rt <- abs(rt.ms1-rt.ms2)
    if (score.mz < mz.tol && score.rt < rt.tol){
      t=nrow(standard.MS1.MS2)+1
      standard.MS1.MS2[t,1] <- 0  
      standard.MS1.MS2[t,2] <- i
      standard.MS1.MS2[t,3] <- j
      standard.MS1.MS2[t,4] <- ms1.standard[i,3]
      standard.MS1.MS2[t,5] <- mz.ms1
      standard.MS1.MS2[t,6] <- mz.ms2
      standard.MS1.MS2[t,7] <- score.mz
      standard.MS1.MS2[t,8] <- rt.ms1
      standard.MS1.MS2[t,9] <- rt.ms2
      standard.MS1.MS2[t,10] <- score.rt
      standard.MS1.MS2[t,11] <- "NA"
      standard.MS1.MS2[t,12] <- spec_t
    }
  }
}
colnames(standard.MS1.MS2)=c("round","node1","node2","name","mz1","mz2","mz.score","rt1","rt2","rt.score","spec.score","intensity")
seed.0 <- remove.node.duplicates(standard.MS1.MS2,prior="intensity")

seed.0.ms2 <- list()
for (i in 1:nrow(seed.0)) {
  seed.0.ms2[[i]] <- standard.spec[[seed.0[i,2]]][[seed.0[i,3]]]
  seed.0.ms2[[i]]@annotation <- c(seed.0[i,4])
}

setwd(data.path)
save(seed.0,file="seed.0.ms1.rda")
save(seed.0.ms2,file="seed.0.ms2.rda")
}