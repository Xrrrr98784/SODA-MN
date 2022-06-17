#sample ms1 data remove redundancy
remove.ms1.feature.duplicates<-function(ms1.data,mz.tol,rt.tol){
  ms1.data<- ms1.data[order(ms1.data$RT), ]
  ms1.data<- ms1.data[order(ms1.data$MZ), ]
  i <- 1
  while (i < nrow(ms1.data)) { 
    j=i+1
    while (j < nrow(ms1.data)+1){
      score.mz <- abs(ms1.data$MZ[j]-ms1.data$MZ[i])*10e5/ms1.data$MZ[i]
      if (score.mz<mz.tol*10){
        score.rt <- abs(ms1.data$RT[j]-ms1.data$RT[i])
        if (score.mz<mz.tol && score.rt<rt.tol ){ 
          if (ms1.data$Int[i] < ms1.data$Int[j]) {
            ms1.data=ms1.data[-c(i),]
            i=i-1
            break
          } else {
            ms1.data=ms1.data[-c(j),]
            j=j-1
          } 
        }
      }else{
        break
      }
      j=j+1
    }
    i=i+1
  }
  return(ms1.data)
}
  