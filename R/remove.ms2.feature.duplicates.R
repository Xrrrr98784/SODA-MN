#find unique sample ms2 data
remove.ms2.feature.duplicates<-function(ms1.data,ms2.spec){
  ms2.data <- list()
  for (i in 1:nrow(ms1.data)) {
    ms2.data[[i]] <- ms2.spec[[ms1.data$Rep[i]]][[ms1.data$No[i]]]
  }
  return(ms2.data)
}