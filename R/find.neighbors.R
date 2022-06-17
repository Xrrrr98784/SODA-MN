#find.neighbors.new
find.neighbors<- function(seed,ms2.data,mz.tol,match.tol,mz.weight,rt.weight,spec.weight,adduct,round,features){
  t=0
  result <- data.frame()
  for (i in 1:nrow(seed)){
    for (j in 1:length(ms2.data)){
      if (j %in% features[,1]==FALSE){
        mz.sd <- seed[i,6]
        mz.ep <- ms2.data[[j]]@precursor
        
        for (k in 1:nrow(adduct)) {
          if (mz.ep-mz.sd>0) {
            score.mz <- abs((mz.ep-mz.sd-adduct[k,2])*10e5/(adduct[k,2]+mz.sd))
          } else{
            # break
            score.mz <- abs((mz.sd-mz.ep-adduct[k,2])*10e5/(mz.sd-adduct[k,2]))
          }
          if (score.mz < mz.tol){
            
            spec.sd <- ms2.data[[seed[i,3]]]@spectrum
            spec.ep <- ms2.data[[j]]@spectrum
            spec.sd.high<-filter(as.data.frame(spec.sd),V2>max(spec.sd[,2])/100)
            spec.ep.high<-filter(as.data.frame(spec.ep),V2>max(spec.ep[,2])/100)
            
            match=0
            if (nrow(spec.sd.high)!=0 && nrow(spec.ep.high)!=0){
              spec.sd.sort=data.frame()
              spec.ep.sort=data.frame()
              for (m in 1:nrow(spec.sd.high)){
                for (n in 1:nrow(spec.ep.high)){
                  if (abs(spec.sd.high[m,1]-spec.ep.high[n,1])*10e5/spec.sd.high[m,1]<mz.tol){
                    match=match+1
                    spec.sd.sort[nrow(spec.sd.sort)+1,1:2]=spec.sd.high[m,1:2]
                    spec.ep.sort[nrow(spec.ep.sort)+1,1:2]=spec.ep.high[n,1:2]
                  }
                }
              }
            }
            
            if (match>match.tol){
              w.sd <- (spec.sd.sort[,1])*spec.sd.sort[,2]/max(spec.sd.sort[,2])#
              w.ep <- (spec.ep.sort[,1])*spec.ep.sort[,2]/max(spec.ep.sort[,2])#
              score.spec <- sum(w.sd*w.ep)/sqrt(sum(w.ep*w.ep)*sum(w.sd*w.sd))
              
              if (score.spec >0.7) {
                rt.sd <- seed[i,9]
                rt.ep <- ms2.data[[j]]@rt
                score.rt <- abs(rt.ep-rt.sd)
                t=nrow(result)+1
                result [t,1] <- round
                result [t,2] <- seed[i,3]
                result [t,3] <- j
                result [t,4] <- adduct[k,1]
                result [t,5] <- mz.sd
                result [t,6] <- mz.ep
                result [t,7] <- score.mz
                result [t,8] <- rt.sd
                result [t,9] <- rt.ep
                result [t,10] <- score.rt
                result [t,11] <- score.spec
                result [t,12] <- (1-score.mz/mz.tol)*mz.weight+(1-score.rt/rt.sd*0.3)*rt.weight+score.spec*spec.weight
                result [t,13] <- max(spec.ep.sort[,2])
                result [t,14] <- match
              }
            }
          }
        }
      }
      if (j%%1000==0) {
        print(paste("It is", r,"round(s) scanning", i,"/",nrow(seed),"seed(s),",j,"/",length(ms2.data),"spectrum,",t,"neighbor(s) are found."))
      }
    }
  }
  colnames(result)=c("round","node1","node2","name","mz1","mz2","mz.score","rt1","rt2","rt.score","spec.score","overall.score","intensity","match")
  return(result)
}