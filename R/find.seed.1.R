find.seed.1<-function(seed.0,seed.0.ms2,ms1.data,ms2.spec,mz.tol,rt.tol,spec.tol,match.tol,filter,filter.ratio,remove.duplicates,prior){
t=0
result.0 <- data.frame()
for (i in 1:nrow(seed.0)){
  for (j in 1:length(ms2.spec)){
    mz.sd <- seed.0[i,6]
    mz.ep <- ms2.spec[[j]]@precursor
    score.mz <- abs(mz.sd-mz.ep)*10e5/mz.sd
    if (score.mz < mz.tol){
      rt.sd <- seed.0[i,9]
      rt.ep <- ms2.spec[[j]]@rt
      score.rt <- abs(rt.sd-rt.ep)
      
      if (score.rt < rt.tol){
        spec.sd <- seed.0.ms2[[i]]@spectrum
        spec.ep <- ms2.spec[[j]]@spectrum
        
        if (filter==TRUE){
        spec.sd.high<-arrange(filter(as.data.frame(spec.sd),V2>max(spec.sd[,2])*filter.ratio),desc(V2))
        spec.ep.high<-arrange(filter(as.data.frame(spec.ep),V2>max(spec.ep[,2])*filter.ratio),desc(V2))
        }else{
          spec.sd.high<-arrange(as.data.frame(spec.sd),desc(V2))
          spec.ep.high<-arrange(as.data.frame(spec.ep),desc(V2))
        }
        
        match=0
        spec.sd.sort=data.frame()
        spec.ep.sort=data.frame()
        if (nrow(spec.sd.high)!=0 && nrow(spec.ep.high)!=0){
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
     
          
          if (score.spec >spec.tol) {
            t=nrow(result.0)+1
            result.0 [t,1] <- 0
            result.0 [t,2] <- i
            result.0 [t,3] <- j
            result.0 [t,4] <- seed.0[i,4]
            result.0 [t,5] <- mz.sd
            result.0 [t,6] <- mz.ep
            result.0 [t,7] <- score.mz
            result.0 [t,8] <- rt.sd
            result.0 [t,9] <- rt.ep
            result.0 [t,10] <- score.rt
            result.0 [t,11] <- score.spec
            result.0 [t,12] <- (1-score.mz/mz.tol)*0.25+(1-score.rt/rt.tol)*0.25+score.spec*0.5
            result.0 [t,13] <- spec.ep.high[1,2]
            result.0 [t,14] <- match
          }
        }
      }
    }
    if (j%%1000==0) {
      print(paste("It is scanning", i,"/",nrow(seed.0),"seed(s),",j,"/",length(ms2.spec),"spectrum,",t,"neighbor(s) are found."))
    }
  }
}
colnames(result.0)=c("round","node1","node2","name","mz1","mz2","mz.score","rt1","rt2","rt.score","spec.score","overall.score","intensity","match")
if (remove.duplicates==TRUE){
  seed.1 <- remove.node.duplicates(result.0,prior=prior)
}else{
  seed.1  <- result.0
}
return(seed.1)
}
