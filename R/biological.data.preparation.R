biological.data.preparation<-function(data.path,code.path,file.name,format){
#read ms2 spectra for biological replicates

library(CluMSID)
path <- paste(data.path,file.name,sep="")
setwd(path)
file.list <- grep(pattern = format, dir(path), value = TRUE)
ms2.spec.single <- list()
for (i in 1:length(file.list)) {
  ms2.spec.single[[i]] <- extractMS2spectra(MSfile=file.list[i], min_peaks = 2, recalibrate_precursor = FALSE,RTlims = NULL)
}

#########################################
print(paste(file.name,"has",length(file.list),"replicates,containing",length(ms2.spec.single[[1]]),length(ms2.spec.single[[2]]),
length(ms2.spec.single[[3]]),length(ms2.spec.single[[4]]),"features respectively",sep=" "))
#########################################

#Extract spectra mz rt and int to a ms1 table

mz.tol <- 10
rt.tol <- 60
ms1.data.single <- list()
for (j in 1:length(ms2.spec.single)){
  ms1.data.single [[j]]=data.frame()
  #ms1.data.single [[j]]=matrix("",length(ms2.spec.single[[j]]),4)
  # using data.frame or matrix give different number of result in next step due to different digit round
  for (i in 1:length(ms2.spec.single[[j]])){
    ms1.data.single [[j]] [i,1] <- i
    ms1.data.single [[j]] [i,2] <- ms2.spec.single[[j]] [[i]]@precursor
    ms1.data.single [[j]] [i,3] <- ms2.spec.single[[j]] [[i]]@rt
    ms1.data.single [[j]] [i,4] <- max(ms2.spec.single[[j]] [[i]]@spectrum[,2])
  }
  ms1.data.single[[j]]=as.data.frame(ms1.data.single[[j]])
  colnames(ms1.data.single[[j]])=c("No","MZ","RT","Int")
}


#remove duplicates for each biological replicates using ms1 table

setwd(code.path)
source('remove.ms1.feature.duplicates.R')
ms1.data.single.unique=list()
for (i in 1:length(ms1.data.single)){
  ms1.data.single.unique[[i]]=remove.ms1.feature.duplicates(ms1.data.single[[i]],mz.tol=10,rt.tol=60)
}

#########################################
print(paste(file.name,"replicates after removing duplicates separately,contain",nrow(ms1.data.single.unique[[1]]),nrow(ms1.data.single.unique[[2]]),
            nrow(ms1.data.single.unique[[3]]),nrow(ms1.data.single.unique[[4]]),"unique features respectively",sep=" "))
#########################################

#find the uinque of all features summed from all biological replicates as alignment reference

ms1.data.all.unique<-data.frame()
for (i in 1:length(ms1.data.single.unique)){
  ms1.data.all.unique=rbind(ms1.data.all.unique,cbind(rep(i,nrow(ms1.data.single.unique[[i]]),1),ms1.data.single.unique[[i]]))
}
colnames(ms1.data.all.unique)=c("Rep","No","MZ","RT","Int")
ms1.data.all.ref<-remove.ms1.feature.duplicates(ms1.data.all.unique,mz.tol=10,rt.tol=60)

#########################################
print(paste(file.name,"replicates after merging,contain",nrow(ms1.data.all.unique),"features; after removing duplicates together,contain",
            nrow(ms1.data.all.ref),"unique features",sep=" "))
#########################################

#align features for biological replicates

ms1.data.align<-cbind(ms1.data.all.ref,matrix(0,nrow(ms1.data.all.ref),length(ms1.data.single.unique)))
colnames(ms1.data.align)=c("Rep","No","MZ.max","RT.max","Int.max",file.list)
ms1.data.align.count<-ms1.data.align
ms1.data.align.mz<-ms1.data.align
ms1.data.align.rt<-ms1.data.align
ms1.data.align.int<-ms1.data.align

for (i in 1:length(ms1.data.single.unique)){
  for (n in 1:nrow(ms1.data.align)){
    for (j in 1:nrow(ms1.data.single.unique[[i]])){
      score.mz=10e5*abs(ms1.data.align$MZ[n]-ms1.data.single.unique[[i]]$MZ[j])/ms1.data.align$MZ[n]
      if (score.mz<mz.tol){
        score.rt=abs(ms1.data.align$RT[n]-ms1.data.single.unique[[i]]$RT[j])
        if (score.rt<rt.tol){
          ms1.data.align.count[n,i+5]<-ms1.data.align[n,i+5]+1
          ms1.data.align.mz[n,i+5]<-ms1.data.single.unique[[i]]$MZ[j]
          ms1.data.align.rt[n,i+5]<-ms1.data.single.unique[[i]]$RT[j]
          ms1.data.align.int[n,i+5]<-ms1.data.single.unique[[i]]$Int[j]
        }
      }
    }
  }
}


#filter out features by cv

sum.count=as.data.frame(apply(as.data.frame(ms1.data.align.count[,6:ncol(ms1.data.align.count)]),1,sum))
colnames(sum.count)="sum.count"
int.cv=as.data.frame(apply(as.data.frame(ms1.data.align.int[,6:ncol(ms1.data.align.int)]),1,sd)/apply(as.data.frame(ms1.data.align.int[,6:ncol(ms1.data.align.int)]),1,mean))
colnames(int.cv)="int.cv"
ms1.data.align.int.clear=ms1.data.align.int[sum.count$sum.count==length(file.list) & int.cv$int.cv<0.3,]
ms1.data.align.rt.clear=ms1.data.align.rt[sum.count$sum.count==length(file.list) & int.cv$int.cv<0.3,]
ms1.data.align.mz.clear=ms1.data.align.mz[sum.count$sum.count==length(file.list) & int.cv$int.cv<0.3,]


ms1.output=ms1.data.align.int.clear
ms1.output$MZ.ave=apply(as.data.frame(ms1.data.align.mz.clear[,6:ncol(ms1.data.align.mz.clear)]),1,mean)
ms1.output$RT.ave=apply(as.data.frame(ms1.data.align.rt.clear[,6:ncol(ms1.data.align.int.clear)]),1,mean)
ms1.output$Int.ave=apply(as.data.frame(ms1.data.align.int.clear[,6:ncol(ms1.data.align.int.clear)]),1,mean)

#########################################
print(paste(file.name,"after after alignment and filtering,contain",nrow(ms1.output),"unique features",sep=" "))
#########################################


setwd(code.path)
source('remove.ms2.feature.duplicates.R')
ms2.output=remove.ms2.feature.duplicates(ms1.output,ms2.spec.single)

Index=1:nrow(ms1.output)
setwd(data.path)
output=list()
output[[1]]=cbind(Index,ms1.output)
output[[2]]=ms2.output
names(output)=c(paste(file.name,".ms1",sep=""),paste(file.name,".ms2",sep=""))
assign(file.name,output)
save(list=file.name,file=paste(file.name,".rda",sep=""))
}