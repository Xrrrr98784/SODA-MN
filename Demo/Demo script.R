setwd("/Users/ruixu/Desktop/SODA-MN/R")
source('biological.data.preparation.R')
source('standard.data.preparation.R')
source('find.seed.1.R')
source('remove.node.duplicates.R')
source("find.neighbors.r")
source("remove.previous.found.r")
source("output.network.files.r")
library(dplyr)

#Read and prepare biological sample data
file.name="Gut_Bac_BRB_1h"
setwd("/Users/ruixu/Desktop/SODA-MN/Demo")
code.path="/Users/ruixu/Desktop/SODA-MN/R/"
data.path="/Users/ruixu/Desktop/SODA-MN/Demo/"

biological.data.preparation(data.path=data.path,code.path=code.path,file.name=file.name,format="mzML")
#biological.data.preparation(data.path=data.path,code.path=code.path,file.name="Gut_Bac_BRB_8h",format="mzML")
#biological.data.preparation(data.path=data.path,code.path=code.path,file.name="Gut_Bac_BRB_24h",format="mzML")

#Read and prepare standard sample data
standard.data.preparation(data.path=data.path,code.path=code.path,file.name="Standards",format="mzXML")

#Read adduct data
adduct <- read.csv('adduct.csv',header = TRUE, sep = ",", quote = "\"",dec = ".", fill = TRUE, comment.char = "")

#load data
load("seed.0.ms1.rda")
load("seed.0.ms2.rda")
load(paste(file.name,".rda",sep=""))
ms1.data=Gut_Bac_BRB_1h[[1]]
ms2.spec=Gut_Bac_BRB_1h[[2]]

#Find seed 1
seed.1=find.seed.1(seed.0=seed.0,             #27 standards ms1 data
                   seed.0.ms2=seed.0.ms2,       #27 standards ms2 data
                   ms1.data=ms1.data,         #QC.no.filter sample data ms1 data
                   ms2.spec=ms2.spec,         #QC.no.filter sample data ms2 data
                   mz.tol=10,                   #mz tolerance
                   rt.tol=600,                   #rt tolerance
                   spec.tol=0.7,                #cosine score tolerance
                   match.tol=3,                 #peaks matched
                   filter=TRUE,                 #filter noise or not
                   filter.ratio=1/100,         #the ratio to highest peak as noise cut-off
                   remove.duplicates=TRUE,
                   prior="intensity")      #remove duplicates or not

#find neighbor compounds
result.individual=list()
seed.individual=list()
seed.individual.pure=list()
seed.individual[[1]]=seed.1
seed.individual.pure[[1]]=seed.1
seed.all.pure=data.frame()
seed.all=data.frame()
for (r in 1:5){
features.has.been.reported <- as.data.frame(matrix(0,1,1))
result.individual[[r]] <- find.neighbors(seed=seed.individual.pure[[r]],       
                                         ms2.data=ms2.spec,
                                         mz.tol=10,
                                         match.tol=3,
                                         mz.weight=0.5,
                                         rt.weight=0,
                                         spec.weight=0.5,
                                         adduct=adduct,
                                         round=r,
                                         features=features.has.been.reported)
seed.individual[[r+1]] <- remove.node.duplicates(result.individual[[r]],prior="match")
if (r==1){
  seed.individual.pure[[r+1]] <- seed.individual[[r+1]] 
}else{
  seed.individual.pure[[r+1]] <- remove.previous.found(pre=seed.all,new=seed.individual[[r+1]])
}
if (nrow(seed.individual.pure[[r+1]])==0) {
  max.round=r
  break
}else{
  max.round=r
}

seed.all.pure <- rbind(seed.all.pure,seed.individual.pure[[r+1]])
seed.all <- rbind(seed.all,seed.individual[[r+1]])
}

output.network.files(seed.all.pure,file.name=file.name)
