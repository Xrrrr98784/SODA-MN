output.network.files<-function(seed.all.pure,file.name){
edge.name <- paste(file.name,".edge.csv",sep="")
edge <- matrix(ncol=10,nrow=nrow(seed.all.pure))
edge[,1] <- seed.all.pure[,2]
edge[,2] <- seed.all.pure[,3]
edge[,3] <- paste(seed.all.pure[,2]," (-) ",seed.all.pure[,3])
edge[,4] <- seed.all.pure[,11]
edge[,5] <- seed.all.pure[,4]
edge[,6] <- c("Cosine")
edge[,7] <- seed.all.pure[,6]-seed.all.pure[,5]
edge[,8] <- c("-")
edge[,9] <- paste(seed.all.pure[,2]," (-) ",seed.all.pure[,3])
edge[,10] <- seed.all.pure[,1]
#edge[,11] <- seed.all.pure[,11]
colnames(edge)=c("node1","node2","name","cosine_score","EdgeAnnotation","EdgeType","mass_difference","shared interaction","shared name","round")
write.csv(edge,edge.name, row.names = FALSE)

node.name <- paste(file.name,".node.csv",sep="")
colnames(seed.1) <- colnames(seed.all.pure)
node.pre <- rbind(seed.1,seed.all.pure)
node <- matrix(ncol=6,nrow=nrow(node.pre))
node[,1] <- node.pre[,3]
node[,2] <- node.pre[,6]
node[,3] <- node.pre[,4]
node[,4] <- node.pre[,13]
node[,5] <- node.pre[,3]
node[,6] <- node.pre[,1]
colnames(node)=c("name","precursor mass","Compound_Name","Intensity","shared name","round")
node[!duplicated(node[,c("name","Compound_Name")])]
write.csv(node,node.name, row.names = FALSE)
}