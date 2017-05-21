options(digits=15)
pair <- read.delim("/home/qzo/Dropbox (ORNL)/project/virus_tree/tree/kmer15cor.txt", header=F,stringsAsFactors=F)
order = pair[order(pair[,1],pair[,2],decreasing=F),]
head(order)

names <-as.character(union(order[,1],order[,2]))
same<-cbind(names,names,rep(0,3905))
colnames(same)<-colnames(order)
kmer15<-rbind(order,same)
order15 = kmer15[order(kmer15[,1],kmer15[,2],decreasing=F),]

rm(pair, order, same)

library(reshape2)

m<-acast(order15,V1~V2,value.var="V3")
library("phytools")


for (i in 1:3905){
  for (j in 1:3905){
    if (is.na(m[i,j])){
      m[i,j]=m[j,i]
    }
  }
}



vir15 <- nj(as.dist(m))


write.tree(vir15,"/home/qzo/Dropbox (ORNL)/project/virus_tree/tree/kmer15.newick")



bacteria <- read.table("~/Dropbox (ORNL)/project/virus_tree/host/bacteria.csv", quote="\"")

index<-rep(NA,1375)

for (i in 1:1375){
  
  index[i]<-which(colnames(m)==bacteria[i,1])
  
} 

n<-m[index,index]
bact<-nj(as.dist(n))

write.tree(bact,"/home/qzo/Dropbox (ORNL)/project/virus_tree/bact_tree/bact_kmer15.newick")



prok_col <- read.delim("~/Dropbox (ORNL)/project/virus_tree/new_host/prok_list.txt", header=T)


index<-rep(NA,1438)

for (i in 1:1438){
  
  index[i]<-which(colnames(m)==prok_col[i,1])
  
} 

pro<-m[index,index]
protree<-nj(as.dist(pro))


write.tree(protree,"/home/qzo/Dropbox (ORNL)/project/virus_tree/new_host/prok_kmer15.newick")



ds_col <- read.delim("~/Dropbox (ORNL)/project/virus_tree/dsDNA_list.txt", header=F)


index<-rep(NA,1826)

for (i in 1:1826){
  
  index[i]<-which(colnames(m)==ds_col[i,1])
  
} 

ds<-m[index,index]
dstree<-nj(as.dist(ds))


write.tree(dstree,"/home/qzo/Dropbox (ORNL)/project/virus_tree/dstree/dskmer15.newick")

