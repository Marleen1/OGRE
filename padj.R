

#dat=read.table("distances.tab",sep="\t",head=F) #too big to read
dat=read.table("p.txt",sep="\t",head=F)
head(dat)
#q=p.adjust(dat[,4],method="BH")
q=p.adjust(dat[,1],method="BH")
head(q)
#dat2=cbind(dat,q)
#write.table(dat2,"distances.padj.tab",quote=F,sep="\t",row.names=F,col.names=F)
write.table(q,"padj.txt",quote=F,sep="\t",row.names=F,col.names=F)
