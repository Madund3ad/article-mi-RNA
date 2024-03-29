#import#################################################
library(rtracklayer)
library(GenomicRanges)
library(GenomicFeatures)
library(plyr)
library(zoo)
library(plyr)
library("ggpubr")
library("ggplot2")
###################################################
ago11=read.table("AGO1  HepG2.bed", sep="\t")
ago12=read.table("AGO1 K562.bed", sep="\t")
ago2=read.table("AGO2 HepG2.bed", sep="\t")

colnames(ago11)=c("chr","start","end")
colnames(ago12)=c("chr","start","end")
colnames(ago2)=c("chr","start","end")

ago11=GenomicRanges::makeGRangesFromDataFrame(fixchrs(ago11))
ago12=GenomicRanges::makeGRangesFromDataFrame(fixchrs(ago12))
ago2=GenomicRanges::makeGRangesFromDataFrame(fixchrs(ago2))

gtf_file_hg <-"Homo_sapiens.GRCh37.87.gtf"
txdb <- makeTxDbFromGFF(gtf_file_hg, format="gtf")
gs=genes(txdb)
gs.tss=resize(gs, width=1, fix='start')
#tss=gs.tss
tss=transcripts(txdb)
tss=resize(tss, width=1, fix='start')

prot=prots[[1]]
makenplot<-function(prot,tss,n,name){
  dat=data.frame(matrix(data=0,nrow=2*n+1,ncol=2))
  dat[,1]=-n:n
  rownames(dat)=dat[,1]
  prec=precede(tss,prot)
  nas=which(is.na(prec))
  ds=distance(tss[-nas],prot[prec[-nas]])
  fol=follow(tss,prot)
  nasfol=which(is.na(fol))
  dsf=distance(tss[-nasfol],prot[fol[-nasfol]])
  ds=-ds
  distc=count(c(ds[ds>-n],dsf[dsf<n]))
  dat[as.character(distc$x),2]=distc$freq[na.omit(match(dat$X1,distc$x))]
  colnames(dat)=c("distance",name)
  return(dat)
}

###################################################
prots=list(ago11,ago12,ago2)
n=10000
dat2=data.frame(matrix(data=0,nrow=2*n+1,ncol=length(prots)+1))
dim(dat2)
for (i in 1:length(prots)){
  dat2[,i+1]=makenplot(prots[[i]],gs.tss,n,"AGO1 HepG2")[,2]
}
dat2[,1]=-n:n
dat2=dat2[-(n+1),]
#norm###################################
dat2[,2]=dat2[,2]/length(ago11)
dat2[,3]=dat2[,3]/length(ago12)
dat2[,4]=dat2[,4]/length(ago2)
    
dat2[,2:4]=dat2[,2:4]*10000
#smooth######################################

dat2[,2]=rollmean(dat2[,2], 7, na.pad=TRUE)
dat2[,3]=rollmean(dat2[,3], 7, na.pad=TRUE)
dat2[,4]=rollmean(dat2[,4], 7, na.pad=TRUE)


colnames(dat2)=c("distance","AGO1_HepG2","AGO1_K562","AGO2_HepG2")
#plot#######################################################################33
p1=ggplot(dat2, aes(x=distance)) + 
  geom_line(aes(y = AGO1_HepG2), color = "darkred")+labs(x = "")+ylim(0,25)
p2=ggplot(dat2, aes(x=distance)) + 
  geom_line(aes(y = AGO1_K562), color="steelblue")+labs(x = "")+ylim(0,25)
p3=ggplot(dat2, aes(x=distance)) + 
  geom_line(aes(y = AGO2_HepG2), color="black")+ylim(0,25)#+labs(x = "Distance to TSS")
#####################################################################
ggarrange(p1, p2, p3 + rremove("x.text"), 
          labels = c("A", "B", "C"),
          common.legend=T,
          #label.x=c("","","Distance to TSS"),
          ncol = 1, nrow = 3)
ggarrange(p1, p3 + rremove("x.text"), 
          labels = c("A", "B"),
          common.legend=T,
          #label.x=c("","","Distance to TSS"),
          ncol = 1, nrow = 2)
