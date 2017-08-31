#Pavana Anur
#CoCA for thymoma

#set path to working directory
setwd("~/Projects/coca_thymoma/fuzzy")

#Methylation, mean calculation for each cluster
#correct NA values

#original file name THYM_methclusterProbeBetas_20150916.txt
methdata=read.table("~/Projects/coca_thymoma/fuzzy/meth_transposed.txt",header=T,sep="\t")
str(methdata)
methclusters=read.table("~/Projects/coca_thymoma/meth_edited.txt",header=T,sep="\t")
str(methclusters)

methtmp=merge(methclusters,methdata,by.x="sample",by.y="sample",all=T)

methcluster1=subset(methtmp,methtmp$type==1)
str(methcluster1)
methmean1=colMeans(methcluster1[,c(-1,-2)],na.rm=T)
#methmean1=unname(methmean1)

methcluster2=subset(methtmp,methtmp$type==2)
str(methcluster2)
methmean2=colMeans(methcluster2[,c(-1,-2)],na.rm=T)
#methmean2=unname(methmean2)

methcluster3=subset(methtmp,methtmp$type==3)
str(methcluster1)
methmean3=colMeans(methcluster3[,c(-1,-2)],na.rm=T)
#methmean1=unname(methmean1)

methmeans=rbind(methmean1,methmean2,methmean3)

write.table(t(methmeans),"methmeans.txt",sep="\t",quote=F)

methdata2=read.table("~/Projects/coca_thymoma/fuzzy/THYM_methclusterProbeBetas_20150916.txt",header=T,sep="\t")
str(methdata2)

methmeans=read.table("~/Projects/coca_thymoma/fuzzy/methmeans.txt",header=T,sep="\t")
str(methmeans)

methmerge=merge(methmeans,methdata2)
methmerge=methmerge[,-1]
methmerge=na.rm(methmerge)

row.names(methmerge)=methmerge$probeID
str(methmerge)

x <- methmerge[1:3]
y <- methmerge[4:120]
methcor=cor(x, y)

write.table(methcor,"methcor.txt",sep="\t",quote=F)
##############################################################
#RPPA, mean calculation for each cluster
rppadata=read.table("~/Projects/coca_thymoma/fuzzy/rppa_transposed.txt",header=T,sep="\t")
str(rppadata)
rppaclusters=read.table("~/Projects/coca_thymoma/THYM_RPPA_clusters.txt",header=T,sep="\t")
str(rppaclusters)

rppatmp=merge(rppaclusters,rppadata)
str(rppatmp)

rppacluster1=subset(rppatmp,rppatmp$RPPAcluster=="C1")
str(rppacluster1)
rppamean1=colMeans(rppacluster1[,c(-1,-2)],na.rm=T)

rppacluster2=subset(rppatmp,rppatmp$RPPAcluster=="C2")
str(rppacluster2)
rppamean2=colMeans(rppacluster2[,c(-1,-2)],na.rm=T)

rppacluster3=subset(rppatmp,rppatmp$RPPAcluster=="C3")
str(rppacluster1)
rppamean3=colMeans(rppacluster3[,c(-1,-2)],na.rm=T)

rppameans=rbind(rppamean1,rppamean2,rppamean3)

write.table(t(rppameans),"rppameans.txt",sep="\t",quote=F)

rppadata2=read.table("~/Projects/coca_thymoma/fuzzy/THYM.protein_exp__mda_rppa_core__mdanderson_org__Level_3__protein_normalization__data.data.txt",header=T,sep="\t")
str(rppadata2)

rppameans=read.table("~/Projects/coca_thymoma/fuzzy/rppameans.txt",header=T,sep="\t")
str(rppameans)

rppamerge=merge(rppameans,rppadata2)

row.names(rppamerge)=rppamerge$protein
rppamerge=rppamerge[,-1]
str(rppamerge)

x <- rppamerge[1:3]
y <- rppamerge[4:93]
rppacor=cor(x, y)

write.table(rppacor,"rppacor.txt",sep="\t",quote=F)

##############################################################
#mirna, mean calculation for each cluster

mirna=read.table("~/Projects/coca_thymoma/THYM.miRNA.RPMabundance.mirbase20names.117tumors.20150908.txt",header=T,sep="\t")

topvar=read.table("~/Projects/coca_thymoma/fuzzy/rowlist_mean303rows.txt",header=T,sep="\t")

str(mirna)
str(topvar)

mirnavar=merge(mirna,topvar)

dim(mirnavar)

mirnavar=mirnavar[,-119]
str(mirnavar)
row.names(mirnavar)=mirnavar$Gene
mirnavar=mirnavar[,-1]

#centre and scale
newmirna=log10(mirnavar+1)
newmirna$mean=rowMeans(newmirna)
newmirna2=apply(newmirna,1,function(x){x[1:117]-x[118]})

write.table(newmirna2,"mirna_varying.txt",sep="\t",quote=F,row.names=T)

mirnadata=read.table("~/Projects/coca_thymoma/fuzzy/mirna_transposed.txt",header=T,sep="\t")
str(mirnadata)
mirnaclusters=read.table("~/Projects/coca_thymoma/THYM_miRNA_hiercluster_k4_117T_20150902.txt",header=T,sep="\t")
str(mirnaclusters)

mirnadata$id=(gsub"\\.","-",mirnadata$id)
mirnatmp=merge(mirnaclusters,mirnadata)
str(mirnatmp)

write.table(mirnatmp,"mirna_varying_tmp.txt",sep="\t",quote=F,row.names=F)

mirnacluster1=subset(mirnatmp,mirnatmp$cluster=="1")
str(mirnacluster1)
mirnamean1=colMeans(mirnacluster1[,c(-1,-2)],na.rm=T)

mirnacluster2=subset(mirnatmp,mirnatmp$cluster=="2")
str(mirnacluster2)
mirnamean2=colMeans(mirnacluster2[,c(-1,-2)],na.rm=T)

mirnacluster3=subset(mirnatmp,mirnatmp$cluster=="3")
str(mirnacluster3)
mirnamean3=colMeans(mirnacluster3[,c(-1,-2)],na.rm=T)

mirnacluster4=subset(mirnatmp,mirnatmp$cluster=="4")
str(mirnacluster4)
mirnamean4=colMeans(mirnacluster4[,c(-1,-2)],na.rm=T)

mirnacluster5=subset(mirnatmp,mirnatmp$cluster=="5")
str(mirnacluster5)
mirnamean5=colMeans(mirnacluster5[,c(-1,-2)],na.rm=T)

mirnameans=rbind(mirnamean1,mirnamean2,mirnamean3,mirnamean4,mirnamean5)

write.table(t(mirnameans),"mirnameans.txt",sep="\t",quote=F)

mirnadata2=read.table("~/Projects/coca_thymoma/fuzzy/mirna_varying.txt",header=T,sep="\t")
str(mirnadata2)

mirnameans=read.table("~/Projects/coca_thymoma/fuzzy/mirnameans.txt",header=T,sep="\t")
str(mirnameans)

mirnamerge=merge(mirnameans,mirnadata2)

row.names(mirnamerge)=mirnamerge$id
mirnamerge=mirnamerge[,-1]
str(mirnamerge)
dim(mirnamerge)
x <- mirnamerge[1:5]
y <- mirnamerge[6:122]

mirnacor=cor(x, y)

write.table(mirnacor,"mirnacor.txt",sep="\t",quote=F)
##############################################################
#CN, mean calculation for each clusters
cndata=read.table("~/Projects/coca_thymoma/fuzzy/cn_data_t.txt",header=T,sep="\t")
str(cndata)

#row.names(cndata)=cndata$chr
#cndata=cndata[,-1]
#cndata$mean=rowMeans(cndata)
#newcndata=apply(cndata,1,function(x){x[1:117]-x[118]})
#write.table(newcndata,"cn_meancentered.txt",sep="\t",quote=F,row.names=T)

cnclusters=read.table("~/Projects/coca_thymoma/data/THYM117_armclusters_ward.D_euclidean_k2-10.txt",header=T,sep="\t")
str(cnclusters)
cnclusters=cnclusters[,c(1,3)]
cntmp=merge(cnclusters,cndata)
str(cntmp)

cncluster1=subset(cntmp,cntmp$k3=="1")
str(cncluster1)
cnmean1=colMeans(cncluster1[,c(-1,-2)],na.rm=T)
#cnmean1=unname(cnmean1)

cncluster2=subset(cntmp,cntmp$k3=="2")
str(cncluster2)
cnmean2=colMeans(cncluster2[,c(-1,-2)],na.rm=T)
#cnmean2=unname(cnmean2)

cncluster3=subset(cntmp,cntmp$k3=="3")
str(cncluster3)
cnmean3=colMeans(cncluster3[,c(-1,-2)],na.rm=T)
#cnmean3=unname(cnmean3)

cnmeans=rbind(cnmean1,cnmean2,cnmean3)

write.table(t(cnmeans),"cnmeans.txt",sep="\t",quote=F)

cndata2=read.table("~/Projects/coca_thymoma/fuzzy/cn_data.txt",header=T,sep="\t")
str(cndata2)

cnmeans=read.table("~/Projects/coca_thymoma/fuzzy/cnmeans.txt",header=T,sep="\t")
str(cnmeans)

cnmerge=merge(cnmeans,cndata2)

row.names(cnmerge)=cnmerge$chr
cnmerge=cnmerge[,-1]
str(cnmerge)
dim(cnmerge)
x <- cnmerge[1:3]
y <- cnmerge[4:120]
cncor=cor(x, y)

write.table(cncor,"cncor.txt",sep="\t",quote=F)
##############################################################
#mRNA, mean calculation for each cluster
rna_top=read.table("~/Projects/coca_thymoma/fuzzy/top1000genes_Consensuscluster.txt",header=T)
rna_rsem=read.table("~/Projects/coca_thymoma/fuzzy/THYM.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",header=T,sep="\t")
rna_rsem=rna_rsem[-1,]
rna_rsem$id=gsub("\\|.*","",rna_rsem$Hybridization.REF)
rna_rsem=rna_rsem[,-1]

#remove 5 samples not in freeze here before other steps

str(rna_rsem)
rna_rsem<-rna_rsem[,!(colnames(rna_rsem) %in% c("TCGA.3G.AB0Q.01A.11R.A42C.07","TCGA.3G.AB0T.01A.11R.A42C.07","TCGA.X7.A8D7.11A.11R.A42C.07","TCGA.X7.A8D6.11A.22R.A42C.07","TCGA.XU.AAXV.01A.11R.A42C.07"))]

top_rsem=merge(rna_top,rna_rsem)
row.names(top_rsem)=top_rsem$id
top_rsem=top_rsem[,-1]

write.table(top_rsem,"top_rsem.txt",sep="\t",quote=F,row.names=T)

top_rsem=read.table("top_rsem.txt",header=T,sep="\t")
row.names(top_rsem)=top_rsem$id
top_rsem=top_rsem[,-1]
top_rsem=log10(top_rsem+1)
top_rsem$mean=rowMeans(top_rsem)
top_rsem=apply(top_rsem,1,function(x){x[1:117]-x[118]})
str(top_rsem)

write.table(top_rsem,"rna_varying.txt",sep="\t",quote=F,col.name=T)

rnadata=read.table("~/Projects/coca_thymoma/fuzzy/rna_varying.txt",header=T,sep="\t")
str(rnadata)
rnaclusters=read.table("~/Projects/coca_thymoma/fuzzy/mRNA_samples_with_clusters_modified.txt",header=T,sep="\t")
str(rnaclusters)

rnatmp=merge(rnaclusters,rnadata)
str(rnatmp)

rnacluster1=subset(rnatmp,rnatmp$rnacluster=="1")
str(rnacluster1)
rnamean1=colMeans(rnacluster1[,c(-1,-2)],na.rm=T)

rnacluster2=subset(rnatmp,rnatmp$rnacluster=="2")
str(rnacluster2)
rnamean2=colMeans(rnacluster2[,c(-1,-2)],na.rm=T)

rnacluster3=subset(rnatmp,rnatmp$rnacluster=="3")
str(rnacluster3)
rnamean3=colMeans(rnacluster3[,c(-1,-2)],na.rm=T)

rnacluster4=subset(rnatmp,rnatmp$rnacluster=="4")
str(rnacluster4)
rnamean4=colMeans(rnacluster4[,c(-1,-2)],na.rm=T)

rnameans=rbind(rnamean1,rnamean2,rnamean3,rnamean4)

write.table(t(rnameans),"rnameans.txt",sep="\t",quote=F)

write.table(t(top_rsem),"rna_rsem_data.txt",sep="\t",quote=F,row.name=T)

rnadata2=read.table("~/Projects/coca_thymoma/fuzzy/rna_rsem_data.txt",header=T,sep="\t")
dim(rnadata2)
head(rnadata2)

rnameans=read.table("~/Projects/coca_thymoma/fuzzy/rnameans.txt",header=T,sep="\t")
dim(rnameans)
#find and replace KRTAP5-10,KRTAP5-9,NKX3-2
rnameans$id=gsub("KRTAP5.10","KRTAP5-10",rnameans$id)
rnameans$id=gsub("KRTAP5.9","KRTAP5-9",rnameans$id)
rnameans$id=gsub("NKX3.2","NKX3-2",rnameans$id)

rnamerge=merge(rnameans,rnadata2)
dim(rnamerge)

row.names(rnamerge)=rnamerge$id
rnamerge=rnamerge[,-1]
str(rnamerge)
dim(rnamerge)

x <- rnamerge[1:4]
y <- rnamerge[5:121]
rnacor=cor(x, y)

write.table(rnacor,"rnacor.txt",sep="\t",quote=F)

##############################################################
#Concat all means from each platform
cn=read.table("~/Projects/coca_thymoma/fuzzy/cn_cor_transposed.txt",header=T,sep="\t")
meth=read.table("~/Projects/coca_thymoma/fuzzy/meth_cor_transposed.txt",header=T,sep="\t")
mir=read.table("~/Projects/coca_thymoma/fuzzy/mir_cor_transposed.txt",header=T,sep="\t")
rna=read.table("~/Projects/coca_thymoma/fuzzy/rna_cor_transposed.txt",header=T,sep="\t")
rppa=read.table("~/Projects/coca_thymoma/fuzzy/rppa_cor_transposed.txt",header=T,sep="\t")


my.list <- list(rna,cn,rppa,mir,meth)

thym_cor_merged=Reduce(function(x, y) merge(x, y,all.x=T,by.x="sample", by.y="sample"), 
my.list, accumulate=F)

write.table(t(thym_cor_merged),"thym_cor_merged.txt",sep="\t",quote=F,row.names=T,col.names=F)

##############################################################
#Input means file for clustering
library(ConsensusClusterPlus)

y=read.table("~/Projects/coca_thymoma/fuzzy/thym_cor_merged.txt",header=T,sep="\t",row.names=1)
str(y)
dim(y)
y[is.na(y)] <- 0

y.rcc<-ConsensusClusterPlus(as.matrix(y),maxK=10,reps=1000,pItem=0.9,pFeature=1,title="THYM_fuzzy",plot="pdf")

write.table(y.rcc[[3]]$consensusClass,"THYM.k3_v2.txt",sep="\t",col.names=NA)
write.table(y.rcc[[4]]$consensusClass,"THYM.k4_v2.txt",sep="\t",col.names=NA)
write.table(y.rcc[[5]]$consensusClass,"THYM.k4_v2.txt",sep="\t",col.names=NA)
write.table(y.rcc[[6]]$consensusClass,"THYM.k6_v2.txt",sep="\t",col.names=NA)
write.table(y.rcc[[7]]$consensusClass,"THYM.k7_v2.txt",sep="\t",col.names=NA)


o3<-colnames(y)[ y.rcc[[3]]$consensusTree$order ]
write.table(o3,"THYM.k3.order_v2.txt",sep="\t",col.names=NA)

o4<-colnames(y)[ y.rcc[[4]]$consensusTree$order ]
write.table(o4,"THYM.k4.order_v2.txt",sep="\t",col.names=NA)

o5<-colnames(y)[ y.rcc[[5]]$consensusTree$order ]
write.table(o5,"THYM.k4.order_v2.txt",sep="\t",col.names=NA)

o6<-colnames(y)[ y.rcc[[6]]$consensusTree$order ]
write.table(o6,"THYM.k6.order_v2.txt",sep="\t",col.names=NA)

o7<-colnames(y)[ y.rcc[[7]]$consensusTree$order ]
write.table(o7,"THYM.k7.order_v2.txt",sep="\t",col.names=NA)


##############################################################
#reorder based on k4 order, order from the file THYM.k4.order_v2.txt 
y=read.table("~/Projects/coca_thymoma/fuzzy/thym_cor_merged.txt",header=T,sep="\t",row.names=1)

yk4=y[,c("TCGA.3G.AB19","TCGA.XU.A92Z","TCGA.ZB.A96I","TCGA.X7.A8M3","TCGA.XU.A931","TCGA.5G.A9ZZ","TCGA.ZC.AAA7","TCGA.ZB.A96V","TCGA.ZB.A966","TCGA.XU.A936","TCGA.XU.A933","TCGA.5U.AB0D","TCGA.3S.A8YW","TCGA.4V.A9QN","TCGA.XM.A8R8","TCGA.ZC.AAAH","TCGA.ZC.AAAA","TCGA.ZB.A96R","TCGA.ZB.A96P","TCGA.ZB.A96O","TCGA.ZB.A96G","TCGA.ZB.A96F","TCGA.ZB.A963","TCGA.YT.A95G","TCGA.XU.AAXZ","TCGA.XU.AAXY","TCGA.XU.AAXX","TCGA.XU.AAXW","TCGA.XU.A932","TCGA.XU.A92Y","TCGA.XU.A92X","TCGA.XU.A92W","TCGA.XM.AAZ1","TCGA.XM.A8RL","TCGA.XM.A8RI","TCGA.XM.A8RD","TCGA.XM.A8RC","TCGA.XM.A8RB","TCGA.XM.A8R9","TCGA.X7.A8M7","TCGA.X7.A8M6","TCGA.X7.A8DG","TCGA.X7.A8DB","TCGA.X7.A8D7","TCGA.5U.AB0F","TCGA.5U.AB0E","TCGA.4X.A9FD","TCGA.4X.A9FB","TCGA.4X.A9F9","TCGA.4V.A9QT","TCGA.4V.A9QS","TCGA.4V.A9QR","TCGA.4V.A9QI","TCGA.3T.AA9L","TCGA.3S.AAYX","TCGA.3Q.A9WF","TCGA.3G.AB0O","TCGA.3G.AB14","TCGA.XU.A930","TCGA.ZB.A961","TCGA.5V.A9RR","TCGA.ZL.A9V6","TCGA.ZB.A96L","TCGA.ZB.A96H","TCGA.ZB.A96C","TCGA.ZB.A969","TCGA.ZB.A965","TCGA.ZB.A964","TCGA.ZB.A962","TCGA.YT.A95D","TCGA.XU.A92V","TCGA.XU.A92T","TCGA.XU.A92Q","TCGA.XU.A92O","TCGA.XM.AAZ3","TCGA.XM.AAZ2","TCGA.XM.A8RG","TCGA.XM.A8RF","TCGA.X7.A8M5","TCGA.4X.A9FA","TCGA.X7.A8D8","TCGA.X7.A8M0","TCGA.4V.A9QJ","TCGA.5K.AAAP","TCGA.4V.A9QL","TCGA.4V.A9QM","TCGA.ZT.A8OM","TCGA.ZC.AAAF","TCGA.ZB.A96Q","TCGA.ZB.A96M","TCGA.ZB.A96K","TCGA.ZB.A96E","TCGA.ZB.A96D","TCGA.ZB.A96B","TCGA.ZB.A96A","TCGA.YT.A95H","TCGA.YT.A95F","TCGA.YT.A95E","TCGA.XU.AAY1","TCGA.XU.AAY0","TCGA.XU.A92U","TCGA.XM.A8RH","TCGA.XM.A8RE","TCGA.XH.A853","TCGA.X7.A8M8","TCGA.X7.A8M4","TCGA.X7.A8M1","TCGA.X7.A8DJ","TCGA.X7.A8DF","TCGA.X7.A8DE","TCGA.X7.A8DD","TCGA.X7.A8D9","TCGA.X7.A8D6","TCGA.4X.A9FC","TCGA.4V.A9QX","TCGA.4V.A9QU","TCGA.4V.A9QW")]

write.table(yk4,"k4_ordered_clusters.txt",sep="\t",quote=F)

library("gplots")
library("devtools")
source('~/Desktop/heatmap.3.R', chdir = TRUE)

yk4=read.table("k4_ordered_clusters.txt",sep="\t",header=T)

rnames=yk4[,1]
yk4mat=data.matrix(yk4[,2:ncol(yk4)])
rownames(yk4mat)=rnames

annot=read.table("~/Projects/coca_thymoma/fuzzy/types_colors.txt",header=T,sep="\t",comment.char="")
str(annot)

hras=read.table("~/Projects/coca_thymoma/hras.txt",header=T,sep="\t")
tp53=read.table("~/Projects/coca_thymoma/tp53.txt",header=T,sep="\t")
gtf2i=read.table("~/Projects/coca_thymoma/gtf2i.txt",header=T,sep="\t")
nras=read.table("~/Projects/coca_thymoma/nras.txt",header=T,sep="\t")

my.list <- list(annot,gtf2i,hras,tp53,nras)

muts_merged=Reduce(function(x, y) merge(x, y,all.x=T,by.x="sample", by.y="sample"), 
                   my.list, accumulate=F)
annot=muts_merged
annot$HRAS=as.character(annot$HRAS)
annot$HRAS[is.na(annot$HRAS)]= ("gray")
annot$TP53=as.character(annot$TP53)
annot$TP53[is.na(annot$TP53)]= ("gray")
annot$GTF2I=as.character(annot$GTF2I)
annot$GTF2I[is.na(annot$GTF2I)]= ("gray")
annot$NRAS=as.character(annot$NRAS)
annot$NRAS[is.na(annot$NRAS)]= ("gray")

rownames(annot)=annot[,1]
annot=annot[,-1]
annotord=annot[c("TCGA-3G-AB19","TCGA-XU-A92Z","TCGA-ZB-A96I","TCGA-X7-A8M3","TCGA-XU-A931","TCGA-5G-A9ZZ","TCGA-ZC-AAA7","TCGA-ZB-A96V","TCGA-ZB-A966","TCGA-XU-A936","TCGA-XU-A933","TCGA-5U-AB0D","TCGA-3S-A8YW","TCGA-4V-A9QN","TCGA-XM-A8R8","TCGA-ZC-AAAH","TCGA-ZC-AAAA","TCGA-ZB-A96R","TCGA-ZB-A96P","TCGA-ZB-A96O","TCGA-ZB-A96G","TCGA-ZB-A96F","TCGA-ZB-A963","TCGA-YT-A95G","TCGA-XU-AAXZ","TCGA-XU-AAXY","TCGA-XU-AAXX","TCGA-XU-AAXW","TCGA-XU-A932","TCGA-XU-A92Y","TCGA-XU-A92X","TCGA-XU-A92W","TCGA-XM-AAZ1","TCGA-XM-A8RL","TCGA-XM-A8RI","TCGA-XM-A8RD","TCGA-XM-A8RC","TCGA-XM-A8RB","TCGA-XM-A8R9","TCGA-X7-A8M7","TCGA-X7-A8M6","TCGA-X7-A8DG","TCGA-X7-A8DB","TCGA-X7-A8D7","TCGA-5U-AB0F","TCGA-5U-AB0E","TCGA-4X-A9FD","TCGA-4X-A9FB","TCGA-4X-A9F9","TCGA-4V-A9QT","TCGA-4V-A9QS","TCGA-4V-A9QR","TCGA-4V-A9QI","TCGA-3T-AA9L","TCGA-3S-AAYX","TCGA-3Q-A9WF","TCGA-3G-AB0O","TCGA-3G-AB14","TCGA-XU-A930","TCGA-ZB-A961","TCGA-5V-A9RR","TCGA-ZL-A9V6","TCGA-ZB-A96L","TCGA-ZB-A96H","TCGA-ZB-A96C","TCGA-ZB-A969","TCGA-ZB-A965","TCGA-ZB-A964","TCGA-ZB-A962","TCGA-YT-A95D","TCGA-XU-A92V","TCGA-XU-A92T","TCGA-XU-A92Q","TCGA-XU-A92O","TCGA-XM-AAZ3","TCGA-XM-AAZ2","TCGA-XM-A8RG","TCGA-XM-A8RF","TCGA-X7-A8M5","TCGA-4X-A9FA","TCGA-X7-A8D8","TCGA-X7-A8M0","TCGA-4V-A9QJ","TCGA-5K-AAAP","TCGA-4V-A9QL","TCGA-4V-A9QM","TCGA-ZT-A8OM","TCGA-ZC-AAAF","TCGA-ZB-A96Q","TCGA-ZB-A96M","TCGA-ZB-A96K","TCGA-ZB-A96E","TCGA-ZB-A96D","TCGA-ZB-A96B","TCGA-ZB-A96A","TCGA-YT-A95H","TCGA-YT-A95F","TCGA-YT-A95E","TCGA-XU-AAY1","TCGA-XU-AAY0","TCGA-XU-A92U","TCGA-XM-A8RH","TCGA-XM-A8RE","TCGA-XH-A853","TCGA-X7-A8M8","TCGA-X7-A8M4","TCGA-X7-A8M1","TCGA-X7-A8DJ","TCGA-X7-A8DF","TCGA-X7-A8DE","TCGA-X7-A8DD","TCGA-X7-A8D9","TCGA-X7-A8D6","TCGA-4X-A9FC","TCGA-4V-A9QX","TCGA-4V-A9QU","TCGA-4V-A9QW"),]

yk4clab=annotord[,c(1,2,3,4,8,7,5,6)]

mydist=function(x) {as.dist(1-cor(na.omit(t(x)),method="pearson"))}
myclust=function(c) {hclust(c,method="average")}

library("RColorBrewer")

setEPS()
postscript("k4_cor_v4.eps")
heatmap.3(yk4mat,Rowv=T,Colv=F,scale="none",trace="none",ColSideColors=as.matrix(yk4clab),hclustfun=myclust,col=colorRampPalette(c('yellow', 'blue')),density.info="none",distfun=mydist,keysize=1,na.color="gray" )
dev.off()

###################################################################
Just for version 5, to correct MG and g2fti muts
yk4clab5=yk4clab[,c(2,7)]
###################################################################

setEPS()
postscript("k4_cor_v5.eps")
heatmap.3(yk4mat,Rowv=T,Colv=F,scale="none",trace="none",ColSideColors=as.matrix(yk4clab5),hclustfun=myclust,col=colorRampPalette(c('yellow', 'blue')),density.info="none",distfun=mydist,keysize=1,na.color="gray" )
dev.off()

plot.new()
setEPS()
postscript("thym_legend.eps")
plot.new()
legend("bottomleft",legend=c("A","AB","B","CA","MN-T","I","IIA","IIB","III","IVA","IVB","yes","no","1","2","3","4"),fill=c("#b2182b","#f4a582","#d6604d","#2166ac","#000000","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#000000","#2166ac","#b2182b","#2166ac","#d1e5f0","#92c5de","#4393c3","#2166ac"),border=FALSE, bty="n", y.intersp = 0.7, cex=0.7)
dev.off()

##################################################################

#heatmap binary input

y=read.table("~/Projects/coca_thymoma/thym_formatted_clusters.txt",header=T,sep="\t",row.names=1)

yk4=y[,c("TCGA.3G.AB19","TCGA.XU.A92Z","TCGA.ZB.A96I","TCGA.X7.A8M3","TCGA.XU.A931","TCGA.5G.A9ZZ","TCGA.ZC.AAA7","TCGA.ZB.A96V","TCGA.ZB.A966","TCGA.XU.A936","TCGA.XU.A933","TCGA.5U.AB0D","TCGA.3S.A8YW","TCGA.4V.A9QN","TCGA.XM.A8R8","TCGA.ZC.AAAH","TCGA.ZC.AAAA","TCGA.ZB.A96R","TCGA.ZB.A96P","TCGA.ZB.A96O","TCGA.ZB.A96G","TCGA.ZB.A96F","TCGA.ZB.A963","TCGA.YT.A95G","TCGA.XU.AAXZ","TCGA.XU.AAXY","TCGA.XU.AAXX","TCGA.XU.AAXW","TCGA.XU.A932","TCGA.XU.A92Y","TCGA.XU.A92X","TCGA.XU.A92W","TCGA.XM.AAZ1","TCGA.XM.A8RL","TCGA.XM.A8RI","TCGA.XM.A8RD","TCGA.XM.A8RC","TCGA.XM.A8RB","TCGA.XM.A8R9","TCGA.X7.A8M7","TCGA.X7.A8M6","TCGA.X7.A8DG","TCGA.X7.A8DB","TCGA.X7.A8D7","TCGA.5U.AB0F","TCGA.5U.AB0E","TCGA.4X.A9FD","TCGA.4X.A9FB","TCGA.4X.A9F9","TCGA.4V.A9QT","TCGA.4V.A9QS","TCGA.4V.A9QR","TCGA.4V.A9QI","TCGA.3T.AA9L","TCGA.3S.AAYX","TCGA.3Q.A9WF","TCGA.3G.AB0O","TCGA.3G.AB14","TCGA.XU.A930","TCGA.ZB.A961","TCGA.5V.A9RR","TCGA.ZL.A9V6","TCGA.ZB.A96L","TCGA.ZB.A96H","TCGA.ZB.A96C","TCGA.ZB.A969","TCGA.ZB.A965","TCGA.ZB.A964","TCGA.ZB.A962","TCGA.YT.A95D","TCGA.XU.A92V","TCGA.XU.A92T","TCGA.XU.A92Q","TCGA.XU.A92O","TCGA.XM.AAZ3","TCGA.XM.AAZ2","TCGA.XM.A8RG","TCGA.XM.A8RF","TCGA.X7.A8M5","TCGA.4X.A9FA","TCGA.X7.A8D8","TCGA.X7.A8M0","TCGA.4V.A9QJ","TCGA.5K.AAAP","TCGA.4V.A9QL","TCGA.4V.A9QM","TCGA.ZT.A8OM","TCGA.ZC.AAAF","TCGA.ZB.A96Q","TCGA.ZB.A96M","TCGA.ZB.A96K","TCGA.ZB.A96E","TCGA.ZB.A96D","TCGA.ZB.A96B","TCGA.ZB.A96A","TCGA.YT.A95H","TCGA.YT.A95F","TCGA.YT.A95E","TCGA.XU.AAY1","TCGA.XU.AAY0","TCGA.XU.A92U","TCGA.XM.A8RH","TCGA.XM.A8RE","TCGA.XH.A853","TCGA.X7.A8M8","TCGA.X7.A8M4","TCGA.X7.A8M1","TCGA.X7.A8DJ","TCGA.X7.A8DF","TCGA.X7.A8DE","TCGA.X7.A8DD","TCGA.X7.A8D9","TCGA.X7.A8D6","TCGA.4X.A9FC","TCGA.4V.A9QX","TCGA.4V.A9QU","TCGA.4V.A9QW")]

yk4=yk4[c(18,8,17,6,7,4,11,12,13,10,3,5,14,15,1,2,9,16),]

yk4mat=data.matrix(yk4[,1:ncol(yk4)])

setEPS()
postscript("k4_types_binary.eps")
heatmap.3(yk4mat,dendrogram="none",Rowv=FALSE,Colv=FALSE,scale="none",trace="none",col=c("white","blue"),density.info="none",na.color="gray")
dev.off()

##############################################################################
#Purity of clusters

cluster_assign=read.table("THYM.k4_v2.txt",header=T,sep="\t")
str(cluster_assign)
colnames(cluster_assign)=c("sample","cluster")
cluster_assign$sample=gsub("\\.","-",cluster_assign$sample)

surv=read.table("~/Projects/coca_thymoma/fuzzy/sample_surv.txt",header=T,sep="\t")
str(surv)

cluster_surv=merge(surv,cluster_assign)
str(cluster_surv)

boxplot(cluster_surv$purity~cluster_surv$cluster,col=c("lightblue","darkgreen","darkblue","lightgreen"),xlab="cluster",ylab="Purity")

############################################################################
#Survival analysis
library(survival)

cluster_surv$days=as.numeric(cluster_surv$days)
cluster_surv$status=as.numeric(cluster_surv$status)

fit1=survfit(Surv(cluster_surv$days,cluster_surv$status)~cluster_surv$cluster)
survdiff(Surv(cluster_surv$days,cluster_surv$status)~cluster_surv$cluster)

plot(fit1, xlab="Survival time in Days", 
     ylab="Survival Probability",lwd=3,  col=c("#92c5de","#2166ac","#f4a582","#b2182b"))
         
legend("bottomright", title="cluster(Total/Events)", c("1(46/3)", "2(14/4)","3(35/0)","4(21/1)","p=0.00859","n=116"),col=c("#92c5de","#2166ac","#f4a582","#b2182b"),lty=c(1,1,1,1),lwd=3)

#############################################################################
#gene mutation enrichment, so one side over exp alt=greater

mut=annot[,c(1,6,8)]
mut_cluster=merge(mut,cluster_assign,by="sample")
str(mut_cluster)

table(mut_cluster$cluster)

write.table(mut_cluster,"mut_cluster.txt",quote=F,sep="\t")

mut_cluster=mut_cluster[,-1]
tbl1=table(mut_cluster$cluster,mut_cluster$HRAS)
tbl1
colnames(tbl1)="HRAS"
tbl2=table(mut_cluster$cluster,mut_cluster$GTF2I )
tbl2
colnames(tbl2)="GTF2I "


a=c(46,14,36,21)
clust_hras=cbind(a,tbl1)
clust_gtf2i=cbind(a,tbl2)
fisher.test(clust_hras,alternative="greater")
fisher.test(clust_gtf2i,alternative="greater")

#####

	Fisher's Exact Test for Count Data

data:  clust_hras
p-value = 0.0003346
alternative hypothesis: greater

> fisher.test(clust_gtf2i,alternative="greater")

	Fisher's Exact Test for Count Data

data:  clust_gtf2i
p-value = 7.403e-07
alternative hypothesis: greater

######




