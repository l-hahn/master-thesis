library(ggplot2)
library(gridExtra)
setwd("~/Schreibtisch/Master/master-thesis/01_Feature-Selection/random-forest/")
Data = read.table("VarImpResults.dat")
Data
colnames(Data) = c("ID","VarImp","Conserv","MAF","Entropy")
Data

ggplot(data=Data,aes(x=VarImp)) + geom_density(aes(fill=ID),alpha=0.5)
plot(Data$VarImp,Data$Entropy,cex=0.4,pch=19)

Data2 = read.table("data_out.attr",header=T)
SortData2 = Data2[order(Data2$VarImp,decreasing = T),]
#Subdata2 = SortData2[1:10000,]
#Subdata2

plot(SortData2$VarImp,SortData2$Conserv,cex=0.4,pch=19,xlab="Variable Importance",ylab="Conservation",main="Simulationdata, Conservation 0.5 to 1 by 0.005")
plot(SortData2$VarImp,SortData2$MAF,cex=0.4,pch=19,xlab="Variable Importance",ylab="MAF",main="Simulationdata, Conservation 0.5 to 1 by 0.005")
plot(SortData2$VarImp,SortData2$Entropy,cex=0.4,pch=19,xlab="Variable Importance",ylab="Entropy",main="Simulationdata, Conservation 0.5 to 1 by 0.005")

plot(y=SortData2$Conserv,x=1:nrow(Data2))
plot(y=SortData2$MAF,x=1:nrow(Data2))
plot(y=SortData2$Entropy,x=1:nrow(Data2))







setwd("~/Schreibtisch/Master/master-thesis/01_Feature-Selection/random-forest/")
DataV = read.table("TopVarImpResults.dat", head=T)
DataV

DataE = read.table("TopEntroResults.dat", head=T)
DataE

DataB = read.table("TopBothLowResults.dat", head=T)
DataB

Data = rbind(DataV,DataE,DataB)
Data$Type = c(rep("HighVarImp",100),rep("HighEntroLowVarImp",100),rep("BothLow",100))
A=ggplot(data=Data,aes(x=Importance)) +
  geom_histogram(aes(fill=Type),alpha=1) +
  scale_fill_discrete(name = "Type", labels = c("BothLow", "HighE-LowV", "HighV")) +
  aes(y=3*stat(count)/sum(stat(count))) +
  scale_y_continuous(labels = scales::percent) +
  labs(y="Percent")
B=ggplot(data=Data,aes(x=MAF)) +
  geom_histogram(aes(fill=Type),alpha=1) +
  scale_fill_discrete(name = "Type", labels = c("BothLow", "HighE-LowV", "HighV")) +
  aes(y=3*stat(count)/sum(stat(count))) +
  scale_y_continuous(labels = scales::percent) +
  labs(y="Percent")
C=ggplot(data=Data,aes(x=Conserv)) +
  geom_histogram(aes(fill=Type),alpha=1) +
  scale_fill_discrete(name = "Type", labels = c("BothLow", "HighE-LowV", "HighV")) +
  aes(y=3*stat(count)/sum(stat(count))) +
  scale_y_continuous(labels = scales::percent) +
  labs(y="Percent")
D=ggplot(data=Data,aes(x=Entro)) +
  geom_histogram(aes(fill=Type),alpha=1) +
  scale_fill_discrete(name = "Type", labels = c("BothLow", "HighE-LowV", "HighV")) +
  aes(y=3*stat(count)/sum(stat(count))) +
  scale_y_continuous(labels = scales::percent) +
  labs(y="Percent")
P=grid.arrange(A,B,C,D)

ggsave("Top100ByType.pdf",P,device = "pdf")

par(mfrow=c(2,2))
DataE2 = read.table("TopEntroResults2.dat",head=T)
plot(density(DataE2$Rank),main="For Maximum Entropy, H = 1.0, 100 Runs",xlab="Rank")
plot(density(DataE2$Importance),main="For Maximum Entropy, H = 1.0, 100 Runs",xlab="Variable Importance")

DataV2 = read.table("TopVarImpResults2.dat",head=T)
plot(density(DataV2$Entro),main="For Maximum Variable Importance, 100 Runs",xlab="Entropy")
plot(density(DataV2$Importance),main="For Maximum Variable Importance, 100 Runs",xlab="Variable Importance")






tapply(DataE2$Rank,DataE2$ID,mean)
tapply(DataE2$Rank,DataE2$ID,median)
tapply(DataE2$Rank,DataE2$ID,sd)

tapply(DataE2$Importance,DataE2$ID,mean)
tapply(DataE2$Importance,DataE2$ID,median)



















######## PCA 

library(ggplot2)
library(gridExtra)
library(hashmap)
library(stringi)
library(FactoMineR)
library(factoextra)
library(randomForest)

MakeSynthData = function(XData){
  ColSample = function(Col){
    sample(Col,length(Col),replace = T)
  }
  SynthData = rbind(XData,apply(DF,2,ColSample))
  SynthData$Label= as.factor(c(rep(0,nrow(XData)),rep(1,nrow(XData))))
  return(SynthData)
}

mypaste = function(Data){
  i = seq(from=1, to=length(Data), by = 2L)
  return(as.character(paste0(Data[i],Data[i+1])))
}

mmnorm = function(Data){
  return((Data - min(Data)) / (max(Data) - min(Data)))
}

SnpBinEncoder = function(Data){
  AlleleSorted = rownames(sort(table(unlist(strsplit(Data,""))), decreasing = T))
  GenotypComb = outer(AlleleSorted,AlleleSorted,"paste0")
  Genotypes = GenotypComb[upper.tri(GenotypComb,diag = T)]
  GenotypesSym = GenotypComb[lower.tri(GenotypComb,diag = F)]
  
  Encoding = 0:(length(Genotypes)-1)+1
  Encoder = hashmap(Genotypes, Encoding)
  Encoding = c(Encoding, Encoder[[stri_reverse(GenotypesSym)]])
  Genotypes = c(Genotypes,GenotypesSym)[order(Encoding)]
  Encoding = Encoding[order(Encoding)]
  Encoder = hashmap(Genotypes, Encoding)
  return(Encoder[[Data]])
}
AlphabetEncoder = function(Data){
  Alleles = c("A","C","G","T","0")
  Genotypes = as.vector(t(outer(Alleles,Alleles,"paste0")))
  Encoding = 1:length(Genotypes)
  Encoder = hashmap(Genotypes, Encoding)
  return(Encoder[[Data]])
}
AlphabetEncoderSym = function(Data){
  Alleles = c("A","C","G","T","0")
  GenotypesComb = t(outer(Alleles,Alleles,"paste0"))
  Genotypes = as.vector(GenotypesComb[lower.tri(t(GenotypesComb),diag = T)])
  Encoding = 1:length(Genotypes)
  Encoder = hashmap(Genotypes, Encoding)
  GenotypesSym = as.vector(GenotypesComb[upper.tri(t(GenotypesComb),diag = F)])
  Encoding = c(Encoding, Encoder[[stri_reverse(GenotypesSym)]])
  Genotypes = c(Genotypes,GenotypesSym)[order(Encoding)]
  Encoding = Encoding[order(Encoding)]
  Encoder = hashmap(Genotypes, Encoding)
  return(Encoder[[Data]])
}


data = read.table("data_out_fine.ped")
MAP = read.table("data_out_fine.map")
PED = as.data.frame(apply(as.data.frame(t(apply(data[-c(1:6)],1,mypaste)), stringsAsFactors = F),2,SnpBinEncoder))
colnames(PED) = MAP$V2


res.pca = PCA(PED, graph = TRUE, scale=T)
summary(res.pca)
fviz_screeplot(res.pca ) 
test = fviz_contrib(res.pca, choice = "var", axes = 1:2)
PcaData = test$data

RfData = read.table("data_out_fine.attr",head=T)
colnames(PcaData) = c("ID","PcaVarImp")

Data = merge(RfData,PcaData,by="ID")
plot(Data$VarImp,Data$PcaVarImp)


R_RF_VarImp_Mean = R_RF_VarImp
for(i in 1:10){
PEDSample = PED[sample(1:ncol(PED))]
SNP_RF_Unsup = randomForest(PEDSample, ntree = 500, mtry=1000 , importance = TRUE)
R_RF_VarImp = data.frame(ID=rownames(SNP_RF_Unsup$importance),Imp=as.vector(SNP_RF_Unsup$importance[,4]), stringsAsFactors = F)[order(as.vector(SNP_RF_Unsup$importance[,4]), decreasing = T),]
R_RF_VarImp[order(R_RF_VarImp$ID),]
colnames(R_RF_VarImp) = c("ID","R_RF_VarImp")
R_RF_VarImp_Mean = merge(R_RF_VarImp_Mean,R_RF_VarImp,by="ID")
}

#R_RF_VarImp
#Data2 = merge(R_RF_VarImp,PcaData,by="ID")
#Data3 = merge(RfData,Data2,by="ID")
#plot(Data3$R_RF_VarImp,Data3$VarImp,cex=0.4,pch=19,xlab="R-RF_VarImp",ylab="Py-RF_VarImp")
#plot(Data3$VarImp,Data3$PcaVarImp,cex=0.4,pch=19,xlab="Py-RF_VarImp",ylab="R-PCA_VarImp")
#plot(Data3$R_RF_VarImp,Data3$PcaVarImp,cex=0.4,pch=19,xlab="R-RF_VarImp",ylab="R-PCA_VarImp")
#cor(Data3$R_RF_VarImp,Data3$VarImp,method="spearman")



IrisUnsup = randomForest(iris[1:4], ntree = 500, importance = TRUE)
R_IrisRF_UnSup_VarImp = data.frame(ID=rownames(IrisUnsup$importance),Imp=as.vector(IrisUnsup$importance[,4]), stringsAsFactors = F)[order(as.vector(IrisUnsup$importance[,4]), decreasing = T),]
colnames(R_IrisRF_UnSup_VarImp) = c("ID","R_RF_ImpUnsup")
R_IrisRF_UnSup_VarImp

IrisSup = randomForest(Species~.,iris, ntree = 500, importance = TRUE)
R_IrisRF_Sup_VarImp = data.frame(ID=rownames(IrisSup$importance),Imp=as.vector(IrisSup$importance[,4]), stringsAsFactors = F)[order(as.vector(IrisSup$importance[,4]), decreasing = T),]
colnames(R_IrisRF_Sup_VarImp) = c("ID","R_RF_ImpSup")
R_IrisRF_Sup_VarImp

res.pca = PCA(iris[1:4], graph = TRUE, scale=T)
#fviz_screeplot(res.pca ) 
test = fviz_contrib(res.pca, choice = "var", axes = 1:2)
PcaDataIris = test$data
colnames(PcaDataIris) = c("ID","R_PCA_Imp")
PcaDataIris

IrisData = merge(R_IrisRF_Sup_VarImp,R_IrisRF_UnSup_VarImp,by="ID")
IrisData = merge(IrisData,PcaDataIris)
IrisData


par(mfrow=c(1,2))
plot(IrisData$R_RF_ImpSup,IrisData$R_PCA_Imp,main="Supervised RF vs PCA in R",xlab="R RF Supervised",ylab="PCA",pch=19,col="red")
plot(IrisData$R_RF_ImpUnsup,IrisData$R_PCA_Imp,main="Unupervised RF vs PCA in R",xlab="R RF Unsupervised",ylab="PCA",pch=19,col="red")














#### Entro Rank of top by sorting
Top = 4
Tmp = Data3[order(Data3$Entropy,decreasing = T),]
Tmp$Idx = 1:nrow(Tmp)
Tmp$Top = c(rep(T,Top),rep(F,nrow(Tmp)-Top))

Entro1SortedPyVarImp = which(Tmp[order(Tmp$VarImp, decreasing = T),]$Top)
Entro1SortedRVarImp = which(Tmp[order(Tmp$R_RF_VarImp, decreasing = T),]$Top)
Entro1SortedPCAVarImp = which(Tmp[order(Tmp$PcaVarImp, decreasing = T),]$Top)

Entro1SortedPCAVarImp
Entro1SortedRVarImp
Entro1SortedPyVarImp







Test = read.table("data_out_fine.attr",head=T)
Test = Test[order(Test$ID),] ##viel Bäume
Test$Group = floor(0:(nrow(Test)-1)/4)+1


Test2 = read.table("data_out_fine.attr",head=T)
Test2 = Test[order(Test2$ID),]
Test2$Group = floor(0:(nrow(Test2)-1)/4)+1


plot(density(tapply(Test$VarImp,Data$Group,var)))
lines(density(tapply(Test2$VarImp,Data$Group,var)),col="red")








Data = read.table("simdata_gt50to100.dat")


### Python3 RF top
Tmp = Data[order(Data$Py_RF_Importance,decreasing = T),]
summary(Tmp[1:50,])

### R RF top
Tmp = Data[order(Data$R_RF_Importance,decreasing = T),]
summary(Tmp[1:50,])

### R PCA top
Tmp = Data[order(Data$R_PCA_Importance,decreasing = T),]
summary(Tmp[1:50,])







































#################### COMPARE TOP 10% PCA, RF, PLINK ########################
setwd("~/Schreibtisch/Master/master-thesis/01_Feature-Selection/random-forest/")
data = read.table("data_out_fine.ped")
MAP = read.table("data_out_fine.map")
PED = as.data.frame(apply(as.data.frame(t(apply(data[-c(1:6)],1,mypaste)), stringsAsFactors = F),2,SnpBinEncoder))
colnames(PED) = MAP$V2


res.pca = PCA(PED, graph = TRUE, scale=T)
summary(res.pca)
fviz_screeplot(res.pca) 
Contrib = fviz_contrib(res.pca, choice = "var", axes = c(1,2,8,9,10))
?fviz_contrib
PcaData = Contrib$data
OrderPcaData = PcaData[order(PcaData$contrib, decreasing = T),]
colnames(OrderPcaData) = c("ID","PcaImportance")
Pca_Important = OrderPcaData[1:50,]


PlinkData = read.table("plink.assoc",head=T)
Plink_Important = PlinkData[order(PlinkData$P,decreasing = F),c("SNP","P")]
colnames(Plink_Important) = c("ID","PlinkImportance")
Plink_Important = Plink_Important[1:50,]

PyRfData = read.table("data_out_fine.attr",head=T)
PyRf_Important = PyRfData[order(PyRfData$VarImp,decreasing = T),]
PyRf_Important = PyRf_Important[1:50,]



TopLists = data.frame(PcaTop = Pca_Important$ID, PlinkTop = Plink_Important$ID, PyRfTop = PyRf_Important$ID, stringsAsFactors = F)
# 
# idx = 1
# NewTmp = c()
# for(Item in strsplit(as.vector(TopLists$PcaTop),"_")){
#   Test = as.character(paste(Item[-3],collapse="_"))
#   NewTmp[idx] = Test
#   idx = idx + 1
# }
# TopLists$PcaTop = NewTmp
# 
# idx = 1
# NewTmp = c()
# for(Item in strsplit(as.vector(TopLists$PlinkTop),"_")){
#   Test = as.character(paste(Item[-3],collapse="_"))
#   NewTmp[idx] = Test
#   idx = idx + 1
# }
# TopLists$PlinkTop = NewTmp
# 
# 
# idx = 1
# NewTmp = c()
# for(Item in strsplit(as.vector(TopLists$PyRfTop),"_")){
#   Test = as.character(paste(Item[-3],collapse="_"))
#   NewTmp[idx] = Test
#   idx = idx + 1
# }
# TopLists$PyRfTop = NewTmp









RfPlink = c()
idx = 1
for(Elem in TopLists$PyRfTop){
  if( Elem %in% TopLists$PlinkTop ){
    RfPlink[idx] = TRUE
  }
  else{
    RfPlink[idx] = FALSE
  }
  idx = idx +1
}
RfPlink
TopLists$PyRfTop[!RfPlink]
TopLists$PyRfTop[RfPlink]
sum(RfPlink)

PlinkRF = c()
idx = 1
for(Elem in TopLists$PlinkTop){
  if( Elem %in% TopLists$PyRfTop ){
    PlinkRF[idx] = TRUE
  }
  else{
    PlinkRF[idx] = FALSE
  }
  idx = idx +1
}
PlinkRF
TopLists$PlinkTop[!PlinkRF]
TopLists$PlinkTop[PlinkRF]
which(PlinkRF)
sum(PlinkRF)


RfPca = c()
idx = 1
for(Elem in TopLists$PyRfTop){
  if( Elem %in% TopLists$PcaTop ){
    RfPca[idx] = TRUE
  }
  else{
    RfPca[idx] = FALSE
  }
  idx = idx +1
}
RfPca
TopLists$PyRfTop[!RfPca]
TopLists$PyRfTop[RfPca]
sum(RfPca)


PcaRf = c()
idx = 1
for(Elem in TopLists$PcaTop){
  if( Elem %in% TopLists$PyRfTop ){
    PcaRf[idx] = TRUE
  }
  else{
    PcaRf[idx] = FALSE
  }
  idx = idx +1
}
PcaRf
TopLists$PcaTop[!PcaRf]
TopLists$PcaTop[PcaRf]
sum(PcaRf)



PcaPlink = c()
idx = 1
for(Elem in TopLists$PcaTop){
  if( Elem %in% TopLists$PlinkTop ){
    PcaPlink[idx] = TRUE
  }
  else{
    PcaPlink[idx] = FALSE
  }
  idx = idx +1
}
PcaPlink
sum(PcaPlink)






pdf("Top10PercPlinkRfOverlap.pdf",paper="special")
data=read.table("simdata-1000_data2_out_2gt_OverlapRfPlink.dat",header = T)
par(mfrow=c(3,1))
plot(data$MTry,data$OverlapSupRfPlink,type="b",pch=20,col="blue",xlab="MTry",ylab="Top 10% Overlap with PLINK",main="SimData with 2 Genotypes (CC,GG)")
legend("bottomright",legend="RF-Supervised",fill="blue")
plot(data$MTry,data$OverlapUnsupRfPlink,type="b",pch=20,col="red",xlab="MTry",ylab="Top 10% Overlap with PLINK",main="SimData with 2 Genotypes (CC,GG)")
legend("bottomright",legend="RF-Unsupervised",fill="red")
plot(data$MTry,data$OverlapSupRfPlink,type="b",pch=20,ylim=c(0,1),col="blue",xlab="MTry",ylab="Top 10% Overlap with PLINK",main="SimData with 2 Genotypes (CC,GG)")
points(data$MTry,data$OverlapUnsupRfPlink,type="b",pch=20,col="red")
legend("bottomright",legend=c("RF-Supervised","RF-Unsupervised"),fill=c("blue","red"))


data=read.table("simdata-1000_data_out_fine_3gt_OverlapRfPlink.dat",header = T)
plot(data$MTry,data$OverlapSupRfPlink,type="b",pch=20,col="blue",xlab="MTry",ylab="Top 10% Overlap with PLINK",main="SimData with 3 Genotypes (CC,CG,GG)")
legend("bottomright",legend="RF-Supervised",fill="blue")
plot(data$MTry,data$OverlapUnsupRfPlink,type="b",pch=20,col="red",xlab="MTry",ylab="Top 10% Overlap with PLINK",main="SimData with 3 Genotypes (CC,CG,GG)")
legend("bottomright",legend="RF-Unsupervised",fill="red")
plot(data$MTry,data$OverlapSupRfPlink,type="b",pch=20,ylim=c(0,1),col="blue",xlab="MTry",ylab="Top 10% Overlap with PLINK",main="SimData with 3 Genotypes (CC,CG,GG)")
points(data$MTry,data$OverlapUnsupRfPlink,type="b",pch=20,col="red")
legend("toplef",legend=c("RF-Supervised","RF-Unsupervised"),fill=c("blue","red"))


data=read.table("bovtubdata_1000-randsamples_withimportantsnps_OverlapRfPlink.dat",header = T)
plot(data$MTry,data$OverlapSupRfPlink,type="b",pch=20,col="blue",xlab="MTry",ylab="Top 10% Overlap with PLINK",main="1000 rand. samples from bovData with important SNPs")
legend("bottomright",legend="RF-Supervised",fill="blue")
plot(data$MTry,data$OverlapUnsupRfPlink,type="b",pch=20,col="red",xlab="MTry",ylab="Top 10% Overlap with PLINK",main="1000 rand. samples from bovData with important SNPs")
legend("bottomright",legend="RF-Unsupervised",fill="red")
plot(data$MTry,data$OverlapSupRfPlink,type="b",pch=20,ylim=c(0,1),col="blue",xlab="MTry",ylab="Top 10% Overlap with PLINK",main="1000 rand. samples from bovData with important SNPs")
points(data$MTry,data$OverlapUnsupRfPlink,type="b",pch=20,col="red")
legend("toplef",legend=c("RF-Supervised","RF-Unsupervised"),fill=c("blue","red"))
dev.off()























setwd("~/Workbench/Python/")
Data = read.table("BreimanOrigShuffleData.dat",head=T)
Data$Mean = apply(Data[c(-1,-12,-13,-14)],1,median)
plot(Data[c(12,13,14,15)])


library(ggplot2)
library(gridExtra)
setwd("~/Schreibtisch/Master/master-thesis/01_Feature-Selection/random-forest/")
Data = read.table("OverlapRfPlink.dat",head=T)
Data2 = read.table("OverlapRfPlink2.dat",head=T)
Data2$PhenoPercent = rep("Breiman",nrow(Data2))

AllData = rbind(Data,Data2)
AllData$PhenoPercent = as.factor(AllData$PhenoPercent)

SubData1 = subset(Data,PhenoPercent < 0.6)
SubData1$PhenoPercent = as.factor(SubData1$PhenoPercent)
A=ggplot(data=SubData1, aes(x=MTry, y=OverlapUnsupRfPlink, group=PhenoPercent)) +
  geom_line(aes(linetype=PhenoPercent,col=PhenoPercent))+
  geom_point(aes(shape=PhenoPercent,col=PhenoPercent))+
  ylim(min(AllData$OverlapUnsupRfPlink), max(AllData$OverlapUnsupRfPlink))


SubData2 = subset(Data,0.6 <= PhenoPercent & PhenoPercent < 0.7)
SubData2$PhenoPercent = as.factor(SubData2$PhenoPercent)
B=ggplot(data=SubData2, aes(x=MTry, y=OverlapUnsupRfPlink, group=PhenoPercent)) +
  geom_line(aes(linetype=PhenoPercent,col=PhenoPercent))+
  geom_point(aes(shape=PhenoPercent,col=PhenoPercent))+
  ylim(min(AllData$OverlapUnsupRfPlink), max(AllData$OverlapUnsupRfPlink))

SubData3 = subset(Data,0.7 <= PhenoPercent & PhenoPercent < 0.8)
SubData3$PhenoPercent = as.factor(SubData3$PhenoPercent)
C=ggplot(data=SubData3, aes(x=MTry, y=OverlapUnsupRfPlink, group=PhenoPercent)) +
  geom_line(aes(linetype=PhenoPercent,col=PhenoPercent))+
  geom_point(aes(shape=PhenoPercent,col=PhenoPercent))+
  ylim(min(AllData$OverlapUnsupRfPlink), max(AllData$OverlapUnsupRfPlink))

SubData4 = subset(Data,0.8 <= PhenoPercent & PhenoPercent < 0.9)
SubData4$PhenoPercent = as.factor(SubData4$PhenoPercent)
D=ggplot(data=SubData4, aes(x=MTry, y=OverlapUnsupRfPlink, group=PhenoPercent)) +
  geom_line(aes(linetype=PhenoPercent,col=PhenoPercent))+
  geom_point(aes(shape=PhenoPercent,col=PhenoPercent))+
  ylim(min(AllData$OverlapUnsupRfPlink), max(AllData$OverlapUnsupRfPlink))


SubData5 = subset(Data,0.9 <= PhenoPercent & PhenoPercent < 1)
SubData5$PhenoPercent = as.factor(SubData5$PhenoPercent)
E=ggplot(data=SubData5, aes(x=MTry, y=OverlapUnsupRfPlink, group=PhenoPercent)) +
  geom_line(aes(linetype=PhenoPercent,col=PhenoPercent))+
  geom_point(aes(shape=PhenoPercent,col=PhenoPercent))+
  ylim(min(AllData$OverlapUnsupRfPlink), max(AllData$OverlapUnsupRfPlink))

SubData6 = Data2
F=ggplot(data=SubData6, aes(x=MTry, y=OverlapUnsupRfPlink, group=PhenoPercent)) +
  geom_line(aes(linetype=PhenoPercent,col=PhenoPercent))+
  geom_point(aes(shape=PhenoPercent,col=PhenoPercent))+
  ylim(min(AllData$OverlapUnsupRfPlink), max(AllData$OverlapUnsupRfPlink))


P=grid.arrange(A,B,C,D,E,F)
ggsave("Phenoswitch",plot = P, device = "pdf")


TabMean = tapply(AllData$OverlapUnsupRfPlink,AllData$PhenoPercent,mean)
TabMedian = tapply(AllData$OverlapUnsupRfPlink,AllData$PhenoPercent,median)
yValuesMean = c(TabMean)
yValuesMedian = c(TabMedian)
xValues=names(TabMean)

plot(1:length(TabMean),yValuesMean,type="b",col="blue",pch=20,xaxt="n")
lines(1:length(TabMean),yValuesMedian,type="b",col="red",pch=20)
axis(1, at=1:length(TabMean), labels=xValues)
legend("bottomleft",legend = c("mean","median"),fill=c("blue","red"))







setwd("~/Schreibtisch/Bioinformatik/")

data = read.table("Test.dat")


CountTable = table(data$Allele1)
Prob = CountTable/sum(CountTable)
Prob[2]





fun1=min
fun2=max

setwd("~/Schreibtisch/Master/master-thesis/01_Feature-Selection/random-forest/")
Data = read.table("OverlapRfPlink_AllFeaturesCompare.dat.dat",head=T)
Data$PhenoPercent = as.factor(Data$PhenoPercent)
TabMean = tapply(Data$OverlapUnsupRfPlink,Data$PhenoPercent,fun1)
TabMedian = tapply(Data$OverlapUnsupRfPlink,Data$PhenoPercent,fun2)

yValuesMean=c(TabMean,fun1(Data$Breiman))
yValuesMedian=c(TabMedian,fun2(Data$Breiman))
xValues=c(names(TabMean),"Breiman")

Min = min(yValuesMean,yValuesMedian)
Max = max(yValuesMean,yValuesMedian)

plot(1:length(yValuesMean),yValuesMean,ylim=c(Min,Max),type="b",col="blue",pch=20,xaxt="n")
lines(1:length(yValuesMean),yValuesMedian,ylim=c(Min,Max),type="b",col="red",pch=20)
axis(1, at=1:length(xValues), labels=xValues)
legend("bottomleft",legend = c("min","max"),fill=c("blue","red"))


















setwd("~/Workbench/Python/")
Data=read.table("BreimanOrigShuffleData_RelativeMin.dat",head=T)
Data$Min = apply(Data[!names(Data) %in% c("ID","Conserv","MAF","entropy")],1,min)
Data$Mean = apply(Data[!names(Data) %in% c("ID","Conserv","MAF","entropy","Min")],1,mean)
Data$Median = apply(Data[!names(Data) %in% c("ID","Conserv","MAF","entropy","Mean","Min")],1,median)
Data$SD = apply(Data[!names(Data) %in% c("ID","Conserv","MAF","entropy","Mean","Median","Min")],1,sd)
Data$ZScore = (Data$Min-mean(Data$Min))/sd(Data$Min)


hist(Data[Data$ZScore>1,]$Conserv)


mypaste = function(D){
  return(paste0(D,sep="", collapse=""))
}




setwd("~/Schreibtisch/Machine_Learning/")
data = read.table("test.ped", stringsAsFactors = F)
data = data[7:10]
SNPs = data.frame(S1=apply(data[1:2],1,mypaste),S2=apply(data[3:4],1,mypaste))
CT=table(SNPs)

H1CT = table(SNPs[1])
H1Prob = H1CT / sum(H1CT) 
H1 = -sum(H1Prob * log2(H1Prob))
H1

H2CT = table(SNPs[2])
H2Prob = H2CT / sum(H2CT) 
H2 = -sum(H2Prob * log2(H2Prob))
H2

H12CT = table(SNPs)
H12Prob = H12CT / sum(H12CT)
H12 = -sum(H12Prob*log2(H12Prob), na.rm=T)
H12

H1+H2-H12
