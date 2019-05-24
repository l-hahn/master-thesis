library(randomForest)
library(ggplot2)
library(gridExtra)
library(hashmap)
library(stringi)

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
  
  Encoding = 0:(length(Genotypes)-1)+0
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

setwd("~/Schreibtisch/Master/master-data/tuberculosis-data/")
data = read.table("data_5000-snps.ped")
MAP = read.table("data_5000-snps.map")
#data = read.table("simdata.ped")
#MAP = read.table("simdata.map")
PED = as.data.frame(apply(as.data.frame(t(apply(data[-c(1:6)],1,mypaste)), stringsAsFactors = F),2,AlphabetEncoderSym))
PEDBin = as.data.frame(apply(as.data.frame(t(apply(data[-c(1:6)],1,mypaste)), stringsAsFactors = F),2,SnpBinEncoder))
colnames(PED) = MAP$V2
colnames(PEDBin) = MAP$V2


SNP_RF_Unsup <- randomForest(PED, ntree = 500)
SnpBin_RF_Unsup <- randomForest(PEDBin, ntree = 500)
R_RF_VarImp = data.frame(ID=rownames(SNP_RF_Unsup$importance),Imp=SNP_RF_Unsup$importance, stringsAsFactors = F)[order(SNP_RF_Unsup$importance, decreasing = T),]
R_RF_Bin_VarImp = data.frame(ID=rownames(SnpBin_RF_Unsup$importance),Imp=SnpBin_RF_Unsup$importance, stringsAsFactors = F)[order(SnpBin_RF_Unsup$importance, decreasing = T),]
rownames(R_RF_VarImp) = 1:nrow(R_RF_VarImp)
rownames(R_RF_Bin_VarImp) = 1:nrow(R_RF_Bin_VarImp)

Py_RF_VarImp = read.table("tuberculosis_sample1000.attr", stringsAsFactors = F)
Py_RF_Bin_VarImp = read.table("tuberculosis_sample1000_bin.attr", stringsAsFactors = F)
#Py_RF_VarImp = read.table("simdata.attr", stringsAsFactors = F)
#Py_RF_Bin_VarImp = read.table("simdata_bin.attr", stringsAsFactors = F)
colnames(Py_RF_VarImp) = c("ID","MeanDecreaseGini")
colnames(Py_RF_Bin_VarImp) = c("ID","MeanDecreaseGini")



R_RF_VarImp$MeanDecreaseGini = mmnorm(R_RF_VarImp$MeanDecreaseGini)
R_RF_VarImp$Important = 1:nrow(R_RF_VarImp) %in% grep(R_RF_VarImp$ID,pattern = "important")
R_RF_Bin_VarImp$MeanDecreaseGini = mmnorm(R_RF_Bin_VarImp$MeanDecreaseGini)
R_RF_Bin_VarImp$Important = 1:nrow(R_RF_Bin_VarImp) %in% grep(R_RF_Bin_VarImp$ID,pattern = "important")
Py_RF_VarImp$MeanDecreaseGini = mmnorm(Py_RF_VarImp$MeanDecreaseGini)
Py_RF_VarImp$Important = 1:nrow(Py_RF_VarImp) %in% grep(Py_RF_VarImp$ID,pattern = "important")
Py_RF_Bin_VarImp$MeanDecreaseGini = mmnorm(Py_RF_Bin_VarImp$MeanDecreaseGini)
Py_RF_Bin_VarImp$Important = 1:nrow(Py_RF_Bin_VarImp) %in% grep(Py_RF_Bin_VarImp$ID,pattern = "important")

R_Pos = grep(R_RF_VarImp$ID,pattern = "important")
R_Bin_Pos = grep(R_RF_Bin_VarImp$ID,pattern = "important")
Py_Pos = grep(Py_RF_VarImp$ID,pattern = "important")
Py_Bin_Pos = grep(Py_RF_Bin_VarImp$ID,pattern = "important")
X_Pos = c(#Py_RF_VarImp$MeanDecreaseGini[Py_Pos],
          #Py_RF_Bin_VarImp$MeanDecreaseGini[Py_Bin_Pos],
          R_RF_VarImp$MeanDecreaseGini[R_Pos],
          R_RF_Bin_VarImp$MeanDecreaseGini[R_Bin_Pos]
          )
ColorScheme = rainbow(4)
Col = rep(ColorScheme,each=length(R_Pos))

PlotData = data.frame(
    ID = c(R_RF_VarImp$ID,
           R_RF_Bin_VarImp$ID
           #Py_RF_VarImp$ID,
           #Py_RF_Bin_VarImp$ID
           ), 
    VarImp = c(R_RF_VarImp$MeanDecreaseGini,
               R_RF_Bin_VarImp$MeanDecreaseGini
               #Py_RF_VarImp$MeanDecreaseGini,
               #Py_RF_Bin_VarImp$MeanDecreaseGini
               ), 
    Type=c(rep("R",nrow(R_RF_VarImp)),
           rep("R-Bin",nrow(R_RF_Bin_VarImp))
           #rep("Python", nrow(Py_RF_VarImp)),
           #rep("Python-Bin", nrow(Py_RF_Bin_VarImp))
           )
  )

ggplot(data=PlotData, aes(x=VarImp, y=..density..))+geom_density(aes(fill = Type), alpha=0.5)+geom_vline(xintercept = X_Pos,colour=Col)+
  scale_fill_manual( values = ColorScheme)


# PyP =ggplot(data=NULL, aes(x=reorder(ID,-MeanDecreaseGini),y=MeanDecreaseGini)) + 
#   geom_bar(stat="identity",aes(fill=Py_RF_VarImp$Important,col=Py_RF_VarImp$Important), data=Py_RF_VarImp)+
#   theme(axis.text.x=element_blank())+
#   scale_fill_discrete(name="Snp type", labels=c("Not Associated", "Associated"))+guides(col=FALSE)+
#   labs(title="SNP Variable Importance via RF\nPython alphabet encoded",x="SNPs sorted by Variable importance", y="Normalised mean decrease Gini")
# 
# PyBP=ggplot(data=NULL, aes(x=reorder(ID,-MeanDecreaseGini),y=MeanDecreaseGini))+ 
#   geom_bar(stat="identity",aes(fill=Py_RF_Bin_VarImp$Important,col=Py_RF_Bin_VarImp$Important), data=Py_RF_Bin_VarImp)+
#   theme(axis.text.x=element_blank())+
#   scale_fill_discrete(name="Snp type", labels=c("Not Associated", "Associated"))+guides(col=FALSE)+
#   labs(title="SNP Variable Importance via RF\nPython additive encoded",x="SNPs sorted by Variable importance", y="Normalised mean decrease Gini")

RP =ggplot(data=NULL, aes(x=reorder(ID,-MeanDecreaseGini),y=MeanDecreaseGini)) + 
  geom_bar(stat="identity",aes(fill=R_RF_VarImp$Important,col=R_RF_VarImp$Important), data=R_RF_VarImp)+
  theme(axis.text.x=element_blank())+
  scale_fill_discrete(name="Snp type", labels=c("Not Associated", "Associated"))+guides(col=FALSE)+
  labs(title="SNP Variable Importance via RF\nR alphabet encoded",x="SNPs sorted by Variable importance", y="Normalised mean decrease Gini")

RBP =ggplot(data=NULL, aes(x=reorder(ID,-MeanDecreaseGini),y=MeanDecreaseGini))+ 
  geom_bar(stat="identity",aes(fill=R_RF_Bin_VarImp$Important,col=R_RF_Bin_VarImp$Important), data=R_RF_Bin_VarImp)+
  theme(axis.text.x=element_blank())+
  scale_fill_discrete(name="Snp type", labels=c("Not Associated", "Associated"))+guides(col=FALSE)+
  labs(title="SNP Variable Importance via RF\nR additive encoded",x="SNPs sorted by Variable importance", y="Normalised mean decrease Gini")

grid.arrange(#PyP,
             #PyBP,
             RP,
             RBP)
