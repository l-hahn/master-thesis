library(randomForest)
library(tictoc)


MakeSynthData = function(XData){
  ColSample = function(Col){
    sample(Col,length(Col),replace = T)
  }
  SynthData = rbind(XData,apply(DF,2,ColSample))
  SynthData$Label= as.factor(c(rep(0,nrow(XData)),rep(1,nrow(XData))))
  return(SynthData)
}





setwd("~/Schreibtisch/Master/")

data = read.table("data.ped")
PED = data[,- c(1:6)]
MAP = read.table("data.map")
colnames(PED) = MAP$V2
rownames(PED) = data$V2

tic()
g <- randomForest(PED[1:1000], keep.forest=FALSE, proximity=TRUE, ntree = 500)
toc()
g$importance


g <- randomForest(iris[,-5], ntree = 5000, importance = T)
g$importance
?randomForest
colSums(g$importance)

