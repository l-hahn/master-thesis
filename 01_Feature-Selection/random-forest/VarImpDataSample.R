MakeSynthData = function(XData){
  ColSample = function(Col){
    sample(Col,length(Col),replace = T)
  }
  SynthData = rbind(XData,apply(DF,2,ColSample))
  SynthData$Label= as.factor(c(rep(0,nrow(XData)),rep(1,nrow(XData))))
  return(SynthData)
}
DF = data.frame(
  A = round(rnorm(100, mean=25, sd=5)),
  B = round(rnorm(100, mean=50, sd=5)),
  C = round(rnorm(100, mean=75, sd=5)),
  D = round(runif(100,min=0,max=100)),
  E = round(runif(100,min=0,max=100))
)
#write.table(DF,"VarImp.dat",quote = F,row.names = F,col.names = F)


g <- randomForest(DF)
g$importance


DFSynth = MakeSynthData(DF)

gTest = randomForest(DFSynth[-ncol(DFSynth)],DFSynth$Label)
gTest$importance


#


g = randomForest(DFSynth[-ncol(DFSynth)],DFSynth$Label)


write.table(DFSynth,"VarImpSynth.dat",quote = F,row.names = F,col.names = F)
write.table(iris[-5],"IrisData",quote = F, row.names = F,col.names = F,sep=" ")
