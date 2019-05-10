library(gaston)
#Sample data creation and write to Bed
data(AGT)
x <- as.bed.matrix(AGT.gen, AGT.fam, AGT.bim)
ld.x = LD(x, c(1,ncol(x)))
write.bed.matrix(x, "LCT")
#Read JSD Matrix
jsdMat = read.table("test.ped.jsd_LD", header = TRUE, sep = " ")
jsdMatrix = as.matrix.data.frame(jsdMat[,-1])
rownames(jsdMatrix) = colnames(jsdMatrix)
#Visualization
LD.plot(ld.x[1:20,1:20])
LD.plot(jsdMatrix[1:20,1:20])
#Export to PDF
LD.plot(ld.x[1:20,1:20], pdf.file = "Normal-LD.pdf", finalize.pdf = TRUE)
LD.plot(jsdMatrix[1:20,1:20], pdf.file = "JSD-LD.pdf", finalize.pdf = TRUE)
LD.plot(ld.x, pdf.file = "Normal-LD.pdf", finalize.pdf = TRUE)
LD.plot(jsdMatrix, pdf.file = "JSD-LD.pdf", finalize.pdf = TRUE)