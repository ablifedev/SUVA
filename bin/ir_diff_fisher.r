#test
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile <- args[2]
outdir <- args[3]

setwd(outdir)
Data <- read.delim(file=infile,header=F)
Pvalue <-c()
chivalue <- c()
Ratiodiff <-c()
for (i in 1:nrow(Data)) {
	compare <- matrix(c(Data[i,5],Data[i,8],Data[i,10],Data[i,13]),nr=2)
	compare <- round(compare)
	result <- fisher.test(compare)		#,alternative = "greater"
	result2 <- chisq.test(compare,simulate.p.value=F,B=40000)
	rdiff <- Data[i,9] - Data[i,14]
	Pvalue <- c(Pvalue,result$p.value)
	chivalue <- c(chivalue,result2$p.value)
	Ratiodiff <- c(Ratiodiff,rdiff)
}
FDR <- p.adjust(chivalue,method = "BH")

Data <- cbind(Data,Ratiodiff)
Data <- cbind(Data,Pvalue)
Data <- cbind(Data,chivalue)
Data <- cbind(Data,FDR)
write.table(Data,file=outfile,row.names=F,append=F,quote=F,sep='\t')