
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile <- args[2]
outdir <- args[3]

exit <- function() {
  print("empty input file")
  .Internal(.invokeRestart(list(NULL, NULL), NULL))
}

if (file.info(infile)$size == 0){
	exit()
}

setwd(outdir)
Data <- read.delim(file=infile,header=F)
Pvalue <-c()
chivalue <- c()
Ratiodiff <-c()
for (i in 1:nrow(Data)) {
	compare <- matrix(c(Data[i,11],Data[i,12],Data[i,16],Data[i,17]),nr=2)
	compare <- round(compare)
	result <- fisher.test(compare)		#,alternative = "greater"
	result2 <- chisq.test(compare,simulate.p.value=F,B=40000)
	# rdiff <- Data[i,13] - Data[i,17]
	Pvalue <- c(Pvalue,result$p.value)
	chivalue <- c(chivalue,result2$p.value)
	# Ratiodiff <- c(Ratiodiff,rdiff)
}

FDR <- p.adjust(chivalue,method = "BH")

# Data <- cbind(Data,Ratiodiff)
Data <- cbind(Data,Pvalue)
Data <- cbind(Data,chivalue)
Data <- cbind(Data,FDR)
write.table(Data,file=outfile,row.names=F,append=F,quote=F,sep='\t')