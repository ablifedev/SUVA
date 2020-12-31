
rm(list=ls())
args <- commandArgs(trailingOnly = TRUE)	#get the arguments from the command line
if (length(args) < 5) {
	stop('Your input arguments is wrong!\n
		args1:\t co_as file\n
		args2:\t test col: 1,2,3\n
		args3:\t ctrl col: 4,5,6\n
		args4:\t out file\n
		args5:\t out directory\n'
	)
} else {
	infile <- args[1]
	test <- args[2]
	ctrl <- args[3]
	outfile <- args[4]
	outdir <- args[5]
}

test <- strsplit(test,",")[[1]]
ctrl <- strsplit(ctrl,",")[[1]]

# print(as.numeric(test)[2])

setwd(outdir)
Data <- read.delim(file=infile,header=T,stringsAsFactors=FALSE)
Pvalue <-c()
Tvalue <-c()
m <-c()
n <-c()

# print(!is.na(as.numeric(Data[2,13])))

for (i in 1:nrow(Data)) {
  x <-c()
  y <-c()
  xnum <- 0
  ynum <- 0
  for (j in as.numeric(test)) { 
    if(!is.na(as.numeric(Data[i,j]))){
      x <- c(x,as.numeric(Data[i,j]))
      xnum <- xnum+1
    }  
  }
  for (j in as.numeric(ctrl)) { 
    if(!is.na(as.numeric(Data[i,j]))){
      y <- c(y,as.numeric(Data[i,j]))
      ynum <- ynum+1
    }   
  }
#  m<-c(m,c(x))
#  n<-c(n,c(y))
  
  xres<-summary(x)
  yres<-summary(y)
  
  
  if(!is.null(x) && !is.null(y) && xnum>1 && ynum>1 && (xres["Min."]!=xres["Max."] || yres["Min."]!=yres["Max."])){
  # if(!is.null(x) && !is.null(y) && xnum>1 && ynum>1 ){
    result <- t.test(x,y,var.equal = TRUE)
    Pvalue <- c(Pvalue,result$p.value)
    Tvalue <- c(Tvalue,result$statistic)
  }else{
    Pvalue <- c(Pvalue,1)
    Tvalue <- c(Tvalue,0)
  }
}

FDR <- p.adjust(Pvalue,method = "BH")

Data <- cbind(Data,Pvalue)
Data <- cbind(Data,Tvalue)
Data <- cbind(Data,FDR)
head(Data)
write.table(Data,file=outfile,row.names=F,append=F,quote=F,sep='\t')