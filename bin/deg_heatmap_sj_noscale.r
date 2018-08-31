#!/usr/bin/env Rscript
####################################################################################
### Copyright (C) 2015-2019 by ABLIFE
####################################################################################
# 名称：HeatMap_vs_sample.r
# 描述：对数据进行log2处理，再做heatmaptu
# 作者：Joseph Wei
# 创建时间：2015-5-20
# 联系方式：yaxunwei@ablife.cc
####################################################################################
### 修改记录
####################################################################################
# Date           Version       Author            ChangeLog
# 2015-5-21      v1.0          Weiyaxun         修改测试版本
#
#####################################################################################

#####################################################################################
#####参数获取
#####################################################################################
#
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stats"))

option_list <- list(
  make_option(c("-f", "--file"),action = "store",type = "character",
              help = "The Input file"),
  make_option(c("-a", "--annofile"),action = "store",type = "character",
              help = "The Anno file"),
  make_option(c("-s","--sample"),action = "store",type = "character",
              help = "The sample file"),
  make_option(c("-b","--start"),action = "store",type = "integer",default = 1,
              help = "appoint the column of start ; default = 1"),
  make_option(c("-e","--end"),action = "store",type = "integer",default = -1,
              help = "appoint the column of end; default = -1"),
  make_option(c("-n","--filename"),action = "store",type = "character",default="heatmap_KCNQ1OT1_negative.pdf",
              help = "The filename of the picture ; default = deg_heatmap.pdf"),
  make_option(c("-o", "--outdir"),action = "store",type = "character",default = "./",
              help = "The outdir;default = ./")
)
opt <- parse_args(OptionParser(option_list = option_list))
start <- as.numeric(opt$start)
end <- as.numeric(opt$end)

# library(cluster)
# library(Biobase)
# library(qvalue)
library("corrplot")
NO_REUSE = F

# opt$file = "PTBP1_negative_target_RPKM.xls"
# opt$file = "XLOC_004865_KCNQ1OT1_mRNA_negative_RPKM_symbol.txt"
# opt$annofile = "anno.txt"
# opt$file = "deg_genes_rpkm_of_all_samples_in_deg.txt"
# opt$file = "data_mRNA.txt"
# opt$file = "data_mRNA.txt"

# # get the filename to use later
# filename <- strsplit(opt$file,"/")[[1]]
# filename <- filename[length(filename)]
# filename <- sub('.txt','',filename)

# # try to reuse earlier-loaded data if possible

# print('Reading matrix file.')
primary_data = read.table(opt$file, header=T, com='', sep="\t", row.names=1, check.names=F)
anno = read.table(opt$annofile, header=T, com='', sep="\t", row.names=1, check.names=F)
# primary_data = read.table("Sample_correlation.dat", header=T, com='', sep="\t", row.names=1, check.names=F)


if(end> 0){
  primary_data <- primary_data[,start:end]
}
primary_data = as.matrix(primary_data)
primary_data <- primary_data[rowSums(primary_data) != 0,]
primary_data = primary_data*100
# primary_data = log2(primary_data+1)   ##10.18修改程序加上log
data = primary_data


# data <- data[7200:7219,]
# data <- subset(data,rowSums(data)<15)
# tail(data)

# sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')
# # write.table(sample_cor, file=paste(opt$filename,".dat",sep=''), quote=F, sep='  ')
# # sample_dist = dist(t(data), method='euclidean')
# # hc_samples = hclust(sample_dist, method='complete')
# col <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
#                           "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
#                           "#4393C3", "#2166AC", "#053061","#67001F", "#B2182B", "#D6604D", 
#                           "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
#                           "#4393C3", "#2166AC", "#053061","#67001F", "#B2182B", "#D6604D", 
#                           "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
#                           "#4393C3", "#2166AC", "#053061","#67001F", "#B2182B", "#D6604D", 
#                           "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
#                           "#4393C3", "#2166AC", "#053061"))(200)
# min <- min(sample_cor)
# corrplot(sample_cor, type="upper", order="hclust", hclust.method="complete",col=col,tl.srt=45,cl.lim=c(min,1))

# library(ALL) #可以使用biocLite("ALL")安装该数据包
# data("ALL")
# library(limma)
# eset<-ALL[,ALL$mol.biol %in% c("BCR/ABL","ALL1/AF4")]
# f<-factor(as.character(eset$mol.biol))
# design<-model.matrix(~f)
# fit<-eBayes(lmFit(eset,design)) #对基因芯片数据进行分析，得到差异表达的数据
# selected  <- p.adjust(fit$p.value[, 2]) <0.001 
# esetSel <- eset[selected,] #选择其中一部分绘制热图
# data<-exprs(esetSel)

library(pheatmap)
library("RColorBrewer")
# pheatmap(data,fontsize=9, fontsize_row=6) #最简单地直接出图
# pheatmap(data, scale = "row", clustering_distance_row = "correlation", fontsize=9, fontsize_row=6) #改变排序算法
# pheatmap(data, color = colorRampPalette(c("navy", "white", "firebrick3"))(50), fontsize=9, fontsize_row=6) #自定义颜色
# pheatmap(data, cluster_row=FALSE, fontsize=9, fontsize_row=6) #关闭按行排序
# pheatmap(data, legend = FALSE, fontsize=9, fontsize_row=6) #关闭图例
# data[lower.tri(data)]=NA



# triangle heatmap
# o = rownames(data)
# sample_dist = dist(t(data), method='euclidean')
# hc = hclust(sample_dist, method='complete')
# # hc = hclust(as.dist(1 - data))
# data = data[hc$order, hc$order]
# data[lower.tri(data)] = NA
# data = data[o, o]
# data
# 
# pheatmap(data, cluster_col = hc, cluster_row = hc,border_color="white",show_colnames=F,cellwidth = 25, cellheight = 25, fontsize=14, width=8,height=8,filename = opt$filename,display_numbers = FALSE)




# correlation
# pheatmap(data, clustering_method = "complete", clustering_distance_rows="euclidean",clustering_distance_cols="euclidean",scale="none",cellwidth = 12, cellheight = 12, fontsize=12, show_rownames = T, filename = "test_cor_chend.pdf") #设定格子的尺寸

# gene heatmap
# pheatmap(data, clustering_method = "complete", clustering_distance_rows="euclidean",clustering_distance_cols="euclidean",scale="row",cellwidth = 35, cellheight = 5, fontsize=12, fontsize_row=12,show_rownames = T, filename = "test_deg.pdf") #设定格子的尺寸

#碱蓬
# pheatmap(data, clustering_method = "complete", clustering_distance_rows="euclidean",cluster_cols=F,scale="row",cellwidth = 35, cellheight = 8, fontsize=12, fontsize_row=8,show_rownames = T, border_color="white",filename = opt$filename,main="Heatmap") #设定格子的尺寸


#原花青素
# data <- ifelse(data>3,3,data)
# data <- ifelse(data < -3,-3,data)

ColorVar        <- c("#1b9e77","#d95f02","#7570b3","#66a61e","#e7298a","#e6ab02","#a6761d","#666666","#92C5DE")

annotype <- unique(anno[,1])
len <- length(annotype)
ColorVar1 <- ColorVar[1:len]
names(ColorVar1) <- annotype
# print(ColorVar1)

annotype <- unique(anno[,2])
len <- length(annotype)
ColorVar2 <- ColorVar[1:len]
names(ColorVar2) <- annotype

annotype <- unique(anno[,3])
annotype <- factor(annotype,levels=c("normal","tissue"))
len <- length(annotype)
len <- len+3
ColorVar3 <- ColorVar[4:len]
names(ColorVar3) <- annotype



anno_colors <- list(cancer = ColorVar1, stage = ColorVar2, type = ColorVar3)


# scaleyellowred <- colorRampPalette(c("#54B447","#54B447","#54B447","#54B447","#36B04A","#0C753B","black","#891619","#C71E24","#E32325","#E32325","#E32325","#E32325"),space = "rgb")(500)


# mihou2
# scaleyellowred <- colorRampPalette(c("#54B447","#36B04A","#0C753B","black","#891619","#C71E24","#E32325"),space = "rgb")(500)
scaleyellowred <- colorRampPalette(c("#08519c","#3182bd","#ffffff","#e6550d","#a63603"),space = "rgb")(500)

# pheatmap(data,color=scaleyellowred, clustering_method = "complete", cluster_rows = F,clustering_distance_cols="euclidean",scale="row",fontsize=9, fontsize_row=8,cellwidth = 7,show_rownames = F, border_color=NA,filename = opt$filename,main="Heatmap",width=8,height=10, annotation = anno, annotation_colors = anno_colors) #猕猴，y不cluster

# pheatmap(data,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",,clustering_distance_cols="euclidean",scale="row",fontsize=10, fontsize_row=10,cellwidth = 9,show_rownames = T, border_color=NA,filename = opt$filename,main="Heatmap",width=12,height=10, annotation = anno, annotation_colors = anno_colors) #白血病

# pheatmap(data,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",,clustering_distance_cols="correlation",scale="row",fontsize=12, fontsize_row=12,cellwidth = 3,show_rownames = F, show_colnames = F, border_color=NA,filename = paste(opt$filename,".png",sep=''),main="",width=8,height=6, annotation = anno, annotation_colors = anno_colors) #白血病



colcount <- ncol(data)

if (colcount<=15){

pheatmap(data,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",clustering_distance_cols="correlation",scale="none",fontsize=12, fontsize_row=12,cellwidth = 10,show_rownames = F, show_colnames = F, border_color=NA,filename = paste(opt$filename,".png",sep=''),main="",width=8,height=6, annotation = anno, annotation_colors = anno_colors) #白血病
}

if (colcount>15){

  pheatmap(data,color=scaleyellowred, clustering_method = "complete", clustering_distance_rows="euclidean",clustering_distance_cols="correlation",scale="none",fontsize=12, fontsize_row=12,cellwidth = 3,show_rownames = F, show_colnames = F, border_color=NA,filename = paste(opt$filename,".png",sep=''),main="",width=8,height=6, annotation = anno, annotation_colors = anno_colors) #白血病
}







# cluster_cols=F,

# ##小梁网
# # data = log2(data+1)
# data[,c(1,2,3,4,5)] <- data[,c(1,2,3,4,5)]/data[,1]
# #deg
# col.pal <- brewer.pal(9,"Blues")
# # colorRampPalette is in the RColorBrewer package.  This creates a colour palette that shades from light yellow to red in RGB space with 100 unique colours
# scaleyellowred <- colorRampPalette(c("lightyellow","pink", "red"),space = "rgb")(500)
# pheatmap(data, color=scaleyellowred,clustering_method = "complete", cluster_cols=F,clustering_distance_rows="euclidean",clustering_distance_cols="euclidean",scale="row",cellwidth = 35,  fontsize=12, fontsize_row=8,show_rownames = T, border_color=NA,filename = "xiaoliangwang4.pdf",main="Heatmap") #设定格子的尺寸


# color.map <- function(mol.biol) { if (mol.biol=="ALL1/AF4") 1 else 2 }
# patientcolors <- unlist(lapply(esetSel$mol.bio, color.map))
# hc<-hclust(dist(t(data)))
# dd.col<-as.dendrogram(hc)
# groups <- cutree(hc,k=7)
# annotation<-data.frame(Var1=factor(patientcolors,labels=c("class1","class2")),Var2=groups)
# pheatmap(data, annotation=annotation, fontsize=9, fontsize_row=6) #为样品分组
# Var1 = c("navy", "skyblue")
# Var2 = c("snow", "steelblue")
# names(Var1) = c("class1", "class2")
# ann_colors = list(Var1 = Var1, Var2 = Var2)
# pheatmap(data, annotation=annotation, annotation_colors = ann_colors, fontsize=9, fontsize_row=6, main="Heatmap Example") #为分组的样品设定颜色

