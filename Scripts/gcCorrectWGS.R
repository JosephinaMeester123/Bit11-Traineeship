#READ DATA

library(gam)
library(ggplot2)
library(mgcv)
# ARGUMENT :
#  1 : input file (counts / bins)
#  2 : output file : normalized counts/bins
#  3 : output file : chromosome based counts (?)
#  4 : output file : basename for median / ff
#  5 : seqFF	   : seqFF value
#  6 : refdir      : reference dir for X & Y correction files

args<-commandArgs(TRUE)
#binned.raw.counts.csv <- args[1]
binned.csv <- args[2]
chr.csv <- args[3]
#basename <- args[4]
basename <- "C:/Users/jeann/Desktop/CBS/"
#seqFF <- as.numeric(args[5])
refdir <- args[6]
SAMPLE <- 108995

gc_content <- read.table("C:/Users/jeann/Desktop/CBS/bins500bp_gc.bed",as.is=T,header=F,sep="\t")
colnames(gc_content) <- c('CHR','BIN.START','BIN.END', 'BIN.GC.CONTENT', 'BIN.N.COUNT')
raw_count <- read.table("C:/Users/jeann/Desktop/CBS/108995.count", as.is=T,header=F,sep="\t")
colnames(raw_count) <- c('CHR','BIN.START','BIN.END', 'COUNT')
binned.raw2 <- merge(gc_content,raw_count,by=c('CHR','BIN.START','BIN.END'))
binned.raw <- cbind(SAMPLE, binned.raw)

#gc_content <- read.table("/home/jmeester/Internship/bins500bp_gc.bed",as.is=T,header=F,sep="\t")
#colnames(gc_content) <- c('CHR','BIN.START','BIN.END', 'BIN.GC.CONTENT', 'BIN.N.COUNT')
#raw_count <- read.table("/home/jmeester/Internship/countfiles/108995.count", as.is=T,header=F,sep=" ")
#colnames(raw_count) <- c('CHR','BIN.START','BIN.END', 'count', 'total_sample', 'total_pop', 'COUNT')
#binned.raw <- merge(gc_content,raw_count,by=c('CHR','BIN.START','BIN.END'))

#summary(binned.raw)


#X_corr_csv <- read.csv(paste(refdir,"/","X_correction.csv",sep=""),header=T)
#Y_corr_csv <- read.csv(paste(refdir,"/","Y_correction.csv",sep=""),header=T)


#X_corr <- merge(x=binned.raw[,c(2:4,7)],y=X_corr_csv,by=c("CHR","BIN.START","BIN.END"),all.x = T)
#Y_corr <- merge(x=binned.raw[,c(2:4,7)],y=Y_corr_csv,by=c("CHR","BIN.START","BIN.END"),all.x = T)

#X_corr[is.na(X_corr$Intercept),"Intercept"] <- 0
#Y_corr[is.na(Y_corr$Intercept),"Intercept"] <- 0
#X_corr[is.na(X_corr$Rico),"Rico"] <- 0
#Y_corr[is.na(Y_corr$Rico),"Rico"] <- 0

#last column
sample.count <- names(binned.raw)[ncol(binned.raw)]
sample.gc <- sub("COUNT","GC",sample.count)

sample.gr10M <- sub("COUNT","GR10M",sample.count)
sample.gr <- sub("COUNT","GR",sample.count)

sample.gc.gr10M <- sub("COUNT","GC.GR10M",sample.count)
sample.gc.gr <- sub("COUNT","GC.GR",sample.count)

sample.auto.gr10M <- sub("COUNT","AUTO.GR10M",sample.count)
sample.auto.gr <- sub("COUNT","AUTO.GR",sample.count)

sample.gc.auto.gr10M <- sub("COUNT","GC.AUTO.GR10M",sample.count)
sample.gc.auto.gr <- sub("COUNT","GC.AUTO.GR",sample.count)

#sample.auto.X.gr10M <- sub("COUNT","AUTO.X.GR10M",sample.count)
#sample.auto.X.gr <- sub("COUNT","AUTO.X.GR",sample.count)
#sample.gc.auto.X.F.gr10M <- sub("COUNT","GC.AUTO.X.F.GR10M",sample.count)
#sample.gc.auto.X.M.gr10M <- sub("COUNT","GC.AUTO.X.M.GR10M",sample.count)
#sample.gc.auto.X.F.gr <- sub("COUNT","GC.AUTO.X.F.GR",sample.count)
#sample.gc.auto.X.M.gr <- sub("COUNT","GC.AUTO.X.M.GR",sample.count)
#sample.auto.Y.gr10M <- sub("COUNT","AUTO.Y.GR10M",sample.count)
#sample.auto.Y.gr <- sub("COUNT","AUTO.Y.GR",sample.count)
#sample.gc.auto.Y.F.gr10M <- sub("COUNT","GC.AUTO.Y.F.GR10M",sample.count)
#sample.gc.auto.Y.M.gr10M <- sub("COUNT","GC.AUTO.Y.M.GR10M",sample.count)
#sample.gc.auto.Y.F.gr <- sub("COUNT","GC.AUTO.Y.F.GR",sample.count)
#sample.gc.auto.Y.M.gr <- sub("COUNT","GC.AUTO.Y.M.GR",sample.count)

#ignore bins that have N's
noNs <- subset(binned.raw,binned.raw$BIN.N.COUNT==0)
noNs.notEmpty <- subset(noNs,noNs[ncol(noNs)]>0)

print("Calculating gam")
noNs.notEmpty.gam <- gam(noNs.notEmpty[,ncol(noNs.notEmpty)]~noNs.notEmpty$BIN.GC.CONTENT)

print("Correcting")
gam.predicted <- data.frame(predict.gam(noNs.notEmpty.gam))

#e <- ggplot(noNs.notEmpty, aes(x=BIN.GC.CONTENT, y=noNs.notEmpty[,ncol(noNs.notEmpty)]))
#e + geom_point() + geom_smooth()


#############################


#print("Calculating loess")
#noNs.notEmpty.loess <- loess(noNs.notEmpty[,ncol(noNs.notEmpty)]~noNs.notEmpty$BIN.GC.CONTENT)

#print("Correcting")
#loess.predicted <- data.frame(predict(noNs.notEmpty.loess,noNs.notEmpty$BIN.GC.CONTENT))


correction.factor <- median(noNs.notEmpty[,ncol(noNs.notEmpty)])/gam.predicted
names(correction.factor)[1] <- "F"
#names(correction.factor)
noNs.notEmpty[[sample.gc]] <- round(100*noNs.notEmpty[,ncol(noNs.notEmpty)]*correction.factor$F)/100

#=====================================================================================================
autoMedian=median(noNs.notEmpty[[sample.gc]][noNs.notEmpty$CHR!="X"&noNs.notEmpty$CHR!="Y"])
#XMedian=median(noNs.notEmpty[[sample.gc]][noNs.notEmpty$CHR=="X"])
#YMedian=median(noNs.notEmpty[[sample.gc]][noNs.notEmpty$CHR=="Y"])
#ffX = 2*(autoMedian - XMedian) / autoMedian
#ffY = 2* YMedian / autoMedian

#print(paste("ffY:",ffY))
#print(paste("ffX:",ffX))
#write(XMedian,file = paste(basename,"XMedian",sep="."),sep="")
#write(YMedian,file = paste(basename,"YMedian",sep="."),sep="")
write(autoMedian,file = paste(basename,"autoMedian",sep="."),sep="")

#write(ffX,file = paste(basename,"fetal.fraction.X",sep="."),sep="")
#write(ffY,file = paste(basename,"fetal.fraction.Y",sep="."),sep="")
#=====================================================================================================

#add corrected GC-count to the original count file
binned.corrected <- merge(binned.raw,noNs.notEmpty,all.x=TRUE,by=names(binned.raw))
binned.corrected[[sample.gc]][is.na(binned.corrected[[sample.gc]])] <- 0

#subset only the autosoma chromosomes
binned.corrected.auto <- subset(binned.corrected,!(binned.corrected$CHR == "X" | binned.corrected$CHR == "Y"))
#binned.corrected.auto.X <- subset(binned.corrected,!(binned.corrected$CHR == "Y"))
#binned.corrected.auto.Y <- subset(binned.corrected,!(binned.corrected$CHR == "X"))

#correct for the total number of counts for that sample (for both raw reads and gc corrected counts)
binned.corrected[[sample.gr10M]] <- 10000000 * binned.corrected[[sample.count]]/sum(binned.corrected[[sample.count]])
binned.corrected[[sample.gc.gr10M]] <-  10000000 * binned.corrected[[sample.gc]]/sum(binned.corrected[[sample.gc]])

#correct for the total number of counts of only the autosomal chromosomes
binned.corrected[[sample.auto.gr10M]] <- 10000000 * binned.corrected[[sample.count]]/sum(binned.corrected.auto[[sample.count]])
#binned.corrected[[sample.auto.X.gr10M]] <- 10000000 * binned.corrected[[sample.count]]/sum(binned.corrected.auto.X[[sample.count]])
#binned.corrected[[sample.auto.Y.gr10M]] <- 10000000 * binned.corrected[[sample.count]]/sum(binned.corrected.auto.Y[[sample.count]])

binned.corrected[[sample.gc.auto.gr10M]] <-  10000000 * binned.corrected[[sample.gc]]/sum(binned.corrected.auto[[sample.gc]])
#binned.corrected[[sample.gc.auto.X.F.gr10M]] <-  10000000 * binned.corrected[[sample.gc]]/sum(binned.corrected.auto.X[[sample.gc]])
#binned.corrected[[sample.gc.auto.X.M.gr10M]] <-  10000000 * binned.corrected[[sample.gc]]/sum(binned.corrected.auto.X[[sample.gc]])
#binned.corrected[[sample.gc.auto.Y.F.gr10M]] <-  10000000 * binned.corrected[[sample.gc]]/sum(binned.corrected.auto.Y[[sample.gc]])
#binned.corrected[[sample.gc.auto.Y.M.gr10M]] <-  10000000 * binned.corrected[[sample.gc]]/sum(binned.corrected.auto.Y[[sample.gc]])

#=====================================================================================================
#if(!is.na(seqFF)){
#	binned.corrected[[sample.gc.auto.X.M.gr10M]] <- binned.corrected[[sample.gc.auto.X.M.gr10M]]-((seqFF*X_corr$Rico)+X_corr$Intercept)
#	binned.corrected[[sample.gc.auto.Y.M.gr10M]] <- binned.corrected[[sample.gc.auto.Y.M.gr10M]]-((seqFF*Y_corr$Rico)+Y_corr$Intercept)
#}
#=====================================================================================================


binned.corrected[[sample.auto.gr10M]][binned.corrected$CHR=="X" | binned.corrected$CHR=="Y"] <- 0
#binned.corrected[[sample.auto.X.gr10M]][binned.corrected$CHR=="Y"] <- 0
#binned.corrected[[sample.auto.Y.gr10M]][binned.corrected$CHR=="X"] <- 0

binned.corrected[[sample.gc.auto.gr10M]][binned.corrected$CHR=="X" | binned.corrected$CHR=="Y"] <- 0
#binned.corrected[[sample.gc.auto.X.F.gr10M]][binned.corrected$CHR=="Y"] <- 0
#binned.corrected[[sample.gc.auto.X.M.gr10M]][binned.corrected$CHR=="Y"] <- 0
#binned.corrected[[sample.gc.auto.Y.F.gr10M]][binned.corrected$CHR=="X"] <- 0
#binned.corrected[[sample.gc.auto.Y.M.gr10M]][binned.corrected$CHR=="X"] <- 0

print("Aggregating per chromosome")

chr.corrected <- aggregate(subset(binned.corrected,select=-c(SAMPLE,CHR,BIN.START,BIN.END,BIN.GC.CONTENT,BIN.N.COUNT)),by=list(binned.corrected$SAMPLE, binned.corrected$CHR),FUN=sum)
colnames(chr.corrected)<-c("SAMPLE","CHR",sample.count,sample.gc,sample.gr,sample.gc.gr,sample.auto.gr,sample.gc.auto.gr)


chr.corrected[[sample.gr]] <- chr.corrected[[sample.gr]] / 10000000
chr.corrected[[sample.gc.gr]] <- chr.corrected[[sample.gc.gr]] / 10000000

chr.corrected[[sample.auto.gr]] <- chr.corrected[[sample.auto.gr]] / 10000000
#chr.corrected[[sample.auto.X.gr]] <- chr.corrected[[sample.auto.X.gr]] / 10000000
#chr.corrected[[sample.auto.Y.gr]] <- chr.corrected[[sample.auto.Y.gr]] / 10000000

chr.corrected[[sample.gc.auto.gr]] <- chr.corrected[[sample.gc.auto.gr]] / 10000000
#chr.corrected[[sample.gc.auto.X.F.gr]] <- chr.corrected[[sample.gc.auto.X.F.gr]] / 10000000
#chr.corrected[[sample.gc.auto.X.M.gr]] <- chr.corrected[[sample.gc.auto.X.M.gr]] / 10000000
#chr.corrected[[sample.gc.auto.Y.F.gr]] <- chr.corrected[[sample.gc.auto.Y.F.gr]] / 10000000
#chr.corrected[[sample.gc.auto.Y.M.gr]] <- chr.corrected[[sample.gc.auto.Y.M.gr]] / 10000000

chr.corrected[[sample.auto.gr]][chr.corrected$CHR=="X" | chr.corrected$CHR=="Y"] <- 0
#chr.corrected[[sample.auto.X.gr]][chr.corrected$CHR=="Y"] <- 0
#chr.corrected[[sample.auto.Y.gr]][chr.corrected$CHR=="X"] <- 0

chr.corrected[[sample.gc.auto.gr]][chr.corrected$CHR=="X" | chr.corrected$CHR=="Y"] <- 0
#chr.corrected[[sample.gc.auto.X.F.gr]][chr.corrected$CHR=="Y"] <- 0
#chr.corrected[[sample.gc.auto.X.M.gr]][chr.corrected$CHR=="Y"] <- 0
#chr.corrected[[sample.gc.auto.Y.F.gr]][chr.corrected$CHR=="X"] <- 0
#chr.corrected[[sample.gc.auto.Y.M.gr]][chr.corrected$CHR=="X"] <- 0

print("Writing")
?write.csv
write.csv(binned.corrected, file=paste(basename,SAMPLE,"binned.csv", sep = "."), row.names=FALSE)
write.csv(chr.corrected, file=paste(basename,SAMPLE,"chr.csv", sep = "."), row.names=FALSE)
