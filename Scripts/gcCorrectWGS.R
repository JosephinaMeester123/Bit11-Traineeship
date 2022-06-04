#READ DATA

#library(gam)
library(mgcv)
# ARGUMENT :
#  1 : input file (counts / bins)
#  2 : sample name
#  3 : output folder

args<-commandArgs(TRUE)
#binned.raw.counts.txt <- args[1]
binned.raw.counts.txt <- "/home/jmeester/Internship/countfiles/108995.count"

#SAMPLE <- args[2]
SAMPLE <- 108995

binned.csv <- paste("/home/jmeester/Internship/GC/", SAMPLE, ".binned.csv", sep="")

chr.csv <- paste("/home/jmeester/Internship/GC/", SAMPLE, ".chr.csv", sep="")

#output folder <- args[3]
basename <- "/home/jmeester/Internship/GC/"


#read gc content file
gc_content <- read.table("/home/jmeester/Internship/bins500bp_gc.bed",as.is=T,header=F,sep="\t")
colnames(gc_content) <- c('CHR','BIN.START','BIN.END', 'BIN.GC.CONTENT', 'BIN.N.COUNT')

#read raw_count file
raw_count <- read.table(binned.raw.counts.txt, as.is=T,header=F,sep="\t")
colnames(raw_count) <- c('CHR','BIN.START','BIN.END', 'COUNT')

#merge gc conent file with raw count file and add sample name
merged <- merge(gc_content,raw_count,by=c('CHR','BIN.START','BIN.END'))
binned.raw <- cbind(SAMPLE, merged)


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


#ignore bins that have N's
noNs <- subset(binned.raw,binned.raw$BIN.N.COUNT==0)
noNs.notEmpty <- subset(noNs,noNs[ncol(noNs)]>0)
colnames(raw_count) <- c('CHR','BIN.START','BIN.END', 'COUNT')

##########################

print("Calculating gam")
noNs.notEmpty.gam <- gam(noNs.notEmpty[,ncol(noNs.notEmpty)]~noNs.notEmpty$BIN.GC.CONTENT)

print("Correcting")
gam.predicted <- data.frame(predict.gam(noNs.notEmpty.gam))

#############################

correction.factor <- median(noNs.notEmpty[,ncol(noNs.notEmpty)])/gam.predicted
names(correction.factor)[1] <- "F"
#names(correction.factor)
noNs.notEmpty[[sample.gc]] <- round(100*noNs.notEmpty[,ncol(noNs.notEmpty)]*correction.factor$F)/100


#=====================================================================================================
autoMedian=median(noNs.notEmpty[[sample.gc]][noNs.notEmpty$CHR!="X"&noNs.notEmpty$CHR!="Y"])
write(autoMedian,file = paste(basename,SAMPLE,".autoMedian",sep=""))
#=====================================================================================================


#add corrected GC-count to the original count file
binned.corrected <- merge(binned.raw,noNs.notEmpty,all.x=TRUE,by=names(binned.raw))
binned.corrected[[sample.gc]][is.na(binned.corrected[[sample.gc]])] <- 0

#subset only the autosomal chromosomes
binned.corrected.auto <- subset(binned.corrected,!(binned.corrected$CHR == "X" | binned.corrected$CHR == "Y"))

#correct for the total number of counts for that sample (for both raw reads and gc corrected counts)
binned.corrected[[sample.gr10M]] <- 10000000 * binned.corrected[[sample.count]]/sum(binned.corrected[[sample.count]])
binned.corrected[[sample.gc.gr10M]] <-  10000000 * binned.corrected[[sample.gc]]/sum(binned.corrected[[sample.gc]])

#correct counts for the total number of counts of only the autosomal chromosomes
binned.corrected[[sample.auto.gr10M]] <- 10000000 * binned.corrected[[sample.count]]/sum(binned.corrected.auto[[sample.count]])

#correct gc corrected counts for the total number of counts of only the autosomal chromosomes
binned.corrected[[sample.gc.auto.gr10M]] <-  10000000 * binned.corrected[[sample.gc]]/sum(binned.corrected.auto[[sample.gc]])

#=====================================================================================================

binned.corrected[[sample.auto.gr10M]][binned.corrected$CHR=="X" | binned.corrected$CHR=="Y"] <- 0

binned.corrected[[sample.gc.auto.gr10M]][binned.corrected$CHR=="X" | binned.corrected$CHR=="Y"] <- 0

print("Aggregating per chromosome")
chr.corrected <- aggregate(subset(binned.corrected,select=-c(SAMPLE,CHR,BIN.START,BIN.END,BIN.GC.CONTENT,BIN.N.COUNT)),by=list(binned.corrected$SAMPLE, binned.corrected$CHR),FUN=sum)
colnames(chr.corrected)<-c("SAMPLE","CHR",sample.count,sample.gc,sample.gr,sample.gc.gr,sample.auto.gr,sample.gc.auto.gr)

chr.corrected[[sample.gr]] <- chr.corrected[[sample.gr]] / 10000000
chr.corrected[[sample.gc.gr]] <- chr.corrected[[sample.gc.gr]] / 10000000

chr.corrected[[sample.auto.gr]] <- chr.corrected[[sample.auto.gr]] / 10000000

chr.corrected[[sample.gc.auto.gr]] <- chr.corrected[[sample.gc.auto.gr]] / 10000000

chr.corrected[[sample.auto.gr]][chr.corrected$CHR=="X" | chr.corrected$CHR=="Y"] <- 0

chr.corrected[[sample.gc.auto.gr]][chr.corrected$CHR=="X" | chr.corrected$CHR=="Y"] <- 0

print("Writing")
write.csv(binned.corrected, file=paste(basename,SAMPLE,"binned.csv", sep = ""), row.names=FALSE)
write.csv(chr.corrected, file=paste(basename,SAMPLE,"chr.csv", sep = ""), row.names=FALSE)
