#load library
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("DNAcopy")

library(grid)
library(gridExtra)
library(lattice)
library(ggplot2)
library(DNAcopy)
options(warn=1)

# read arguments
#args <- commandArgs(TRUE)
#bin_file <- args[1]
#samplename <- args[2]
samplename <- 127869
#gender <- args[3]
gender <- "MALE"
#ref_map <- args[4]

#output_map <- args[5]
output_map <- "C:/Users/jeann/Desktop/CBS/"

#tresholdadel <- as.numeric(args[6])
tresholdadel <- as.numeric(-2)

#tresholdadup <- as.numeric(args[7])
tresholdadup <- as.numeric(2)

#tresholdxdel <- as.numeric(args[8])
tresholdxdel <- as.numeric(-3)

#tresholdxdup <- as.numeric(args[9])
tresholdxdup <- as.numeric(3)

#adel <- as.integer(args[10])
adel <- as.integer(8)

#adup <- as.integer(args[11])
adup <- as.integer(8)

#xdel <- as.integer(args[12])
xdel <- as.integer(8)

#xdup <- as.integer(args[13])
xdup <- as.integer(8)

#diagBandcsv <- args[14]

#cytoBandcsv <- args[15]

bintreshold <- matrix(c(adel,adup,xdel,xdup),ncol=2,nrow=2,byrow=T)
#           ||  del  |  dup
#==========================
# Autosomen ||  adel | adup
#           ||-------------
# Xchrom    ||  xdel | xdup

# set parameters
genders = c('jongen' = 'M', 'meisje' = 'F', 'male' = 'M', 'female' = 'F', 'm' = 'M', 'f' = 'F')
min.width.bins = 2
alpha=0.001

panel.fill<-function(x,y,mychr,CBSsegments){
  normalIndex<-which(abs(y)<10)
  outliersIndexhigh <- which(y >=  10)
  outliersIndexlow  <- which(y <= -10)
  
  panel.xyplot(x[normalIndex],y[normalIndex],pch=20,cex=0.3,col="black") # normal points
  panel.points(x[outliersIndexhigh],10,pch=3,cex=0.3,col="black")        # outliers high
  panel.points(x[outliersIndexlow],-10,pch=3,cex=0.3,col="black")        # outliers low

  chrSegments = subset(CBSsegments,CBSsegments$chrom==mychr)
  
  myLwd = 0.4
  
  diagnosticBand <- read.csv("C:/Users/jeann/Desktop/diagnosticBand.b37.csv",header=TRUE)
  cytoBand <- read.csv("C:/Users/jeann/Desktop/cytoBand.b37.csv",header=TRUE)
  
  cyto.H = 0.10
  cyto.Y = 0.07 #0.85
  
  cyto.chr <-  subset(cytoBand,cytoBand$CHR==as.character(mychr))
  diag.chr <-  subset(diagnosticBand,diagnosticBand$CHR==as.character(mychr))
  cyto.max <- max(cyto.chr$END)
  
  chrSegments$loc.start = chrSegments$loc.start/cyto.max
  chrSegments$loc.end = chrSegments$loc.end/cyto.max 
  
  for(i in seq(1,dim(chrSegments)[1])){
    if(chrSegments[i,"Z"] >= 10){
      panel.segments(x0=chrSegments[i,"loc.start"],x1=chrSegments[i,"loc.end"],y0=10,y1=10,col="orange",lwd=3)
    } else if(chrSegments[i,"Z"] <= -10){
      panel.segments(x0=chrSegments[i,"loc.start"],x1=chrSegments[i,"loc.end"],y0=-10,y1=-10,col="orange",lwd=3)
    } else {
      panel.segments(x0=chrSegments[i,"loc.start"],x1=chrSegments[i,"loc.end"],y0=chrSegments[i,"Z"],y1=chrSegments[i,"Z"],col="red",lwd=3)
    }
  }
  
  panel.abline(h=0,col=rgb(0,0,0,0.5),lwd=2) #zero
  if(mychr == "X"){
    panel.abline(h=tresholdxdup,col=rgb(0,0,1,0.5),lwd=2,lty=3)   #CNV treshold
    panel.abline(h=tresholdxdel,col=rgb(0,0,1,0.5),lwd=2,lty=3)
  } else {
    panel.abline(h=tresholdadup,col=rgb(0,0,1,0.5),lwd=2,lty=3)   #CNV treshold
    panel.abline(h=tresholdadel,col=rgb(0,0,1,0.5),lwd=2,lty=3)
  }

  for (i in 1:nrow(cyto.chr)) {
    row <- cyto.chr[i,]
    from <- row$START/cyto.max
    to <- row$END/cyto.max
    
    if (row$TYPE!="acen") { #acen done at the end
      alpha <- 0
      if (row$TYPE=="gpos25") {alpha <- 0.25}
      else if (row$TYPE=="gpos50" || row$TYPE=="stalk") {alpha <- 0.5}
      else if (row$TYPE=="gpos75") {alpha <- 0.75}
      else if (row$TYPE=="gpos100" || row$TYPE=="gvar") {alpha <- 1}
      grid.rect(x=unit(from,"native"),y=cyto.Y,width = unit(to-from,"native"), height = cyto.H, just=c("left","bottom"), gp=gpar(fill=rgb(0,0,0,alpha),col="transparent",lwd=myLwd))
    }
  }
  
  myY = 0.20
  myH = 0.75
  myLwd = 2
  
  if (nrow(diag.chr) > 0) {
    for (i in 1:nrow(diag.chr)) {
      row <- diag.chr[i,]
      from <- row$START/cyto.max
      to <- row$END/cyto.max
      
      myLty <- "dotted"
      
      grid.lines(x=unit(from,"native"),y=c(myY,myY+myH),gp=gpar(col=uzcol.light2,lty=myLty,lwd=myLwd))
      grid.lines(x=unit(to,"native"),y=c(myY,myY+myH),gp=gpar(col=uzcol.light2,lty=myLty,lwd=myLwd))
    }
  }
  
  #finalize cytoband 
  myLwd = 0.4
  grid.rect(x=unit(0,"native"),y=cyto.Y,width = unit(1,"native"), height = cyto.H, just=c("left","bottom"), gp=gpar(fill="transparent",col="black",lwd=myLwd))  #cytogen band frame
  
  for (i in 1:nrow(cyto.chr)) {
    row <- cyto.chr[i,]
    
    if (row$TYPE=="acen") {
      from <- row$START/cyto.max
      to <- row$END/cyto.max
      grid.rect(x=unit(from,"native"),y=cyto.Y-0.01,width = unit(to-from,"native"), height = cyto.H+0.02, just=c("left","bottom"), gp=gpar(fill="white",col="white",lwd=myLwd))
      myXList <- c(from,from,to)
      myYList <- c(cyto.Y,cyto.Y+cyto.H,cyto.Y+cyto.H/2)
      if (grepl("^q",row$NAME)) {
        myXList <- c(from,to,to)
        myYList <- c(cyto.Y+cyto.H/2,cyto.Y,cyto.Y+cyto.H)
      }
      grid.polygon(x=unit(myXList,"native"),y=myYList, gp=gpar(fill="indianred3",col="indianred3",lwd=myLwd))
    }
  }
}

# read data
#bin_data = read.table(bin_file,as.is=T,header=T,sep=",")
bin_data <- read.table("C:/Users/jeann/Desktop/127869.count.stdout.5890115.pbs.biomina.be",as.is=T,header=F,sep="\t")
names(bin_data)[1:4] <- list("CHR", "BIN.START", "BIN.END", "COUNT")

sample.count <- names(bin_data)[7]
#sample.gc.auto.X.F.gr10M <- sub("COUNT","GC.AUTO.X.F.GR10M",sample.count)
#sample.gc.auto.X.M.gr10M <- sub("COUNT","GC.AUTO.X.M.GR10M",sample.count)
#sample.gc.auto.gr10M <- sub("COUNT","GC.AUTO.GR10M",sample.count)


# ref data
## input file is loaded here --> contains header, csv file
#ref_bin_file = paste(ref_map,"referenceNormals.stats.binned.csv",sep="/");
ref_data <- read.table("C:/Users/jeann/Desktop/referenceNormals.stats.binned.csv",as.is=T,header=F,sep=",")
names(ref_data)[1:5] <- list("CHR", "BIN.START", "BIN.END", "COUNT.MEAN", "COUNT.SD")


# get gender specific ref data.
#gender_ref_data <- ref_data;
#if(!(is.na(genders[tolower(gender)]))){
#    gender_ref_bin_file <- sub("stats.binned.csv",paste('stats.',genders[tolower(gender)],'.binned.csv',sep=''),ref_bin_file,fixed=TRUE)
#    gender_ref_data <- read.table(gender_ref_bin_file,as.is=T,header=T,sep=",")
#} else {

#    message(paste("Using standard reference data for gender specific reference data (unknown gender: ",gender,")",sep=""));
#    gender_ref_bin_file <- sub("stats.binned.csv",paste('stats.','M','.binned.csv',sep=''),ref_bin_file,fixed=TRUE)
#    gender_ref_data <- read.table(gender_ref_bin_file,as.is=T,header=T,sep=",")
#    message(paste("Using standard male data for gender specific reference data for unknown gender: ",gender,sep=""));
    # reason male ref: 
    # - figures chr X and Y on report are sample vs male when undetermined
    # - more likely to have male < gendertreshold than XXY
#}

#autosomes: with normal ref
#x: with gender ref

achrom = bin_data[which(bin_data$CHR != "chrX" & bin_data$CHR != "chrY"),]	#autosomes
achrom_full <- merge(achrom,ref_data,by=c('CHR','BIN.START','BIN.END'))
achrom_full$CHR <- gsub("[a-zA-Z ]", "", achrom_full$CHR)


# Y: no CNVcalling
# sd_full$ratio = sd_full[,16] / sd_full[,35]
#	sd_full$z = (sd_full[,16] - sd_full[35]) / sd_full[,36]

# no gender data
achrom_full$ratio = achrom_full[,4] / achrom_full[,5]
achrom_full$z = (achrom_full[,4] - achrom_full[,5]) / achrom_full[,6]
achrom_full$logR = log2(achrom_full$ratio)

#browseVignettes("DNAcopy")
aCNA <- CNA(as.matrix(achrom_full$z),achrom_full$CHR,as.numeric(achrom_full$BIN.START),data.type="logratio",sampleid=samplename)
aObj <- segment(aCNA,verbose=0,min.width = min.width.bins,alpha=alpha)
aSeg <- aObj$output
aSeg$Sample <- samplename

CNVtable = aSeg

CNVtable$loc.end = CNVtable$loc.end + 49999

#add Ratio to table
tmp = matrix(nrow=0,ncol=8)
for(i in seq(1,dim(CNVtable)[1])){
  chr = as.character(CNVtable[i,"chrom"])
  start = as.integer(CNVtable[i,"loc.start"])
  end = as.integer(CNVtable[i,"loc.end"])
  if(chr=="X"){
    cnv_full = xchrom_full[which(xchrom_full$CHR == chr & xchrom_full$BIN.START >= start & xchrom_full$BIN.END <= end),]
  } else {
    cnv_full = achrom_full[which(achrom_full$CHR == chr & achrom_full$BIN.START >= start & achrom_full$BIN.END <= end),]
  }
  tmp = rbind(tmp,c(as.character(samplename),as.character(chr),start,end,
                    as.integer(dim(cnv_full)[1]),as.integer(dim(cnv_full)[1])-sum(!is.finite(cnv_full$ratio)),
                    as.numeric(round(mean(cnv_full[is.finite(cnv_full$z),"z"],na.rm=T),2)),
                    as.numeric(round(mean(cnv_full[is.finite(cnv_full$ratio),"ratio"],na.rm=T),2))))
}
CNVtable = data.frame(tmp,stringsAsFactors = F)
colnames(CNVtable)<-c("Sample","chrom","loc.start","loc.end","bins","NoNA.bins","Z","Ratio")
CNVtable$Z = as.numeric(CNVtable$Z)
CNVtable$Ratio = as.numeric(CNVtable$Ratio)
CNVtable$loc.start = as.integer(CNVtable$loc.start)
CNVtable$loc.end = as.integer(CNVtable$loc.end)
CNVtable$bins = as.integer(CNVtable$bins)
CNVtable$NoNA.bins = as.integer(CNVtable$NoNA.bins)

CNVtable=CNVtable[sort.list(suppressWarnings(as.integer(CNVtable$chrom)),na.last=T),] 
#as.integer: X => NA: needs suppressWarnings + X gets sorted last by na.last=T

SegmentListoutput <- paste(output_map,paste(samplename,'CBS.SegmentList.csv',sep='.'),sep="/");
write.table(CNVtable,file=SegmentListoutput,sep="\t",row.names=FALSE,quote=F);

nopadding.heights <-
  list(top.padding = 1,
       main.key.padding = 0,
       key.axis.padding = 0,
       axis.xlab.padding = 0,
       xlab.key.padding = 0,
       key.sub.padding = 0,
       bottom.padding = 1)
nopadding.widths <-
  list(left.padding = 1,
       key.ylab.padding = 0,
       ylab.axis.padding = 0,
       axis.key.padding = 0,
       right.padding = 0)


#CHR plots
for(chromosome in seq(1,22)){
  pdf(paste(output_map,paste(samplename,".","chr",chromosome,".CBSplot.pdf",sep=""),sep="/"))
  achrom_chr = subset(achrom_full,achrom_full$CHR==as.character(chromosome))
  achrom_chr$CHR = factor(achrom_chr$CHR)
  achrom_chr$BIN.START.REL = achrom_chr$BIN.START / max(achrom_chr$BIN.START)
  print(xyplot(achrom_chr$z~achrom_chr$BIN.START.REL,
         xlab="Position (Mbp)",ylab=list(label="Z-score",y=unit(0.6,"native")),panel=panel.fill,mychr=as.character(chromosome),
         ylim=c(-15,10.2),CBSsegments=CNVtable,
         strip=F,scales=list(col=1,tck=c(1,0),
                             y=list(draw=T,at=seq(-10,10,2)),
                             x=list(draw=T,at=seq(0,1,0.1),labels = round(seq(0,10,1)*max(achrom_chr$BIN.START)/10e6,1))),
         par.settings=list(strip.background=list(col="grey"), 
                           layout.heights = nopadding.heights, 
                           layout.widths = nopadding.widths,
                           axis.line = list(col=0)),
         main=paste("CBS-Plot ",samplename,":chr",chromosome,sep="")))
  dev.off()
}

if(exists("xObj")){
  pdf(paste(output_map,paste(samplename,".","chr","X",".CBSplot.pdf",sep=""),sep="/"))
  xchrom_full$BIN.START.REL = xchrom_full$BIN.START / max(xchrom_full$BIN.START)
  print(xyplot(xchrom_full$z~xchrom_full$BIN.START.REL,
         xlab="Position (Mbp)",ylab=list(label="Z-score",y=unit(0.6,"native")),panel=panel.fill,mychr="X",
         ylim=c(-15,10.2),CBSsegments=CNVtable,
         strip=F,scales=list(col=1,tck=c(1,0),
                             y=list(draw=T,at=seq(-10,10,2)),
                             x=list(draw=T,at=seq(0,1,0.1),labels = round(seq(0,10,1)*max(xchrom_full$BIN.START)/10e6,1))),
         par.settings=list(strip.background=list(col="grey"), 
                           layout.heights = nopadding.heights, 
                           layout.widths = nopadding.widths,
                           axis.line = list(col=0)),
         main=paste("CBS-Plot ",samplename,":chrX",sep="")))
  dev.off()
}

#FILTER SCORE
trueCNVs = CNVtable[which(
                    CNVtable$chrom != "X" & (CNVtable$Z >= tresholdadup | CNVtable$Z <= tresholdadel) |
                    CNVtable$chrom == "X" & (CNVtable$Z >= tresholdxdup | CNVtable$Z <= tresholdxdel)
                    )
                  ,]

#Combine consecutive CNVs: 
#if start = end of previous row => combine and remove current row + adjust index
i = 2
while(i <= dim(trueCNVs)[1]){
  if(trueCNVs[i,"chrom"]==trueCNVs[i-1,"chrom"] & trueCNVs[i,"loc.start"]==(trueCNVs[i-1,"loc.end"]+1) & sign(trueCNVs[i,"Z"])*sign(trueCNVs[i-1,"Z"]) > 0 ){
    trueCNVs[i-1,"loc.end"] = trueCNVs[i,"loc.end"]
    trueCNVs[i-1,"bins"] = trueCNVs[i,"bins"] + trueCNVs[i-1,"bins"] #sum bin size
    trueCNVs[i-1,"NoNA.bins"] = trueCNVs[i,"NoNA.bins"] + trueCNVs[i-1,"NoNA.bins"] #sum bin size
    trueCNVs = trueCNVs[-i,]
    i = i - 1
  }
  i = i + 1
}

#FILTER BINS
GJB6_start = 20750000
GJB6_end = 21100000

#FILTER BINS
trueCNVs = trueCNVs[ which(
  (as.character(trueCNVs$chrom) != "X" & as.numeric(trueCNVs$Z) < 0 & as.numeric(trueCNVs$bins) >= adel) |
  (as.character(trueCNVs$chrom) != "X" & as.numeric(trueCNVs$Z) > 0 & as.numeric(trueCNVs$bins) >= adup) |
  (as.character(trueCNVs$chrom) == "X" & as.numeric(trueCNVs$Z) < 0 & as.numeric(trueCNVs$bins) >= xdel) |
  (as.character(trueCNVs$chrom) == "X" & as.numeric(trueCNVs$Z) > 0 & as.numeric(trueCNVs$bins) >= xdup) |
  (
    as.character(trueCNVs$chrom) == "13" &
    as.numeric(trueCNVs$loc.start) <= GJB6_end &
    as.numeric(trueCNVs$loc.end) >= GJB6_start &
    as.numeric(trueCNVs$Z) < 0 &
    as.numeric(trueCNVs$bins) >= 6
  )
),]
list=matrix(ncol=7,nrow=0)
if(dim(trueCNVs)[1]>0){
  for(i in seq(1,dim(trueCNVs)[1])){
    chr = as.character(trueCNVs[i,"chrom"])
    start = as.integer(trueCNVs[i,"loc.start"])
    end = as.integer(trueCNVs[i,"loc.end"])
    if(chr=="X"){
      sd_full = xchrom_full[which(xchrom_full$CHR == chr & xchrom_full$BIN.START >= start - 10*50000  & xchrom_full$BIN.END <= end + 1 + 10*50000),]
      cnv_full = xchrom_full[which(xchrom_full$CHR == chr & xchrom_full$BIN.START >= start & xchrom_full$BIN.END <= end),]
    } else {
      sd_full = achrom_full[which(achrom_full$CHR == chr & achrom_full$BIN.START >= start - 10*50000  & achrom_full$BIN.END <= end + 1 + 10*50000),]
      cnv_full = achrom_full[which(achrom_full$CHR == chr & achrom_full$BIN.START >= start & achrom_full$BIN.END <= end),]
    }
    
    list = rbind(list,c(samplename,chr,paste(start,end,sep="-"),dim(cnv_full)[1],dim(cnv_full)[1]-sum(!is.finite(cnv_full$ratio)),
                        round(mean(cnv_full[is.finite(cnv_full$z),"z"],na.rm=T),2),round(mean(cnv_full[is.finite(cnv_full$ratio),"ratio"],na.rm=T),2)))
    
    sd_full$pos = sd_full$BIN.START + 25000
    
    dataoutput <- paste(output_map,paste(samplename, paste("chr",chr,sep=""),paste(start,end,sep="-"),'CNVData.csv',sep='.'),sep="/");
    write.table(sd_full,file=dataoutput,sep="\t",row.names=FALSE,quote=F);
    
    sd_full$logR[which(sd_full$logR < -4)] <- -3.9
    sd_full$logR[which(sd_full$logR > 4)] <- 3.9
    sd_full$z[which(sd_full$z < -4)] <- -3.9
    sd_full$z[which(sd_full$z > 4)] <- 3.9
    sd_full$ratio[which(sd_full$ratio < -4)] <- -3.9
    sd_full$ratio[which(sd_full$ratio > 4)] <- 3.9
    
    # plot.
    f <- ggplot(sd_full,aes(x=pos, y=value , color=Variable)) + geom_point(aes(y=logR,col="log2R")) + geom_point(aes(y=z,col='Z-score'))
    # legendlabels
    f <- f + scale_color_hue(labels=c(bquote(paste("",Log[2],"R")),"Z-score"))
    # set scale
    f <- f + scale_y_continuous(limits=c(-4,4),breaks=c(-3,-1.5,0,1.5,3),labels=c(-3,-1.5,0,1.5,3),minor_breaks=seq(-4,4))
    # add title
    f <- f + labs(x='Position',y='Value', title=paste(samplename,":",'chr',chr,':',start,'-',end,sep=""))
    f <- f + theme(plot.title = element_text(hjust=0.5))
    # add guiding-lines
    f <- f + annotate("rect",xmin=start,xmax=end,ymin=-4,ymax=4,color='red',alpha=.1)
    #f <- f + geom_segment(x=start,xend=start,y=-4,yend=4,color='red');
    #f <- f + geom_segment(x=end,xend=end,y=-4,yend=4,color='red');
    f <- f + geom_hline(yintercept=1.5, color = "black", size = 0.2);
    f <- f + geom_hline(yintercept=-1.5, color = "black", size = 0.2);
    file <- paste(output_map,paste(samplename, paste("chr",chr,sep=""),paste(start,end,sep="-"),'CNVPlot.pdf',sep='.'),sep="/")
    # save
    ggsave(file,f,height=8,width=20,units="cm")
  
  }
}
CNVlist = data.frame(list)
colnames(CNVlist) = c("Sample","Chromosome","Position","Number.Bins","Informative.Bins","Z.Mean","Read.Ratio")
CNVlist=CNVlist[sort.list(suppressWarnings(as.integer(as.character(CNVlist$Chromosome))),na.last=T),] 
#suppres warnings: as.integer changes X to NA, which are then set as last lines

listoutput <- paste(output_map,paste(samplename, "CNVs.txt", sep="."),sep="/")
write.table(CNVlist,file=listoutput,sep=" ",row.names=FALSE,quote=F)
