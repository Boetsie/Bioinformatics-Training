library(Biostrings)
library(Hmisc)
library(plyr)
library(ggplot2)
library(seqinr)
options(scipen=20)

# Parse an attribute field from the GFF attributes field by name
getAttributeField <- function (x, field, attrsep = ";") {
  s = strsplit(x, split = attrsep, fixed = TRUE)
  sapply(s, function(atts) {
    a = strsplit(atts, split = "=", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv = a[[m]][2]
    }
    else {
      rv = as.character(NA)
    }
    return(rv)
  })
}

# Load a GFF file into a data.frame.  The attributes column is left alone
contigDict <- function(gffFile) {
  gff = scan(gffFile, what="raw", sep="\n")
  gff <- gff[grep("##sequence-header", gff)]
  gff <- ldply(gff, function(x) unlist(strsplit(x, split=" ")))
  seqid <- gff[,2]
  contig <- gff[,3]
  if(length(unlist(strsplit(contig[1],"ref......\\|"))) > 1){
    contig <- ldply(gff[,3], function(x) unlist(strsplit(x, split="ref......\\|")))
    contig <- contig[,length(contig)]
  }
  df <- data.frame(seqid, contig)
  colnames(df) <- c('seqid','contig')
  return(df)
}
getSeqid <- function(df, lookup){
  Contig = levels(factor(df$contig))
  df$seqid <- subset(lookup, contig==Contig)$seqid
  return(df)
}
getContig <- function(df, lookup){
  Seqid = levels(factor(df$seqid))
  df$contig <- subset(lookup, seqid==Seqid)$contig
  return(df)
}

gffRead <- function(gffFile, nrows = -1) {
  cat("Reading ", gffFile, sep="")
  gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
                   header=FALSE, comment.char="#", nrows = nrows,
                   colClasses=c("character", "character", "character", "integer",  
                                "integer",
                                "integer", "character", "character", "character"))
  colnames(gff) = c("seqid", "source", "type", "start", "end",
                    "score", "strand", "phase", "attributes")
  return(gff)
}


# Helper to load a PacBio modifications.gff file.  Parse some custom attribute fields relevant
# for modification detection.
readModificationsGff <- function(gffFile, nrows=-1)
{
  gffData <- gffRead(gffFile, nrows)
  gffData$coverage <- as.numeric(getAttributeField(gffData$attributes, 'coverage'))
  gffData$context <- as.character(getAttributeField(gffData$attributes, 'context'))
  gffData$CognateBase <- substr(gffData$context, 21, 21)
  gffData$IPDRatio <- as.numeric(getAttributeField(gffData$attributes, 'IPDRatio'))  
  gffData$attributes <- NULL
  ctgdict <- contigDict(gffFile)
  gffData <- ddply(gffData, "seqid", function(x) getContig(x, ctgdict))
  print(head(gffData))
  print("Summary of Coverage:")
  print(summary(gffData$coverage))
  return(gffData)
}
formatEquationText <- function(cutoffs){
  EqA <- paste("score > ",format(cutoffs$slope[1], digits=3), '*cvg AND cvg > ', cutoffs$x[1])
  EqC <- paste("score > ",format(cutoffs$slope[2], digits=3), '*cvg AND cvg > ', cutoffs$x[1])
  EqG <- paste("score > ",format(cutoffs$slope[3], digits=3), '*cvg AND cvg > ', cutoffs$x[1])
  EqT <- paste("score > ",format(cutoffs$slope[4], digits=3), '*cvg AND cvg > ', cutoffs$x[1])
  equation=c(EqA, EqC, EqG, EqT)
  return(equation)
}
aboveThreshholdHits <- function(hits, cutoffs){
  goodHits <- subset(hits, (CognateBase=="A" & coverage > cutoffs$x[cutoffs$CognateBase=="A"] & score > cutoffs$slope[cutoffs$CognateBase=="A"] * coverage + cutoffs$intercept[cutoffs$CognateBase=="A"]) | (CognateBase=="C" & coverage > cutoffs$x[cutoffs$CognateBase=="C"] & score > cutoffs$slope[cutoffs$CognateBase=="C"] * coverage + cutoffs$intercept[cutoffs$CognateBase=="C"]) | (CognateBase=="G" & coverage > cutoffs$x[cutoffs$CognateBase=="G"] & score > cutoffs$slope[cutoffs$CognateBase=="G"] * coverage + cutoffs$intercept[cutoffs$CognateBase=="G"]) | (CognateBase=="T" & coverage > cutoffs$x[cutoffs$CognateBase=="T"] & score > cutoffs$slope[cutoffs$CognateBase=="T"] * coverage + cutoffs$intercept[cutoffs$CognateBase=="T"]))
  return(goodHits)
}

# example set of motifs to label

# Sequence of motifs
# motifs <- c('CGACCAG', 'gantc', 'CGACNNNNNNNTRGG', 'CCYANNNNNNNGTCG')

# Modified position in motif
# motifSites <- c(6,2,3,4)

# Label each element of contextString as containing a motifs, or 'None'
# motifs is a list of motif strings, possibly including ambiguous bases
# motifSites is the index of the modified position in the motif
# labels is an optional label to apply - by default the string of the matching 
# motif will be used

labelContexts <- function(contextStrings, motifs, motifSites, labels=motifs, contextStringCenter=21){
  labelVect <- character(length(contextStrings))
  labelVect[1:length(contextStrings)] <- 'None'

  for (m in 1:length(motifs))
  {
    isHit <- as.vector(isMatchingAt(motifs[m], 
                                    DNAStringSet(as.character(contextStrings)), 
                                    at=contextStringCenter + 1 -motifSites[m],  
                                    fixed=F))
    labelVect[isHit] <- motifs[m]
  }
  factor(labelVect)
}

# Generate a data.frame with one row for each example of motif in the reference
# meth pos is the position of the m6A in the motif, offset=0 means to expect a 
# peak at the metylated position, offset=-5 is the normal position of the 
# secondary peak for m6A.  If there are Ns in the motif tb.end must be set 
# to the end of the first block of non-Ns. Also takes a lookup table from the header
# of the modifications.gff file to add the seqid to the df, so you can write a circos track.

genomeAnnotation <- function(refPath, motifs, positions, gff, offsets=c(0)){
  dnaSeq <- read.DNAStringSet(refPath)
  df <- data.frame()
  dnaSeqRc <- reverseComplement(dnaSeq)
  
  for (ref in 1:length(dnaSeq)){
    for(m in 1:length(motifs)){  
      motif <- motifs[m]
      methPos <- positions[m]
      
      maskFwd <- maskMotif(dnaSeq[[ref]], 'N')      
      maskRev <- maskMotif(dnaSeqRc[[ref]], 'N')
      
      matchFwd <- matchPattern(motif, maskFwd, fixed=F)
      matchRev <- matchPattern(motif, maskRev, fixed=F)  
      
      for (o in offsets){
        modName <- motif
        if(o > 0) modName = paste(motif, "p", o, sep='')
        if(o < 0) modName = paste(motif, "n", o, sep='')
        
        onTarget <- 'On'
        if(abs(o) > 0) onTarget <- 'Off'
        
        fwdDf <- data.frame(start=start(matchFwd))
        if(nrow(fwdDf) > 0){
          fwdDf$seqid <- ref
          fwdDf$contig <- names(dnaSeq)[[ref]]
          
          fwdDf$strand <- '+'
          fwdDf$start <- fwdDf$start + methPos - 1 + o
          fwdDf$motif <- modName
          fwdDf$onTarget <- onTarget  
          df <- rbind(df, fwdDf)
        }
        
       revDf <- data.frame(start=start(matchRev))        
       if(nrow(revDf) > 0){
          revDf$seqid <- ref
          revDf$contig <- names(dnaSeq)[[ref]]
          
          revDf$strand <- '-'
          revDf$start <- length(dnaSeq[[ref]]) + 1 - (revDf$start + methPos - 1 + o)
          revDf$motif <- modName
          revDf$onTarget <- onTarget
          df <- rbind(df, revDf)
        }
      }
    }
  }

  df <- df[,c('strand', 'start', 'motif',  'onTarget', 'contig','seqid')]
  contig <- df$contig
  if(length(unlist(strsplit(df$contig[1],"ref......\\|"))) > 1){
    contig <- ldply(df$contig, function(x) unlist(strsplit(x, split="ref......\\|")))
    contig <- contig[,length(contig)]
  }
  if(length(unlist(strsplit(df$contig[1]," "))) > 1){
    contig <- ldply(df$contig, function(x) unlist(strsplit(x, split=" ")))
    contig <- contig[,1]
  }
  df$contig <- contig
  ctgdict <- contigDict(gff)
  df <- ddply(df, 'contig', function(x) getSeqid(x,ctgdict))
  return(df)
}

# Write the contextStrings to a FASTA file
writeContextFasta <- function(dat, filename){
  stringSet <- DNAStringSet(dat$context)
  names(stringSet) <- paste(as.character(dat$start),dat$strand,sep='')
  write.XStringSet(stringSet, filename)
}
getUnmodifiedMotifKin <- function(mm, allkin){
  allkin <- allkin[c("tpl","strand","refId","ipdRatio","score","coverage")]
  names(allkin) <- c('start','strand','seqid','IPDRatio','score','coverage')
  allkin$strand[allkin$strand == 1] <- '-'
  allkin$strand[allkin$strand == 0] <- '+'
  mmsub <- mm[c("start","strand","seqid","contig","end","type","CognateBase","motif")]
  mga <- merge(mmsub, allkin, by=c('start','strand','seqid'))
  return(mga)
}
# Plot kinetics = c('ipdRatio','score','coverage')
plotKinetics <- function(dat,ref,statMode,plotRange,rangePerPlot){
  dat <- subset(dat,seqid=ref)
  tableNames <- names(dat)
  idx <- match(statMode,tableNames)
  posIdx <- match('tpl',tableNames)
  baseIdx <- match('base',tableNames)
  plotTitle <- statMode
  fontSize <- floor(800/rangePerPlot)
  baseline <- median(dat[,idx],na.rm=TRUE)
  if (statMode =='coverage') baseline <- 0
  yFwdStr <- 'Forward Strand'
  yRevStr <- 'Reverse Strand'
 
  posRange <- plotRange
  nPlots <- ceiling((posRange[2]-posRange[1])/rangePerPlot)
  
  for (n in 1:nPlots)  {
    sWin <- (n-1)*rangePerPlot+posRange[1]
    eWin <- sWin+rangePerPlot
    xbin <- 10
    subDat <- subset(dat,tpl >= sWin & tpl <= eWin)
    fwdDat <- subDat[subDat$strand=='+',c(posIdx,idx,baseIdx)]
    revDat <- subDat[subDat$strand=='-',c(posIdx,idx,baseIdx)]
    yMin <- -0.5
    yMax <- max(subDat[,idx],na.rm=TRUE)*1.05
    if (yMax < 5) yMax = 5
    plotDat <- data.frame(tpl = sWin:eWin, statFwd = baseline, statRev = baseline)
    plotDat$baseRev = "x"
    plotDat$baseFwd = "x"
    
    posMatch <- match(plotDat$tpl,fwdDat$tpl)
    fwdIdx <- posMatch[which(is.na(posMatch)==FALSE)]
    plotDat$statFwd[!is.na(posMatch)] = fwdDat[fwdIdx,2]
    plotDat$baseFwd[!is.na(posMatch)] = as.character(fwdDat$base[fwdIdx])
    
    posMatch <- match(plotDat$tpl,revDat$tpl)
    revIdx <- posMatch[which(is.na(posMatch)==FALSE)]
    plotDat$statRev[!is.na(posMatch)] = revDat[revIdx,2]
    plotDat$baseRev[!is.na(posMatch)] = as.character(revDat$base[revIdx])
    plotDat$baseline <- baseline
    
    p1 <- ggplot(data=plotDat, aes(xmin=tpl-0.4,xmax=tpl+0.4,ymin=baseline,ymax=statFwd))+geom_rect(fill='salmon',label='')+ylab(yFwdStr)+ylim(c(yMin,yMax))+ scale_x_continuous(breaks=seq(sWin,eWin,xbin))+opts(axis.text.x=theme_blank(),legend.position = "none",title=plotTitle,plot.margin= unit(c(0, 2, -.5, 0), "lines"))+coord_cartesian(xlim = c((sWin-1),(eWin+1)))
    pp1 <- p1 + geom_text(aes(x=tpl, y=0, vjust=1, label=toupper(baseFwd)), size=fontSize)+opts(axis.text.x = theme_blank(),axis.label.x = theme_blank())
    
    p2 <- ggplot(data=plotDat, aes(xmin=tpl-0.4,xmax=tpl+0.4,ymin=baseline,ymax=statRev))+geom_rect(fill='salmon',label='')+xlab('')+ylab(yRevStr)+ylim(c(yMax,yMin))+ scale_x_continuous(breaks=seq(sWin,eWin,xbin))+opts(legend.position = "none",plot.margin= unit(c(-0.5, 2, 0, 0), "lines"))+xlab('Ref Position')+coord_cartesian(xlim = c((sWin-1),(eWin+1)))
    pp2 <- p2 + geom_text(aes(x=tpl, y=0, vjust=0, label=toupper(baseRev)), size=fontSize)
    
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1)))
    print(pp1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
    print(pp2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  }
}

makeCircosTracks <- function(mm, circosPrefix){
  mm <- mm[,c("motif","seqid","start","end","strand","IPDRatio")]
  mm <- subset(mm, end!='NA')
  motifs <- levels(factor(mm$motif))

  for (word in motifs) {
    df <- subset(mm, motif==word)
    fwd <- subset(df, strand=="+")
    fwd$IPDRatio <- laply(fwd$IPDRatio, function(x) x=-x)
    fwdHits <- fwd[,c("seqid","start","end","IPDRatio")]

    rev <- subset(df, strand=="-")
    revHits <- rev[,c("seqid","start","end","IPDRatio")]

  #Note, if you see a read aligned to the + strand of the reference,
  #it was really sequenced from the - strand, and vice versa. 

    fwdFile = paste(circosPrefix, word,'.spikes.fwd.txt',sep="")
    revFile = paste(circosPrefix, word,'.spikes.rev.txt',sep="")
    cat("Writing forward tracks file.")
    cat("\n")
    write.table(file=fwdFile, revHits, col.names=FALSE, row.names=FALSE, quote=FALSE)
    cat("Writing reverse tracks file.")
    cat("\n")
    write.table(file=revFile, fwdHits, col.names=FALSE, row.names=FALSE, quote=FALSE)
  }
}
makeKaryotypeFile <- function(reference, gff, circosPrefix){
  sequences <- read.fasta(reference)
  contig <- names(sequences)
  end <- laply(sequences, length)
  karyotype <- data.frame(contig,end)
  ctgdict <- contigDict(gff)
  karyotype <- merge(karyotype, ctgdict, sort=FALSE)
  karyotype$start <- 0
  karyotype$label <- "chr - "
  karyotype <- karyotype[,c("label","seqid","contig","start","end")]
  write.table(file=paste(circosPrefix, "karyotype.txt", sep=""), karyotype, col.names=FALSE,row.names=FALSE, quote=FALSE)

  indent <- c("r0=0.95r,r1=1.0r", "")
  while (length(indent)<length(contig)){
    indent <- append(indent, "r0=0.95r,r1=1.0r")
    indent <- append(indent, "")
  }
  if (length(indent) > length(contig)) {
    indent <- indent[1:(length(indent)-1)]
  }
  highlights <- data.frame(seqid=karyotype$seqid, start=karyotype$start, end=karyotype$end, radius=indent)
  write.table(file=paste(circosPrefix, "chromosomes.txt", sep=""), highlights ,col.names=FALSE, row.names=FALSE, quote=FALSE)
}

makeCircosMotifHighlights <- function(mm, circosPrefix){
  motifs <- levels(factor(mm$motif))
  mm <- mm[,c("motif", "seqid", "start", "strand")]
  mm$end <- as.integer(mm$start) + nchar(mm$motif)
  for (word in motifs) {
    df <- subset(mm, motif==word)
    fwd <- subset(df, strand=="+")
    outFwd <- fwd[,c("seqid","start","end")]

    rev <- subset(df, strand=="-")
    outRev <- rev[,c("seqid","start","end")]

    fwdFile = paste(circosPrefix, word, '.motifPositions.fwd.txt', sep='')
    revFile = paste(circosPrefix, word, '.motifPositions.rev.txt', sep='')
    cat("Writing forward tracks file.")
    cat("\n")
    write.table(file=fwdFile,outRev,col.names=FALSE,row.names=FALSE,quote=FALSE)
    cat("Writing reverse tracks file.")
    cat("\n")
    write.table(file=revFile,outFwd,col.names=FALSE,row.names=FALSE,quote=FALSE)
  }
}
makeCircosCoverageTrack <- function(bedFile, gffFile, outFile){
  bed = scan(bedFile, what="raw", sep="\n")
  bed <- bed[2:length(bed)]
  bed <- ldply(bed, function(x) unlist(strsplit(x, split="\t")))
  bed$seqid <- getSeqid(bed, contigDict(gffFile))
  cvg <- data.frame(seqid=bed$seqid, start=bed$V2, end=bed$V3, coverage=bed$V5)

  cvgFile = paste(outFile,'.cvg.txt',sep='')
  cat("Writing coverage track file.")
  write.table(file=cvgFile, cvg, col.names=FALSE, row.names=FALSE, quote=FALSE)
}
makeCircosBaseline <- function(baselineMotif, refPath, gff, allkin){
  baselinePos <- c(1)
  baselineMotifAnnot <- genomeAnnotation(refPath, baselineMotif, baselinePos, gff) #find the positions of the baseline motif
  baselineMotifAnnot$end <- baselineMotifAnnot$start + nchar(baselineMotif)
  baselineMotifAnnot$type <- "baseline_motif"
  baselineMotifAnnot$CognateBase <- substr(baselineMotif[1], baselinePos,baselinePos)
  baselineKin <- getUnmodifiedMotifKin(baselineMotifAnnot, allkin)
  makeCircosTracks(baselineKin, circosPrefix)
}
mergeMotifAndBaselineData <- function(circosPrefix, motifs, baselineMotif){
  fwdBaselineFile = paste(circosPrefix, baselineMotif, '.spikes.fwd.txt', sep='')
  fwdBaseline = scan(fwdBaselineFile, what="raw", sep="\n")
  revBaselineFile = paste(circosPrefix, baselineMotif, '.spikes.rev.txt', sep='')
  revBaseline = scan(revBaselineFile, what="raw", sep="\n")

  
  for (word in motifs) {
    fwdMotifAndBaselineFile <- paste(circosPrefix, word, '.wBaseline.spikes.fwd.txt', sep='')
    revMotifAndBaselineFile <- paste(circosPrefix, word, '.wBaseline.spikes.rev.txt', sep='')

    fwdMotifFile = paste(circosPrefix, word, '.spikes.fwd.txt', sep='')
    fwdMotif <- scan(fwdMotifFile, what="raw", sep="\n")
    cat("Merging baseline and forward tracks file.")
    cat("\n")
    fwdMotifAndBaseline <- append(fwdMotif, fwdBaseline)
    write.table(file=fwdMotifAndBaselineFile, fwdMotifAndBaseline, col.names=FALSE, row.names=FALSE, quote=FALSE)

    revMotifFile = paste(circosPrefix, word, '.spikes.rev.txt', sep='')
    revMotif <- scan(revMotifFile, what="raw", sep="\n")
    cat("Merging baseline and reverse tracks file.")
    cat("\n")
    revMotifAndBaseline <- append(revMotif, revBaseline)
    write.table(file=revMotifAndBaselineFile, revMotifAndBaseline, col.names=FALSE, row.names=FALSE, quote=FALSE)
  }
}
makeCircosConfFiles <- function(motifs, circosPrefix, confModel, baselineMotif){ #fix to add replacement of karyotype and chr files w/circos prefix
  colors = c('vdorange', 'green', 'vdblue','vvdpurple')
  while (length(colors) < length(motifs)) colors <- append(colors,colors)
  colors <- colors[1:length(motifs)]
  colorMap <- data.frame("motif"= motifs, "color"= colors)
  
  for (m in motifs){
    thisMotif = subset(colorMap, motif==m)
    color = as.character(thisMotif$color)
    configInfo <- readLines(con=confModel)
    configInfo[which(configInfo=="karyotype = KARYOTYPE.txt")] <- paste("karyotype = ", circosPrefix,"karyotype.txt", sep="")
    configInfo[which(configInfo=="file = CHROMOSOMES.txt")] <- paste("file = ", circosPrefix,"chromosomes.txt", sep="")
    configInfo[which(configInfo=="fill_color=COLORNAME")] <- paste("fill_color =", color)
    configInfo[which(configInfo=="stroke_color=COLORNAME")] <- paste("stroke_color =",color)
    configInfo[which(configInfo=="file = WORD.circosPlot.png")] <- paste("file = ", circosPrefix, m, ".circosPlot.png", sep="")
    configInfo[which(configInfo=="file = WORD.motifPositions.fwd.txt")] <- paste("file = ", circosPrefix, m, ".motifPositions.fwd.txt", sep="")
    configInfo[which(configInfo=="file = WORD.motifPositions.rev.txt")] <- paste("file = ", circosPrefix, m, ".motifPositions.rev.txt", sep="")
    configInfo[which(configInfo=="file = WORD.spikes.fwd.txt")] <- paste("file = ", circosPrefix, m, ".wBaseline.spikes.fwd.txt", sep="")
    configInfo[which(configInfo=="file = WORD.spikes.rev.txt")] <- paste("file = ", circosPrefix, m, ".wBaseline.spikes.rev.txt", sep="")
    writeLines(configInfo, con = paste(circosPrefix, m, ".conf", sep=""), sep = "\n")
  }
}
splitChromosomes <- function(xsomelist, circosTrack){
  circos <- scan(circosTrack, what="raw", sep="\n")
  circos <- ldply(circos, function(x) unlist(strsplit(x, split=" "))) #seqid is always first column
  for (xsome in xsomelist) {
      karyotype = scan(xsome, what="raw", sep="\n")
      karyotype <- ldply(karyotype, function(x) unlist(strsplit(x, split=" ")))
      seqids <- karyotype$V4
      newTrack <- ddply(circos, "V1", function(df) if (levels(factor(df$V1)) %in% seqids) return(df))
      suffix <- strsplit(xsome, split="\\.")[[1]]
      suffix <- suffix[length(suffix)-1]
      newTrackFile <- paste(circosTrack, suffix, sep="")
      write.table(file=newTrackFile, newTrack, col.names=FALSE, row.names=FALSE, quote=FALSE)
    }
}

