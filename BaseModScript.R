# This script is for use with SMRTPortal 1.3.3 modifications.csv and modifications.gff (no motif annotation) output files.
# It reads gff3 compliant files.
# Meredith Ashby, Khai Luong, Jonas Korlach 10/2012

source("BaseModFunctions.R")

#### USER ENTERED VARIABLES HERE
# NOTE: The modification detection should be done on the reference you wish to use for circos, or the
# coordinates will not match up.

refPath <- 'ecoli_reference.fasta'
gff <- 'modifications.gff' 
csv <-  'modificationsSubset.csv' 
pdf_out_path <- 'ecoli.'
circosPrefix <- 'circos.'
confModel <- 'oneMotif.conf'

#### READ IN THE GFF FILE
hits <- readModificationsGff(gff)
workHits <- hits # That took forever! Let's keep a record of the whole gff file by making a copy we can play with.

#### IDENTIFY AND EXAMINE THE HIGH CONFIDENCE HITS

###PAGE 1
# Distribution of base modification scores broken out by base.
pdf(paste(pdf_out_path, 'ScoreDist.pdf', sep=""), width=11, height=8.5, onefile=T)
d <- qplot(score, colour=CognateBase, geom='freqpoly', data=subset(hits, CognateBase %in% c('A','T','G','C')), binwidth=1) +
  coord_cartesian(xlim=c(0,400)) +
  scale_x_continuous(breaks=seq(50,350,50)) +
  xlab("Modification QV")
show(d)
dev.off()

###PAGE 2
# Scatter plot of base modification score vs coverage on one set of axes, colored by base
pdf(paste(pdf_out_path, 'ScoreScatter.pdf', sep=""), width=11, height=8.5, onefile=T)
p <- qplot(coverage, score, colour=CognateBase, size=I(2.0), data=subset(hits, CognateBase %in% c('A','T','G','C') & score > 25))+
  xlab("Per Strand Coverage")+
  ylab("Modification QV")
show(p)
dev.off()


##################### OPTIONAL ADDITIONAL FILTERING OF HITS FROM THE GFF FILE
###PAGE 3
# Scatter plot, allowing USER TO CHOOSE THE CUTOFFS which will feed forward into downstream analysis, after iteratively adjusting the slope for each base

cutoffs <- data.frame(CognateBase=c('A','C','G','T'), slope=c(140/150,300/150,300/150,300/150), intercept=c(0,0,0,0), x=c(20,20,20,20))
cutoffs

# change the values for slope, intercept, and x for each base. X is a vertical line on the plot that
# specifies a minimum coverage, independent of score. 'Slope' and 'intercept' define an additional
# minimum score which depends on coverage. For a straight line cut independent of coverage, use
# slope = 0 and intercept = desired minimum score.

equations <- data.frame(CognateBase=c('A','C','G','T'), equation=formatEquationText(cutoffs)) # print formatting, for the plot.

###PAGE 4
# This plot shows how the cutoff choices made above will filter hits.
pdf(paste(pdf_out_path, 'HitFiltering.pdf', sep=""), width=11, height=8.5, onefile=T)
p <- qplot(coverage, score, colour=CognateBase, size=I(0.5), data=subset(hits, CognateBase %in% c('A','T','G','C')), xlim=c(0,100), ylim=c(0,200)) +
  facet_wrap(~CognateBase) +
  xlab("Per Strand Coverage")+
  ylab("Modification QV")+
  geom_abline(data=cutoffs, aes(slope=slope, intercept=intercept)) +
  geom_vline(data=cutoffs, aes(xintercept=x)) +
  geom_text(data=equations, aes(label=equation, x=30, y=5), vjust=0, hjust=0)
show(p)
dev.off()

### Select and order a subset of hits using the above cutoffs.
workHits <- aboveThreshholdHits(hits, cutoffs)
workHits <- workHits[order(workHits$score, decreasing=T),]

###### INPUT MOTIF INFORMATION FROM SMRTPORTAL OR MOTIFFINDER

# Beginning in v1.3.3, SMRTPortal performs automated modification identificaiton and motif finding on
# genomic positions with significantly altered sequencing kinetics.
# You can use the functions below to further characterize motifs of interest - modified or unmodified.
# Input the motifs you wish to examine further below to use R for additional analysis, or to generate
# circos plots.

motifs <- c('GATC','GCACNNNNNNGTT','AACNNNNNNGTGC')
positions <- c(2,3,2,1)
motifLabels <- labelContexts(workHits$context, motifs, positions)
table(motifLabels)  # reports how many motifs are found within the 41 base contexts in the gff file.

#### GENOME ANNOTATION

# Once you've found the motifs of interest among the hits, calculate how many instances of each motif within
# the entire genome were detected as modified. Merge all-genome motifs (identified using the reference fasta)
# with the results of the gff/contexts search to see what fraction of motifs in the genome are modified. 

genomeAnnotations <- genomeAnnotation(refPath, motifs, positions, gff) #find the positions of all motifs
mm <- merge(workHits, genomeAnnotations, all = TRUE)  
mm$motif[is.na(mm$motif)] <- 'no_motif'
mm$type[is.na(mm$type)] <- 'not_detected'
table(mm$type, mm$motif)

###PAGE 5
# Compare the score distributions for motif cognate-base positions to scores in the 'no_motif' set
pdf(paste(pdf_out_path, 'motif.scores.pdf', sep=""), width=11, height=8.5, onefile=T)
d=qplot(data=subset(mm, type!='not_detected'), score, colour=motif, geom='density', xlim=c(0,200))
show(d)
dev.off()

###PAGE 6
# Compare the IPD ratio distributions of the found motifs.
pdf(paste(pdf_out_path, 'motif.IPDratios.pdf', sep=""), width=11, height=8.5, onefile=T)
d=qplot(data=subset(mm, type!='Not_Detected'), IPDRatio, colour=motif, geom='density')
show(d)
dev.off()

### PLOT WHAT *ISN'T* IN THE GFF FILE
# Plot the score vs coverage information for all motif positions found in the reference.
# To do this, we have to get kinetic information on unmethylated motifs from the modifications.csv file.

allKinData <- csv.get(csv)
mga <- getUnmodifiedMotifKin(mm, allKinData)

# You can also subset the csv file to only the lines you need for this task + circos plotting
# by temporarily adding 'GA' to your motif list. This file will load faster if you need to revisit the analysis.
# allMotifKin <- subset(allKinData, tpl %in% mm$start)
# write.table(allMotifKin, file = 'my.modificationsSubset.csv', quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)

###PAGE 7
pdf(paste(pdf_out_path, 'motif.scoresVScvg.pdf', sep=""), width=11, height=8.5, onefile=T)
d <- qplot(data=mga, x=coverage, y=score, theme="bw", geom="point", colour=CognateBase)
d + ylab("Score") + xlab("Coverage") + facet_grid(motif~type)
show(d)
dev.off()

###PAGE 8
pdf(paste(pdf_out_path, 'motif.HitvsMissCvg.pdf', sep=""), width=11, height=8.5, onefile=T)
d <- qplot(data=mga, x=coverage, y=IPDRatio, theme="bw", geom="point", colour=type, ylim=c(0,15))
d + ylab("IPD Ratio") + xlab("Coverage") + facet_grid(motif~.) 
show(d)
dev.off()

# R has very useful tools for making nearly any type of plot.  The 'names' command reveals
# all the information held in the dataframe 'mga' that you may be interested in further examining.
# Check out the ggplot online manual for documentation and additional plotting options:
# docs.ggplot2.org/current/
names(mga)

#### MAKING CIRCOS PLOTS TO VISUALIZE YOUR HITS AND MOTIFS

# This sequence of commands will create the input files and a conf file needed to generate the red
# spiked plots seen in PacBio basemod publications. To actually make the plot you will need to install
# and run circos.

# Create a karyotype track using information from the reference fasta and the gff header.
makeKaryotypeFile(refPath, gff, circosPrefix)   #This preserves the chromosome order found in the the reference

# Make circos highlights tracks showing the genome positions of each motif,
# both modified and unmodified.
makeCircosMotifHighlights(mm, circosPrefix)

# Make circos basemod score tracks (the red spikes) for each set of motif-associated hits.
makeCircosTracks(mm, circosPrefix)

# Take data points from the base modification csv file to draw a baseline in the circos plot. Choose a 'motif'
# that will be common enough to define a baseline but not so common as to take forever to run. For a plasmid,
# one might choose 'A'; for a bacterial genome 'GA' would be appropriate.
baselineMotif <- c("GA")
makeCircosBaseline(baselineMotif, refPath, gff, allKinData)
mergeMotifAndBaselineData(circosPrefix, motifs, baselineMotif)

# From here, you have all the input files needed to assemble the circos config file. 
makeCircosConfFiles(motifs, circosPrefix, confModel, baselineMotif)

# To generate the plot, from the bash prompt (not in R!):
# $ circos -conf (word).conf

