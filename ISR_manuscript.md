ISR Manuscript Markdown
================
Troy Ellison
4/13/2019

Analysis Script for comparing pipelines & ISR-type Analysis Updated: 05/2018

#### Experiment Notes

Notes: Here we will compare the 3 pipelines and do ISR-type analysis

#### Environment Notes

R version 3.4.0 (2017-04-21) -- "You Stupid Darkness" Platform: x86\_64-apple-darwin15.6.0 (64-bit) RStudio Version 1.0.143

 

### Load Required Libraries

``` r
library(plotrix)
library(GUniFrac) #Rarefy
library(gplots) #heatmap
library(reshape2) #decostand
library(data.table) #setDT
library(Biostrings) #readDNAstringset
```

Set Seed = 12345

#### Define Custom Function

``` r
'%!in%' <- function(x,y)!('%in%'(x,y))
```

#### Run AD\_settings Script for Loading Functions

``` r
source("/Volumes/GriffenLeysLab/Troy/AD_16S/Markdown/ISR/AD_settings.R")
```

    ## Loading required package: magrittr

    ## 
    ## Attaching package: 'ggpubr'

    ## The following object is masked from 'package:ape':
    ## 
    ##     rotate

 

### Load R Objects

``` r
ISR_long_R.atab <- loadER("/Volumes/GriffenLeysLab/Ran/ISR_2018/ISR/.RData","ISR_long_R.atab") # Load ASVs from DADA2 processing
meta <- loadER("/Volumes/GriffenLeysLab/Ran/ISR_2018/ISR/.RData","meta") # Load meta file
meta$level <- as.factor("level")

ISR_long_R_seqs <-loadER("/Volumes/GriffenLeysLab/Ran/ISR_2018/ISR/.RData","ISR_long_R_seqs") # Load sequences
```

### Load 16S Samples and Sampling Times

``` r
# Input actual sampling time differences
time_diff <- read.table(file="/Volumes/GriffenLeysLab/Ran/ISR_2018/ISR/time_diff.txt", header = T) # input time differences from clinical sampling schedule (based on meta file :: sampling dates)
time_diff$mo <- signif(time_diff$mo, 2)

# Loading 16S DADA2 Processed Samples:
V12_long_R.atab <- loadER("/Volumes/GriffenLeysLab/Ran/ISR_2018/v12_dada/.RData","V12_long_R.atab")

# Loading 16S CORE matched samples:
V13_long_R.core <- loadER("/Volumes/GriffenLeysLab/Troy/ISR_manuscript_markdown/v13/.RData", "V13_long_R.core")
```

 

#### Random Accuracy Rates

``` r
ISR_random_forest(V13.plots$rar_prab)

ISR_random_forest(V12.plots$rar_prab)

ISR_random_forest(ISR.plots$rar_prab)
```

 

### Plotting Barplots of richness and Boxplot of mean centroid distances

#### Statistical Calculations

``` r
# Calculating betadispersion
mod.v13 <- betadisper(V13.plots$prab_bcdist, factor(meta$level))
mod.v12 <- betadisper(V12.plots$prab_bcdist, factor(meta$level))
mod.ISR <- betadisper(ISR.plots$prab_bcdist, factor(meta$level))

# Statistical Testing
wilcox.test(mod.v13$distances, mod.v12$distances) 
wilcox.test(mod.v12$distances, mod.ISR$distances) 
wilcox.test(mod.v13$distances, mod.ISR$distances) 

# Calculating Shared OTUs/ASVs with custom function otu_stats:
V13_sub <- otu_stats(V13.plots$rar_prab)
V12_sub <- otu_stats(V12.plots$rar_prab)
ISR_sub <- otu_stats(ISR.plots$rar_prab)
```

#### Combining Stats

``` r
pipeline_stats <- rbind.data.frame(V13_sub$stats,V12_sub$stats, ISR_sub$stats)
colnames(pipeline_stats) <- c("Unique.OTUs.ASVs","Shared.OTUs.ASVs","total")
pipeline_stats$Pipeline <- c("16S OTUs","16S ASVs ","ISR ASVs")

# Set par for plots
par(mar = c(13,5,4,4), mgp=c(3.5, 0.5, 0)) # barplot
par(mfrow=c(1,4))
```

 

### Generate Barplots

 

``` r
par(mar = c(8,6,3,1), mgp=c(4, 0.5, 0.5)) # barplot
barplot(pipeline_stats$Unique.OTUs.ASVs, names.arg = pipeline_stats$Pipeline,
        col=c("#424E5C","#365BB7","#5B1A8E"),las=2, 
        cex.lab=1.75,cex.names = 1.5, cex.axis = 1.2,
        ylab = "% Variants Unique to 1 Subject")
```

 

``` r
par(mar = c(8,6,3,1), mgp=c(4, 0.5, 0.5)) # barplot
barplot(pipeline_stats$Shared.OTUs.ASVs, names.arg = pipeline_stats$Pipeline,
        col=c("#424E5C","#365BB7","#5B1A8E"),las=2, 
        cex.lab=1.75,cex.names = 1.5,cex.axis = 1.2,
        ylab = "% Variants Shared by 5 Subjects")
```

 

### Generate Boxplots

``` r
par(mar = c(7,6,1,1), mgp=c(4, 0.5, 0.5)) # barplot
boxplot(mod.v13$distances, mod.v12$distances, mod.ISR$distances,
        ylim=c(0.2,0.7), col=c("#424E5C","#365BB7","#5B1A8E"),
        xaxt="n",
        boxwex = 0.25,
        cex.lab=1.75, cex.axis=1.2,cex.names = 1.5,
        frame=F,
        #main="Comparison of Pipelines",
        ylab="Centroid Distances")
axis(side = 1, at = c(1,2,3), tck=-0.01,
     labels = c("16S Ref.","16S HR ","ISR HR"), 
     cex.axis=1.5,
     las=2)


# 16S Sp vs 16S dada2
# create the y coordinate of the line
y <- 0.68
# set an offset for tick lengths
offset <- 0.005
# draw first horizontal line
lines(c(1,2),c(y, y))
# draw ticks
lines(c(1,1),c(y, y-offset))
lines(c(2,2),c(y, y-offset))
# draw asterics
text(1.5,y+offset,"****")

# 16S DADA2 vs ISR
# create the y coordinate of the line
y <- 0.69
# draw first horizontal line
lines(c(2,3),c(y, y))
# draw ticks
lines(c(2,2),c(y, y-offset))
lines(c(3,3),c(y, y-offset))
# draw asterics
text(2.5,y+offset,"**")

# 16S Species vs ISR
# create the y coordinate of the line
y <- 0.70
# draw first horizontal line
lines(c(1,3),c(y, y))
# draw ticks
lines(c(1,1),c(y, y-offset))
lines(c(3,3),c(y, y-offset))
# draw asterics
text(2,y+offset,"****")
```

 

### Plot NMDS, Dendogram, and Stability Figures

``` r
V13.plots <- ISR_dist_plots(V13_long_R.core, meta, "16S V1-V3", 1.4,0.8,0.7,0.8)
```

``` r
V12.plots <- ISR_dist_plots(V12_long_R.atab, meta,"16S DADA2", 1.6,1.8,1.8,1.9)
```

``` r
ISR.plots <- ISR_dist_plots(ISR_long_R.atab, meta,  "ISR", 1.8,1,0.9,0.9)
```

 

### Subspecies Level Analysis

#### Step1: Select ISR ASVs present in &gt;1 sample

``` r
# Select ISR ASVs present in >1 sample:
ISR_long_R.atab.prab <- data.frame((ISR_long_R.atab > 0)*1)
ISR_long_R.atab.prab.sel <- ISR_long_R.atab.prab[,colSums(ISR_long_R.atab.prab) > 1] # 1000 sequences selected
ISR_blast_accnos <- colnames(ISR_long_R.atab.prab.sel)

ISR_long_R.atab.sel <- ISR_long_R.atab[,colnames(ISR_long_R.atab) %in% colnames(ISR_long_R.atab.prab.sel)] # 30 x 1000, raw counts
ISR_long_R.atab.scrap <- ISR_long_R.atab[,colnames(ISR_long_R.atab) %!in% colnames(ISR_long_R.atab.prab.sel)]

ISR_long_R.fasta = readDNAStringSet("/Volumes/GriffenLeysLab/Ran/ISR_2018/ISR/ISR_long_R.fa")
seq_name = names(ISR_long_R.fasta)
sequence = paste(ISR_long_R.fasta)
ISR_long_R.seqs.df <- data.frame(seq_name, sequence)

# Selecting rarefied sequences:

ISR_seqs.fa <- ISR_long_R.fasta[names(ISR_long_R.fasta) %in% ISR_blast_accnos]
writeXStringSet(ISR_seqs.fa, "ISR_blast.seqs.fa")
```

#### Step2: Blast selected ISR ASV sequences with ISR Database

``` r
# BLAST with latest ISR Database, using shell script, from UNIX command line:
# /scripts/./ISR_blast_nosql.sh ISR_blast ISR_db_jan2018.fasta  12
# Total reads:1000  Matches:999 No Hits:  1

# ISR Strain Analysis

# Load BLAST results: 

ISR_blastres <- read.table(file="/Volumes/GriffenLeysLab/Ran/ISR_project/ISR_manuscript_04_2018/Final_analysis/species_analysis/ISR_blast_blast_top2.txt", header=F, sep="\t") 
# Load Taxonomy file
ISR_db_tax <- read.table(file="/Volumes/GriffenLeysLab/Ran/scripts_dbs/blast_db/ISR_db_jan2018_taxonomy.txt", header=T, sep="\t")

# Processing Blast results:
row.names(ISR_blastres) <- ISR_blastres$V1
ISR_blastres$V1 <- NULL
colnames(ISR_blastres) <- c("db1","score1","m_size1","db2","score2","m_size2")

# Assign species level taxonomy
ISR_blastres$otu1 <- ISR_db_tax$CORE_otu[match(ISR_blastres$db1,ISR_db_tax$isr_db_id)]
ISR_blastres$otu2 <- ISR_db_tax$CORE_otu[match(ISR_blastres$db2,ISR_db_tax$isr_db_id)]

# Filter Blast matches with score > 0.1
ISR_blast_matches <- ISR_blastres[ISR_blastres$score1 > 0.1,]

# Verifying ISR database BLAST validity
same_blast_scores <- ISR_blast_matches[(ISR_blast_matches$score1 == ISR_blast_matches$score2),] 
blast_mismatches <- same_blast_scores[same_blast_scores$otu1 != same_blast_scores$otu2,]
```

 

### Selected 90% as cutoff for BLAST Matches

``` r
# Subsetting BLAST matches for scores >= 90%
ISR_matches_over90 <- ISR_blast_matches[ISR_blast_matches$score1 >= 0.90,] 
ISR_matches_over90 <- ISR_matches_over90[,c(1:3, 7)]

# Checking unmatched sequences
ISR_nonmatches.fa <- ISR_seqs.fa[names(ISR_seqs.fa) %!in% rownames(ISR_matches_over90)] 
writeXStringSet(ISR_nonmatches.fa, "ISR_nonmatches.fasta")
# Will be matched against larger NCBI databases

# Subsetting ISR pivot for BLAST matches at 90%
ISR_long_R.matches <- ISR_long_R.atab[,colnames(ISR_long_R.atab) %in% row.names(ISR_matches_over90)] 

# NMDS with only matches over 90
ISR_blast_matches.plots <- ISR_dist_plots(ISR_long_R.matches,meta,"ISR BLAST Matches",2.9,1.9,1.9,1.9)
```

``` r
# NMDS and Hclust looks similar to the one made with all ASVs

write.table(ISR_long_R.matches, file="ISR_long_R.over90.txt", sep="\t", quote=F, col.names = NA)

# Making CORE OTU based Counts table
ISR_long_R.matches.melt <- melt(t(ISR_long_R.matches))
ISR_long_R.matches.melt$otu <- ISR_matches_over90$otu1[match(ISR_long_R.matches.melt$Var1,row.names(ISR_matches_over90))]
ISR_long_R.338.counts <- dcast(ISR_long_R.matches.melt, Var2 ~ otu, fun.aggregate= sum) 
row.names(ISR_long_R.338.counts) <- ISR_long_R.338.counts$Var2
ISR_long_R.338.counts$Var2 <- NULL


# Making ISR strain table
ISR_strain_table <- as.data.frame(sort(table(droplevels(ISR_matches_over90$otu1)), decreasing=T))
colnames(ISR_strain_table) <- c("CORE OTU","No. of Strains")

ISR_strain_table$Sequences <- colSums(ISR_long_R.338.counts)[match(ISR_strain_table$`CORE OTU`,colnames(ISR_long_R.338.counts))]
write.table(ISR_strain_table, file="ISR_strains_table.txt", sep="\t", col.names = NA, quote = F)

# Sorting ISR strain table
ISR_strain_table_sorted <- ISR_strain_table[order(-ISR_strain_table$Sequences),]
strain_table <- ISR_strain_table_sorted[c(1:15),]
mean(strain_table$`No. of Strains`) # 15.26
write.table(strain_table, file="strain_table_top15.txt", sep="\t", col.names = NA, quote = F)
```

 

### ISR Type Plots

#### Haemophilus parainfluenzae

``` r
# Define color palette (source: stack overflow)
c25 <- c("dodgerblue2","#E31A1C", # red
         "green4",
         "#6A3D9A", # purple
         "#FF7F00", # orange
         "black","gold1",
         "skyblue2","#FB9A99", # lt pink
         "palegreen2",
         "#CAB2D6", # lt purple
         "#FDBF6F", # lt orange
         "gray70", "khaki2",
         "maroon","orchid1","deeppink1","blue1","steelblue4",
         "darkturquoise","green1","yellow4","yellow3",
         "darkorange4","brown")

# ISR-type plots for 3 most diverse species-groups
hp.plots <- strain_plots("Haemophilus parainfluenzae", "HP", 0.2)
```

 

#### Granulicatella adiacens

``` r
ga.plots <- strain_plots("Granulicatella adiacens", "GA", 0.07)
```

 

#### Streptococcus mitis pneumoniae infantis oralis

``` r
sm.plots <- strain_plots("Streptococcus mitis pneumoniae infantis oralis", "SM", 0.7)
```

``` r
# Plot key for ISR-types barplot
x <- 1:25
my_labels <- x

# default par
par(mar=c(0, 0, 0, 0))
plot(c(0,1),type="n", axes=F, xlab="", ylab="")
legend("center", legend=my_labels, pch=22, pt.bg=c25, cex=1, box.col = "white", pt.cex=3, ncol=10)
```

 

Functions
=========

### ISR\_Dist\_Plots

Function for plotting distance-based beta diversity analysis results This function plots NMDS, H-C dendograms, and stability plots Needs function for stability plot to be loaded first

``` r
# Define function
ISR_dist_plots <- function(counts,meta,name,xlm1, xlm2, ylm1, ylm2) {
  
  #set seed
  set.seed(999)
  
  # Load required libraries
  library(vegan)
  library(GUniFrac)
  library(ecodist)
  library(pvclust)
  library(dendextend)
  library(data.table)
  
  # Subset Counts table based on meta
  counts.sub <- counts[row.names(counts) %in% meta$sample, ]
  
  # Drop empty columns
  counts.sub <- counts.sub[,colSums(counts.sub) > 0]
  
  # Rarefy Counts for Pr/Ab
  counts.sub.rfy <- rrarefy(counts.sub, min(rowSums(counts.sub)))
  counts.sub.rar <- data.frame(counts.sub.rfy)
  #counts.sub.rar <- data.frame(counts.sub.rfy$otu.tab.rff)
  
  # Convert to Presence/Abscence
  counts.sub.prab <- data.frame((counts.sub.rar > 0)*1)
  
  # Compute NMDS
  set.seed(99999);counts.sub.prab.mds <- metaMDS(counts.sub.prab, wascores=F, autotransform = F, trymax=1000)
  
  # NMDS on relative abundance, output at PDF"
  pdf(paste(name,"_Prab_NMDS.pdf",sep=""), width = 8, height = 6)
  
  # Set default par for MDS
  par(mar=c(5.1, 4.1, 4.1, 4.1), mgp=c(3, 1, 0), las=0)
  par(bg = "white")
  # Plot NMDS
  plot(counts.sub.prab.mds$points, 
       pch=21,
       bg=c("red","blue","green","purple","orange")[meta$subject], 
       main=paste("Presence-Abscence ",name), 
       xlim = c(-xlm1,xlm2), ylim = c(-ylm1,ylm2), 
       cex=2.5, yaxt="n", xaxt="n") 
  
  # Add ellipse
  prab.ord <- ordiellipse(counts.sub.prab.mds$points, 
                          meta$subject, conf =  0.95, 
                          col=c("red","blue","green","purple","orange"))
  
  # Plot axes
  par(bg = "white")
  par(las=2)
  axis(2,cex.axis=1.5)
  par(las=1)
  axis(1,cex.axis=1.5)
  
  dev.off()
  
  # Set default par for MDS
  par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
  par(bg = "white")
  # Plot NMDS
  plot(counts.sub.prab.mds$points, 
       pch=21,
       bg=c("red","blue","green","purple","orange")[meta$subject], 
       main=paste("Presence-Abscence ",name), 
       xlim = c(-xlm1,xlm2), ylim = c(-ylm1,ylm2), 
       cex=2.5, yaxt="n", xaxt="n") 
  
  # Add ellipse
  prab.ord <- ordiellipse(counts.sub.prab.mds$points, 
                          meta$subject, conf =  0.95, 
                          col=c("red","blue","green","purple","orange"))
  
  # Plot axes
  par(bg = "white")
  par(las=2)
  axis(2,cex.axis=1.5)
  par(las=1)
  axis(1,cex.axis=1.5)
  
  # Legend plot
  pdf(paste(name,"_sub_legend.pdf",sep=""), width = 4, height = 6)
  # Empty Plot for legend
  plot(c(0,1),type="n", axes=F, xlab="", ylab="")
  
  # Add legends
  legend("bottomleft", pch = 21, pt.bg = c("red","blue","green","purple","orange"), 
         legend=levels(meta$subject), title="subjects", cex=2, text.width=0.4, box.col = "white")
  
  dev.off()
  
  # Calculate Bray-Curtis Distances
  counts.sub.prab.bcdist <- vegdist(counts.sub.prab)
  
  # Stability Plot
  ISR_stability_base(counts.sub.prab.bcdist, paste(name, "_PrAb",sep=""))
  
  # Calculate Hierarchical Clustering P-values
  counts.sub.prab.pvclust <- pvclust(t(counts.sub.prab),
                                     method.dist = "cor",
                                     method.hclust="complete", nboot=1000)
  
  # Convert to class 'dendogram'
  counts.sub.prab.dend <- as.dendrogram(counts.sub.prab.pvclust)
  
  # Add relevant color labels
  labels_colors(counts.sub.prab.dend) <- c("red","blue","green3","purple","darkorange")[meta$subject][order.dendrogram(counts.sub.prab.dend)]
  
  # Set hang height
  counts.sub.prab.dend <- hang.dendrogram(counts.sub.prab.dend, hang_height=0.1)
  
  # Plot Dendogram
  pdf(paste(name,"_dend.pdf",sep=""), width = 8, height = 6)

  plot(counts.sub.prab.dend, 
       main=paste("Presence-Abscence ",name),
       ylim = c(0,1.2), cex=1, yaxt="n")
  par(las=2)
  axis(2,cex.axis=1)
  #text(counts.sub.prab.pvclust, col="black", cex=0.7)
  pvrect(counts.sub.prab.pvclust, alpha = 0.95, pv = "bp", border="red4", lwd=1)
  
  dev.off() 
  
  plot(counts.sub.prab.dend, 
       main=paste("Presence-Abscence ",name),
       ylim = c(0,1.2), cex=1, yaxt="n")
  par(las=2)
  axis(2,cex.axis=1)
  #text(counts.sub.prab.pvclust, col="black", cex=0.7)
  pvrect(counts.sub.prab.pvclust, alpha = 0.95, pv = "bp", border="red4", lwd=1)
  
  # Output 
  list(
    counts = counts.sub,
    counts_rar = counts.sub.rar,
    rar_prab = counts.sub.prab,
    prab_bcdist = counts.sub.prab.bcdist,
    prab_mds = counts.sub.prab.mds,
    prab_dend = counts.sub.prab.dend,
    pvclust = counts.sub.prab.pvclust,
    prab_ord = prab.ord,
    prab_perm = adonis(formula = counts.sub.prab ~ meta$subject, strata = meta$run)
  )
}
```

 

### ISR\_Stability\_Base

Function to plot distance from base timepoint

``` r
ISR_stability_base <- function(counts.sub.prab.bcdist, name){
  
  #set seed
  set.seed(12345)
  
  #Load libraries
  library(lme4)
  
  # Melt the distance matrix to extract sample specific values
  prab_dist <- melt.dist(counts.sub.prab.bcdist)
  
  # Compute specific fields
  prab_dist$sub1 <- sapply(strsplit(as.character(prab_dist$Var1),"_"), '[', 1)
  prab_dist$tp1 <- sapply(strsplit(as.character(prab_dist$Var1),"_"), '[', 2)
  prab_dist$sub2 <- sapply(strsplit(as.character(prab_dist$Var2),"_"), '[', 1)
  prab_dist$tp2 <- sapply(strsplit(as.character(prab_dist$Var2),"_"), '[', 2)
  prab_dist$tpdiff <- abs(as.numeric(sapply(strsplit(as.character(prab_dist$tp2),""), 
                                            '[', 2))-as.numeric(sapply(strsplit(as.character(prab_dist$tp1),""), '[', 2)))
  
  # Adding same/diff markers
  prab_dist$subdiff <- "same"
  prab_dist[prab_dist$sub1!=prab_dist$sub2,]$subdiff <- "diff"
  
  # Adding correct time duration information
  prab_dist$tp1m <- time_diff$mo[match(prab_dist$tp1,time_diff$tp)]
  prab_dist$tp2m <- time_diff$mo[match(prab_dist$tp2,time_diff$tp)]
  
  # Computing absolute values for time differences
  prab_dist$td <- abs(prab_dist$tp2m - prab_dist$tp1m)
  
  # Divide into same and diff subject tables
  prab_dist.same <- prab_dist[prab_dist$subdiff=="same",]
  prab_dist.diff <- prab_dist[prab_dist$subdiff=="diff",]
  
  # Divergence from baseline
  # Same subject
  prab_dist.same.base <- prab_dist.same[prab_dist.same$tp1=="T1",]
  
  # Different subject
  # Different time points only
  prab_dist.diff.tpdiff <- prab_dist.diff[prab_dist.diff$tp1!=prab_dist.diff$tp2,]
  
  # Difference from base only
  prab_dist.diff.tpdiff.base <- prab_dist.diff.tpdiff[abs(prab_dist.diff.tpdiff$td) %in% time_diff$mo,]
  
  prab_dist.diff.tpdiff.base <- prab_dist.diff.tpdiff.base[order(prab_dist.diff.tpdiff.base$tp1),]
  prab_dist.diff.tpdiff.base <- prab_dist.diff.tpdiff.base[order(prab_dist.diff.tpdiff.base$tp2),]
  
  prab_dist.same.base <<-  prab_dist.same.base
  
  # Calculate linear model
  
  # Linear Model
  # Diff subjects
  prab_dist.dif.lm <- lm(prab_dist.diff$value ~ abs(prab_dist.diff$td))
  
  # Same subjects
  #prab_dist.same.base.lm <- lm(prab_dist.same.base$value ~ prab_dist.same.base$tp2m)
  
  # Linear Mixed effects model using 
  prab_dist.same.base.lm <- lmer(value ~ tp2m + ( 1 | sub1), data=prab_dist.same.base)
  
  # Calculate Summary stats
  prab_dist.same.base.lm.sum <- summary(prab_dist.same.base.lm)
  
  prab_dist.same.base.lm.sum <<- prab_dist.same.base.lm.sum
  
  # Calculate P value with ANOVA test:
  
  # anova(prab_dist.diff.tpdiff.base.lm) # P=0.05847
  prab_dist.same.base.lm.an <- anova(prab_dist.same.base.lm)
  prab_dist.same.base.lm.an <<- prab_dist.same.base.lm.an
  print(prab_dist.same.base.lm.an)
  
  # Plot stability as PDF output
  pdf(paste(name,"_Stability.pdf",sep=""), width = 8, height = 6)
  plot(prab_dist.same.base$value ~ jitter(prab_dist.same.base$tp2m,0.4),
       xlim=c(0,12), ylim=c(0,1),
       bg=c("red","blue","green","purple","orange")[factor(prab_dist.same.base$sub1)],
       pch=21, cex=1.4, cex.lab=1,
       xlab="Time in Months", ylab="Dissimilarity Index",
       main=paste(name,"Stability From Baseline"))
  #axis(1, at =c(3.1,4.6,8.1,9.0,10.7))
  
  # Plot trend lines
  abline(prab_dist.same.base.lm.sum$coefficients[1],  prab_dist.same.base.lm.sum$coefficients[2], col="red", lty=2, lwd=2)
  abline(prab_dist.dif.lm, lty=2, lwd=2) # for all inter-sub distances
  
  # Print slope on plot
  text(0, y=(signif(prab_dist.same.base.lm.sum$coefficients[1],2)+0.1),
       paste("Slope:",signif(prab_dist.same.base.lm.sum$coefficients[2],2)), pos=4, col="red")
  
  text(0, y=(signif(prab_dist.dif.lm$coefficients[1],2))+0.05, "Inter-Subject Distances", pos=4)
  
  dev.off()
  
  plot(prab_dist.same.base$value ~ jitter(prab_dist.same.base$tp2m,0.4),
       xlim=c(0,12), ylim=c(0,1),
       bg=c("red","blue","green","purple","orange")[factor(prab_dist.same.base$sub1)],
       pch=21, cex=1.4, cex.lab=1,
       xlab="Time in Months", ylab="Dissimilarity Index",
       main=paste(name,"Stability From Baseline"))
  #axis(1, at =c(3.1,4.6,8.1,9.0,10.7))
  
  # Plot trend lines
  abline(prab_dist.same.base.lm.sum$coefficients[1],  prab_dist.same.base.lm.sum$coefficients[2], col="red", lty=2, lwd=2)
  abline(prab_dist.dif.lm, lty=2, lwd=2) # for all inter-sub distances
  
  # Print slope on plot
  text(0, y=(signif(prab_dist.same.base.lm.sum$coefficients[1],2)+0.1),
       paste("Slope:",signif(prab_dist.same.base.lm.sum$coefficients[2],2)), pos=4, col="red")
  
  text(0, y=(signif(prab_dist.dif.lm$coefficients[1],2))+0.05, "Inter-Subject Distances", pos=4)
  
}  
```

 

#### ISR\_Random\_Forest

Function to calculate accuracy from Random Forest Classifier

``` r
ISR_random_forest <- function(prabs){
  
  # load required libraries
  library(randomForest)
  library(caret)
  
  # set seed
  set.seed(12345)
  
  # Clean up community membership matrix
  prabs <- data.frame(prabs)
  prabs <- prabs[,colSums(prabs) > 0]
  
  #Merge community membership matrix with subject info from meta file
  prabs.merged <- cbind.data.frame(prabs,meta$subject)
  prabs.merged$subject <- prabs.merged$`meta$subject`
  prabs.merged$`meta$subject` <- NULL
  
  numRows = nrow(prabs.merged)
  gp <- runif(numRows)
  mixed_data <- prabs.merged[order(gp),] #Mix up the rows of data frame
  
  m <- c()
  for(n in seq(from=1, to=30, by=1)) {
    set.seed(12345)
    data_train <- mixed_data[-n, -ncol(mixed_data)]
    
    data_test  <- mixed_data[n, -ncol(mixed_data)]
    
    labels_train <- mixed_data[-n,ncol(mixed_data)]
    labels_target <- mixed_data[n,ncol(mixed_data)]
    
    #Using Random Forest classifier
    rf <- randomForest(labels_train ~ ., data=data_train, ntree=100,proximity=T)
    rf_tab <- table(predict(rf, data_test), labels_target)
    conf <- confusionMatrix(rf_tab)
    c <- conf$overall[1]
    m <- c(m, c)
  }
  
  print(mean(m)*100)
  
}
```

 

#### OTU\_Stats

Function to generate statistics for ASV/OTU counts

``` r
otu_stats <- function(rar_prab) {
    
  # Combine presence/absence for all time points for each subject
  asv_prab.sub <- setDT(rar_prab)[, as.list(colSums(.SD)), 
                                 by = gl(ceiling(nrow(rar_prab)/6), 6, 
                                         nrow(rar_prab), labels = levels(meta$subject))]
  asv_prab.sub <- as.data.frame(asv_prab.sub)
  row.names(asv_prab.sub) <- asv_prab.sub$gl
  asv_prab.sub$gl <- NULL
  
  # Convert to presence/absence by subject
  asv_prab.sub.prab <- data.frame((asv_prab.sub > 0)*1)
  
  # Remove OTU/ASVs not present in any subject
  asv_prab.sub.prab <- asv_prab.sub.prab[,colSums(asv_prab.sub.prab) > 0]
  
  # Which OTU/ASVs are present in only 1 subject:
  asv_prab.sub.prab.unique <- asv_prab.sub.prab[,colSums(asv_prab.sub.prab) == 1]
  asv_prab.unique.stat <- signif(ncol(asv_prab.sub.prab.unique)/ncol(asv_prab.sub.prab)*100,2)
  
  # Which OTU/ASVs are present in all 5 subjects:
  asv_prab.sub.prab.all <- asv_prab.sub.prab[,colSums(asv_prab.sub.prab) == 5]
  asv_prab.shared.stat <- signif(ncol(asv_prab.sub.prab.all)/ncol(asv_prab.sub.prab)*100,2)
  
  # Output 
  list(
    stats = c(asv_prab.unique.stat,asv_prab.shared.stat, ncol(asv_prab.sub.prab)),
    sub_prab = asv_prab.sub.prab
  )
   
}
```

 

#### Strain\_Plots

Function to generate ISR-type figures

``` r
strain_plots <- function(otu, code, yl) {
  
  # Set seed
  set.seed(12345)
  
  # Retreive SVs for the species:
  code.strains <- ISR_long_R.matches[,colnames(ISR_long_R.matches) %in% row.names(ISR_matches_over90[ISR_matches_over90$otu1==otu,])]
  
  # Add "others" sequences column
  code.strains$others  <- rowSums(ISR_long_R.atab) - rowSums(code.strains)
  
  # Standardize
  code.strains.pct <- decostand(code.strains, method="total")
  
  # Remove others columns
  code.strains.pct$others <- NULL
  code.strains$others <- NULL
  
  # Order by column sums
  code.strains.pct.ord <- code.strains.pct[,order(colSums(code.strains.pct), decreasing=T)]
  
  # Save standardized, ordered table
  ord <- code.strains.pct.ord
  
  # Name the ISR-types
  colnames(ord) <- paste(code,1:ncol(ord),sep="")
  
  # Make directory for the species
  dir.create(code)
  
  # Output strain table (abundance ordered, rarified):
  write.table(ord, sep="\t", col.names = NA, file=paste(code,"/strains_pct_ord.txt",sep=""))
  
  # Subset sequences
  code.strains.seqs <- ISR_long_R_seqs[as.numeric(substring(colnames(code.strains.pct.ord),first = 4))]
  
  # Export as fasta file
  for (i in 1:length(code.strains.seqs)){
    sink(paste(code,"/",code,"_strains_seqs.fasta",sep=""),append = T)  # Append is set to true, delete previous file if any of same name
    cat(paste(">",colnames(ord)[i],sep=""),code.strains.seqs[i],sep="\n")
    sink()
  }
  
  # Fig 1: Relative Abundance Bar Plot [sequence count of each ISR-type by percentage of total sequences]
  pdf(paste(code,"/",code,"_strains_relab.pdf",sep=""), width = 9, height = 4)
  par(mar = c(5,6,1,1), mgp=c(4, 0.5, 0.5)) # barplot
  barplot(colSums(ord)/sum(colSums(ord))*100,
          ylab = "% of Sequences",
          ylim=c(0,20),
          las=2, col="mediumseagreen", border = "mediumseagreen", cex.axis = 1, cex.names = 1, cex.lab=1.5)
  dev.off()
  
  barplot(colSums(ord)/sum(colSums(ord))*100,
          ylab = "% of Sequences",
          ylim=c(0,20),
          las=2, col="mediumseagreen", border = "mediumseagreen", cex.axis = 1, cex.names = 1, cex.lab=1.5)

  
  
  # Split by subjects
  ord.split <- data.frame(t(data.frame(lapply(split(ord, ceiling((1:nrow(ord))/nlevels(meta$run))), colSums), check.names = F)))
  
  # Convert to pr/ab
  ord.split.prab <- data.frame((ord.split > 0 )*1)
  
  # Fig 2: Prevalence in subjects (sorted by abundance)
  
  pdf(paste(code,"/",code,"_strains_prev_sub_sorted.pdf",sep=""), width = 9, height = 4)
  
  par(mar = c(5,6,1,1), mgp=c(4, 0.5, 0.5)) # barplot
  barplot(colSums(ord.split.prab)/5*100,
          ylab = "% of Subjects", 
          ylim=c(0,100),
          las=2, col="darkslategrey", border = "darkslategrey", cex.axis = 1, cex.names = 1, cex.lab=1.5)
  
  dev.off()
  
  par(mar = c(5,6,1,1), mgp=c(4, 0.5, 0.5)) # barplot
  barplot(colSums(ord.split.prab)/5*100,
          ylab = "% of Subjects", 
          ylim=c(0,100),
          las=2, col="darkslategrey", border = "darkslategrey", cex.axis = 1, cex.names = 1, cex.lab=1.5)
  
  # Reset par
  par(mar = c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0)) 
  
  
  ord.prab <- data.frame((ord > 0 )*1)
  ord.prab.split <- data.frame(t(data.frame(lapply(split(ord.prab, 
                                                         ceiling((1:nrow(ord.prab))/nlevels(meta$run))), colSums), check.names = F)))
  
  # Figure 3: Heatmap
  
  pdf(paste(code,"/",code,"_strains_heatmap.pdf",sep=""), width = 8, height = 5)
  
  heatmap.2(as.matrix(ord.prab.split), Rowv = NA, scale="none",
            rowsep=c(1:4), 
            colsep = c(0:(ncol(ord.prab.split)-1)),
            labRow = FALSE,
            Colv = NA,
            density.info="none",sepwidth = c(0.25,0.25),
            trace="none", dendrogram = "none",
            col=c("ghostwhite","#C7E9C0","#A1D99B","#74C476","#41AB5D","#238B45","#005A32")) 
  
  dev.off()
  
  heatmap.2(as.matrix(ord.prab.split), Rowv = NA, scale="none",
            rowsep=c(1:4), 
            colsep = c(0:(ncol(ord.prab.split)-1)),
            labRow = FALSE,
            Colv = NA,
            density.info="none",sepwidth = c(0.25,0.25),
            trace="none", dendrogram = "none",
            col=c("ghostwhite","#C7E9C0","#A1D99B","#74C476","#41AB5D","#238B45","#005A32")) 
  
  # Levels of most abundant strains
  
  stacked <- strain_stacked_barplots(code, ord, yl)
  
  # Figure S: Rarefaction
  
  pdf(paste(code,"/",code,"_rarefaction_curve.pdf",sep=""), width = 8, height = 6)
  rarecurve(code.strains, step = 20, 
            col = c("red","blue","green","purple","orange")[meta$subject],
            main=otu,
            ylab="ISR-Types",
            xlab="Sequence Counts",
            cex = 0.2)
  dev.off()
  
  rarecurve(code.strains, step = 20, 
            col = c("red","blue","green","purple","orange")[meta$subject],
            main=otu,
            ylab="ISR-Types",
            xlab="Sequence Counts",
            cex = 0.2)
  
  # Output
  list(
    counts = code.strains,
    stacked <- stacked,
    pct_ord = code.strains.pct.ord,
    pct_ord_named = ord,
    ord_split_prab = ord.split.prab
  )
}
```

 

#### Strain\_Stacked\_Barplots

Function for plotting stacked bar graphs for ISR-types

``` r
strain_stacked_barplots <- function(code, std_ord, yl) {
  
  # Set seed
  set.seed(12345)
  
  tbl <- data.frame(std_ord)
  
  if (ncol(tbl) > 25) {
    numb = 25
  } else {
    numb = ncol(tbl)
  }
  
  tbl.top25 <-  tbl[,c(1:numb)]
  
  row.names(tbl.top25) <- gsub("_", "  ", row.names(tbl.top25))
  
  # Set par for bar plot
  
  pdf(paste(code,"/",code,"_stacked_barplots.pdf",sep=""), width = 9, height = 4)
  par(mar = c(5,6,1,1), mgp=c(4, 0.5, 0.5)) # barplot
  
  barplot(as.matrix(t(tbl.top25)),
          ylab = "% of Sequences",
          space=c(0.1,0.1,0.1,0.1,0.1,0.1,1,
                  0.1,0.1,0.1,0.1,0.1,1,
                  0.1,0.1,0.1,0.1,0.1,1,
                  0.1,0.1,0.1,0.1,0.1,1),
          xlab=paste(code, "ISR-types"), 
          ylim=c(0,yl),
          las=2, col=c25[as.numeric(substring(colnames(tbl.top25), 3))])
  
  dev.off()
  
  par(mar = c(5,6,1,1), mgp=c(4, 0.5, 0.5)) # barplot
  barplot(as.matrix(t(tbl.top25)),
          ylab = "% of Sequences",
          space=c(0.1,0.1,0.1,0.1,0.1,0.1,1,
                  0.1,0.1,0.1,0.1,0.1,1,
                  0.1,0.1,0.1,0.1,0.1,1,
                  0.1,0.1,0.1,0.1,0.1,1),
          xlab=paste(code, "ISR-types"), 
          ylim=c(0,yl),
          las=2, col=c25[as.numeric(substring(colnames(tbl.top25), 3))])
  
  hp.top.25 <<- tbl.top25
  
  # Output
  list(
    tbl = tbl,
    tbl_top25 = tbl.top25
  )
}
```
