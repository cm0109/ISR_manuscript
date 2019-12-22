# This script loads the custom functions required for generating the analysis for the ISR manuscript (Mukherjee et al. Microbiome 2018)

# Define 'not include' function
'%!in%' <- function(x,y)!('%in%'(x,y))


# Function to plot distance from base timepoint
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
  #pdf(paste(name,"_Stability.pdf",sep=""), width = 8, height = 6)
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
  
  #dev.off()
  
}  




# R function for plotting distance-based beta diversity analysis results
# This function plots NMDS, H-C dendograms, and stability plots
# Needs function for stability plot to be loaded first

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
  counts.sub.prab.mds <- invisible(metaMDS(counts.sub.prab, wascores=F, autotransform = F, trymax=1000))
  
  # NMDS on relative abundance, output at PDF"
  #pdf(paste(name,"_Prab_NMDS.pdf",sep=""), width = 8, height = 6)
  
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
  
  #dev.off()
  
  # Legend plot
  #pdf(paste(name,"_sub_legend.pdf",sep=""), width = 4, height = 6)
  # Empty Plot for legend
  plot(c(0,1),type="n", axes=F, xlab="", ylab="")
  
  # Add legends
  legend("bottomleft", pch = 21, pt.bg = c("red","blue","green","purple","orange"), 
         legend=levels(meta$subject), title="subjects", cex=2, text.width=0.4, box.col = "white")
  
  #dev.off()
  
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
  #pdf(paste(name,"_dend.pdf",sep=""), width = 8, height = 6)
  
  plot(counts.sub.prab.dend, 
       main=paste("Presence-Abscence ",name),
       ylim = c(0,1.2), cex=1, yaxt="n")
  par(las=2)
  axis(2,cex.axis=1)
  #text(counts.sub.prab.pvclust, col="black", cex=0.7)
  pvrect(counts.sub.prab.pvclust, alpha = 0.95, pv = "bp", border="red4", lwd=1)
  
  #dev.off() 
  
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





# Function to calculate accuracy from Random Forest Classifier
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


#  R script to generate statistics for ASV/OTU counts

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

# R script to generate ISR-type figures

# Define ISR strain plots function
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
  #write.table(ord, sep="\t", col.names = NA, file=paste(code,"/strains_pct_ord.txt",sep=""))
  
  # Fig 1: Relative Abundance Bar Plot [sequence count of each ISR-type by percentage of total sequences]
  #pdf(paste(code,"/",code,"_strains_relab.pdf",sep=""), width = 9, height = 4)
  par(mar = c(5,6,1,1), mgp=c(4, 0.5, 0.5)) # barplot
  barplot(colSums(ord)/sum(colSums(ord))*100,
          ylab = "% of Sequences",
          ylim=c(0,20),
          las=2, col="mediumseagreen", border = "mediumseagreen", cex.axis = 1, cex.names = 1, cex.lab=1.5)
  #dev.off()
  
  
  
  # Split by subjects
  ord.split <- data.frame(t(data.frame(lapply(split(ord, ceiling((1:nrow(ord))/nlevels(meta$run))), colSums), check.names = F)))
  
  # Convert to pr/ab
  ord.split.prab <- data.frame((ord.split > 0 )*1)
  
  # Fig 2: Prevalence in subjects (sorted by abundance)
  
  #pdf(paste(code,"/",code,"_strains_prev_sub_sorted.pdf",sep=""), width = 9, height = 4)
  
  par(mar = c(5,6,1,1), mgp=c(4, 0.5, 0.5)) # barplot
  barplot(colSums(ord.split.prab)/5*100,
          ylab = "% of Subjects", 
          ylim=c(0,100),
          las=2, col="darkslategrey", border = "darkslategrey", cex.axis = 1, cex.names = 1, cex.lab=1.5)
  
  #dev.off()
  
  # Reset par
  par(mar = c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0)) 
  
  
  ord.prab <- data.frame((ord > 0 )*1)
  ord.prab.split <- data.frame(t(data.frame(lapply(split(ord.prab, 
                                                         ceiling((1:nrow(ord.prab))/nlevels(meta$run))), colSums), check.names = F)))
  
  # Figure 3: Heatmap
  
  #pdf(paste(code,"/",code,"_strains_heatmap.pdf",sep=""), width = 8, height = 5)
  
  heatmap.2(as.matrix(ord.prab.split), Rowv = NA, scale="none",
            rowsep=c(1:4), 
            colsep = c(0:(ncol(ord.prab.split)-1)),
            labRow = FALSE,
            Colv = NA,
            density.info="none",sepwidth = c(0.25,0.25),
            trace="none", dendrogram = "none",
            col=c("ghostwhite","#C7E9C0","#A1D99B","#74C476","#41AB5D","#238B45","#005A32")) 
  
  #dev.off()
  
  # Levels of most abundant strains
  
  stacked <- strain_stacked_barplots(code, ord, yl)
  
  # Figure S: Rarefaction
  
  #pdf(paste(code,"/",code,"_rarefaction_curve.pdf",sep=""), width = 8, height = 6)
  rarecurve(code.strains, step = 20, 
            col = c("red","blue","green","purple","orange")[meta$subject],
            main=otu,
            ylab="ISR-Types",
            xlab="Sequence Counts",
            cex = 0.2)
  #dev.off()
  
  # Output
  list(
    counts = code.strains,
    stacked <- stacked,
    pct_ord = code.strains.pct.ord,
    pct_ord_named = ord,
    ord_split_prab = ord.split.prab
  )
}


# R script for plotting stacked bar graphs for ISR-types

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
  
  #pdf(paste(code,"/",code,"_stacked_barplots.pdf",sep=""), width = 9, height = 4)
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
  
  #dev.off()
  
  hp.top.25 <<- tbl.top25
  
  # Output
  list(
    tbl = tbl,
    tbl_top25 = tbl.top25
  )
}

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
