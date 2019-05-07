# Overall settings for analysis: Adoption Study
# Updated: 03/19/2019

# Load required libraries
library(reshape2)
library(vegan)
library(ggplot2)
library(GUniFrac)
library(ggpubr)
library(readxl)
library(RColorBrewer)
library(seqRFLP)

# Set options
options(scipen=999)

# Color Schemes
# Feed colors
feed_cols <- c("mediumseagreen","khaki3","coral1")
age_colors <- colorRampPalette(brewer.pal(9, "YlGn")[2:9])

# Custom functions used in this analysis

# For loading R objects from other environments
loadER <- function (rdata, object) {
  local({
    load(rdata)
    eval(parse(text = object))
  })
}



# Function mellt.dist to melt distance objects
melt.dist <- function(x,only_identity = FALSE,omit_identity=TRUE) {
  library(reshape2)
  if(!is.matrix(x)) {
    x <- as.matrix(x)
  }
  y <- melt(x)
  if(omit_identity) {
    l <- list();
    nr <- nrow(x);
    nin <- c();
    for(i in 0:(nr - 1)){
      rs <- (i * nr)+1;
      if(only_identity==FALSE) {
        l[[i+1]] <- (rs+i):(rs+nr-1)
      } else {
        l[[i+1]] <- (rs+i)
      }
    }
    l <- -c(unlist(l))
    y <- y[l,]
  }
  return(y)
}

# Same fam plot function
AD.melt.dist.sf <- function(asv_table,meta){
  
  # Calculate Distance matrix
  dist <- vegdist(asv_table, method = "bray")
  dist1 <- vegdist(asv_table, method = "jaccard")
  
  # Melt dist
  dist_melt <- melt.dist(dist) # BC
  dist1_melt <- melt.dist(dist1) # JC
  
  # Add relevant meta fields
  dist_melt$fam1 <- meta$family_id[match(dist_melt$Var1, meta$sample)] # BC
  dist_melt$fam2 <- meta$family_id[match(dist_melt$Var2, meta$sample)] # BC
  
  dist1_melt$fam1 <- meta$family_id[match(dist1_melt$Var1, meta$sample)] # JC
  dist1_melt$fam2 <- meta$family_id[match(dist1_melt$Var2, meta$sample)] # JC
  
  # Adding same/diff markers
  dist_melt$family <- "Same Family" # BC
  dist_melt[dist_melt$fam1!=dist_melt$fam2,]$family <- "Diff Family" # BC
  
  dist1_melt$family <- "Same Family" # JC
  dist1_melt[dist1_melt$fam1!=dist1_melt$fam2,]$family <- "Diff Family" # JC
  
  # Subset for samefam
  dist_melt_sf <- dist_melt[dist_melt$family=="Same Family",] # BC
  dist_melt_sf$status <- meta$status[match(dist_melt_sf$Var2, meta$sample)] # BC
  
  dist1_melt_sf <- dist1_melt[dist1_melt$family=="Same Family",] # JC
  dist1_melt_sf$status <- meta$status[match(dist1_melt_sf$Var2, meta$sample)] # JC
  
  # Subset for different family
  dist_melt_df <- dist_melt[dist_melt$family=="Diff Family",] # BC
  
  
  # Output
  list(
    dm = dist_melt,
    dmj = dist1_melt,
    dmsf = dist_melt_sf,
    dmdf = dist_melt_df,
    jm = dist1_melt,
    jmsf = dist1_melt_sf
  )
}


# Updated unrelated distance randomization script
unrelated_dist1 <- function(dfdm){
  mylist <- list()
  mylist1 <- list()
  
  # Split on 1st sample
  dfdm_split <- split(dfdm, droplevels(dfdm$Var1))
  
  # Loop to extact random row from split dfs
  for (i in 1:length(dfdm_split)){
    df <- dfdm_split[[i]]
    mylist[[i]] <- df[sample(nrow(df), 1),]
  }
  
  # Make new df
  mydf <- do.call("rbind",mylist)
  
  # Split on 2nd sample
  mydf_split <- split(mydf, droplevels(mydf$Var2))
  
  # Loop to extact random row from split dfs
  for (i in 1:length(mydf_split)){
    df1 <- mydf_split[[i]]
    mylist1[[i]] <- df1[sample(nrow(df1), 1),]
  }
  
  # Make new df
  mydf1 <- do.call("rbind",mylist1)
  
  return(list(
    df = mydf1
  ))
  
}


# Overall Analysis Strategy

# Main Input
# 1. Clinical Metadata File (hiseq_all_samples_clinical_030419.xlsx)
# 2. Sequencing statistics (AD_DADA2_stats_all.xlsx, AD_16S_D2_stats.xlsx, [16S counts table from sql])
# 3. Count tables (ISR seqtabs, 16Sd2 seqtabs, 16Score counts)

# Step1: Clean sequencing stats to remove low read samples, repeats, non-clinical samples etc

# Step2: Prepare specific meta files

# Generate relevant meta file for I. ISR, II. 16Sd2, III. 16Score
# For each, generate meta file for a. Mother-child comparions, b. Family comparisons, c. Sibling comparisons

# Step3: Subset counts by meta for analysis

# Analysis Strategy

# Combined Analysis
# Separate Analysis



