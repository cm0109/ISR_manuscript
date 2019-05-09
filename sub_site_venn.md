Subgingival Mother-Child Matching
================
Troy Ellison
5/6/2019

Script Notes
============

This script takes Mother-Child Pair Combinations of subgingival samples and output
1. How many sequences the Mother-Child pair share
2. How many sequences are present in the mother's sample only
3. How many sequences are present in the child's sequence only
These statistics are then separated based on the Mother-Child relations (Biological, Adopted, or differrent family), then displayed in venn diagrams

 

Environment Notes
=================

Run in R Studio 1.0.153
Using R version 3.5.2 (2018-12-20) -- "Eggshell Igloo"
Copyright (C) 2018 The R Foundation for Statistical Computing
Platform: x86\_64-apple-darwin15.6.0 (64-bit)

 

### Load Required Libraries

``` r
library(dplyr)
library(data.table)
library(vegan)
library(VennDiagram)
library(ggplot2)
```

 

### Load Required R Objects

``` r
load("/Volumes/GriffenLeysLab/Troy/site_venn/ISR_sub_mch_prab.rds")
load("/Volumes/GriffenLeysLab/Troy/ISR_comb_meta_mch.rds")
```

 

Mother-Child Matching
---------------------

 

### Create subset tables of mother sequences and child sequences

``` r
ISR_comb_meta_mch_m <- ISR_comb_meta_mch[ISR_comb_meta_mch$subject_type == "Biological Mother" | ISR_comb_meta_mch$subject_type == "Adopted Mother",]
ISR_comb_meta_mch_c <- ISR_comb_meta_mch[ISR_comb_meta_mch$subject_type == "Biological Child" | ISR_comb_meta_mch$subject_type == "Adopted Child",]
ISR.mother.sub <- ISR_sub_mch_prab[rownames(ISR_sub_mch_prab) %in% ISR_comb_meta_mch_m$sample, ]
ISR.child.sub <- ISR_sub_mch_prab[rownames(ISR_sub_mch_prab) %in% ISR_comb_meta_mch_c$sample, ]
```

 

### Construct Dataframe for Mother-Child combinations

``` r
ISR.matches.sub <- data.frame(Mother = character(),
                          Child = character(),
                          Matches = integer(),
                          MotherOnly = integer(),
                          ChildOnly = integer())
```

 

### Fill Dataframe with Mother-Child combinations

``` r
#Loop through the mother and child subtables
for(mother in 1:nrow(ISR.mother.sub)){
  for(child in 1:nrow(ISR.child.sub)){
    
    #Set sequence counters for matched, mother only and child only    
    matched <- 0
    motherOnly <- 0
    childOnly <- 0
    
    #Loop through each row in the mother subtable    
    for(seq in 1:ncol(ISR.mother.sub)){
      
      #Check if the mother and/or child is positive for the sequence
      if(ISR.child.sub[child, seq] == ISR.mother.sub[mother,seq] && ISR.mother.sub[mother,seq]==1 ){
        
        #If both positive, increment the matched counter
        matched <- matched+1
        
        #If mother only positive, increment the mother only counter
      } else if(ISR.mother.sub[mother,seq] == 1){
        motherOnly <- motherOnly+1
        
        #If child only positive, increment the child only counter
      } else if(ISR.child.sub[child, seq] == 1){
        childOnly <- childOnly+1
      }
    }
    
    #Create a new row of the mother-child combination with the matched, mother-only, and child-only counters,
    #and bind it to the matches dataframe
    ISR.row <- data.frame(rownames(ISR.mother.sub[mother,]), rownames(ISR.child.sub[child,]), matched, motherOnly, childOnly)
    ISR.matches.sub <- rbind(ISR.matches.sub, ISR.row)
  }
  
}
```

 

Venn Diagram Plotting
---------------------

 

### Adding Family Info and Subset Sequences of the Same Family

``` r
# Combining loop results
colnames(ISR.matches.sub) <- c("Sub1","Sub2","Shared","Sub1_only","Sub2_only")
View(ISR.matches.sub)

# Add family info
ISR.matches.sub$fam_id1 <- ISR_comb_meta_mch$family_id[match(ISR.matches.sub$Sub1, ISR_comb_meta_mch$sample)]
ISR.matches.sub$fam_id2 <- ISR_comb_meta_mch$family_id[match(ISR.matches.sub$Sub2, ISR_comb_meta_mch$sample)]
ISR.matches.sub$fam <- "Same Family"
ISR.matches.sub[ISR.matches.sub$fam_id1 != ISR.matches.sub$fam_id2,]$fam <- "Different Family"

# Subset same fam
ISR.matches_sf.sub <- ISR.matches.sub[ISR.matches.sub$fam == "Same Family",]
colnames(ISR.matches_sf.sub) <- c("Sub1","Sub2","Shared","Mother_only","Child_only", "fam_id1","fam_id2","fam")
View(ISR.matches_sf.sub)

# Add adoption status 
ISR.matches_sf.sub$status <- ISR_comb_meta_mch$status[match(ISR.matches_sf.sub$fam_id1, ISR_comb_meta_mch$family_id)]
```

 

### Create Dataframe for Biological Mother-Child Pairs

``` r
sub.ISR.matches_sf_bio <- ISR.matches_sf.sub[ISR.matches_sf.sub$status == "Biological",]
View(sub.ISR.matches_sf_bio)
nrow(sub.ISR.matches_sf_bio) # 54
```

    ## [1] 54

``` r
# Extract columns for plotting
sub.ISR_strncnt_sf_bio <- sub.ISR.matches_sf_bio[,c("fam_id1", "Shared","Mother_only","Child_only")]

# Clean up table for computing %
row.names(sub.ISR_strncnt_sf_bio) <- sub.ISR_strncnt_sf_bio$fam_id1
sub.ISR_strncnt_sf_bio$fam_id1 <- NULL
View(sub.ISR_strncnt_sf_bio)

# Convert to numeric
sub.ISR_strncnt_sf_bio$Shared <- as.numeric(as.character(sub.ISR_strncnt_sf_bio$Shared))
sub.ISR_strncnt_sf_bio$Child_only <- as.numeric(as.character(sub.ISR_strncnt_sf_bio$Child_only))
sub.ISR_strncnt_sf_bio$Mother_only <- as.numeric(as.character(sub.ISR_strncnt_sf_bio$Mother_only))

# Conver to %
sub.ISR_strncnt_sf_bio.pct <- decostand(sub.ISR_strncnt_sf_bio, method = "total")

# Create df for plotting
sub.bio_pie_df <- data.frame(
  group = c("Shared", "Child Only", "Mother Only"),
  value = c(mean(sub.ISR_strncnt_sf_bio.pct$Shared)*100, mean(sub.ISR_strncnt_sf_bio.pct$Child_only)*100, mean(sub.ISR_strncnt_sf_bio.pct$Mother_only)*100)
)
```

 

### Create Dataframe for Adopted Mother-Child Pairs

``` r
# Subset for Adopted
sub.ISR.matches_sf_adp <- ISR.matches_sf.sub[ISR.matches_sf.sub$status == "Adopted",]
View(sub.ISR.matches_sf_adp)
nrow(sub.ISR.matches_sf_adp) # 43
```

    ## [1] 43

``` r
# Extract columns for plotting
sub.ISR_strncnt_sf_adp <- sub.ISR.matches_sf_adp[,c("fam_id1", "Shared","Mother_only","Child_only")]

# Clean up table for computing %
row.names(sub.ISR_strncnt_sf_adp) <- sub.ISR_strncnt_sf_adp$fam_id1
sub.ISR_strncnt_sf_adp$fam_id1 <- NULL
View(sub.ISR_strncnt_sf_adp)

# Convert to numeric
sub.ISR_strncnt_sf_adp$Shared <- as.numeric(as.character(sub.ISR_strncnt_sf_adp$Shared))
sub.ISR_strncnt_sf_adp$Child_only <- as.numeric(as.character(sub.ISR_strncnt_sf_adp$Child_only))
sub.ISR_strncnt_sf_adp$Mother_only <- as.numeric(as.character(sub.ISR_strncnt_sf_adp$Mother_only))

# Conver to %
sub.ISR_strncnt_sf_adp.pct <- decostand(sub.ISR_strncnt_sf_adp, method = "total")

# Create df for plotting
sub.adp_pie_df <- data.frame(
  group = c("Shared", "Child Only", "Mother Only"),
  value = c(mean(sub.ISR_strncnt_sf_adp.pct$Shared)*100, mean(sub.ISR_strncnt_sf_adp.pct$Child_only)*100, mean(sub.ISR_strncnt_sf_adp.pct$Mother_only)*100)
)
```

 

### Create Dataframe for Different Family Mother-Child Pairs

``` r
# Subset diff fam
sub.ISR.matches_df <- ISR.matches.sub[ISR.matches.sub$fam == "Different Family",]
View(sub.ISR.matches_df)

# Add adoption status 
sub.ISR.matches_df$sub_type1 <- ISR_comb_meta_mch$subject_type[match(sub.ISR.matches_df$Sub1, ISR_comb_meta_mch$sample)]
sub.ISR.matches_df$sub_type2 <- ISR_comb_meta_mch$subject_type[match(sub.ISR.matches_df$Sub2, ISR_comb_meta_mch$sample)]

# Remove bio/adp from sub type
sub.ISR.matches_df$sub_type1 <- gsub("Biological ", "", sub.ISR.matches_df$sub_type1)
sub.ISR.matches_df$sub_type1 <- gsub("Adopted ", "", sub.ISR.matches_df$sub_type1)
sub.ISR.matches_df$sub_type2 <- gsub("Biological ", "", sub.ISR.matches_df$sub_type2)
sub.ISR.matches_df$sub_type2 <- gsub("Adopted ", "", sub.ISR.matches_df$sub_type2)

# Add rel field
sub.ISR.matches_df$rel <- paste(sub.ISR.matches_df$sub_type1, sub.ISR.matches_df$sub_type2, sep="_")

# Select mother-child only
sub.ISR.matches_df_mc <- sub.ISR.matches_df[sub.ISR.matches_df$rel != "Child_Child" & sub.ISR.matches_df$rel != "Mother_Mother",]
View(sub.ISR.matches_df_mc)
sub.ISR.matches_df_mc$rel %>% unique # "Child_Mother" "Mother_Child"
```

    ## [1] "Mother_Child"

``` r
# Separate "Child_Mother" and "Mother_Child" for re-ordering fields
sub.ISR.matches_df_mc1 <- sub.ISR.matches_df_mc[sub.ISR.matches_df_mc$rel == "Child_Mother",]
sub.ISR.matches_df_mc2 <- sub.ISR.matches_df_mc[sub.ISR.matches_df_mc$rel == "Mother_Child",]

# Subset columns
sub.ISR.matches_df_mc1 <- sub.ISR.matches_df_mc1[,c("Sub1", "Sub2", "Shared","Sub1_only","Sub2_only")]
sub.ISR.matches_df_mc2 <- sub.ISR.matches_df_mc2[,c("Sub1", "Sub2", "Shared","Sub1_only","Sub2_only")]

# Rename columns
colnames(sub.ISR.matches_df_mc1) <- c("Sub1","Sub2","Shared","Child_only","Mother_only")
colnames(sub.ISR.matches_df_mc2) <- c("Sub1","Sub2","Shared","Mother_only","Child_only")

# Re-order columns for mc2
sub.ISR.matches_df_mc2 <- sub.ISR.matches_df_mc2[,c(1,2,3,5,4)]

# Combine mc1 and mc2
sub.ISR.matches_df_mc_comb <- rbind.data.frame(sub.ISR.matches_df_mc1, sub.ISR.matches_df_mc2)
sub.ISR_strncnt_diffam <- sub.ISR.matches_df_mc_comb[,c(3:5)]

# Convert to numeric
sub.ISR_strncnt_diffam$Shared <- as.numeric(as.character(sub.ISR_strncnt_diffam$Shared))
sub.ISR_strncnt_diffam$Child_only <- as.numeric(as.character(sub.ISR_strncnt_diffam$Child_only))
sub.ISR_strncnt_diffam$Mother_only <- as.numeric(as.character(sub.ISR_strncnt_diffam$Mother_only))

# Conver to %
sub.ISR_strncnt_diffam.pct <- decostand(sub.ISR_strncnt_diffam, method = "total")

# Create df for plotting
sub.diffam_pie_df <- data.frame(
  group = c("Shared", "Child Only", "Mother Only"),
  value = c(round(mean(sub.ISR_strncnt_diffam.pct$Shared)*100,1), mean(sub.ISR_strncnt_diffam.pct$Child_only)*100, mean(sub.ISR_strncnt_diffam.pct$Mother_only)*100)
)
```

 

Venn Diagrams
-------------

 

### Biological Mother-Child Pairs

<img src="sub_site_venn_files/figure-markdown_github/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

 

### Adopted Mother-Child Pairs

<img src="sub_site_venn_files/figure-markdown_github/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

 

### Different Family Mother-Child Pairs

<img src="sub_site_venn_files/figure-markdown_github/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

 

### Legend

<img src="sub_site_venn_files/figure-markdown_github/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />
