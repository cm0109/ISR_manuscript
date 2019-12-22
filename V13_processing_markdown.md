V13 Processing
================
Chiranjit Mukherjee
5/10/2019

Notes
-----

R Script for processing the V13 sequences

 

### Query for V13\_long\_R\_results.txt

    echo "select query_id, score, match_size, otu, genus_group, family, order_, class, phylum, sample, 
    location, run from 16S_all_tax  where score >= 0.98 and longitudinal='T' and location='R';" 
    | mysql -u root ISR_2017_long > V13_long_R_results.txt

 

### Libraries

``` r
library(reshape2)
```

 

### Load BLAST result file

``` r
V13_long_R_core.results <- read.table(file="/Volumes/GriffenLeysLab/Troy/ISR_manuscript_markdown/V13_long_R_results.txt", header=T, sep="\t")
```

 

### Construct OTU table

``` r
# Make OTU table:
V13_long_R_core.results$value <- 1
V13_long_R_core.pivot <- dcast(V13_long_R_core.results, sample ~ otu, fun.aggregate = sum)

# Assign Row Names
row.names(V13_long_R_core.pivot) <- V13_long_R_core.pivot$sample
V13_long_R_core.pivot$sample <- NULL

V13_long_R.core <- V13_long_R_core.pivot[,colSums(V13_long_R_core.pivot) > 0]
```

 

### Clean-Up

``` r
# Clean up names to match meta file:

rownames(V13_long_R.core) <- sub("_16S", "", rownames(V13_long_R.core))
rownames(V13_long_R.core) <- sub("MP_", "", rownames(V13_long_R.core))

# Rename for analysis:
rownames(V13_long_R.core) <- sub("_R_run3", "_T1", rownames(V13_long_R.core))
rownames(V13_long_R.core) <- sub("_R_run4", "_T2", rownames(V13_long_R.core))
rownames(V13_long_R.core) <- sub("_R_run5", "_T3", rownames(V13_long_R.core))
rownames(V13_long_R.core) <- sub("_R_run6", "_T4", rownames(V13_long_R.core))
rownames(V13_long_R.core) <- sub("_R_run7", "_T5", rownames(V13_long_R.core))
rownames(V13_long_R.core) <- sub("_R_run8", "_T6", rownames(V13_long_R.core))
rownames(V13_long_R.core) <- sub("1003_", "S1_", rownames(V13_long_R.core))
rownames(V13_long_R.core) <- sub("1004_", "S2_", rownames(V13_long_R.core))
rownames(V13_long_R.core) <- sub("1006_", "S3_", rownames(V13_long_R.core))
rownames(V13_long_R.core) <- sub("1007_", "S4_", rownames(V13_long_R.core))
rownames(V13_long_R.core) <- sub("2002_", "S5_", rownames(V13_long_R.core))
```

 

### Output Counts Table

``` r
write.table(V13_long_R.core, "V13_long_R.core.txt", quote=FALSE, sep="\t", col.names = NA)
```
