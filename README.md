# ISR_manuscript


### AD_settings.R:  
Functions necessary to run the ISR_manuscript R Script

### DADA2 Markdown:  
R Script for processing ISR sequences with DADA2 Final ISR Analysis for Publication Processing 30 ISR Longitudinal "Right" side samples with DADA2 to generate ASVs  
**REQUIRES:** Directory of zipped forward FASTQ files

### V12 DADA2 Processing:  
R Script for processing 16S V1-V2 Longitudinal "Right" side samples with DADA2  
**REQUIRES:** Directory of zipped forward FASTQ files

### V13 Processing:  
Script for Processing V13 sequences  
**REQUIRES:** Database query results table (see script for query)

### ISR Manuscript:  
Analysis Script for comparing pipelines & ISR-type Analysis  
**REQUIRES:** Functions from AD_Settings.R, DADA2 .RData (ASV tables, meta file sequences), time differences table from clinical sampling schedules, V12 DADA2 .RDATA (ASV table), V13 Processing .RDATA (Counts Table), sequence data "ISR_long_R.fa", ISR Blast Results Table, ISR Database taxonomy table

#### Libraries Used:  
    reshape2  
    vegan  
    ggplot2  
    GUniFrac  
    ggpubr  
    readxl  
    RColorBrewer  
    seqRFLP  
    plotrix  
    gplots  
    reshape2
    data.table  
    Biostrings  
    dada2 (packageVersion"dada2")
    phyloseq
    
