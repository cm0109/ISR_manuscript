# Repository for Analysis Related to Mukherjee Microbiome 2018 Publication

## See: https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0535-z


### DADA2 Markdown:   (Demo only)
R Script for processing ISR sequences with DADA2 Final ISR Analysis for Publication Processing 30 ISR Longitudinal "Right" side samples with DADA2 to generate ASVs  
**REQUIRES:** Directory of zipped forward FASTQ files

### V12 DADA2 Processing:  (Demo only)
R Script for processing 16S V1-V2 Longitudinal "Right" side samples with DADA2  
**REQUIRES:** Directory of zipped forward FASTQ files

### V13 Processing:  (Demo only)
Script for Processing V13 sequences  
**REQUIRES:** Database query results table (see script for query)



## Analysis Scripts

### ISR_custom_functions.R:  
Custom R functions to produce the visualizations and analysis


### ISR MS Demo:  
Analysis Script for comparing pipelines & strain level Analysis  
**REQUIRES:** Functions from ISR_custom_functions.R, Input files

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
    
