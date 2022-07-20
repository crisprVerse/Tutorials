Building a gRNA database for the human transcriptome using CasRx
(RfxCas13d)
================

-   [Introduction](#introduction)
-   [Loading necessary packages](#loading-necessary-packages)
    -   [Specifying the genome](#specifying-the-genome)
    -   [Specifying the genome index](#specifying-the-genome-index)
    -   [Specifying the nuclease](#specifying-the-nuclease)
    -   [Specifying on-target scoring
        methods](#specifying-on-target-scoring-methods)
    -   [Specifying gene models and TSS
        annotations](#specifying-gene-models-and-tss-annotations)
-   [Building a complete annotation for a given
    gene](#building-a-complete-annotation-for-a-given-gene)
    -   [Converting the `GuideSet` object to a list of
        data.frames](#converting-the-guideset-object-to-a-list-of-dataframes)
-   [Building a complete annotation for all
    genes](#building-a-complete-annotation-for-all-genes)
-   [Reproducibility](#reproducibility)

Authors: Jean-Philippe Fortin

Date: July 16, 2022

# Introduction

In this tutorial, we provide reproducible code to design and annotates
all gRNAs in a transcriptome, for a given organism, and for a
DNA-targeting nuclease such as Cas13d. As an example, we design CRISPRkd
gRNAs for the human transcriptome, in GRCh38 coordinates, for the
nuclease RfxCas13d (CasRx).

# Loading necessary packages

We first load the necessary packages:

``` r
library(crisprBase)
library(crisprScore)
library(crisprDesign)
library(crisprDesignData)
library(BSgenome.Hsapiens.UCSC.hg38)
```

### Specifying the genome

We specify a `BSGenome` object that contains the reference genome of the
target organism:

``` r
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
```

### Specifying the genome index

We specify the file path of a Bowtie index build on a transcriptome
(different from an index built on a reference genome). Briefly, the
index was build on a fasta file storing mRNA sequences for a given gene
model (here Ensembl release 104). See the Genome index tutorial to learn
how to build such an index.

``` r
bowtie_index <- "/Users/fortinj2/crisprIndices/bowtie/ensembl_human_104/ensembl_human_104"
```

### Specifying the nuclease

The `crisprBase` provides functionalities for creating custom
`crisprNuclease` objects, and also provides already-available
`crisprNuclease` objects such as the commonly-used `SpCas9`,
`enAsCas12a`, `CasRx` nucleases.

Here, we will use the DNA-targeting nuclease `CasRx`:

``` r
data(CasRx, package="crisprBase")
crisprNuclease=CasRx
```

### Specifying on-target scoring methods

We specify which on-target scoring methods we want to use:

``` r
scoring_methods=c("casrxrf")
```

More information about which scoring methods are available for a given
nuclease can be obtained using the following:

``` r
crisprScore::scoringMethodsInfo
```

    ##        method   nuclease left right       type      label len
    ## 1    ruleset1     SpCas9  -24     5  On-target   RuleSet1  30
    ## 2     azimuth     SpCas9  -24     5  On-target    Azimuth  30
    ## 3      deephf     SpCas9  -20     2  On-target     DeepHF  23
    ## 4      lindel     SpCas9  -33    31  On-target     Lindel  65
    ## 5         mit     SpCas9  -20     2 Off-target        MIT  23
    ## 6         cfd     SpCas9  -20     2 Off-target        CFD  23
    ## 7    deepcpf1   AsCas12a   -4    29  On-target   DeepCpf1  34
    ## 8     enpamgb enAsCas12a   -4    29  On-target    EnPAMGB  34
    ## 9  crisprscan     SpCas9  -26     8  On-target CRISPRscan  35
    ## 10    casrxrf      CasRx   NA    NA  On-target   CasRx-RF  NA
    ## 11   crisprai     SpCas9  -19     2  On-target   CRISPRai  22
    ## 12 crisprater     SpCas9  -20    -1  On-target CRISPRater  20
    ## 13 deepspcas9     SpCas9  -24     5  On-target DeepSpCas9  30
    ## 14   ruleset3     SpCas9  -24     5  On-target   RuleSet3  30

### Specifying gene models and TSS annotations

We also need a gene model formatted as a `TxDbObject`. The
`crisprDesignData` contains such objects for both the human and mouse
genomes, in GRCh38 (hg38) and GRCm38 (mm10) coordinates, respectively.
Ensembl gene models were used to generate such objects.

``` r
data(txdb_human, package="crisprDesignData")
txObject <- txdb_human
```

The `crisprDesignData` vignette shows how to create such objects from
other transcriptomes.

# Building a complete annotation for a given gene

The `designCompleteAnnotation` function in `crisprDesign` provides a
one-step workflow to design and annotate gRNAs targeting a coding gene
for a user-specific combination of parameters.

Here, we design all gRNAs targeting the main KRAS isoform
(ENST00000311936) for CRISPRkd applications using CasRx:

``` r
gs <- designCompleteAnnotation(queryValue="ENST00000311936",
                               queryColumn="tx_id",
                               modality="CRISPRkd",
                               bsgenome=bsgenome,
                               bowtie_index=bowtie_index,
                               crisprNuclease=CasRx,
                               txObject=txObject,
                               n_mismatches=1,
                               scoring_methods=scoring_methods)
```

    ## [designCompleteAnnotation] Adding sequence statistics 
    ## [designCompleteAnnotation] Adding spacer alignments

    ## Loading required namespace: crisprBwa

    ## [runCrisprBowtie] Searching for CasRx protospacers 
    ## [designCompleteAnnotation] Adding gene annotation 
    ## [designCompleteAnnotation] Adding on-target scores

    ## [addOnTargetScores] Adding casrxrf scores.

    ## [addCasRxScores] Calculating MFE features 
    ## [addCasRxScores] Calculating accessibility features 
    ## [addCasRxScores] Calculating hybridization features 
    ## [addCasRxScores] Calculating nucleotide density features 
    ## [getCasRxRFScores] Calculating scores

    ## snapshotDate(): 2022-04-26

    ## see ?crisprScoreData and browseVignettes('crisprScoreData') for documentation

    ## loading from cache

The resulting object is a `GuideSet` object containing the fully
annotated gRNAs, as described in the `crisprDesign` vignette. See the
Tutorial `CRISPRkd` to learn more about `GuideSet` objects and how to
access their rich annotations for RNA-targeting nucleases such as CasRx.

``` r
gs
```

    ## GuideSet object with 5283 ranges and 19 metadata columns:
    ##                               seqnames    ranges strand |
    ##                                  <Rle> <IRanges>  <Rle> |
    ##      ENST00000311936_1 ENST00000311936        24      + |
    ##      ENST00000311936_2 ENST00000311936        25      + |
    ##      ENST00000311936_3 ENST00000311936        26      + |
    ##      ENST00000311936_4 ENST00000311936        27      + |
    ##      ENST00000311936_5 ENST00000311936        28      + |
    ##                    ...             ...       ...    ... .
    ##   ENST00000311936_5279 ENST00000311936      5302      + |
    ##   ENST00000311936_5280 ENST00000311936      5303      + |
    ##   ENST00000311936_5281 ENST00000311936      5304      + |
    ##   ENST00000311936_5282 ENST00000311936      5305      + |
    ##   ENST00000311936_5283 ENST00000311936      5306      + |
    ##                                    protospacer            pam  pam_site
    ##                                 <DNAStringSet> <DNAStringSet> <numeric>
    ##      ENST00000311936_1 CTAGGCGGCGGCCGCGGCGGCGG              A        24
    ##      ENST00000311936_2 TAGGCGGCGGCCGCGGCGGCGGA              G        25
    ##      ENST00000311936_3 AGGCGGCGGCCGCGGCGGCGGAG              G        26
    ##      ENST00000311936_4 GGCGGCGGCCGCGGCGGCGGAGG              C        27
    ##      ENST00000311936_5 GCGGCGGCCGCGGCGGCGGAGGC              A        28
    ##                    ...                     ...            ...       ...
    ##   ENST00000311936_5279 GTAATGTAATAAAAATAGTTACA              G      5302
    ##   ENST00000311936_5280 TAATGTAATAAAAATAGTTACAG              T      5303
    ##   ENST00000311936_5281 AATGTAATAAAAATAGTTACAGT              G      5304
    ##   ENST00000311936_5282 ATGTAATAAAAATAGTTACAGTG              A      5305
    ##   ENST00000311936_5283 TGTAATAAAAATAGTTACAGTGA              C      5306
    ##                         cut_site          region percentGC     polyA     polyC
    ##                        <numeric>     <character> <numeric> <logical> <logical>
    ##      ENST00000311936_1        NA ENST00000311936      91.3     FALSE     FALSE
    ##      ENST00000311936_2        NA ENST00000311936      87.0     FALSE     FALSE
    ##      ENST00000311936_3        NA ENST00000311936      91.3     FALSE     FALSE
    ##      ENST00000311936_4        NA ENST00000311936      95.7     FALSE     FALSE
    ##      ENST00000311936_5        NA ENST00000311936      95.7     FALSE     FALSE
    ##                    ...       ...             ...       ...       ...       ...
    ##   ENST00000311936_5279        NA ENST00000311936      17.4     FALSE     FALSE
    ##   ENST00000311936_5280        NA ENST00000311936      17.4     FALSE     FALSE
    ##   ENST00000311936_5281        NA ENST00000311936      17.4     FALSE     FALSE
    ##   ENST00000311936_5282        NA ENST00000311936      21.7     FALSE     FALSE
    ##   ENST00000311936_5283        NA ENST00000311936      21.7     FALSE     FALSE
    ##                            polyG     polyT startingGGGGG     n0_tx     n1_tx
    ##                        <logical> <logical>     <logical> <numeric> <numeric>
    ##      ENST00000311936_1     FALSE     FALSE         FALSE         4         0
    ##      ENST00000311936_2     FALSE     FALSE         FALSE         4         0
    ##      ENST00000311936_3     FALSE     FALSE         FALSE         4         9
    ##      ENST00000311936_4     FALSE     FALSE         FALSE         4        81
    ##      ENST00000311936_5     FALSE     FALSE         FALSE         4        67
    ##                    ...       ...       ...           ...       ...       ...
    ##   ENST00000311936_5279     FALSE      TRUE         FALSE         2         0
    ##   ENST00000311936_5280     FALSE      TRUE         FALSE         2         0
    ##   ENST00000311936_5281     FALSE      TRUE         FALSE         2         0
    ##   ENST00000311936_5282     FALSE      TRUE         FALSE         2         0
    ##   ENST00000311936_5283     FALSE      TRUE         FALSE         2         0
    ##                          n0_gene   n1_gene
    ##                        <numeric> <numeric>
    ##      ENST00000311936_1         1         0
    ##      ENST00000311936_2         1         0
    ##      ENST00000311936_3         1         6
    ##      ENST00000311936_4         1        26
    ##      ENST00000311936_5         1        18
    ##                    ...       ...       ...
    ##   ENST00000311936_5279         1         0
    ##   ENST00000311936_5280         1         0
    ##   ENST00000311936_5281         1         0
    ##   ENST00000311936_5282         1         0
    ##   ENST00000311936_5283         1         0
    ##                                                                                     alignments
    ##                                                                                  <GRangesList>
    ##      ENST00000311936_1  ENST00000256078:24:+,ENST00000256078:1579:+,ENST00000356971:3858:+,...
    ##      ENST00000311936_2  ENST00000311936:26:+,ENST00000311936:1459:+,ENST00000359466:1898:+,...
    ##      ENST00000311936_3 ENST00000640969:224:+,ENST00000311936:1504:+,ENST00000409806:4289:+,...
    ##      ENST00000311936_4  ENST00000539137:76:+,ENST00000311936:1954:+,ENST00000357816:7677:+,...
    ##      ENST00000311936_5  ENST00000629930:60:+,ENST00000256078:2079:+,ENST00000357992:2883:+,...
    ##                    ...                                                                     ...
    ##   ENST00000311936_5279  ENST00000311936:57:+,ENST00000256078:2076:+,ENST00000356971:3859:+,...
    ##   ENST00000311936_5280  ENST00000556131:57:+,ENST00000311936:1952:+,ENST00000357039:1661:+,...
    ##   ENST00000311936_5281  ENST00000557334:64:+,ENST00000256078:2077:+,ENST00000357137:5313:+,...
    ##   ENST00000311936_5282  ENST00000307877:86:+,ENST00000311936:1953:+,ENST00000357310:5122:+,...
    ##   ENST00000311936_5283 ENST00000331327:214:+,ENST00000256078:2078:+,ENST00000357635:2042:+,...
    ##                                                                  geneAnnotation
    ##                                                            <SplitDataFrameList>
    ##      ENST00000311936_1    ENST00000311936_1:ENST00000311936:ENSG00000133703:...
    ##      ENST00000311936_2    ENST00000311936_2:ENST00000311936:ENSG00000133703:...
    ##      ENST00000311936_3    ENST00000311936_3:ENST00000311936:ENSG00000133703:...
    ##      ENST00000311936_4    ENST00000311936_4:ENST00000311936:ENSG00000133703:...
    ##      ENST00000311936_5    ENST00000311936_5:ENST00000311936:ENSG00000133703:...
    ##                    ...                                                      ...
    ##   ENST00000311936_5279 ENST00000311936_5279:ENST00000311936:ENSG00000133703:...
    ##   ENST00000311936_5280 ENST00000311936_5280:ENST00000311936:ENSG00000133703:...
    ##   ENST00000311936_5281 ENST00000311936_5281:ENST00000311936:ENSG00000133703:...
    ##   ENST00000311936_5282 ENST00000311936_5282:ENST00000311936:ENSG00000133703:...
    ##   ENST00000311936_5283 ENST00000311936_5283:ENST00000311936:ENSG00000133703:...
    ##                             enzymeAnnotation score_casrxrf
    ##                         <SplitDataFrameList>     <numeric>
    ##      ENST00000311936_1 FALSE:FALSE:FALSE:...            NA
    ##      ENST00000311936_2 FALSE:FALSE:FALSE:...      0.387198
    ##      ENST00000311936_3 FALSE:FALSE:FALSE:...            NA
    ##      ENST00000311936_4 FALSE:FALSE:FALSE:...            NA
    ##      ENST00000311936_5 FALSE:FALSE:FALSE:...            NA
    ##                    ...                   ...           ...
    ##   ENST00000311936_5279 FALSE:FALSE:FALSE:...            NA
    ##   ENST00000311936_5280 FALSE:FALSE:FALSE:...            NA
    ##   ENST00000311936_5281 FALSE:FALSE:FALSE:...            NA
    ##   ENST00000311936_5282 FALSE:FALSE:FALSE:...            NA
    ##   ENST00000311936_5283 FALSE:FALSE:FALSE:...            NA
    ##   -------
    ##   seqinfo: 1 sequence from custom genome
    ##   crisprNuclease: CasRx

### Converting the `GuideSet` object to a list of data.frames

The `flattenGuideSet` function in `crisprDesign` is a convenience
function to convert the complex `GuideSet` object into a set of 5
data.frames that can be saved as plain text files:

``` r
dfs <- flattenGuideSet(gs)
```

We can look at the names of the data.frames:

``` r
names(dfs)
```

    ## [1] "primary"          "alignments"       "geneAnnotation"   "enzymeAnnotation"

As an example, let’s look at the first rows of the primary data.frame:

``` r
head(dfs$primary)
```

    ##                  ID                  spacer             protospacer
    ## 1 ENST00000311936_1 CCGCCGCCGCGGCCGCCGCCTAG CTAGGCGGCGGCCGCGGCGGCGG
    ## 2 ENST00000311936_2 TCCGCCGCCGCGGCCGCCGCCTA TAGGCGGCGGCCGCGGCGGCGGA
    ## 3 ENST00000311936_3 CTCCGCCGCCGCGGCCGCCGCCT AGGCGGCGGCCGCGGCGGCGGAG
    ## 4 ENST00000311936_4 CCTCCGCCGCCGCGGCCGCCGCC GGCGGCGGCCGCGGCGGCGGAGG
    ## 5 ENST00000311936_5 GCCTCCGCCGCCGCGGCCGCCGC GCGGCGGCCGCGGCGGCGGAGGC
    ## 6 ENST00000311936_6 TGCCTCCGCCGCCGCGGCCGCCG CGGCGGCCGCGGCGGCGGAGGCA
    ##               chr start end strand pam pam_site cut_site          region
    ## 1 ENST00000311936     1  23      +   A       24       NA ENST00000311936
    ## 2 ENST00000311936     2  24      +   G       25       NA ENST00000311936
    ## 3 ENST00000311936     3  25      +   G       26       NA ENST00000311936
    ## 4 ENST00000311936     4  26      +   C       27       NA ENST00000311936
    ## 5 ENST00000311936     5  27      +   A       28       NA ENST00000311936
    ## 6 ENST00000311936     6  28      +   G       29       NA ENST00000311936
    ##   percentGC polyA polyC polyG polyT startingGGGGG n0_tx n1_tx n0_gene n1_gene
    ## 1      91.3 FALSE FALSE FALSE FALSE         FALSE     4     0       1       0
    ## 2      87.0 FALSE FALSE FALSE FALSE         FALSE     4     0       1       0
    ## 3      91.3 FALSE FALSE FALSE FALSE         FALSE     4     9       1       6
    ## 4      95.7 FALSE FALSE FALSE FALSE         FALSE     4    81       1      26
    ## 5      95.7 FALSE FALSE FALSE FALSE         FALSE     4    67       1      18
    ## 6      91.3 FALSE FALSE FALSE FALSE         FALSE     4    44       1       5
    ##   score_casrxrf
    ## 1            NA
    ## 2     0.3871984
    ## 3            NA
    ## 4            NA
    ## 5            NA
    ## 6            NA

# Building a complete annotation for all genes

We first get all possibles transcripts from our gene model:

``` r
tx_ids <- unique(txObject$cds$tx_id)
head(tx_ids)
```

    ## [1] "ENST00000641515" "ENST00000616016" "ENST00000618323" "ENST00000342066"
    ## [5] "ENST00000338591" "ENST00000622660"

and build an index to loop over for generating the GuideSet objects:

``` r
tx_index <- seq_along(tx_ids)
head(tx_index)
```

    ## [1] 1 2 3 4 5 6

We specify where to save the GuideSets:

``` r
dir <- "./crisprko_casrx_hg38"
if (!dir.exists(dir)){
  dir.create(dir, recursive=TRUE)
}
```

We are now ready to generate and save all data with the function
`designCompleteAnnotation` from `crisprDesign`. The function was
designed to be as comprehensive as possible to design and annotate gRNAs
in one step. It does the following:

-   Extract the DNA/RNA sequences with `queryTss`/`queryTxDB`
-   Design gRNAs with `findSpacers`
-   Remove gRNAs targeting repeat elements with `removeRepeats`
-   Characterize spacer sequences with `addSequenceFeatures`
-   Find on- and off-targets with `addSpacerAlignmentsIterative`
-   Add gene annotation with `addGeneAnnotation`
-   Add TSS annotation with `addTssAnnotation`
-   Add on-target efficiency scores with `addOnTargetScores`
-   Add off-target specificity scores with `addOffTargetScores`
-   Add SNP annotation with `addSNPAnnotation`
-   Add restriction enzymes information with `addRestrictionEnzymes`

We are now looping over all transcripts to generate the data:

``` r
lapply(tx_index, function(tx){
    gs <- designCompleteAnnotation(queryValue=tx,
                                   queryColumn="tx_id",
                                   modality="CRISPRkd",
                                   bsgenome=bsgenome,
                                   bowtie_index=bowtie_index,
                                   crisprNuclease=CasRx,
                                   txObject=txObject,
                                   n_mismatches=3,
                                   scoring_methods=scoring_methods)
    write.rds(gs, file=file.path(dir, paste0(tx, ".rds")))
})
```

This loop can be modified by the user to use an embarrassingly-parallel
approach to save time using the `BiocParallel` package, for instance.

# Reproducibility

``` r
sessionInfo()
```

    ## R Under development (unstable) (2022-03-21 r81954)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Catalina 10.15.7
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] BSgenome.Hsapiens.UCSC.hg38_1.4.4 BSgenome_1.63.5                  
    ##  [3] rtracklayer_1.55.4                Biostrings_2.63.2                
    ##  [5] XVector_0.35.0                    GenomicRanges_1.47.6             
    ##  [7] GenomeInfoDb_1.31.6               IRanges_2.29.1                   
    ##  [9] S4Vectors_0.33.11                 crisprDesignData_0.99.8          
    ## [11] crisprDesign_0.99.102             crisprScore_1.1.9                
    ## [13] crisprScoreData_1.1.3             ExperimentHub_2.3.5              
    ## [15] AnnotationHub_3.3.9               BiocFileCache_2.3.4              
    ## [17] dbplyr_2.1.1                      BiocGenerics_0.41.2              
    ## [19] crisprBase_1.1.2                 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] rjson_0.2.21                  ellipsis_0.3.2               
    ##  [3] Rbowtie_1.35.0                rstudioapi_0.13              
    ##  [5] bit64_4.0.5                   interactiveDisplayBase_1.33.0
    ##  [7] AnnotationDbi_1.57.1          fansi_1.0.2                  
    ##  [9] xml2_1.3.3                    cachem_1.0.6                 
    ## [11] knitr_1.37                    jsonlite_1.8.0               
    ## [13] Rsamtools_2.11.0              png_0.1-7                    
    ## [15] shiny_1.7.1                   BiocManager_1.30.16          
    ## [17] readr_2.1.2                   compiler_4.2.0               
    ## [19] httr_1.4.2                    basilisk_1.9.2               
    ## [21] assertthat_0.2.1              Matrix_1.4-0                 
    ## [23] fastmap_1.1.0                 cli_3.3.0                    
    ## [25] later_1.3.0                   htmltools_0.5.2              
    ## [27] prettyunits_1.1.1             tools_4.2.0                  
    ## [29] glue_1.6.2                    GenomeInfoDbData_1.2.7       
    ## [31] crisprBowtie_1.1.1            dplyr_1.0.8                  
    ## [33] rappdirs_0.3.3                Rcpp_1.0.8.3                 
    ## [35] Biobase_2.55.0                vctrs_0.3.8                  
    ## [37] crisprBwa_1.1.2               xfun_0.30                    
    ## [39] stringr_1.4.0                 mime_0.12                    
    ## [41] lifecycle_1.0.1               restfulr_0.0.13              
    ## [43] XML_3.99-0.9                  zlibbioc_1.41.0              
    ## [45] basilisk.utils_1.9.1          vroom_1.5.7                  
    ## [47] VariantAnnotation_1.41.3      hms_1.1.1                    
    ## [49] promises_1.2.0.1              MatrixGenerics_1.7.0         
    ## [51] parallel_4.2.0                SummarizedExperiment_1.25.3  
    ## [53] yaml_2.3.5                    curl_4.3.2                   
    ## [55] memoise_2.0.1                 reticulate_1.24              
    ## [57] biomaRt_2.51.3                stringi_1.7.6                
    ## [59] RSQLite_2.2.12                BiocVersion_3.15.0           
    ## [61] BiocIO_1.5.0                  randomForest_4.7-1           
    ## [63] GenomicFeatures_1.47.13       filelock_1.0.2               
    ## [65] BiocParallel_1.29.18          rlang_1.0.2                  
    ## [67] pkgconfig_2.0.3               matrixStats_0.61.0           
    ## [69] bitops_1.0-7                  evaluate_0.15                
    ## [71] lattice_0.20-45               purrr_0.3.4                  
    ## [73] GenomicAlignments_1.31.2      bit_4.0.4                    
    ## [75] tidyselect_1.1.2              magrittr_2.0.2               
    ## [77] R6_2.5.1                      generics_0.1.2               
    ## [79] DelayedArray_0.21.2           DBI_1.1.2                    
    ## [81] withr_2.5.0                   pillar_1.7.0                 
    ## [83] KEGGREST_1.35.0               RCurl_1.98-1.6               
    ## [85] tibble_3.1.6                  dir.expiry_1.3.0             
    ## [87] crayon_1.5.0                  utf8_1.2.2                   
    ## [89] tzdb_0.2.0                    rmarkdown_2.13               
    ## [91] progress_1.2.2                grid_4.2.0                   
    ## [93] blob_1.2.2                    digest_0.6.29                
    ## [95] xtable_1.8-4                  httpuv_1.6.5                 
    ## [97] Rbwa_1.1.0
