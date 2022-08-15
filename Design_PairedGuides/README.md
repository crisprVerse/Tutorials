Using crisprDesign to design paired gRNAs
================

-   [Introduction](#introduction)
-   [Getting started](#getting-started)
    -   [Installation](#installation)
    -   [Terminology](#terminology)
    -   [Paired gRNA design overview](#paired-grna-design-overview)
-   [A simple example: paired gRNAs flanking an
    exon](#a-simple-example-paired-grnas-flanking-an-exon)
-   [Use cases](#use-cases)
    -   [Double nicking with
        CRISPR/Cas9](#double-nicking-with-crisprcas9)
    -   [Dual-promoter systems](#dual-promoter-systems)
    -   [Multiplexing gRNAs with arrayed spacers
        (enAsCas12a)](#multiplexing-grnas-with-arrayed-spacers-enascas12a)
    -   [Nanopore Cas9-targeted sequencing
        (nCATS)](#nanopore-cas9-targeted-sequencing-ncats)
-   [Session Info](#session-info)
-   [References](#references)

Authors: Jean-Philippe Fortin, Luke Hoberecht

Date: 15 August, 2022

# Introduction

In this tutorial, we illustrate the main functionalities of
`crisprDesign` for designing pairs of gRNAs.

# Getting started

## Installation

First, we install the necessary packages for this tutorial from
Bioconductor using the following commands:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("crisprBase")
BiocManager::install("crisprDesign")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
```

as well as the data package `crisprDesignData` from GitHub:

``` r
install.packages("devtools")
devtools::install_github("Jfortin1/crisprDesignData")
```

## Terminology

See the [CRISPRko Cas9 design
tutorial](https://github.com/crisprVerse/Tutorials/tree/master/Design_CRISPRko_Cas9)
to get familiar with the terminology used throughout this tutorial.

## Paired gRNA design overview

There is a number of applications that require the design of pairs of
gRNAs. Here is a list of cases that will be covered in this tutorial:

1.  Double nicking with CRISPR/Cas9 (Ran et al. et al. 2013)
2.  Dual-promoter screening systems (Han et al. 2017)
3.  Multiplexing gRNAs with arrayed spacers (enAsCas12a) (DeWeirdt et
    al. et al. 2021)
4.  Nanopore Cas9-targeted sequencing (nCATS) (Gilpatrick et al. 2020)

Before we dive into the different applications, we will describe a
general use case to go over the main paired gRNA design features.

# A simple example: paired gRNAs flanking an exon

To illustrate the general concept behind paired gRNA design, we will
show here how to design pairs of gRNAs flanking the second exon of the
canonical transcript (ENST00000311936) of the human gene KRAS
(ENSG00000133703).

We first start by loading the necessary packages:

``` r
library(crisprDesign)
library(crisprDesignData)
library(crisprBase)
library(BSgenome.Hsapiens.UCSC.hg38)
```

We will be designing gRNAs for the SpCas9 nuclease. We load the `SpCas9`
nuclease object from the `crisprBase` package (see the `crisprBase`
[vignette](https://github.com/crisprVerse/crisprBase) for instructions
on how to create or load alternative nucleases):

``` r
data(SpCas9, package="crisprBase")
```

Let’s get the genomic coordinates of the second exon. Let’s first obtain
from `crisprDesignData` a `GRangesList` object that defines the genomic
coordinates (in hg38 coordinates) of coding genes in the human genome:

``` r
data(txdb_human, package="crisprDesignData")
```

We can get the exon coodinates using the function `queryTxObject`:

``` r
exons <- queryTxObject(txObject=txdb_human,
                       featureType="exons",
                       queryColumn="tx_id",
                       queryValue="ENST00000311936")
exons
```

    ## GRanges object with 5 ranges and 14 metadata columns:
    ##            seqnames            ranges strand |           tx_id         gene_id
    ##               <Rle>         <IRanges>  <Rle> |     <character>     <character>
    ##   region_1    chr12 25250751-25250929      - | ENST00000311936 ENSG00000133703
    ##   region_2    chr12 25245274-25245395      - | ENST00000311936 ENSG00000133703
    ##   region_3    chr12 25227234-25227412      - | ENST00000311936 ENSG00000133703
    ##   region_4    chr12 25225614-25225773      - | ENST00000311936 ENSG00000133703
    ##   region_5    chr12 25205246-25209911      - | ENST00000311936 ENSG00000133703
    ##                 protein_id        tx_type gene_symbol         exon_id exon_rank
    ##                <character>    <character> <character>     <character> <integer>
    ##   region_1            <NA> protein_coding        KRAS ENSE00003903543         1
    ##   region_2 ENSP00000256078 protein_coding        KRAS ENSE00000936617         2
    ##   region_3 ENSP00000256078 protein_coding        KRAS ENSE00001719809         3
    ##   region_4 ENSP00000256078 protein_coding        KRAS ENSE00001644818         4
    ##   region_5 ENSP00000308495 protein_coding        KRAS ENSE00002456976         5
    ##            cds_start   cds_end  tx_start    tx_end   cds_len exon_start
    ##            <integer> <integer> <integer> <integer> <integer>  <integer>
    ##   region_1      <NA>      <NA>  25205246  25250929       567       <NA>
    ##   region_2  25245274  25245384  25205246  25250929       567       <NA>
    ##   region_3  25227234  25227412  25205246  25250929       567       <NA>
    ##   region_4  25225614  25225773  25205246  25250929       567       <NA>
    ##   region_5  25209795  25209911  25205246  25250929       567       <NA>
    ##             exon_end
    ##            <integer>
    ##   region_1      <NA>
    ##   region_2      <NA>
    ##   region_3      <NA>
    ##   region_4      <NA>
    ##   region_5      <NA>
    ##   -------
    ##   seqinfo: 25 sequences (1 circular) from hg38 genome

and we select the second exon:

``` r
exon <- exons[exons$exon_rank==2]
names(exon) <- "exon_kras"
exon
```

    ## GRanges object with 1 range and 14 metadata columns:
    ##             seqnames            ranges strand |           tx_id         gene_id
    ##                <Rle>         <IRanges>  <Rle> |     <character>     <character>
    ##   exon_kras    chr12 25245274-25245395      - | ENST00000311936 ENSG00000133703
    ##                  protein_id        tx_type gene_symbol         exon_id
    ##                 <character>    <character> <character>     <character>
    ##   exon_kras ENSP00000256078 protein_coding        KRAS ENSE00000936617
    ##             exon_rank cds_start   cds_end  tx_start    tx_end   cds_len
    ##             <integer> <integer> <integer> <integer> <integer> <integer>
    ##   exon_kras         2  25245274  25245384  25205246  25250929       567
    ##             exon_start  exon_end
    ##              <integer> <integer>
    ##   exon_kras       <NA>      <NA>
    ##   -------
    ##   seqinfo: 25 sequences (1 circular) from hg38 genome

The exon is on chr12, and spans the region 25245274-25245395 (122
nucleotides in length). We aim to design gRNAs pairs for which one gRNA
is located upstream of the exon, and another located downstream of the
exon. Let’s define those regions to be 500 nucleotides on each side:

``` r
library(IRanges)
regionUpstream   <- IRanges::flank(exon, width=100, start=FALSE)
regionDownstream <- IRanges::flank(exon, width=100, start=TRUE)
names(regionUpstream) <- "upstreamTarget"
names(regionDownstream) <- "downstreamTarget"
```

We specify the BSgenome for hg38 coordinates

``` r
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
```

and we are now ready to find all spacer pairs:

``` r
pairs <- findSpacerPairs(x1=regionUpstream,
                         x2=regionDownstream,
                         bsgenome=bsgenome,
                         crisprNuclease=SpCas9)
```

The `x1` and `x2` arguments specify the genomic regions in which gRNAs
at position 1 and position 2 should be targeting, respectively. The
function finds all possible pair combinations between spacers found in
the region specified by `x1` and spacers found in the region specified
by `x2`. Let’ first name our pairs:

``` r
names(pairs) <- paste0("pair_", seq_along(pairs))
```

Let’s see what the results look like:

``` r
head(pairs, n=3)
```

    ## PairedGuideSet object with 3 pairs and 4 metadata columns:
    ##                     first           second | pamOrientation pamDistance
    ##                <GuideSet>       <GuideSet> |    <character>   <numeric>
    ##   pair_1 chr12:25245201:+ chr12:25245397:- |             in         196
    ##   pair_2 chr12:25245215:- chr12:25245397:- |            rev         182
    ##   pair_3 chr12:25245233:+ chr12:25245397:- |             in         164
    ##          spacerDistance cutLength
    ##               <integer> <numeric>
    ##   pair_1            198       202
    ##   pair_2            163       182
    ##   pair_3            166       170

The returned object is a `PairedGuideSet`, which can be though of a list
of two `GuideSet` objects. The first and second `GuideSet` store
information about gRNAs at position 1 and position 2, respectively. They
can be accessed using the `first` and `second` functions:

``` r
grnas1 <- first(pairs)
head(grnas1, n=3)
```

    ## GuideSet object with 3 ranges and 5 metadata columns:
    ##            seqnames    ranges strand |          protospacer            pam
    ##               <Rle> <IRanges>  <Rle> |       <DNAStringSet> <DNAStringSet>
    ##   spacer_1    chr12  25245201      + | GTAATAAGTACTCATGAAAA            TGG
    ##   spacer_2    chr12  25245215      - | CCATTCTTTGATACAGATAA            AGG
    ##   spacer_3    chr12  25245233      + | CCTTTATCTGTATCAAAGAA            TGG
    ##             pam_site  cut_site         region
    ##            <numeric> <numeric>    <character>
    ##   spacer_1  25245201  25245198 upstreamTarget
    ##   spacer_2  25245215  25245218 upstreamTarget
    ##   spacer_3  25245233  25245230 upstreamTarget
    ##   -------
    ##   seqinfo: 640 sequences (1 circular) from hg38 genome
    ##   crisprNuclease: SpCas9

and

``` r
grnas2 <- second(pairs)
head(grnas2, n=3)
```

    ## GuideSet object with 3 ranges and 5 metadata columns:
    ##            seqnames    ranges strand |          protospacer            pam
    ##               <Rle> <IRanges>  <Rle> |       <DNAStringSet> <DNAStringSet>
    ##   spacer_1    chr12  25245397      - | TTTTCATTATTTTTATTATA            AGG
    ##   spacer_1    chr12  25245397      - | TTTTCATTATTTTTATTATA            AGG
    ##   spacer_1    chr12  25245397      - | TTTTCATTATTTTTATTATA            AGG
    ##             pam_site  cut_site           region
    ##            <numeric> <numeric>      <character>
    ##   spacer_1  25245397  25245400 downstreamTarget
    ##   spacer_1  25245397  25245400 downstreamTarget
    ##   spacer_1  25245397  25245400 downstreamTarget
    ##   -------
    ##   seqinfo: 640 sequences (1 circular) from hg38 genome
    ##   crisprNuclease: SpCas9

The `pamOrientation` function returns the PAM orientation of the pairs:

``` r
head(pamOrientation(pairs))
```

    ## [1] "in"  "rev" "in"  "rev" "rev" "fwd"

and takes 4 different values: `in` (for PAM-in configuration), `out`
(for PAM-out configuration), `fwd` (both gRNAs target the forward
strand), and `rev` (both gRNAs target the reverse strand); see figure
below for an illustration of the PAM orientations for the SpCas9
nuclease. The importance of the PAM orientation is application-specific,
and will be discussed in the relevant sections below.

<img src="./figures/paired_simplified.svg" title="Different PAM orientations for Cas9 paired gRNAs" alt="Different PAM orientations for Cas9 paired gRNAs" width="75%" style="display: block; margin: auto;" />

The function `pamDistance` returns the distance between the PAM sites of
the two gRNAs. The function `cutLength` returns the distance between the
cut sites of the two gRNAs, and the function `spacerDistance` returns
the distance between the two spacer sequences of the gRNAs.

Most functionalities available for designing single gRNAs (`GuideSet`
annotation functions described in [this
tutorial](https://github.com/crisprVerse/Tutorials/tree/master/Design_CRISPRko_Cas9))
work similarly for `PairedGuideSet` objects. This includes:

-   `addSequenceFeatures`
-   `addSpacerAlignments`
-   `addGeneAnnotation`
-   `addTssAnnotation`
-   `addOnTargetScores`
-   `addOffTargetScores`
-   `addPamScores`
-   `addSNPAnnotation`
-   `addRestrictionEnzymes`
-   `addCompositeScores`
-   `addConservationScores`

Each function adds an annotation to the first and second `GuideSet`
objects stored in the `PairedGuideSet`. Let’s look at an example using
`addSequenceFeatures`:

``` r
pairs <- addSequenceFeatures(pairs)
```

and let’s look at the `GuideSet` in the first position:

``` r
head(first(pairs), n=3)
```

    ## GuideSet object with 3 ranges and 12 metadata columns:
    ##            seqnames    ranges strand |          protospacer            pam
    ##               <Rle> <IRanges>  <Rle> |       <DNAStringSet> <DNAStringSet>
    ##   spacer_1    chr12  25245201      + | GTAATAAGTACTCATGAAAA            TGG
    ##   spacer_2    chr12  25245215      - | CCATTCTTTGATACAGATAA            AGG
    ##   spacer_3    chr12  25245233      + | CCTTTATCTGTATCAAAGAA            TGG
    ##             pam_site  cut_site         region          coordID percentGC
    ##            <numeric> <numeric>    <character>      <character> <numeric>
    ##   spacer_1  25245201  25245198 upstreamTarget chr12_25245201_+        25
    ##   spacer_2  25245215  25245218 upstreamTarget chr12_25245215_-        30
    ##   spacer_3  25245233  25245230 upstreamTarget chr12_25245233_+        30
    ##                polyA     polyC     polyG     polyT startingGGGGG
    ##            <logical> <logical> <logical> <logical>     <logical>
    ##   spacer_1      TRUE     FALSE     FALSE     FALSE         FALSE
    ##   spacer_2     FALSE     FALSE     FALSE     FALSE         FALSE
    ##   spacer_3     FALSE     FALSE     FALSE     FALSE         FALSE
    ##   -------
    ##   seqinfo: 640 sequences (1 circular) from hg38 genome
    ##   crisprNuclease: SpCas9

This comes in handy to filter out pairs with unwanted sgRNA
characteristics, e.g. sgRNA with polyT stretches:

``` r
good1 <- !first(pairs)$polyT
good2 <- !second(pairs)$polyT
pairs <- pairs[good1 & good2]
```

To select the final candidate pairs, one could filter out pairs with low
predicted on-target activity. Let’s add the DeepHF on-target activity
score:

``` r
pairs <- addOnTargetScores(pairs, methods="deephf")
```

    ## [addOnTargetScores] Adding deephf scores.

    ## snapshotDate(): 2022-04-26

    ## see ?crisprScoreData and browseVignettes('crisprScoreData') for documentation

    ## loading from cache

and only keep pairs for which both gRNAs have a score greater than 0.5:

``` r
good1 <- first(pairs)$score_deephf>=0.5
good2 <- second(pairs)$score_deephf>=0.5
pairs <- pairs[good1 & good2]
```

This leaves us with 2 candidate pairs:

``` r
pairs
```

    ## PairedGuideSet object with 2 pairs and 4 metadata columns:
    ##                      first           second | pamOrientation pamDistance
    ##                 <GuideSet>       <GuideSet> |    <character>   <numeric>
    ##   pair_14 chr12:25245239:- chr12:25245472:- |            rev         233
    ##   pair_19 chr12:25245239:- chr12:25245475:- |            rev         236
    ##           spacerDistance cutLength
    ##                <integer> <numeric>
    ##   pair_14            214       233
    ##   pair_19            217       236

One can get the spacer sequences using the `spacers` accessor function
as usual:

``` r
spacers(pairs)
```

    ## DataFrame with 2 rows and 2 columns
    ##                  first               second
    ##         <DNAStringSet>       <DNAStringSet>
    ## 1 AATATGCATATTACTGGTGC TTTGTATTAAAAGGTACTGG
    ## 2 AATATGCATATTACTGGTGC GAGTTTGTATTAAAAGGTAC

Finally, one might want to check for off-targets. Let’s specify the path
of the bowtie index for the human reference genome:

``` r
bowtie_index <- "/Users/fortinj2/crisprIndices/bowtie/hg38/hg38"
```

For instructions on how to build a Bowtie index from a given reference
genome, see the [genome index
tutorial](https://github.com/crisprVerse/Tutorials/tree/master/Building_Genome_Indices)
or the [crisprBowtie page](https://github.com/crisprVerse/crisprBowtie).

We are now ready to search for off-targets with up to 3 mismatches:

``` r
pairs <- addSpacerAlignments(pairs,
                             txObject=txdb_human,
                             aligner_index=bowtie_index,
                             bsgenome=bsgenome,
                             n_mismatches=3)
```

    ## Loading required namespace: crisprBwa

    ## [runCrisprBowtie] Using BSgenome.Hsapiens.UCSC.hg38 
    ## [runCrisprBowtie] Searching for SpCas9 protospacers

Let’s only keep pairs that have no off-targets in other coding
sequences:

``` r
good1 <- first(pairs)$n1_c==0 & first(pairs)$n2_c==0 & first(pairs)$n3_c==0
good2 <- second(pairs)$n1_c==0 & second(pairs)$n2_c==0 & second(pairs)$n3_c==0
pairs <- pairs[good1 & good2]
pairs
```

    ## PairedGuideSet object with 2 pairs and 4 metadata columns:
    ##                      first           second | pamOrientation pamDistance
    ##                 <GuideSet>       <GuideSet> |    <character>   <numeric>
    ##   pair_14 chr12:25245239:- chr12:25245472:- |            rev         233
    ##   pair_19 chr12:25245239:- chr12:25245475:- |            rev         236
    ##           spacerDistance cutLength
    ##                <integer> <numeric>
    ##   pair_14            214       233
    ##   pair_19            217       236

# Use cases

## Double nicking with CRISPR/Cas9

Discussed in (Ran et al. et al. 2013)

## Dual-promoter systems

Discussed in (Han et al. 2017)

## Multiplexing gRNAs with arrayed spacers (enAsCas12a)

Discussed in (DeWeirdt et al. et al. 2021)

## Nanopore Cas9-targeted sequencing (nCATS)

Discussed in (Han et al. 2017)

# Session Info

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
    ##  [1] crisprScoreData_1.1.3             ExperimentHub_2.3.5              
    ##  [3] AnnotationHub_3.3.9               BiocFileCache_2.3.4              
    ##  [5] dbplyr_2.1.1                      BSgenome.Hsapiens.UCSC.hg38_1.4.4
    ##  [7] BSgenome_1.64.0                   rtracklayer_1.55.4               
    ##  [9] Biostrings_2.64.0                 XVector_0.35.0                   
    ## [11] GenomicRanges_1.48.0              GenomeInfoDb_1.32.2              
    ## [13] IRanges_2.30.0                    S4Vectors_0.33.11                
    ## [15] BiocGenerics_0.42.0               crisprDesignData_0.99.14         
    ## [17] crisprDesign_0.99.127             crisprBase_1.1.5                 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] rjson_0.2.21                  ellipsis_0.3.2               
    ##  [3] Rbowtie_1.36.0                rstudioapi_0.13              
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
    ## [37] crisprBwa_1.1.2               crisprScore_1.1.13           
    ## [39] xfun_0.30                     stringr_1.4.0                
    ## [41] mime_0.12                     lifecycle_1.0.1              
    ## [43] restfulr_0.0.13               XML_3.99-0.9                 
    ## [45] zlibbioc_1.41.0               basilisk.utils_1.9.1         
    ## [47] vroom_1.5.7                   VariantAnnotation_1.41.3     
    ## [49] hms_1.1.1                     promises_1.2.0.1             
    ## [51] MatrixGenerics_1.7.0          parallel_4.2.0               
    ## [53] SummarizedExperiment_1.25.3   yaml_2.3.5                   
    ## [55] curl_4.3.2                    memoise_2.0.1                
    ## [57] reticulate_1.25               biomaRt_2.51.3               
    ## [59] stringi_1.7.6                 RSQLite_2.2.12               
    ## [61] BiocVersion_3.15.0            highr_0.9                    
    ## [63] BiocIO_1.5.0                  randomForest_4.7-1           
    ## [65] GenomicFeatures_1.47.13       filelock_1.0.2               
    ## [67] BiocParallel_1.29.18          rlang_1.0.4                  
    ## [69] pkgconfig_2.0.3               matrixStats_0.61.0           
    ## [71] bitops_1.0-7                  evaluate_0.15                
    ## [73] lattice_0.20-45               purrr_0.3.4                  
    ## [75] GenomicAlignments_1.31.2      bit_4.0.4                    
    ## [77] tidyselect_1.1.2              magrittr_2.0.2               
    ## [79] R6_2.5.1                      generics_0.1.2               
    ## [81] DelayedArray_0.21.2           DBI_1.1.2                    
    ## [83] pillar_1.7.0                  withr_2.5.0                  
    ## [85] KEGGREST_1.35.0               RCurl_1.98-1.6               
    ## [87] tibble_3.1.6                  dir.expiry_1.3.0             
    ## [89] crayon_1.5.0                  utf8_1.2.2                   
    ## [91] tzdb_0.2.0                    rmarkdown_2.13               
    ## [93] progress_1.2.2                grid_4.2.0                   
    ## [95] blob_1.2.2                    digest_0.6.29                
    ## [97] xtable_1.8-4                  httpuv_1.6.5                 
    ## [99] Rbwa_1.1.0

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-deweirdt2021optimization" class="csl-entry">

DeWeirdt, Peter C, Kendall R Sanson, Annabel K Sangree, Mudra Hegde,
Ruth E Hanna, Marissa N Feeley, Audrey L Griffith, et al., et al. 2021.
“Optimization of AsCas12a for Combinatorial Genetic Screens in Human
Cells.” *Nature Biotechnology* 39 (1): 94–104.

</div>

<div id="ref-gilpatrick2020targeted" class="csl-entry">

Gilpatrick, Timothy, Isac Lee, James E Graham, Etienne Raimondeau,
Rebecca Bowen, Andrew Heron, Bradley Downs, Saraswati Sukumar, Fritz J
Sedlazeck, and Winston Timp. 2020. “Targeted Nanopore Sequencing with
Cas9-Guided Adapter Ligation.” *Nature Biotechnology* 38 (4): 433–38.

</div>

<div id="ref-han2017synergistic" class="csl-entry">

Han, Kyuho, Edwin E Jeng, Gaelen T Hess, David W Morgens, Amy Li, and
Michael C Bassik. 2017. “Synergistic Drug Combinations for Cancer
Identified in a CRISPR Screen for Pairwise Genetic Interactions.”
*Nature Biotechnology* 35 (5): 463–74.

</div>

<div id="ref-ran2013double" class="csl-entry">

Ran, F Ann, Patrick D Hsu, Chie-Yu Lin, Jonathan S Gootenberg, Silvana
Konermann, Alexandro E Trevino, David A Scott, et al., et al. 2013.
“Double Nicking by RNA-Guided CRISPR Cas9 for Enhanced Genome Editing
Specificity.” *Cell* 154 (6): 1380–89.

</div>

</div>
