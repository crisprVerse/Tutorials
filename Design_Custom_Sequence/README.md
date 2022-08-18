Using crisprDesign to design gRNAs for custom sequences
================

-   [Introduction](#introduction)
-   [Getting started](#getting-started)
    -   [Installation](#installation)
-   [Use case: designing gRNAs against
    EGFP](#use-case-designing-grnas-against-egfp)
    -   [Loading necessary packages](#loading-necessary-packages)
    -   [Obtaining the DNA sequence](#obtaining-the-dna-sequence)
    -   [Constructing the `GuideSet`
        object:](#constructing-the-guideset-object)
    -   [Finding off-targets in the human genome to find gRNAs specific
        to
        EGFP](#finding-off-targets-in-the-human-genome-to-find-grnas-specific-to-egfp)
    -   [Predicting on-target activity](#predicting-on-target-activity)
    -   [Final selection](#final-selection)
-   [Session Info](#session-info)
-   [References](#references)

Authors: Jean-Philippe Fortin, Luke Hoberecht

Date: 17 August, 2022

# Introduction

In this tutorial, we illustrate the main functionalities of
`crisprDesign` for designing gRNAs for custom sequences. To design gRNAs
for targets located in an organism genome, see the [introductory
CRISPRko
tutorial](https://github.com/crisprVerse/Tutorials/tree/master/Design_CRISPRko_Cas9).

# Getting started

## Installation

First, we install the necessary packages from Bioconductor using the
following commands:

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

# Use case: designing gRNAs against EGFP

Suppose we are engineering a human cell line to express the enhanced
green fluorescent protein (EGFP) marker, and that we want to design
gRNAs that knockout EGFP as experimental controls. Such control gRNAs
should target EGFP with (1) high efficiency, and (2) should be specific
to EGFP, that is, should not target the cell genome (human genome in
this case). Supposed also that the cell line is also stably expressing
SpCas9.

## Loading necessary packages

We first start by loading the necessary packages:

``` r
library(Biostrings)
library(crisprBase)
library(crisprDesign)
library(crisprDesignData)
library(BSgenome.Hsapiens.UCSC.hg38)
```

## Obtaining the DNA sequence

In the folder `data`, we have included a fasta file containing the DNA
sequence of the EGFP marker. The sequence was obtained from the
[SnapGene
website](https://www.snapgene.com/resources/plasmid-files/?set=fluorescent_protein_genes_and_plasmids&plasmid=EGFP)

We can read in the fasta file using the `readDNAStringSet` function from
the package `Biostrings`:

``` r
dna <- Biostrings::readDNAStringSet("data/egfp.fa")
names(dna) <- "EGFP"
dna
```

    ## DNAStringSet object of length 1:
    ##     width seq                                               names               
    ## [1]   720 ATGGTGAGCAAGGGCGAGGAGCT...GCATGGACGAGCTGTACAAGTAA EGFP

This could also be simply constructed from a regular string:

``` r
dna <- "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA"
dna <- DNAStringSet(dna)
names(dna) <- "EGFP"
```

(Note that the function also accepts a simple string, which would be
internally converted into a `DNAStringSet`). This is the custom sequence
input that we will use to design gRNAs.

## Constructing the `GuideSet` object:

Next, we design all possible SpCas0 gRNAs targeting EGFP. First, we load
the SpCas9 object from the `crisprBase` package:

``` r
data(SpCas9, package="crisprBase")
```

and we design gRNAs using the function `findSpacers` from
`crisprDesign`:

``` r
gs <- findSpacers(dna, 
                  crisprNuclease=SpCas9)
head(gs)
```

    ## GuideSet object with 6 ranges and 5 metadata columns:
    ##            seqnames    ranges strand |          protospacer            pam
    ##               <Rle> <IRanges>  <Rle> |       <DNAStringSet> <DNAStringSet>
    ##   spacer_1     EGFP        30      + | AAGGGCGAGGAGCTGTTCAC            CGG
    ##   spacer_2     EGFP        31      + | AGGGCGAGGAGCTGTTCACC            GGG
    ##   spacer_3     EGFP        31      - | GACCAGGATGGGCACCACCC            CGG
    ##   spacer_4     EGFP        32      + | GGGCGAGGAGCTGTTCACCG            GGG
    ##   spacer_5     EGFP        35      + | CGAGGAGCTGTTCACCGGGG            TGG
    ##   spacer_6     EGFP        42      - | CCGTCCAGCTCGACCAGGAT            GGG
    ##             pam_site  cut_site      region
    ##            <numeric> <numeric> <character>
    ##   spacer_1        30        27        EGFP
    ##   spacer_2        31        28        EGFP
    ##   spacer_3        31        34        EGFP
    ##   spacer_4        32        29        EGFP
    ##   spacer_5        35        32        EGFP
    ##   spacer_6        42        45        EGFP
    ##   -------
    ##   seqinfo: 1 sequence from custom genome
    ##   crisprNuclease: SpCas9

The resulting output is a regular `GuideSet` object, and all
functionalities described in the [introductory CRISPRko
tutorial](https://github.com/crisprVerse/Tutorials/tree/master/Design_CRISPRko_Cas9)
can be applied here as well.

There are a few key differences to note with respect to a `GuideSet`
object constructed using a reference genome. First, the name of the
input DNA sequence (EGFP) is used as the chromosome name stored in the
`seqnames` field. Second, the `pam_site` and `cut_site` coordinates are
all relative to the first nucleotide of the custom DNA sequence.
Finally, the `GuideSet` object stores the input sequence, which can be
accessed using the function `customSequences`:

``` r
customSequences(gs)
```

    ## DNAStringSet object of length 1:
    ##     width seq                                               names               
    ## [1]   720 ATGGTGAGCAAGGGCGAGGAGCT...GCATGGACGAGCTGTACAAGTAA EGFP

## Finding off-targets in the human genome to find gRNAs specific to EGFP

Now that we have designed all possible gRNAs targeting EGFP, we will
filter out gRNAs that have on- and off-targets located in the human
genome. We will use the bowtie aligner to find targets, so we need to
first specify the path of a bowtie index constructed on the human
genome:

``` r
# Path of the hg38 bowtie index on my personal laptop:
bowtie_index <- "/Users/fortinj2/crisprIndices/bowtie/hg38/hg38"
```

For instructions on how to build a Bowtie index from a given reference
genome, see the [genome index
tutorial](https://github.com/crisprVerse/Tutorials/tree/master/Building_Genome_Indices)
or the [crisprBowtie page](https://github.com/Jfortin1/crisprBowtie) .

To annotate off-targets with genomic context, for instance to know
whether or not they are located in coding regions, we will also need a
gene model object. We will use the gene model object `txdb_human` from
`crisprDesignData`, which contains genomic coordinates of all human
protein-coding genes. See the [crisprDesignData
package](https://github.com/crisprVerse/crisprDesignData) for more
details.

``` r
data(txdb_human, package="crisprDesignData")
```

We are now ready to find all on- and off-targets using the
`addSpacerAlignments` function from `crisprDesign`:

``` r
gs <- addSpacerAlignments(gs,
                          aligner="bowtie",
                          aligner_index=bowtie_index,
                          bsgenome=BSgenome.Hsapiens.UCSC.hg38,
                          n_mismatches=3,
                          txObject=txdb_human)
```

    ## [runCrisprBowtie] Using BSgenome.Hsapiens.UCSC.hg38 
    ## [runCrisprBowtie] Searching for SpCas9 protospacers

``` r
gs
```

    ## GuideSet object with 119 ranges and 14 metadata columns:
    ##              seqnames    ranges strand |          protospacer            pam
    ##                 <Rle> <IRanges>  <Rle> |       <DNAStringSet> <DNAStringSet>
    ##     spacer_1     EGFP        30      + | AAGGGCGAGGAGCTGTTCAC            CGG
    ##     spacer_2     EGFP        31      + | AGGGCGAGGAGCTGTTCACC            GGG
    ##     spacer_3     EGFP        31      - | GACCAGGATGGGCACCACCC            CGG
    ##     spacer_4     EGFP        32      + | GGGCGAGGAGCTGTTCACCG            GGG
    ##     spacer_5     EGFP        35      + | CGAGGAGCTGTTCACCGGGG            TGG
    ##          ...      ...       ...    ... .                  ...            ...
    ##   spacer_115     EGFP       684      + | CTGGAGTTCGTGACCGCCGC            CGG
    ##   spacer_116     EGFP       685      + | TGGAGTTCGTGACCGCCGCC            GGG
    ##   spacer_117     EGFP       685      - | GTCCATGCCGAGAGTGATCC            CGG
    ##   spacer_118     EGFP       696      + | ACCGCCGCCGGGATCACTCT            CGG
    ##   spacer_119     EGFP       701      + | CGCCGGGATCACTCTCGGCA            TGG
    ##               pam_site  cut_site      region        n0        n1        n2
    ##              <numeric> <numeric> <character> <numeric> <numeric> <numeric>
    ##     spacer_1        30        27        EGFP         0         0         1
    ##     spacer_2        31        28        EGFP         0         0         2
    ##     spacer_3        31        34        EGFP         0         0         1
    ##     spacer_4        32        29        EGFP         0         0         0
    ##     spacer_5        35        32        EGFP         0         0         0
    ##          ...       ...       ...         ...       ...       ...       ...
    ##   spacer_115       684       681        EGFP         0         0         1
    ##   spacer_116       685       682        EGFP         0         0         0
    ##   spacer_117       685       688        EGFP         0         0         0
    ##   spacer_118       696       693        EGFP         0         0         0
    ##   spacer_119       701       698        EGFP         0         0         0
    ##                     n3      n0_c      n1_c      n2_c      n3_c
    ##              <numeric> <numeric> <numeric> <numeric> <numeric>
    ##     spacer_1         7         0         0         0         2
    ##     spacer_2         6         0         0         0         0
    ##     spacer_3        23         0         0         0         3
    ##     spacer_4         6         0         0         0         1
    ##     spacer_5         4         0         0         0         2
    ##          ...       ...       ...       ...       ...       ...
    ##   spacer_115         5         0         0         0         0
    ##   spacer_116         2         0         0         0         0
    ##   spacer_117         2         0         0         0         0
    ##   spacer_118         2         0         0         0         1
    ##   spacer_119         0         0         0         0         0
    ##                                                          alignments
    ##                                                       <GRangesList>
    ##     spacer_1   chr4:151052480:-,chr1:29350153:-,chrX:44232571:-,...
    ##     spacer_2  chr6:115711229:+,chr6:52987132:+,chr3:186748384:+,...
    ##     spacer_3 chr6:149167095:+,chr17:37844933:-,chr17:82019964:-,...
    ##     spacer_4 chr17:82484917:+,chr11:35064010:-,chr18:48539310:-,...
    ##     spacer_5    chr4:426282:-,chr19:16847037:+,chr19:14471687:+,...
    ##          ...                                                    ...
    ##   spacer_115  chr2:142494626:+,chr1:54281906:+,chr2:117131280:-,...
    ##   spacer_116                      chr4:139453357:-,chr17:82078529:+
    ##   spacer_117                        chr1:53536150:+,chr7:24540197:+
    ##   spacer_118                       chr18:77345044:+,chr3:51662810:-
    ##   spacer_119                                                       
    ##   -------
    ##   seqinfo: 1 sequence from custom genome
    ##   crisprNuclease: SpCas9

## Predicting on-target activity

We also want to make sure to filter out gRNAs that are predicted to have
poor on-target activity. To do so, we annotate gRNAs with the DeepHF
on-target activity score:

``` r
gs <- addOnTargetScores(gs, methods="deephf")
```

    ## [addOnTargetScores] Adding deephf scores.

    ## snapshotDate(): 2022-04-26

    ## see ?crisprScoreData and browseVignettes('crisprScoreData') for documentation

    ## loading from cache

Finally, we characterize the spacer sequences using the
`addSequenceFeatures` function from `crisprDesign`:

``` r
gs <- addSequenceFeatures(gs)
```

## Final selection

For our use case, we will only retain gRNAs that do not map to the human
genome (`n0=0`), don’t have any 1 or 2-mismatch off-targets (`n1=0` and
`n2=0`), and do not have 3-mismatch off-targets located in coding
regions (`n3_c=0`):

``` r
gs <- gs[gs$n0==0 & gs$n1==0 & gs$n2==0 & gs$n3_c==0]
```

We also remove gRNAs that contain polyT sequences

``` r
gs <- gs[!gs$polyT,]
```

and only keep gRNAs that don’t have extreme GC content:

``` r
gs <- gs[gs$percentGC>=20 & gs$percentGC<=80]
```

Finally, we rank gRNAs from the highest to the lowest on-target activity
score:

``` r
gs <- gs[order(-gs$score_deephf)]
head(gs)
```

    ## GuideSet object with 6 ranges and 21 metadata columns:
    ##              seqnames    ranges strand |          protospacer            pam
    ##                 <Rle> <IRanges>  <Rle> |       <DNAStringSet> <DNAStringSet>
    ##    spacer_64     EGFP       359      + | GAAGTTCGAGGGCGACACCC            TGG
    ##   spacer_114     EGFP       682      - | CATGCCGAGAGTGATCCCGG            CGG
    ##    spacer_77     EGFP       446      - | CGGCCATGATATAGACGTTG            TGG
    ##    spacer_16     EGFP        98      + | CAAGTTCAGCGTGTCCGGCG            AGG
    ##    spacer_45     EGFP       229      - | GTCGTGCTGCTTCATGTGGT            CGG
    ##   spacer_107     EGFP       635      - | TGTGATCGCGCTTCTCGTTG            GGG
    ##               pam_site  cut_site      region        n0        n1        n2
    ##              <numeric> <numeric> <character> <numeric> <numeric> <numeric>
    ##    spacer_64       359       356        EGFP         0         0         0
    ##   spacer_114       682       685        EGFP         0         0         0
    ##    spacer_77       446       449        EGFP         0         0         0
    ##    spacer_16        98        95        EGFP         0         0         0
    ##    spacer_45       229       232        EGFP         0         0         0
    ##   spacer_107       635       638        EGFP         0         0         0
    ##                     n3      n0_c      n1_c      n2_c      n3_c
    ##              <numeric> <numeric> <numeric> <numeric> <numeric>
    ##    spacer_64         0         0         0         0         0
    ##   spacer_114         0         0         0         0         0
    ##    spacer_77         1         0         0         0         0
    ##    spacer_16         3         0         0         0         0
    ##    spacer_45         9         0         0         0         0
    ##   spacer_107         2         0         0         0         0
    ##                                                       alignments score_deephf
    ##                                                    <GRangesList>    <numeric>
    ##    spacer_64                                                         0.716188
    ##   spacer_114                                                         0.700199
    ##    spacer_77                                    chr3:140247159:+     0.686111
    ##    spacer_16    chr5:132092368:+,chr6:42782483:+,chrX:97074134:+     0.670439
    ##    spacer_45 chr7:87066799:+,chr4:89369669:-,chrX:82007103:+,...     0.664180
    ##   spacer_107                    chr2:128980117:+,chr1:53729319:-     0.654151
    ##              percentGC     polyA     polyC     polyG     polyT startingGGGGG
    ##              <numeric> <logical> <logical> <logical> <logical>     <logical>
    ##    spacer_64        65     FALSE     FALSE     FALSE     FALSE         FALSE
    ##   spacer_114        65     FALSE     FALSE     FALSE     FALSE         FALSE
    ##    spacer_77        50     FALSE     FALSE     FALSE     FALSE         FALSE
    ##    spacer_16        65     FALSE     FALSE     FALSE     FALSE         FALSE
    ##    spacer_45        55     FALSE     FALSE     FALSE     FALSE         FALSE
    ##   spacer_107        55     FALSE     FALSE     FALSE     FALSE         FALSE
    ##   -------
    ##   seqinfo: 1 sequence from custom genome
    ##   crisprNuclease: SpCas9

Users can select the top gRNAs as their control gRNAs.

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
    ##  [9] GenomicRanges_1.48.0              crisprDesignData_0.99.14         
    ## [11] crisprDesign_0.99.130             crisprBase_1.1.5                 
    ## [13] Biostrings_2.64.0                 GenomeInfoDb_1.32.2              
    ## [15] XVector_0.35.0                    IRanges_2.30.0                   
    ## [17] S4Vectors_0.33.11                 BiocGenerics_0.42.0              
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
    ## [61] BiocVersion_3.15.0            BiocIO_1.5.0                 
    ## [63] randomForest_4.7-1            GenomicFeatures_1.47.13      
    ## [65] filelock_1.0.2                BiocParallel_1.29.18         
    ## [67] rlang_1.0.4                   pkgconfig_2.0.3              
    ## [69] matrixStats_0.61.0            bitops_1.0-7                 
    ## [71] evaluate_0.15                 lattice_0.20-45              
    ## [73] purrr_0.3.4                   GenomicAlignments_1.31.2     
    ## [75] bit_4.0.4                     tidyselect_1.1.2             
    ## [77] magrittr_2.0.2                R6_2.5.1                     
    ## [79] generics_0.1.2                DelayedArray_0.21.2          
    ## [81] DBI_1.1.2                     pillar_1.7.0                 
    ## [83] withr_2.5.0                   KEGGREST_1.35.0              
    ## [85] RCurl_1.98-1.6                tibble_3.1.6                 
    ## [87] dir.expiry_1.3.0              crayon_1.5.0                 
    ## [89] utf8_1.2.2                    tzdb_0.2.0                   
    ## [91] rmarkdown_2.13                progress_1.2.2               
    ## [93] grid_4.2.0                    blob_1.2.2                   
    ## [95] digest_0.6.29                 xtable_1.8-4                 
    ## [97] httpuv_1.6.5                  Rbwa_1.1.0

# References
