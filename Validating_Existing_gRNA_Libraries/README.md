Validating existing gRNA libraries
================

-   [Introduction](#introduction)
-   [Loading necessary packages](#loading-necessary-packages)
-   [Reading in data](#reading-in-data)
-   [Building a `GuideSet` object](#building-a-guideset-object)
-   [Off-target characterization](#off-target-characterization)
-   [Session Info](#session-info)

Authors: Jean-Philippe Fortin

Date: July 8, 2022

# Introduction

In this vignette, we characterize a small mouse CRISPR knockout
(CRISPRko) library that was designed to target tumor suppressors. The
library was obtained from Addgene, and is stored in `.\extdata`.

# Loading necessary packages

``` r
library(crisprDesign)
library(crisprBowtie)
library(crisprBase)
library(readxl)
library(BSgenome.Mmusculus.UCSC.mm10)
bsgenome <- BSgenome.Mmusculus.UCSC.mm10
```

We also load `crisprDesignData`, which is a data package containing
already-processed Ensembl objects for gene annotation of human and mouse
gRNAs:

``` r
library(crisprDesignData)
```

# Reading in data

``` r
data <- read_excel("extdata/mtsg-grnas-readcounts.xlsx")
data <- as.data.frame(data)[,1:2]
colnames(data) <- c("ID", "spacer_20mer")

# Getting genes names:
data$gene_symbol <- sapply(strsplit(data$ID, split="_"), function(x)x[[1]])
head(data)
```

    ##           ID         spacer_20mer gene_symbol
    ## 1   Fat1_sg1 GGGCAGTGTTTCAAAATCCA        Fat1
    ## 2   Fat1_sg2 GGAACACGAGCCGTCAGCGG        Fat1
    ## 3   Fat1_sg3 GGATTTCTGTTCTGCATCAA        Fat1
    ## 4   Fat1_sg4 GGTCCCATCTGTTGCCTCCA        Fat1
    ## 5   Fat1_sg5 GTTTGGAGATCCACTCGATA        Fat1
    ## 6 Arid1b_sg1 GTACCCAGTGCAAGCTACAG      Arid1b

# Building a `GuideSet` object

We first define the nuclease for the analysis. We here use the standard
wildtype Cas9 (SpCas9) from the `crisprBase` package:

``` r
data(SpCas9, package="crisprBase")
crisprNuclease <- SpCas9
crisprNuclease
```

    ## Class: CrisprNuclease
    ##   Name: SpCas9
    ##   Target type: DNA
    ##   Metadata: list of length 1
    ##   PAMs: NGG, NAG, NGA
    ##   Weights: 1, 0.2593, 0.0694
    ##   Spacer length: 20
    ##   PAM side: 3prime
    ##     Distance from PAM: 0
    ##   Prototype protospacers: 5'--SSSSSSSSSSSSSSSSSSSS[NGG]--3', 5'--SSSSSSSSSSSSSSSSSSSS[NAG]--3', 5'--SSSSSSSSSSSSSSSSSSSS[NGA]--3'

The default length of the spacer sequences is 20nt. This can be changed
to a different length if needed, for instance 19nt:

``` r
# Not sunr
spacerLength(SpCas9) <- 19
```

We next need to define a bowtie index that we will use for alignment:

``` r
bowtie_index <- "/Users/fortinj2/crisprIndices/bowtie/mm10/mm10"
```

We first map the gRNAs to the reference genome with perfect match to
obtain genomic coordinates:

``` r
spacers <- unique(data$spacer_20mer)
aln <- runCrisprBowtie(spacers,
                       crisprNuclease=crisprNuclease,
                       bowtie_index=bowtie_index,
                       n_mismatches=0)
```

    ## [runCrisprBowtie] Searching for SpCas9 protospacers

``` r
head(aln)
```

    ##                 spacer          protospacer pam                chr  pam_site
    ## 1 GAAAACAGCCAAGGTTTGTA GAAAACAGCCAAGGTTTGTA CGG               chr8 106659780
    ## 2 GAAAACCCTGAAGTGCCCAC GAAAACCCTGAAGTGCCCAC GGG              chr17  29099403
    ## 3 GAAAACCCTGAAGTGCCCAC GAAAACCCTGAAGTGCCCAC GGG chr17_JH584267_alt   1798300
    ## 4 GAAAAGGGAAGACCAGCCCC GAAAAGGGAAGACCAGCCCC TGG               chr3 152219735
    ## 5 GAACCGACAAACAGTCCTGG GAACCGACAAACAGTCCTGG AGG              chr14  31255196
    ## 6 GAACCTAGATTTTGAGACAG GAACCTAGATTTTGAGACAG GGG              chr17  33952329
    ##   strand n_mismatches canonical
    ## 1      +            0      TRUE
    ## 2      +            0      TRUE
    ## 3      +            0      TRUE
    ## 4      +            0      TRUE
    ## 5      +            0      TRUE
    ## 6      +            0      TRUE

`n_mismatches=0` specifies that we require a perfect match between
spacer and protospacer sequences (on-targets).

Non-targeting controls should not have any alignments to the genome, and
some guides might have multiple alignments if they were not designed
carefully. For such guides, that’s OK, we can pick up pick one genomic
coordinate for now, and the multiple alignments annotation will be
handled later on.

We keep only alignments to the standard chromosomes:

``` r
chrs <- paste0("chr",c(1:22, "X", "Y"))
aln <- aln[aln$chr %in% chrs,,drop=FALSE]
```

We add the genomic coordinates to the data.frame:

``` r
wh <- match(data$spacer_20mer, aln$spacer)
data$chr <- aln$chr[wh]
data$pam_site <- aln$pam_site[wh]
data$pam <- aln$pam[wh]
data$strand <- aln$strand[wh]
head(data)
```

    ##           ID         spacer_20mer gene_symbol   chr pam_site pam strand
    ## 1   Fat1_sg1 GGGCAGTGTTTCAAAATCCA        Fat1  chr8 45023166 AGG      -
    ## 2   Fat1_sg2 GGAACACGAGCCGTCAGCGG        Fat1  chr8 45024178 TGG      +
    ## 3   Fat1_sg3 GGATTTCTGTTCTGCATCAA        Fat1  chr8 45013029 GGG      -
    ## 4   Fat1_sg4 GGTCCCATCTGTTGCCTCCA        Fat1  chr8 45013065 CGG      -
    ## 5   Fat1_sg5 GTTTGGAGATCCACTCGATA        Fat1  chr8 45010453 TGG      -
    ## 6 Arid1b_sg1 GTACCCAGTGCAAGCTACAG      Arid1b chr17  5291008 CGG      +

We can now build a proper `GuideSet` object in `crisprDesign` that will
allow us to do (more) sophisticated analyses.

We need to filter out first guides that don’t have a match to the
genome:

``` r
data <- data[!is.na(data$pam_site),,drop=FALSE]
```

Ready to build the object with the constructor function in crisprDesign:

``` r
gs <- GuideSet(protospacers=data$spacer_20mer,
               pams=data$pam,
               pam_site=data$pam_site,
               seqnames=data$chr,
               strand=data$strand,
               CrisprNuclease=crisprNuclease,
               bsgenome=bsgenome)
gs$gene_symbol <- data$gene_symbol
```

While it is not necessary, we will unique names to the spacers:

``` r
names(gs) <- paste0("gRNA_", seq_along(gs))
```

The `GuideSet` object, and `crisprDesign`, provide rich functionalities
to annotate and manipulate gRNAs; see the main `crisprDesign` vignette
for more details. For the rest of this vignette, we only focus on
characterizing the off-targets.

# Off-target characterization

Having a `GuideSet` object, it is now a piece of cake to characterize
the off-targets. We characterize off-targets using the bowtie aligner,
with up to 3 mismatches between the spacer (gRNA) and protospacer
(target DNA) sequences. The function `addSpacerAlignments` accomplishes
that.

It has an optional argument `txObject` that can be used to provide gene
model data to put the off-targets in a gene model context. We made such
objects available for human and mouse in the `crisprDesignData` package
(see `txdb_human` and \`txdb_mouse).

``` r
data(txdb_mouse, package="crisprDesignData")
txObject <- txdb_mouse
gs <- addSpacerAlignments(gs,
                          txObject=txObject,
                          aligner="bowtie",
                          aligner_index=bowtie_index,
                          bsgenome=bsgenome,
                          n_mismatches=2)
```

    ## Loading required namespace: crisprBwa

    ## [runCrisprBowtie] Using BSgenome.Mmusculus.UCSC.mm10 
    ## [runCrisprBowtie] Searching for SpCas9 protospacers

The alignments are stored in a metadata column called `alignments`. See
?getSpacerAlignments for more details about what the different columns
are.

As an example, we can access the on- and off-target alignments of the
first gRNA using the following:

``` r
aln <- gs$alignments[[1]]
aln
```

    ## GRanges object with 3 ranges and 14 metadata columns:
    ##            seqnames    ranges strand |               spacer
    ##               <Rle> <IRanges>  <Rle> |       <DNAStringSet>
    ##     gRNA_1     chr8  45023166      - | GGGCAGTGTTTCAAAATCCA
    ##    gRNA_97     chr9 104311568      + | GGCAAAACAGAATAGGAGAG
    ##   gRNA_218     chr1 113261571      + | GAACCTAGATTTTGAGACAG
    ##                     protospacer            pam  pam_site n_mismatches canonical
    ##                  <DNAStringSet> <DNAStringSet> <numeric>    <integer> <logical>
    ##     gRNA_1 GGGCAGTGTTTCAAAATCCA            AGG  45023166            0      TRUE
    ##    gRNA_97 GGCAAAACAGAATATGAGTG            GGG 104311568            2      TRUE
    ##   gRNA_218 GAACCTAACTTTTGAGACAG            AGG 113261571            2      TRUE
    ##             cut_site         cds    fiveUTRs   threeUTRs       exons
    ##            <numeric> <character> <character> <character> <character>
    ##     gRNA_1  45023169        Fat1        <NA>        <NA>        Fat1
    ##    gRNA_97 104311565        <NA>        <NA>        <NA>        <NA>
    ##   gRNA_218 113261568        <NA>        <NA>        <NA>        <NA>
    ##                introns  intergenic intergenic_distance
    ##            <character> <character>           <integer>
    ##     gRNA_1        <NA>        <NA>                <NA>
    ##    gRNA_97        Acpp        <NA>                <NA>
    ##   gRNA_218        <NA>     Gm28189               27771
    ##   -------
    ##   seqinfo: 22 sequences (1 circular) from mm10 genome

We can also add CFD and MIT scores to the off-targets to characterize
the likelihood of SpCas9 to cut at the off-targets:

``` r
gs <- addOffTargetScores(gs)
```

The scores range from 0 to 1, and a higher score indicates a higher
probability of the off-target to occur.

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
    ##  [1] crisprDesignData_0.99.2            BSgenome.Mmusculus.UCSC.mm10_1.4.3
    ##  [3] BSgenome_1.63.5                    rtracklayer_1.55.4                
    ##  [5] Biostrings_2.63.2                  XVector_0.35.0                    
    ##  [7] GenomicRanges_1.47.6               GenomeInfoDb_1.31.6               
    ##  [9] IRanges_2.29.1                     S4Vectors_0.33.11                 
    ## [11] BiocGenerics_0.41.2                readxl_1.3.1                      
    ## [13] crisprBowtie_1.1.1                 crisprDesign_0.99.93              
    ## [15] crisprBase_1.1.2                  
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] rjson_0.2.21                  ellipsis_0.3.2               
    ##   [3] Rbowtie_1.35.0                rstudioapi_0.13              
    ##   [5] bit64_4.0.5                   interactiveDisplayBase_1.33.0
    ##   [7] AnnotationDbi_1.57.1          fansi_1.0.2                  
    ##   [9] xml2_1.3.3                    cachem_1.0.6                 
    ##  [11] knitr_1.37                    jsonlite_1.8.0               
    ##  [13] Rsamtools_2.11.0              dbplyr_2.1.1                 
    ##  [15] png_0.1-7                     shiny_1.7.1                  
    ##  [17] BiocManager_1.30.16           readr_2.1.2                  
    ##  [19] compiler_4.2.0                httr_1.4.2                   
    ##  [21] basilisk_1.9.1                assertthat_0.2.1             
    ##  [23] Matrix_1.4-0                  fastmap_1.1.0                
    ##  [25] cli_3.3.0                     later_1.3.0                  
    ##  [27] htmltools_0.5.2               prettyunits_1.1.1            
    ##  [29] tools_4.2.0                   glue_1.6.2                   
    ##  [31] GenomeInfoDbData_1.2.7        dplyr_1.0.8                  
    ##  [33] rappdirs_0.3.3                Rcpp_1.0.8.3                 
    ##  [35] Biobase_2.55.0                cellranger_1.1.0             
    ##  [37] vctrs_0.3.8                   ExperimentHub_2.3.5          
    ##  [39] crisprBwa_1.1.2               crisprScore_1.1.6            
    ##  [41] xfun_0.30                     stringr_1.4.0                
    ##  [43] mime_0.12                     lifecycle_1.0.1              
    ##  [45] restfulr_0.0.13               XML_3.99-0.9                 
    ##  [47] AnnotationHub_3.3.9           zlibbioc_1.41.0              
    ##  [49] basilisk.utils_1.5.0          vroom_1.5.7                  
    ##  [51] VariantAnnotation_1.41.3      hms_1.1.1                    
    ##  [53] promises_1.2.0.1              MatrixGenerics_1.7.0         
    ##  [55] parallel_4.2.0                SummarizedExperiment_1.25.3  
    ##  [57] yaml_2.3.5                    curl_4.3.2                   
    ##  [59] memoise_2.0.1                 reticulate_1.24              
    ##  [61] biomaRt_2.51.3                stringi_1.7.6                
    ##  [63] RSQLite_2.2.12                BiocVersion_3.15.0           
    ##  [65] BiocIO_1.5.0                  randomForest_4.7-1           
    ##  [67] crisprScoreData_1.1.3         GenomicFeatures_1.47.13      
    ##  [69] filelock_1.0.2                BiocParallel_1.29.18         
    ##  [71] rlang_1.0.2                   pkgconfig_2.0.3              
    ##  [73] matrixStats_0.61.0            bitops_1.0-7                 
    ##  [75] evaluate_0.15                 lattice_0.20-45              
    ##  [77] purrr_0.3.4                   GenomicAlignments_1.31.2     
    ##  [79] bit_4.0.4                     tidyselect_1.1.2             
    ##  [81] magrittr_2.0.2                R6_2.5.1                     
    ##  [83] generics_0.1.2                DelayedArray_0.21.2          
    ##  [85] DBI_1.1.2                     pillar_1.7.0                 
    ##  [87] KEGGREST_1.35.0               RCurl_1.98-1.6               
    ##  [89] tibble_3.1.6                  dir.expiry_1.3.0             
    ##  [91] crayon_1.5.0                  utf8_1.2.2                   
    ##  [93] BiocFileCache_2.3.4           tzdb_0.2.0                   
    ##  [95] rmarkdown_2.13                progress_1.2.2               
    ##  [97] grid_4.2.0                    blob_1.2.2                   
    ##  [99] digest_0.6.29                 xtable_1.8-4                 
    ## [101] httpuv_1.6.5                  Rbwa_1.1.0
