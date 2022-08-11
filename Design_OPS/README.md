Using crisprDesign to design gRNAs for optical pooled screening (OPS)
================

-   [Introduction](#introduction)
-   [Installation](#installation)
-   [Terminology](#terminology)
-   [Design for optical pooled screening
    (OPS)](#design-for-optical-pooled-screening-ops)
    -   [Loading packages](#loading-packages)
    -   [Creating the GuideSet](#creating-the-guideset)
    -   [Adding OPS barcodes](#adding-ops-barcodes)
    -   [Barcode distance matrix](#barcode-distance-matrix)
    -   [Designing OPS libraries](#designing-ops-libraries)
    -   [Adding gRNAs to an existing OPS
        library](#adding-grnas-to-an-existing-ops-library)
-   [Session Info](#session-info)
-   [References](#references)

Authors: Jean-Philippe Fortin, Luke Hoberecht

Date: 11 August, 2022

# Introduction

Optical pooled screening (OPS) combines image-based sequencing (in situ
sequencing) of gRNAs and optical phenotyping on the same physical wells
(Feldman et al. 2019). In such experiments, guide RNA (gRNA) spacer
sequences are partially sequenced from the 5-prime end; the length of
these truncated sequences, or barcodes, which corresponds to the number
of sequencing cycles, is fixed and chosen by the experimentalist. From a
gRNA design perspective, additional constraints are needed to ensure
sufficient dissimilarity between the rruncated barcodes for their
identification during the analysis.

This tutorial will demonstrate how to design gRNAs for use in optical
pooled screens, with emphasis on the constraints described above. Common
gRNA design steps that are not specific to OPS are omitted in this
tutorial (e.g. off-target search, or on-target activity prediction)
here. Users can peruse through the list of [available
tutorials](https://github.com/crisprVerse/Tutorials) for more
information regarind application-specific gRNA design rules.

# Installation

The Bioconductor packages needed in this tutorial can be installed using
the following commands:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("crisprDesign")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
```

The GitHub packages needed in this tutorial can be installed using the
following commands:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("crisprVerse/crisprDesignData")
```

# Terminology

See the [CRISPRko design
vignette](https://github.com/crisprVerse/Tutorials/tree/master/Design_CRISPRko_Cas9)
to get familiar with the terminology used throughout this tutorial.

# Design for optical pooled screening (OPS)

To illustrate the functionalities of `crisprDesign` for designing OPS
libraries, we will design a small CRISPRko OPS library targeting 3 genes
of the human RAS family: KRAS, HRAS, and NRAS. We will use the SpCas9
nuclease.

We will design gRNAs for an experiment that uses 8 in situ sequencing
cycles:

``` r
n_cycles=8
```

## Loading packages

Before we start, we first load the necessary packages for this tutorial:

``` r
library(crisprBase)
library(crisprDesign)
library(crisprDesignData)
library(BSgenome.Hsapiens.UCSC.hg38)
```

## Creating the GuideSet

We begin by loading the SpCas9 `CrisprNuclease` object from the
`crisprBase` package

``` r
data(SpCas9, package="crisprBase")
```

as well as data containing gene regions for the human genome:

``` r
data(txdb_human, package="crisprDesignData")
```

For more information on `txdb_human` and how to create similar gene
annotation objects, see the [Building a gene annotation
object](https://github.com/crisprVerse/Tutorials/tree/master/Building_Gene_Annotation)
tutorial.

Next, we find the CDS coordinates for our genes using the
`queryTxObject` function:

``` r
target_genes <- c("KRAS", "HRAS", "NRAS")
target_regions <- queryTxObject(txdb_human,
                                featureType="cds",
                                queryColumn="gene_symbol",
                                queryValue=target_genes)
```

then build our `GuideSet` with the `findSpacers` function:

``` r
gs <- findSpacers(target_regions,
                  crisprNuclease=SpCas9,
                  bsgenome=BSgenome.Hsapiens.UCSC.hg38)
```

As we will want to distinguish which gene each spacer targets, we will
add `gene_symbol` and `gene_id` columns from `target_regions`.

``` r
gene_info <- target_regions[gs$region]
gs$gene_symbol <- gene_info$gene_symbol
gs$gene_id <- gene_info$gene_id
```

## Adding OPS barcodes

We can add our OPS barcodes to the GuideSet with the `addOpsBarcodes`
function. This function extracts the `n_cycles` nucleotides from the
5-prime end of our spacers and stores them in the `opsBarcode` column:

``` r
gs <- addOpsBarcodes(gs,
                     n_cycles=n_cycles)
head(gs$opsBarcode)
```

    ## DNAStringSet object of length 6:
    ##     width seq                                               names               
    ## [1]     8 CCATGTGT                                          spacer_1
    ## [2]     8 TTGTATGG                                          spacer_2
    ## [3]     8 CCTTGTTA                                          spacer_3
    ## [4]     8 GATGGGAC                                          spacer_4
    ## [5]     8 TGATGGGA                                          spacer_5
    ## [6]     8 AGCAGTGA                                          spacer_6

## Barcode distance matrix

We can pass our barcodes to the function `getBarcodeDistanceMatrix` to
calculate the nucleotide distance between them. The `dist_method`
argument determines the type of distance to calculate: `"hamming"`,
which only considers substitutions (default) or `"levenstein"`, which
also allows for insertions and deletions.

As a brief demonstration, let’s look at the distances between the first
few barcodes in our `GuideSet`. We set the `binarize` argument (more on
this parameter later) to `FALSE` to show distances:

``` r
barcodes <- gs$opsBarcode
dist <- getBarcodeDistanceMatrix(barcodes[1:5],
                                 binnarize=FALSE)
dist
```

    ## 5 x 5 sparse Matrix of class "dsCMatrix"
    ##          CCATGTGT TTGTATGG CCTTGTTA GATGGGAC TGATGGGA
    ## CCATGTGT        .        5        3        7        4
    ## TTGTATGG        5        .        6        8        5
    ## CCTTGTTA        3        6        .        6        5
    ## GATGGGAC        7        8        6        .        6
    ## TGATGGGA        4        5        5        6        .

Note that the output is a sparse matrix, so the barcodes along the
diagonal (i.e., compared against themselves) return `.`, or a distance
of zero. To compare one set of barcodes against another, we can pass the
other set to the `targetBarcodes` argument (the former barcode set being
passed to the `queryBarcodes` argument, which is compared against itself
when `targetBarcodes` is `NULL`):

``` r
dist <- getBarcodeDistanceMatrix(barcodes[1:5],
                                 targetBarcodes=barcodes[6:10],
                                 binnarize=FALSE)
dist
```

    ## 5 x 5 sparse Matrix of class "dgCMatrix"
    ##          AGCAGTGA CAGCAGTG AACTCAAC AAACTCAA TGCTGTTG
    ## CCATGTGT        5        7        7        7        5
    ## TTGTATGG        6        5        7        8        4
    ## CCTTGTTA        5        6        7        7        4
    ## GATGGGAC        7        6        5        6        7
    ## TGATGGGA        4        7        7        6        4

The question we are interested in with respect to barcode distances is
whether this distance is sufficiently dissimilar for accurate
identification of spacers during sequencing. This minimum distance edit
(`min_dist_edit`) relies on the accuracy of various steps in the
experiment. Suppose, as a conservative estimate, that we can expect no
more than two edits per barcode in our example. A `min_dist_edit` of `3`
should suffice. Setting the `binnarize` argument to `TRUE`, and passing
our minimum distance edit value to `min_dist_edit` will binarize the
output, flagging barcodes (with a value of `1`) that are too similar and
should not both be included in our library:

``` r
dist <- getBarcodeDistanceMatrix(barcodes[1:5],
                                 barcodes[6:10],
                                 binnarize=TRUE,
                                 min_dist_edit=3)
dist
```

    ## 5 x 5 sparse Matrix of class "dtCMatrix"
    ##          AGCAGTGA CAGCAGTG AACTCAAC AAACTCAA TGCTGTTG
    ## CCATGTGT        .        .        .        .        .
    ## TTGTATGG        .        .        .        .        .
    ## CCTTGTTA        .        .        .        .        .
    ## GATGGGAC        .        .        .        .        .
    ## TGATGGGA        .        .        .        .        .

Using this function with large sets of barcodes can be taxing on memory.
To manage this, it is recommended to set `splitByChunks=TRUE` and
specify the number of chunks with `n_chunks` (see
`?getBarcodeDistanceMatrix`).

## Designing OPS libraries

The `designOpsLibrary` function allows users to perform a complete
end-to-end OPS library design. We will design our library with 4 gRNAs
per gene using the `n_guides` and `gene_field` (to identify gRNAs by
gene target) parameters. We will also use the same distance method and
minimum distance edit parameters as in the example above.

NOTE: it is advised to first complete other steps in gRNA design
(annotating, filtering, and ranking gRNAs in the `GuideSet`) prior to
using this function; this will ensure the library contains the best
gRNAs. As this example did not rank gRNAs, we are notified that rankings
are assigned by the order in which gRNAs appear in our input.

``` r
df <- data.frame(ID=names(gs),
                 spacer=gs$protospacer,
                 opsBarcode=gs$opsBarcode,
                 gene_symbol=gs$gene_symbol)
opsLibrary <- designOpsLibrary(df,
                               n_guides=4,
                               gene_field="gene_symbol",
                               min_dist_edit=5,
                               dist_method="hamming")
```

    ## Since 'rank' column is not provided, using default order has ranking.

``` r
opsLibrary
```

    ##                  ID               spacer opsBarcode gene_symbol rank
    ## HRAS      spacer_73 ACTTGCAGCTCATGCAGCCG   ACTTGCAG        HRAS   10
    ## HRAS1     spacer_76 CTGAACCCTCCTGATGAGAG   CTGAACCC        HRAS   13
    ## HRAS2     spacer_79 CAGCCGGGGCCACTCTCATC   CAGCCGGG        HRAS   16
    ## HRAS3    spacer_131 TGGGTCACATGGGTCCCGGG   TGGGTCAC        HRAS   68
    ## KRAS     spacer_531 AAAGAAAAGATGAGCAAAGA   AAAGAAAA        KRAS    1
    ## KRAS1    spacer_533 TTCTCGAACTAATGTATAGA   TTCTCGAA        KRAS    3
    ## KRAS2    spacer_539 GGAGGATGCTTTTTATACAT   GGAGGATG        KRAS    9
    ## KRAS3    spacer_564 AACTCTTTTAATTTGTTCTC   AACTCTTT        KRAS   34
    ## NRAS       spacer_1 CCATGTGTGGTGATGTAACA   CCATGTGT        NRAS    1
    ## spacer_2   spacer_2 TTGTATGGGATTGCCATGTG   TTGTATGG        NRAS    2
    ## spacer_4   spacer_4 GATGGGACTCAGGGTTGTAT   GATGGGAC        NRAS    4
    ## NRAS1      spacer_6 AGCAGTGATGATGGGACTCA   AGCAGTGA        NRAS    6

## Adding gRNAs to an existing OPS library

Suppose we later wish to add another gene target to our library, but
also want to retain the gRNAs that are currently in our library. We can
append these additional gRNAs with the `updateOpsLibrary` function. This
function has the same parameters as `designOpsLibrary`, with an
additional `opsLibrary` argument to which we pass our original OPS
library.

To demonstrate, we will add the MRAS gene to our library. We first
construct the `GuideSet` for MRAS:

``` r
target_region <- queryTxObject(txdb_human,
                               featureType="cds",
                               queryColumn="gene_symbol",
                               queryValue="MRAS")
gs_mras <- findSpacers(target_region,
                       crisprNuclease=SpCas9,
                       bsgenome=BSgenome.Hsapiens.UCSC.hg38)
gs_mras$gene_symbol <- "MRAS"
gs_mras$gene_id <- "ENSG00000158186"
```

then add barcodes and construct the `data.frame`:

``` r
## add OPS barcodes
gs_mras <- addOpsBarcodes(gs_mras,
                          n_cycles=n_cycles)

## construct data.frame
df_mras <- data.frame(ID=names(gs_mras),
                      spacer=gs_mras$protospacer,
                      opsBarcode=gs_mras$opsBarcode,
                      gene_symbol=gs_mras$gene_symbol)
```

which we then pass with our other parameters to `updateOpsLibrary`:

``` r
opsLibrary <- updateOpsLibrary(opsLibrary,
                               df_mras,
                               n_guides=4,
                               gene_field="gene_symbol",
                               min_dist_edit=5,
                               dist_method="hamming")
```

    ## Since 'rank' column is not provided, using default order has ranking.

``` r
opsLibrary
```

    ##                  ID               spacer opsBarcode gene_symbol rank
    ## HRAS      spacer_73 ACTTGCAGCTCATGCAGCCG   ACTTGCAG        HRAS   10
    ## HRAS1     spacer_76 CTGAACCCTCCTGATGAGAG   CTGAACCC        HRAS   13
    ## HRAS2     spacer_79 CAGCCGGGGCCACTCTCATC   CAGCCGGG        HRAS   16
    ## HRAS3    spacer_131 TGGGTCACATGGGTCCCGGG   TGGGTCAC        HRAS   68
    ## KRAS     spacer_531 AAAGAAAAGATGAGCAAAGA   AAAGAAAA        KRAS    1
    ## KRAS1    spacer_533 TTCTCGAACTAATGTATAGA   TTCTCGAA        KRAS    3
    ## KRAS2    spacer_539 GGAGGATGCTTTTTATACAT   GGAGGATG        KRAS    9
    ## KRAS3    spacer_564 AACTCTTTTAATTTGTTCTC   AACTCTTT        KRAS   34
    ## MRAS       spacer_4 GGGGAGGTTGTCACTGGGGA   GGGGAGGT        MRAS    4
    ## MRAS1     spacer_19 CCACCACCAGCTTGTATGTG   CCACCACC        MRAS   19
    ## MRAS2     spacer_34 CCCCACATACAAGCTGGTGG   CCCCACAT        MRAS   34
    ## MRAS3     spacer_37 CACATACAAGCTGGTGGTGG   CACATACA        MRAS   37
    ## NRAS       spacer_1 CCATGTGTGGTGATGTAACA   CCATGTGT        NRAS    1
    ## spacer_2   spacer_2 TTGTATGGGATTGCCATGTG   TTGTATGG        NRAS    2
    ## spacer_4   spacer_4 GATGGGACTCAGGGTTGTAT   GATGGGAC        NRAS    4
    ## NRAS1      spacer_6 AGCAGTGATGATGGGACTCA   AGCAGTGA        NRAS    6

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
    ##  [1] BSgenome.Hsapiens.UCSC.hg38_1.4.4 BSgenome_1.64.0                  
    ##  [3] rtracklayer_1.55.4                Biostrings_2.64.0                
    ##  [5] XVector_0.35.0                    GenomicRanges_1.48.0             
    ##  [7] GenomeInfoDb_1.32.2               IRanges_2.30.0                   
    ##  [9] S4Vectors_0.33.11                 BiocGenerics_0.42.0              
    ## [11] crisprDesignData_0.99.13          crisprDesign_0.99.117            
    ## [13] crisprBase_1.1.5                 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bitops_1.0-7                  matrixStats_0.61.0           
    ##  [3] bit64_4.0.5                   filelock_1.0.2               
    ##  [5] progress_1.2.2                httr_1.4.2                   
    ##  [7] tools_4.2.0                   utf8_1.2.2                   
    ##  [9] R6_2.5.1                      DBI_1.1.2                    
    ## [11] tidyselect_1.1.2              prettyunits_1.1.1            
    ## [13] bit_4.0.4                     curl_4.3.2                   
    ## [15] compiler_4.2.0                crisprBowtie_1.1.1           
    ## [17] cli_3.3.0                     Biobase_2.55.0               
    ## [19] basilisk.utils_1.9.1          crisprScoreData_1.1.3        
    ## [21] xml2_1.3.3                    DelayedArray_0.21.2          
    ## [23] randomForest_4.7-1            readr_2.1.2                  
    ## [25] rappdirs_0.3.3                stringr_1.4.0                
    ## [27] digest_0.6.29                 Rsamtools_2.11.0             
    ## [29] rmarkdown_2.13                crisprScore_1.1.13           
    ## [31] basilisk_1.9.2                pkgconfig_2.0.3              
    ## [33] htmltools_0.5.2               MatrixGenerics_1.7.0         
    ## [35] dbplyr_2.1.1                  fastmap_1.1.0                
    ## [37] rlang_1.0.4                   rstudioapi_0.13              
    ## [39] RSQLite_2.2.12                shiny_1.7.1                  
    ## [41] BiocIO_1.5.0                  generics_0.1.2               
    ## [43] jsonlite_1.8.0                BiocParallel_1.29.18         
    ## [45] dplyr_1.0.8                   VariantAnnotation_1.41.3     
    ## [47] RCurl_1.98-1.6                magrittr_2.0.2               
    ## [49] GenomeInfoDbData_1.2.7        Matrix_1.4-0                 
    ## [51] Rcpp_1.0.8.3                  fansi_1.0.2                  
    ## [53] reticulate_1.25               Rbowtie_1.36.0               
    ## [55] lifecycle_1.0.1               stringi_1.7.6                
    ## [57] yaml_2.3.5                    SummarizedExperiment_1.25.3  
    ## [59] zlibbioc_1.41.0               BiocFileCache_2.3.4          
    ## [61] AnnotationHub_3.3.9           grid_4.2.0                   
    ## [63] blob_1.2.2                    promises_1.2.0.1             
    ## [65] parallel_4.2.0                ExperimentHub_2.3.5          
    ## [67] crayon_1.5.0                  dir.expiry_1.3.0             
    ## [69] lattice_0.20-45               GenomicFeatures_1.47.13      
    ## [71] hms_1.1.1                     KEGGREST_1.35.0              
    ## [73] knitr_1.37                    pillar_1.7.0                 
    ## [75] rjson_0.2.21                  biomaRt_2.51.3               
    ## [77] BiocVersion_3.15.0            XML_3.99-0.9                 
    ## [79] glue_1.6.2                    evaluate_0.15                
    ## [81] BiocManager_1.30.16           httpuv_1.6.5                 
    ## [83] png_0.1-7                     vctrs_0.3.8                  
    ## [85] tzdb_0.2.0                    purrr_0.3.4                  
    ## [87] assertthat_0.2.1              cachem_1.0.6                 
    ## [89] xfun_0.30                     mime_0.12                    
    ## [91] xtable_1.8-4                  restfulr_0.0.13              
    ## [93] later_1.3.0                   tibble_3.1.6                 
    ## [95] GenomicAlignments_1.31.2      AnnotationDbi_1.57.1         
    ## [97] memoise_2.0.1                 interactiveDisplayBase_1.33.0
    ## [99] ellipsis_0.3.2

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-ops" class="csl-entry">

Feldman, David, Avtar Singh, Jonathan L Schmid-Burgk, Rebecca J Carlson,
Anja Mezger, Anthony J Garrity, Feng Zhang, and Paul C Blainey. 2019.
“Optical Pooled Screens in Human Cells.” *Cell* 179 (3): 787–99.

</div>

</div>
