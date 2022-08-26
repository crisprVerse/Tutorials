Validating existing gRNA libraries
================
Jean-Philippe Fortin, Luke Hoberecht

-   <a href="#introduction" id="toc-introduction">Introduction</a>
-   <a href="#loading-necessary-packages"
    id="toc-loading-necessary-packages">Loading necessary packages</a>
-   <a href="#reading-in-data" id="toc-reading-in-data">Reading in data</a>
-   <a href="#building-a-guideset-object"
    id="toc-building-a-guideset-object">Building a <code>GuideSet</code>
    object</a>
-   <a href="#off-target-characterization"
    id="toc-off-target-characterization">Off-target characterization</a>
-   <a href="#session-info" id="toc-session-info">Session Info</a>

# Introduction

In this vignette, we characterize a small mouse CRISPR knockout
(CRISPRko) library that was designed to target tumor suppressors. The
library was obtained from Addgene, and is stored in the folder `extdata`
in the current directory.

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
# Not run
spacerLength(SpCas9) <- 19
```

We next need to define a bowtie index that we will use for alignment:

``` r
bowtie_index <- "/Users/fortinj2/crisprIndices/bowtie/mm10/mm10"
```

For instructions on how to build a Bowtie index from a given reference
genome, see the [genome index
tutorial](https://github.com/crisprVerse/Tutorials/tree/master/Building_Genome_Indices)
or the [crisprBowtie page](https://github.com/crisprVerse/crisprBowtie)
.

We first map the gRNAs to the reference genome with perfect match to
obtain genomic coordinates of those gRNAs:

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

Finally, we create unique ids to identify the spacer sequences:

``` r
ids <- paste0("gRNA_", seq_len(nrow(data)))
head(ids)
```

    ## [1] "gRNA_1" "gRNA_2" "gRNA_3" "gRNA_4" "gRNA_5" "gRNA_6"

We are now ready to build the `GuideSet` object using the constructor
function `GuideSet` from `crisprDesign`:

``` r
gs <- GuideSet(ids=ids,
               protospacers=data$spacer_20mer,
               pams=data$pam,
               pam_site=data$pam_site,
               seqnames=data$chr,
               strand=data$strand,
               CrisprNuclease=crisprNuclease,
               bsgenome=bsgenome)
gs$gene_symbol <- data$gene_symbol
```

The `GuideSet` object, and
[crisprDesign](https://github.com/crisprVerse/crisprDesign), provide
rich functionalities to annotate and manipulate gRNAs. See the [CRISPRko
design
tutorial](https://github.com/crisprVerse/Tutorials/tree/master/Design_CRISPRko_Cas9)
to get an overview of the functionalities. For the rest of this
tutorial, we only focus on characterizing the off-targets.

# Off-target characterization

Having a `GuideSet` object, it is now a piece of cake to characterize
the off-targets. We characterize off-targets using the bowtie aligner,
with up to 3 mismatches between the spacer (gRNA) and protospacer
(target DNA) sequences. The function `addSpacerAlignments` accomplishes
that.

It has an optional argument `txObject` that can be used to provide gene
model data to put the off-targets in a gene model context. We made such
objects available for human and mouse in the
[crisprDesignData](https://github.com/crisprVerse/crisprDesignData)
package (see `txdb_human` and `txdb_mouse`).

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
`?getSpacerAlignments` for more details about what the different columns
are.

As an example, we can access the on- and off-target alignments of the
first gRNA using the following commands:

``` r
aln <- gs$alignments[[1]]
aln
```

    ## GRanges object with 2 ranges and 14 metadata columns:
    ##          seqnames    ranges strand |               spacer          protospacer
    ##             <Rle> <IRanges>  <Rle> |       <DNAStringSet>       <DNAStringSet>
    ##   gRNA_1     chr8  45023166      - | GGGCAGTGTTTCAAAATCCA GGGCAGTGTTTCAAAATCCA
    ##   gRNA_1    chr12  80656576      - | GGGCAGTGTTTCAAAATCCA AGGCAGGGTTTCAAAATCCA
    ##                     pam  pam_site n_mismatches canonical  cut_site         cds
    ##          <DNAStringSet> <numeric>    <integer> <logical> <numeric> <character>
    ##   gRNA_1            AGG  45023166            0      TRUE  45023169        Fat1
    ##   gRNA_1            AGG  80656576            2      TRUE  80656579        <NA>
    ##             fiveUTRs   threeUTRs       exons     introns  intergenic
    ##          <character> <character> <character> <character> <character>
    ##   gRNA_1        <NA>        <NA>        Fat1        <NA>        <NA>
    ##   gRNA_1        <NA>        <NA>        <NA>     Slc39a9        <NA>
    ##          intergenic_distance
    ##                    <integer>
    ##   gRNA_1                <NA>
    ##   gRNA_1                <NA>
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

    ## R version 4.2.1 (2022-06-23)
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
    ##  [1] crisprDesignData_0.99.17           BSgenome.Mmusculus.UCSC.mm10_1.4.3
    ##  [3] BSgenome_1.65.2                    rtracklayer_1.57.0                
    ##  [5] Biostrings_2.65.2                  XVector_0.37.0                    
    ##  [7] GenomicRanges_1.49.1               GenomeInfoDb_1.33.5               
    ##  [9] IRanges_2.31.2                     S4Vectors_0.35.1                  
    ## [11] BiocGenerics_0.43.1                readxl_1.4.1                      
    ## [13] crisprBowtie_1.1.1                 crisprDesign_0.99.133             
    ## [15] crisprBase_1.1.5                  
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] rjson_0.2.21                  ellipsis_0.3.2               
    ##   [3] Rbowtie_1.37.0                rstudioapi_0.14              
    ##   [5] bit64_4.0.5                   interactiveDisplayBase_1.35.0
    ##   [7] AnnotationDbi_1.59.1          fansi_1.0.3                  
    ##   [9] xml2_1.3.3                    codetools_0.2-18             
    ##  [11] cachem_1.0.6                  knitr_1.40                   
    ##  [13] jsonlite_1.8.0                Rsamtools_2.13.4             
    ##  [15] dbplyr_2.2.1                  png_0.1-7                    
    ##  [17] shiny_1.7.2                   BiocManager_1.30.18          
    ##  [19] readr_2.1.2                   compiler_4.2.1               
    ##  [21] httr_1.4.4                    basilisk_1.9.2               
    ##  [23] assertthat_0.2.1              Matrix_1.4-1                 
    ##  [25] fastmap_1.1.0                 cli_3.3.0                    
    ##  [27] later_1.3.0                   htmltools_0.5.3              
    ##  [29] prettyunits_1.1.1             tools_4.2.1                  
    ##  [31] glue_1.6.2                    GenomeInfoDbData_1.2.8       
    ##  [33] dplyr_1.0.9                   rappdirs_0.3.3               
    ##  [35] Rcpp_1.0.9                    Biobase_2.57.1               
    ##  [37] cellranger_1.1.0              vctrs_0.4.1                  
    ##  [39] ExperimentHub_2.5.0           crisprBwa_1.1.3              
    ##  [41] crisprScore_1.1.14            xfun_0.32                    
    ##  [43] stringr_1.4.1                 mime_0.12                    
    ##  [45] lifecycle_1.0.1               restfulr_0.0.15              
    ##  [47] XML_3.99-0.10                 AnnotationHub_3.5.0          
    ##  [49] zlibbioc_1.43.0               basilisk.utils_1.9.1         
    ##  [51] vroom_1.5.7                   VariantAnnotation_1.43.3     
    ##  [53] hms_1.1.2                     promises_1.2.0.1             
    ##  [55] MatrixGenerics_1.9.1          parallel_4.2.1               
    ##  [57] SummarizedExperiment_1.27.1   yaml_2.3.5                   
    ##  [59] curl_4.3.2                    memoise_2.0.1                
    ##  [61] reticulate_1.25               biomaRt_2.53.2               
    ##  [63] stringi_1.7.8                 RSQLite_2.2.16               
    ##  [65] BiocVersion_3.16.0            BiocIO_1.7.1                 
    ##  [67] randomForest_4.7-1.1          crisprScoreData_1.1.3        
    ##  [69] GenomicFeatures_1.49.6        filelock_1.0.2               
    ##  [71] BiocParallel_1.31.12          rlang_1.0.4                  
    ##  [73] pkgconfig_2.0.3               matrixStats_0.62.0           
    ##  [75] bitops_1.0-7                  evaluate_0.16                
    ##  [77] lattice_0.20-45               purrr_0.3.4                  
    ##  [79] GenomicAlignments_1.33.1      bit_4.0.4                    
    ##  [81] tidyselect_1.1.2              magrittr_2.0.3               
    ##  [83] R6_2.5.1                      generics_0.1.3               
    ##  [85] DelayedArray_0.23.1           DBI_1.1.3                    
    ##  [87] pillar_1.8.1                  KEGGREST_1.37.3              
    ##  [89] RCurl_1.98-1.8                tibble_3.1.8                 
    ##  [91] dir.expiry_1.5.0              crayon_1.5.1                 
    ##  [93] utf8_1.2.2                    BiocFileCache_2.5.0          
    ##  [95] tzdb_0.3.0                    rmarkdown_2.15.2             
    ##  [97] progress_1.2.2                grid_4.2.1                   
    ##  [99] blob_1.2.3                    digest_0.6.29                
    ## [101] xtable_1.8-4                  httpuv_1.6.5                 
    ## [103] Rbwa_1.1.0
