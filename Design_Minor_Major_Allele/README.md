Using crisprDesign to design gRNAs with minor and major alleles
================
Jean-Philippe Fortin, Luke Hoberecht

-   <a href="#introduction" id="toc-introduction">Introduction</a>
-   <a href="#installation" id="toc-installation">Installation</a>
-   <a href="#terminology" id="toc-terminology">Terminology</a>
-   <a href="#grna-design" id="toc-grna-design">gRNA design</a>
    -   <a href="#loading-packages" id="toc-loading-packages">Loading
        packages</a>
    -   <a href="#designing-grnas-for-human-hg38-with-injected-major-alleles"
        id="toc-designing-grnas-for-human-hg38-with-injected-major-alleles">Designing
        gRNAs for human (hg38) with injected major alleles</a>
    -   <a href="#designing-grnas-for-human-hg38-with-injected-minor-alleles"
        id="toc-designing-grnas-for-human-hg38-with-injected-minor-alleles">Designing
        gRNAs for human (hg38) with injected minor alleles</a>
-   <a href="#session-info" id="toc-session-info">Session Info</a>

# Introduction

Genetic variants such as single nucleotide polymorphisms (SNPs) can be
problematic in guide RNA (gRNA) design, as different alleles can result
in unintended gRNA:DNA mismatches for on-targets that reduce gRNA
efficacy. To circumvent this, it is advisable to generally avoid
targeting sequences that contain variants. However, this may not always
be possible, due to a small target window and/or few target options, or
desirable, if, for example, a CRISPR application intends to target a
pathogenic variant.

Functions in `crisprDesign` are well equipped to handle these cases.
gRNAs overlapping SNPs can be identified with the `addSNPAnnotation`
function, as documented in the [CRISPRko design with
Cas9](https://github.com/crisprVerse/Tutorials/tree/master/Design_CRISPRko_Cas9)
tutorial. Should the user wish to target region despite (or because of)
the presence of variants, the user only needs to take care in the choice
of `BSgenome` when constructing the `GuideSet` (an alternative option is
to supply a custom target sequence; see the [Working with a custom
sequence](https://github.com/crisprVerse/Tutorials/tree/master/Design_Custom_Sequence)
tutorial for more information).

This tutorial covers use cases for `BSgenome` objects that store
variants of the reference human genome (hg38) injected with major and
minor alleles. It assumes the reader is familiar with constructing and
using gene annotation objects (see the [Building a gene annotation
object](https://github.com/crisprVerse/Tutorials/tree/master/Building_Gene_Annotation)
tutorial) and `GuideSet` objects (see the [CRISPRko design with
Cas9](https://github.com/crisprVerse/Tutorials/tree/master/Design_CRISPRko_Cas9)
tutorial) so that the content may focus on the utility of the `BSgenome`
variants discussed herein. Please consult the applicable tutorials if
necessary.

Finally, it goes without saying that the user should be knowledgeable of
the sequence(s), including possible variations in such, he or she is
designing gRNAs for.

# Installation

See the [Installation
tutorial](https://github.com/crisprVerse/Tutorials/tree/master/Installation)
to learn how to install the crisprVerse packages necessary for this
tutorial.

# Terminology

See the [CRISPRko design
tutorial](https://github.com/crisprVerse/Tutorials/tree/master/Design_CRISPRko_Cas9)
to get familiar with the terminology used throughout this tutorial.

# gRNA design

## Loading packages

We first load the necessary packages for this tutorial:

``` r
library(crisprBase)
library(crisprDesign)
library(crisprDesignData)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg38.dbSNP151.major)
library(BSgenome.Hsapiens.UCSC.hg38.dbSNP151.minor)
```

We will also load the gene annotation model `txdb_human` from the
`crisprDesignData` package:

``` r
data(txdb_human, package="crisprDesignData")
```

The `BSgenome.Hsapiens.UCSC.hg38.dbSNP151.major` and
`BSgenome.Hsapiens.UCSC.hg38.dbSNP151.minor` packages are `BSgenome`
packages that contain the major and minor alleles of the human reference
genome hg38 based dbSNP151.For more information, type the following:

``` r
help(BSgenome.Hsapiens.UCSC.hg38.dbSNP151.major)
help(BSgenome.Hsapiens.UCSC.hg38.dbSNP151.minor)
```

## Designing gRNAs for human (hg38) with injected major alleles

It is worth noting that the human reference genome sequence
(GRCh38.p12), does not give the major allele (i.e. most common allele in
a population) at all nucleotide locations. Indeed, given that it was
historically constructed from a small set of human genomes, it contains
minor alleles that were common across this set of human genomes.

For example, in the coding region (CDS) of the human gene SMC3, the
reference genome contains the minor allele of the SNP rs2419565; the
reference allele (A) frequency is 0.00577, and the alternative allele
(G) frequency is 0.99423 as indicated in the table here
[here](https://www.ncbi.nlm.nih.gov/snp/rs2419565)).

Designing gRNAs targeting the CDS of SMC3 with the SpCas9 nuclease
returns one gRNA that overlaps this SNP. Below, we first construct a
`GuideSet` object using the reference `BSgenome` object:

``` r
smc3 <- queryTxObject(txdb_human,
                      featureType="cds",
                      queryColumn="gene_symbol",
                      queryValue="SMC3")
gs_reference <- findSpacers(smc3,
                            crisprNuclease=SpCas9,
                            bsgenome=BSgenome.Hsapiens.UCSC.hg38)
```

and a `GuideSet` object with the `BSgenome` object that contains the
major alleles:

``` r
gs_major <- findSpacers(smc3,
                        crisprNuclease=SpCas9,
                        bsgenome=BSgenome.Hsapiens.UCSC.hg38.dbSNP151.major)
```

Let’s compare the protospacer sequences from both objects:

``` r

protospacers(gs_reference["spacer_199"])
## DNAStringSet object of length 1:
##     width seq                                               names               
## [1]    20 TAGGGGTTACAAATCAATCA                              spacer_199
protospacers(gs_major["spacer_199"])
## DNAStringSet object of length 1:
##     width seq                                               names               
## [1]    20 TAGGGGTTACAAATCGATCA                              spacer_199
```

The variant occurs in the seed sequence of this gRNA, 5 bases upstream
of the `pam_site`, so a gRNA:DNA mismatch at this location is likely
detrimental to its efficacy. Also, as this major allele occurs at \>99%
frequency, it may be more beneficial to design gRNAs in this example
using the major allele genome contained in
`BSgenome.Hsapiens.UCSC.hg38.dbSNP151.major`.

## Designing gRNAs for human (hg38) with injected minor alleles

It may be desirable, in some applications, to target a genic sequence
that contains a minor allele (i.e. less common allele) rather than the
major or reference allele. For example, if a particular minor allele is
pathogenic and the host cell has a single copy of that allele, the user
may want to target that pathogenic variant and disrupt its behavior
while leaving the other copy undisturbed.

As an example, using `BSgenome.Hsapiens.UCSC.hg38.dbSNP151.minor`, we
can target the pathogenic minor allele
([rs398122995](https://www.ncbi.nlm.nih.gov/clinvar/variation/92240/?oq=rs398122995&m=NM_001378454.1(ALMS1):c.1897C%3ET%20(p.Gln633Ter)))
located in the human gene ALMS1:

``` r
alms1 <- queryTxObject(txdb_human, 'cds', 'gene_symbol', 'ALMS1')
gs_minor <- findSpacers(alms1,
                        crisprNuclease=SpCas9,
                        bsgenome=BSgenome.Hsapiens.UCSC.hg38.dbSNP151.minor)
gs_minor <- unique(gs_minor)
```

We also include, for comparison, the resulting `GuideSet` using the
reference genome sequence:

``` r
gs_reference <- findSpacers(alms1,
                            crisprNuclease=SpCas9,
                            bsgenome=BSgenome.Hsapiens.UCSC.hg38)
gs_reference <- unique(gs_reference)
```

and compare the two versions of the gRNA:

``` r
gs_reference["spacer_615"]
```

    ## GuideSet object with 1 range and 5 metadata columns:
    ##              seqnames    ranges strand |          protospacer            pam
    ##                 <Rle> <IRanges>  <Rle> |       <DNAStringSet> <DNAStringSet>
    ##   spacer_615     chr2  73448425      - | TACTCTCTGGTAACTCTTGT            TGG
    ##               pam_site  cut_site      region
    ##              <numeric> <numeric> <character>
    ##   spacer_615  73448425  73448428    region_7
    ##   -------
    ##   seqinfo: 640 sequences (1 circular) from hg38 genome
    ##   crisprNuclease: SpCas9

``` r
gs_minor["spacer_612"]
```

    ## GuideSet object with 1 range and 5 metadata columns:
    ##              seqnames    ranges strand |          protospacer            pam
    ##                 <Rle> <IRanges>  <Rle> |       <DNAStringSet> <DNAStringSet>
    ##   spacer_612     chr2  73448425      - | TACTCTCTGGTAACTCTTGC            TGG
    ##               pam_site  cut_site      region
    ##              <numeric> <numeric> <character>
    ##   spacer_612  73448425  73448428    region_7
    ##   -------
    ##   seqinfo: 595 sequences (1 circular) from hg38 genome
    ##   crisprNuclease: SpCas9

The variant occurs 1 base upstream of the `pam_site`, and likely
influences gRNA activity, that is, we can design a gRNA that targets the
minor allele and has a much lower affinity for the reference, or major
allele.

Note that while the two `GuideSet` objects differ only by their
`BSgenome` object, we need to provide different indices to access
protospacers at equivalent `pam_site`s. This is due to variants in one
`BSgenome` (in this case theone with minor alleles) eliminating PAM
sequences, that is, one of the Gs in NGG is changed to another base such
that SpCas9 does not recognize it. This, where permissible, is also an
effective way of ensuring gRNAs only target a specific sequence if that
sequence contains the desired variant.

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
    ##  [1] BSgenome.Hsapiens.UCSC.hg38.dbSNP151.minor_0.0.9999
    ##  [2] BSgenome.Hsapiens.UCSC.hg38.dbSNP151.major_0.0.9999
    ##  [3] BSgenome.Hsapiens.UCSC.hg38_1.4.4                  
    ##  [4] BSgenome_1.65.2                                    
    ##  [5] rtracklayer_1.57.0                                 
    ##  [6] Biostrings_2.65.2                                  
    ##  [7] XVector_0.37.0                                     
    ##  [8] GenomicRanges_1.49.1                               
    ##  [9] GenomeInfoDb_1.33.5                                
    ## [10] IRanges_2.31.2                                     
    ## [11] S4Vectors_0.35.1                                   
    ## [12] BiocGenerics_0.43.1                                
    ## [13] crisprDesignData_0.99.17                           
    ## [14] crisprDesign_0.99.133                              
    ## [15] crisprBase_1.1.5                                   
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] bitops_1.0-7                  matrixStats_0.62.0           
    ##   [3] bit64_4.0.5                   filelock_1.0.2               
    ##   [5] progress_1.2.2                httr_1.4.4                   
    ##   [7] tools_4.2.1                   utf8_1.2.2                   
    ##   [9] R6_2.5.1                      DBI_1.1.3                    
    ##  [11] tidyselect_1.1.2              prettyunits_1.1.1            
    ##  [13] bit_4.0.4                     curl_4.3.2                   
    ##  [15] compiler_4.2.1                crisprBowtie_1.1.1           
    ##  [17] cli_3.3.0                     Biobase_2.57.1               
    ##  [19] basilisk.utils_1.9.1          crisprScoreData_1.1.3        
    ##  [21] xml2_1.3.3                    DelayedArray_0.23.1          
    ##  [23] randomForest_4.7-1.1          readr_2.1.2                  
    ##  [25] rappdirs_0.3.3                stringr_1.4.1                
    ##  [27] digest_0.6.29                 Rsamtools_2.13.4             
    ##  [29] rmarkdown_2.15.2              crisprScore_1.1.14           
    ##  [31] basilisk_1.9.2                pkgconfig_2.0.3              
    ##  [33] htmltools_0.5.3               MatrixGenerics_1.9.1         
    ##  [35] dbplyr_2.2.1                  fastmap_1.1.0                
    ##  [37] rlang_1.0.4                   rstudioapi_0.14              
    ##  [39] RSQLite_2.2.16                shiny_1.7.2                  
    ##  [41] BiocIO_1.7.1                  generics_0.1.3               
    ##  [43] jsonlite_1.8.0                BiocParallel_1.31.12         
    ##  [45] dplyr_1.0.9                   VariantAnnotation_1.43.3     
    ##  [47] RCurl_1.98-1.8                magrittr_2.0.3               
    ##  [49] GenomeInfoDbData_1.2.8        Matrix_1.4-1                 
    ##  [51] Rcpp_1.0.9                    fansi_1.0.3                  
    ##  [53] reticulate_1.25               Rbowtie_1.37.0               
    ##  [55] lifecycle_1.0.1               stringi_1.7.8                
    ##  [57] yaml_2.3.5                    SummarizedExperiment_1.27.1  
    ##  [59] zlibbioc_1.43.0               BiocFileCache_2.5.0          
    ##  [61] AnnotationHub_3.5.0           grid_4.2.1                   
    ##  [63] blob_1.2.3                    promises_1.2.0.1             
    ##  [65] parallel_4.2.1                ExperimentHub_2.5.0          
    ##  [67] crayon_1.5.1                  dir.expiry_1.5.0             
    ##  [69] lattice_0.20-45               GenomicFeatures_1.49.6       
    ##  [71] hms_1.1.2                     KEGGREST_1.37.3              
    ##  [73] knitr_1.40                    pillar_1.8.1                 
    ##  [75] rjson_0.2.21                  codetools_0.2-18             
    ##  [77] biomaRt_2.53.2                BiocVersion_3.16.0           
    ##  [79] XML_3.99-0.10                 glue_1.6.2                   
    ##  [81] evaluate_0.16                 BiocManager_1.30.18          
    ##  [83] httpuv_1.6.5                  png_0.1-7                    
    ##  [85] vctrs_0.4.1                   tzdb_0.3.0                   
    ##  [87] purrr_0.3.4                   assertthat_0.2.1             
    ##  [89] cachem_1.0.6                  xfun_0.32                    
    ##  [91] mime_0.12                     xtable_1.8-4                 
    ##  [93] restfulr_0.0.15               later_1.3.0                  
    ##  [95] tibble_3.1.8                  GenomicAlignments_1.33.1     
    ##  [97] AnnotationDbi_1.59.1          memoise_2.0.1                
    ##  [99] interactiveDisplayBase_1.35.0 ellipsis_0.3.2
