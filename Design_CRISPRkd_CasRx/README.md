Using crisprDesign to design gRNAs for CRISPRkd with CasRx
================

-   <a href="#introduction" id="toc-introduction">Introduction</a>
-   <a href="#installation" id="toc-installation">Installation</a>
-   <a href="#terminology" id="toc-terminology">Terminology</a>
-   <a href="#crispr-knockdown-with-cas13d"
    id="toc-crispr-knockdown-with-cas13d">CRISPR knockdown with Cas13d</a>
    -   <a href="#loading-packages" id="toc-loading-packages">Loading
        packages</a>
    -   <a href="#creating-the-guideset" id="toc-creating-the-guideset">Creating
        the GuideSet</a>
    -   <a href="#annotating-the-guideset"
        id="toc-annotating-the-guideset">Annotating the GuideSet</a>
        -   <a href="#adding-spacer-alignments"
            id="toc-adding-spacer-alignments">Adding spacer alignments</a>
-   <a href="#session-info" id="toc-session-info">Session Info</a>
-   <a href="#references" id="toc-references">References</a>

Authors: Jean-Philippe Fortin, Luke Hoberecht

Date: 05 August, 2022

# Introduction

`crisprDesign` is a comprehensive software package for designing and
annotating CRISPR guide RNA (gRNA) sequences, including the
characterization of on-targets and off-targets, gene context annotation,
and SNP annotation (human only). The software was developed to be as
applicable and generalizable as possible. It currently support four
types of CRISPR modalities (modes of perturbations): CRISPR knockout
(CRISPRko), CRISPR activation (CRISPRa), CRISPR inhibition (CRISPRi) and
CRISPR base editing (CRISPRbe) (see Kampmann (2018) for a review of
CRISPR modalities).

This package utilizes the `crisprBase` package to enable gRNA design for
any CRISPR nuclease via the `CrisprNuclease` class. Nucleases that are
commonly used in the field are provided, including DNA-targeting
nucleases (e.g. SpCas9, AsCas12a) and RNA-targeting nuclease (e.g. CasRx
(RfxCas13d)).

`crisprDesign` is fully developed to work with the genome of any
organism, and can also be used to design gRNAs targeting custom DNA
sequences.

# Installation

`crisprDesign` can be installed from Bioconductor using the following
commands in an R session:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("crisprDesign")
```

Users interested in contributing to `crisprDesign` might want to look at
the following CRISPR-related package dependencies:

-   `crisprBase`: core CRISPR functions and S4 objects
-   `crisprBowtie`: aligns gRNA spacers to genomes using the ungapped
    aligner `bowtie`
-   `crisprBwa`: aligns gRNA spacers to genomes using the ungapped
    aligner `BWA`
-   `crisprScore`: implements state-of-the-art on- and off-target
    scoring algorithms
-   `crisprScoreData`: pre-trained models necessary for `crisprScore`

You can contribute to the package by submitting pull requests to our
[GitHub repo](https://github.com/Jfortin1/crisprDesign).

# Terminology

See the [CRISPRko design
vignette](https://github.com/crisprVerse/Tutorials/tree/master/Design_CRISPRko_Cas9)
to get familiar with the terminology used throughout this tutorial.

# CRISPR knockdown with Cas13d

## Loading packages

The following examples use the `crisprBase`, `crisprDesign`,
`crisprDesignData`, and `BSgenome.Hsapiens.UCSC.hg38` packages. Before
we begin, let’s install and load the necessary packages (`crisprDesign`
(and tacitly its dependency `crisprBase`) are installed in the
**Installation** section).

``` r
library(crisprBase)
library(crisprDesign)

install.packages("devtools")
devtools::install_github("Jfortin1/crisprDesignData")
library(crisprDesignData)

BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
```

## Creating the GuideSet

In this example, we will design gRNAs for an RNA-targeting nuclease
using `crisprDesign`, specifically, we will use the CasRx (Cas13d)
nuclease (Konermann et al. 2018) to target the primary isoform of the
human KRAS gene. In contrast to DNA-targeting nucleases, the target
spacers are composed of mRNA sequences rather than DNA genomic
sequences. \[disclosure…\]

We begin by loading the CasRx `CrisprNuclease` object from `crisprBase`:

``` r
data(CasRx, package="crisprBase")
CasRx
```

    ## Class: CrisprNuclease
    ##   Name: CasRx
    ##   Target type: RNA
    ##   Metadata: list of length 2
    ##   PFS: N
    ##   Weights: 1
    ##   Spacer length: 23
    ##   PFS side: 3prime
    ##     Distance from PFS: 0
    ##   Prototype protospacers: 5'--SSSSSSSSSSSSSSSSSSSSSSS[N]--3'

The PFS sequence (the equivalent of a PAM sequence for RNA-targeting
nucleases) for CasRx is `N`, meaning there is no specific PFS sequences
preferred by CasRx.

Next, we extract the mRNA sequence for our target transcript (Ensembl
transcript ID: ENST00000311936). For this we need the coordinates for
our transcript, which are conveniently stored in the `crisprDesignData`
package as `txdb_human` (for more information on `txdb_human` and how to
create similar gene annotation objects, see the [Building a gene
annotation
object](https://github.com/crisprVerse/Tutorials/tree/master/Building_Gene_Annotation)
tutorial), as well as a `BSgenome` object containing the exonic
sequences.

``` r
data("txdb_human", package="crisprDesignData")
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
```

We obtain our mRNA sequence with `getMrnaSequences`:

``` r
txid <- "ENST00000311936"
mrna <- getMrnaSequences(txids=txid,
                         bsgenome=bsgenome,
                         txObject=txdb_human)
mrna
```

    ## DNAStringSet object of length 1:
    ##     width seq                                               names               
    ## [1]  5306 CTAGGCGGCGGCCGCGGCGGCGG...GTAATAAAAATAGTTACAGTGAC ENST00000311936

We then use the function `findSpacers` to design our gRNAs. For the sake
of brevity we will only consider a subset of 100 gRNAs:

``` r
gs <- findSpacers(mrna[["ENST00000311936"]],
                  crisprNuclease=CasRx)
gs <- gs[1000:1100]
head(gs)
```

    ## GuideSet object with 6 ranges and 5 metadata columns:
    ##               seqnames    ranges strand |             protospacer
    ##                  <Rle> <IRanges>  <Rle> |          <DNAStringSet>
    ##   spacer_1000 region_1      1023      + | CCTAAGATTTCTGTCTTGGGGCT
    ##   spacer_1001 region_1      1024      + | CTAAGATTTCTGTCTTGGGGCTT
    ##   spacer_1002 region_1      1025      + | TAAGATTTCTGTCTTGGGGCTTT
    ##   spacer_1003 region_1      1026      + | AAGATTTCTGTCTTGGGGCTTTT
    ##   spacer_1004 region_1      1027      + | AGATTTCTGTCTTGGGGCTTTTG
    ##   spacer_1005 region_1      1028      + | GATTTCTGTCTTGGGGCTTTTGG
    ##                          pam  pam_site  cut_site      region
    ##               <DNAStringSet> <numeric> <numeric> <character>
    ##   spacer_1000              T      1023        NA    region_1
    ##   spacer_1001              T      1024        NA    region_1
    ##   spacer_1002              T      1025        NA    region_1
    ##   spacer_1003              G      1026        NA    region_1
    ##   spacer_1004              G      1027        NA    region_1
    ##   spacer_1005              T      1028        NA    region_1
    ##   -------
    ##   seqinfo: 1 sequence from custom genome
    ##   crisprNuclease: CasRx

Note that all protospacer sequences are located on the original strand
of the mRNA sequence. For RNA-targeting nucleases, the spacer and
protospacer sequences are the reverse complement of each other. (Compare
the output of the code below with a `GuideSet` that uses a DNA-targeting
nuclease–for such `GuideSet`s, the output of `spacers` and
`protospacers` are identical.)

``` r
head(spacers(gs))
```

    ## DNAStringSet object of length 6:
    ##     width seq                                               names               
    ## [1]    23 AGCCCCAAGACAGAAATCTTAGG                           spacer_1000
    ## [2]    23 AAGCCCCAAGACAGAAATCTTAG                           spacer_1001
    ## [3]    23 AAAGCCCCAAGACAGAAATCTTA                           spacer_1002
    ## [4]    23 AAAAGCCCCAAGACAGAAATCTT                           spacer_1003
    ## [5]    23 CAAAAGCCCCAAGACAGAAATCT                           spacer_1004
    ## [6]    23 CCAAAAGCCCCAAGACAGAAATC                           spacer_1005

``` r
head(protospacers(gs))
```

    ## DNAStringSet object of length 6:
    ##     width seq                                               names               
    ## [1]    23 CCTAAGATTTCTGTCTTGGGGCT                           spacer_1000
    ## [2]    23 CTAAGATTTCTGTCTTGGGGCTT                           spacer_1001
    ## [3]    23 TAAGATTTCTGTCTTGGGGCTTT                           spacer_1002
    ## [4]    23 AAGATTTCTGTCTTGGGGCTTTT                           spacer_1003
    ## [5]    23 AGATTTCTGTCTTGGGGCTTTTG                           spacer_1004
    ## [6]    23 GATTTCTGTCTTGGGGCTTTTGG                           spacer_1005

## Annotating the GuideSet

Next, we annotate our candidate gRNAs to assess quality. There are
several functions in `crisprDesign` that provide annotation for features
that are nonspecific to CRISPRkd, for which we refer the reader to the
[CRISPRko design with
Cas9](https://github.com/crisprVerse/Tutorials/tree/master/Design_CRISPRko_Cas9)
tutorial for more information. The sections below will cover annotation
functions that are of particular interest to, or deserve extra care for
CRISPRkd applications.

### Adding spacer alignments

Since our CRISPR nuclease targets RNA rather than DNA, off-target
searches should be restricted to the transcriptome. We can perform such
a search using one of two methods. For the first method, we set the
`aligner` argument to `"biostrings"` and pass a `DNAStringSet`
representation of the transcriptome to the argument `custom_seq`. We can
create this representation with `getMrnaSequences` and all transcript
IDs found in `txdb_human`. For the sake brevity and demonstration,
however, we will limit our search to mRNAs belonging to the KRAS gene.
We will also search for off-targets having up to one mismatch and pass
`txdb_human` to the `txObject` argument so that our alignments will be
accompanied with gene annotation.

``` r
tx_ids <- c("ENST00000256078", "ENST00000311936",
            "ENST00000557334", "ENST00000556131")
kras_mrnas <- getMrnaSequences(txids=tx_ids,
                               bsgenome=bsgenome,
                               txObject=txdb_human)
gs <- addSpacerAlignments(gs,
                          aligner="biostrings",
                          txObject=txdb_human,
                          n_mismatches=1,
                          custom_seq=kras_mrnas)
```

``` r
tail(gs)
```

    ## GuideSet object with 6 ranges and 10 metadata columns:
    ##               seqnames    ranges strand |             protospacer
    ##                  <Rle> <IRanges>  <Rle> |          <DNAStringSet>
    ##   spacer_1095 region_1      1118      + | AGCTTTTGAATCATCCCTATTCT
    ##   spacer_1096 region_1      1119      + | GCTTTTGAATCATCCCTATTCTG
    ##   spacer_1097 region_1      1120      + | CTTTTGAATCATCCCTATTCTGT
    ##   spacer_1098 region_1      1121      + | TTTTGAATCATCCCTATTCTGTG
    ##   spacer_1099 region_1      1122      + | TTTGAATCATCCCTATTCTGTGT
    ##   spacer_1100 region_1      1123      + | TTGAATCATCCCTATTCTGTGTT
    ##                          pam  pam_site  cut_site      region     n0_tx
    ##               <DNAStringSet> <numeric> <numeric> <character> <numeric>
    ##   spacer_1095              G      1118        NA    region_1         3
    ##   spacer_1096              T      1119        NA    region_1         3
    ##   spacer_1097              G      1120        NA    region_1         3
    ##   spacer_1098              T      1121        NA    region_1         3
    ##   spacer_1099              T      1122        NA    region_1         3
    ##   spacer_1100              T      1123        NA    region_1         3
    ##                   n1_tx   n0_gene   n1_gene
    ##               <numeric> <numeric> <numeric>
    ##   spacer_1095         0         1         0
    ##   spacer_1096         0         1         0
    ##   spacer_1097         0         1         0
    ##   spacer_1098         0         1         0
    ##   spacer_1099         0         1         0
    ##   spacer_1100         0         1         0
    ##                                                                        alignments
    ##                                                                     <GRangesList>
    ##   spacer_1095 ENST00000557334:722:+,ENST00000311936:1088:+,ENST00000256078:1246:+
    ##   spacer_1096 ENST00000256078:1179:+,ENST00000557334:756:+,ENST00000311936:1122:+
    ##   spacer_1097 ENST00000311936:1055:+,ENST00000256078:1213:+,ENST00000557334:790:+
    ##   spacer_1098 ENST00000557334:723:+,ENST00000311936:1089:+,ENST00000256078:1247:+
    ##   spacer_1099 ENST00000256078:1180:+,ENST00000557334:757:+,ENST00000311936:1123:+
    ##   spacer_1100 ENST00000311936:1056:+,ENST00000256078:1214:+,ENST00000557334:791:+
    ##   -------
    ##   seqinfo: 1 sequence from custom genome
    ##   crisprNuclease: CasRx

The columns `n0_gene` and `n0_tx` report the number of on-targets at the
gene- and transcript-level, respectively. For instance, for all spacers
shown above `n0_tx` is equal to 3, meaning they all map to three
isoforms of KRAS (since our `custom_seq` is made up of KRAS isoforms
only). We can use the `onTargets` accessor function to look at each
alignment for the first spacer:

``` r
# onTargets(gs["spacer_1095"]) # bug: incorrect spacer alignments selected
onTargets(gs)[names(onTargets(gs)) == "spacer_1095"]
```

    ## GRanges object with 3 ranges and 9 metadata columns:
    ##                      seqnames    ranges strand |                 spacer
    ##                         <Rle> <IRanges>  <Rle> |            <character>
    ##   spacer_1095 ENST00000256078      1242      + | AGAATAGGGATGATTCAAAA..
    ##   spacer_1095 ENST00000311936      1118      + | AGAATAGGGATGATTCAAAA..
    ##   spacer_1095 ENST00000557334       786      + | AGAATAGGGATGATTCAAAA..
    ##                           protospacer            pam  pam_site n_mismatches
    ##                        <DNAStringSet> <DNAStringSet> <numeric>    <numeric>
    ##   spacer_1095 AGCTTTTGAATCATCCCTATTCT              G      1242            0
    ##   spacer_1095 AGCTTTTGAATCATCCCTATTCT              G      1118            0
    ##   spacer_1095 AGCTTTTGAATCATCCCTATTCT              G       786            0
    ##               canonical  cut_site         gene_id gene_symbol
    ##               <logical> <numeric>     <character> <character>
    ##   spacer_1095      TRUE        NA ENSG00000133703        KRAS
    ##   spacer_1095      TRUE        NA ENSG00000133703        KRAS
    ##   spacer_1095      TRUE        NA ENSG00000133703        KRAS
    ##   -------
    ##   seqinfo: 4 sequences from custom genome

The second method uses the `bowtie` (or `bwa`) aligner. This requires
building a transcriptome bowtie (or BWA) index file first. See the
[Building genome indices for short read
aligners](https://github.com/crisprVerse/Tutorials/tree/master/Building_Genome_Indices)
tutorial for more information.

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
    ##  [3] rtracklayer_1.56.1                Biostrings_2.64.0                
    ##  [5] XVector_0.36.0                    GenomicRanges_1.48.0             
    ##  [7] GenomeInfoDb_1.32.2               IRanges_2.30.0                   
    ##  [9] S4Vectors_0.34.0                  BiocGenerics_0.42.0              
    ## [11] crisprDesignData_0.99.12          crisprDesign_0.99.110            
    ## [13] crisprBase_1.1.2                 
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] rjson_0.2.21                  ellipsis_0.3.2               
    ##   [3] Rbowtie_1.36.0                fs_1.5.2                     
    ##   [5] rstudioapi_0.13               remotes_2.4.2                
    ##   [7] bit64_4.0.5                   interactiveDisplayBase_1.34.0
    ##   [9] AnnotationDbi_1.58.0          fansi_1.0.3                  
    ##  [11] xml2_1.3.3                    codetools_0.2-18             
    ##  [13] cachem_1.0.6                  knitr_1.39                   
    ##  [15] pkgload_1.3.0                 jsonlite_1.8.0               
    ##  [17] Rsamtools_2.12.0              dbplyr_2.2.1                 
    ##  [19] png_0.1-7                     shiny_1.7.2                  
    ##  [21] BiocManager_1.30.18           readr_2.1.2                  
    ##  [23] compiler_4.2.0                httr_1.4.3                   
    ##  [25] basilisk_1.9.2                assertthat_0.2.1             
    ##  [27] Matrix_1.4-1                  fastmap_1.1.0                
    ##  [29] cli_3.3.0                     later_1.3.0                  
    ##  [31] htmltools_0.5.3               prettyunits_1.1.1            
    ##  [33] tools_4.2.0                   glue_1.6.2                   
    ##  [35] GenomeInfoDbData_1.2.8        crisprBowtie_1.0.0           
    ##  [37] dplyr_1.0.9                   rappdirs_0.3.3               
    ##  [39] Rcpp_1.0.9                    Biobase_2.56.0               
    ##  [41] vctrs_0.4.1                   ExperimentHub_2.4.0          
    ##  [43] crisprBwa_1.0.0               crisprScore_1.1.14           
    ##  [45] xfun_0.31                     stringr_1.4.0                
    ##  [47] ps_1.7.1                      mime_0.12                    
    ##  [49] miniUI_0.1.1.1                lifecycle_1.0.1              
    ##  [51] restfulr_0.0.15               devtools_2.4.4               
    ##  [53] XML_3.99-0.10                 AnnotationHub_3.4.0          
    ##  [55] zlibbioc_1.42.0               basilisk.utils_1.9.1         
    ##  [57] VariantAnnotation_1.42.1      hms_1.1.1                    
    ##  [59] promises_1.2.0.1              MatrixGenerics_1.8.1         
    ##  [61] parallel_4.2.0                SummarizedExperiment_1.26.1  
    ##  [63] yaml_2.3.5                    curl_4.3.2                   
    ##  [65] memoise_2.0.1                 reticulate_1.25              
    ##  [67] biomaRt_2.52.0                stringi_1.7.8                
    ##  [69] RSQLite_2.2.15                BiocVersion_3.15.2           
    ##  [71] BiocIO_1.6.0                  randomForest_4.7-1.1         
    ##  [73] crisprScoreData_1.1.3         GenomicFeatures_1.48.3       
    ##  [75] filelock_1.0.2                pkgbuild_1.3.1               
    ##  [77] BiocParallel_1.30.3           rlang_1.0.4                  
    ##  [79] pkgconfig_2.0.3               matrixStats_0.62.0           
    ##  [81] bitops_1.0-7                  evaluate_0.15                
    ##  [83] lattice_0.20-45               purrr_0.3.4                  
    ##  [85] htmlwidgets_1.5.4             GenomicAlignments_1.32.1     
    ##  [87] bit_4.0.4                     tidyselect_1.1.2             
    ##  [89] processx_3.7.0                magrittr_2.0.3               
    ##  [91] R6_2.5.1                      profvis_0.3.7                
    ##  [93] generics_0.1.3                DelayedArray_0.22.0          
    ##  [95] DBI_1.1.3                     pillar_1.8.0                 
    ##  [97] KEGGREST_1.36.3               RCurl_1.98-1.8               
    ##  [99] tibble_3.1.8                  dir.expiry_1.4.0             
    ## [101] crayon_1.5.1                  utf8_1.2.2                   
    ## [103] BiocFileCache_2.4.0           urlchecker_1.0.1             
    ## [105] tzdb_0.3.0                    rmarkdown_2.14               
    ## [107] usethis_2.1.6                 progress_1.2.2               
    ## [109] grid_4.2.0                    blob_1.2.3                   
    ## [111] callr_3.7.1                   digest_0.6.29                
    ## [113] xtable_1.8-4                  httpuv_1.6.5                 
    ## [115] Rbwa_1.0.0                    sessioninfo_1.2.2

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-crispracrisprireview" class="csl-entry">

Kampmann, Martin. 2018. “CRISPRi and CRISPRa Screens in Mammalian Cells
for Precision Biology and Medicine.” *ACS Chemical Biology* 13 (2):
406–16.

</div>

<div id="ref-cas13d" class="csl-entry">

Konermann, Silvana, Peter Lotfy, Nicholas J Brideau, Jennifer Oki, Maxim
N Shokhirev, and Patrick D Hsu. 2018. “Transcriptome Engineering with
RNA-Targeting Type VI-d CRISPR Effectors.” *Cell* 173 (3): 665–76.

</div>

</div>
