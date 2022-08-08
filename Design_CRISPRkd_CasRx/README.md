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

Date: 08 August, 2022

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

The Bioconductor packages needed in this vignette can be downloaded
using the following commands:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("crisprDesign")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
```

The GitHub packages needed in this vignette can be downloaded using the
following commands:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("Jfortin1/crisprDesignData")
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
we begin, let’s load the necessary packages:

``` r
library(crisprBase)
library(crisprDesign)
library(crisprDesignData)
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

We then use the function `findSpacers` to design our gRNAs.

``` r
gs <- findSpacers(mrna[["ENST00000311936"]],
                  crisprNuclease=CasRx)
head(gs)
```

    ## GuideSet object with 6 ranges and 5 metadata columns:
    ##            seqnames    ranges strand |             protospacer            pam
    ##               <Rle> <IRanges>  <Rle> |          <DNAStringSet> <DNAStringSet>
    ##   spacer_1 region_1        24      + | CTAGGCGGCGGCCGCGGCGGCGG              A
    ##   spacer_2 region_1        25      + | TAGGCGGCGGCCGCGGCGGCGGA              G
    ##   spacer_3 region_1        26      + | AGGCGGCGGCCGCGGCGGCGGAG              G
    ##   spacer_4 region_1        27      + | GGCGGCGGCCGCGGCGGCGGAGG              C
    ##   spacer_5 region_1        28      + | GCGGCGGCCGCGGCGGCGGAGGC              A
    ##   spacer_6 region_1        29      + | CGGCGGCCGCGGCGGCGGAGGCA              G
    ##             pam_site  cut_site      region
    ##            <numeric> <numeric> <character>
    ##   spacer_1        24        NA    region_1
    ##   spacer_2        25        NA    region_1
    ##   spacer_3        26        NA    region_1
    ##   spacer_4        27        NA    region_1
    ##   spacer_5        28        NA    region_1
    ##   spacer_6        29        NA    region_1
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
    ## [1]    23 CCGCCGCCGCGGCCGCCGCCTAG                           spacer_1
    ## [2]    23 TCCGCCGCCGCGGCCGCCGCCTA                           spacer_2
    ## [3]    23 CTCCGCCGCCGCGGCCGCCGCCT                           spacer_3
    ## [4]    23 CCTCCGCCGCCGCGGCCGCCGCC                           spacer_4
    ## [5]    23 GCCTCCGCCGCCGCGGCCGCCGC                           spacer_5
    ## [6]    23 TGCCTCCGCCGCCGCGGCCGCCG                           spacer_6

``` r
head(protospacers(gs))
```

    ## DNAStringSet object of length 6:
    ##     width seq                                               names               
    ## [1]    23 CTAGGCGGCGGCCGCGGCGGCGG                           spacer_1
    ## [2]    23 TAGGCGGCGGCCGCGGCGGCGGA                           spacer_2
    ## [3]    23 AGGCGGCGGCCGCGGCGGCGGAG                           spacer_3
    ## [4]    23 GGCGGCGGCCGCGGCGGCGGAGG                           spacer_4
    ## [5]    23 GCGGCGGCCGCGGCGGCGGAGGC                           spacer_5
    ## [6]    23 CGGCGGCCGCGGCGGCGGAGGCA                           spacer_6

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
a search using one of two methods.

#### Adding spacer alignments with Biostrings

For the first method, we set the `aligner` argument to `"biostrings"`
and pass a `DNAStringSet` representation of the transcriptome to the
argument `custom_seq`. We can create this representation with
`getMrnaSequences` and all transcript IDs found in `txdb_human`. The
code below uses this method to search for off-targets having up to one
mismatch and passes `txdb_human` to the `txObject` argument so that the
alignments will be accompanied with gene annotation.

``` r
exon_ids <- unique(txdb_human$exons$tx_id)
mrnasHuman <- getMrnaSequences(exon_ids,
                               bsgenome=BSgenome.Hsapiens.UCSC.hg38,
                               txObject=txdb_human)
## long run time
results <- addSpacerAlignments(gs,
                               aligner="biostrings",
                               txObject=txdb_human,
                               n_mismatches=1,
                               custom_seq=mrnasHuman)
```

NOTE: since `mrnasHuman` contains many sequences (\>100k), this method
has a very long run time; for transcriptome-wide searches, or for
searches against a large number of sequences, we recommend the following
method instead.

#### Adding spacer alignments with bowtie or BWA

The second method uses the `bowtie` (or `bwa`) aligner. This requires
building a transcriptome bowtie (or BWA) index file first. See the
[Building genome indices for short read
aligners](https://github.com/crisprVerse/Tutorials/tree/master/Building_Genome_Indices)
tutorial for more information.

Here we set `aligner` to `"bowtie"` and pass a precomputed transcriptome
bowtie index to `aligner_index` to find off-targets:

``` r
bowtie_index <- "/Users/hoberecl/crisprIndices/bowtie/ensembl_human_104/ensembl_human_104"
results <- addSpacerAlignments(gs,
                               aligner="bowtie",
                               aligner_index=bowtie_index,
                               txObject=txdb_human,
                               n_mismatches=1)
```

``` r
head(results)
```

    ## GuideSet object with 6 ranges and 10 metadata columns:
    ##            seqnames    ranges strand |             protospacer            pam
    ##               <Rle> <IRanges>  <Rle> |          <DNAStringSet> <DNAStringSet>
    ##   spacer_1 region_1        24      + | CTAGGCGGCGGCCGCGGCGGCGG              A
    ##   spacer_2 region_1        25      + | TAGGCGGCGGCCGCGGCGGCGGA              G
    ##   spacer_3 region_1        26      + | AGGCGGCGGCCGCGGCGGCGGAG              G
    ##   spacer_4 region_1        27      + | GGCGGCGGCCGCGGCGGCGGAGG              C
    ##   spacer_5 region_1        28      + | GCGGCGGCCGCGGCGGCGGAGGC              A
    ##   spacer_6 region_1        29      + | CGGCGGCCGCGGCGGCGGAGGCA              G
    ##             pam_site  cut_site      region     n0_tx     n1_tx   n0_gene
    ##            <numeric> <numeric> <character> <numeric> <numeric> <numeric>
    ##   spacer_1        24        NA    region_1         4         0         1
    ##   spacer_2        25        NA    region_1         4         0         1
    ##   spacer_3        26        NA    region_1         4         9         1
    ##   spacer_4        27        NA    region_1         4        81         1
    ##   spacer_5        28        NA    region_1         4        67         1
    ##   spacer_6        29        NA    region_1         4        44         1
    ##              n1_gene
    ##            <numeric>
    ##   spacer_1         0
    ##   spacer_2         0
    ##   spacer_3         6
    ##   spacer_4        26
    ##   spacer_5        18
    ##   spacer_6         5
    ##                                                                    alignments
    ##                                                                 <GRangesList>
    ##   spacer_1 ENST00000256078:24:+,ENST00000311936:24:+,ENST00000556131:24:+,...
    ##   spacer_2 ENST00000256078:25:+,ENST00000311936:25:+,ENST00000556131:25:+,...
    ##   spacer_3 ENST00000256078:26:+,ENST00000311936:26:+,ENST00000556131:26:+,...
    ##   spacer_4 ENST00000256078:27:+,ENST00000311936:27:+,ENST00000556131:27:+,...
    ##   spacer_5 ENST00000256078:28:+,ENST00000311936:28:+,ENST00000556131:28:+,...
    ##   spacer_6 ENST00000256078:29:+,ENST00000311936:29:+,ENST00000556131:29:+,...
    ##   -------
    ##   seqinfo: 1 sequence from custom genome
    ##   crisprNuclease: CasRx

The columns `n0_gene` and `n0_tx` report the number of on-targets at the
gene- and transcript-level, respectively. For instance, each spacer
shown above shows `n0_gene` equal to 1 and `n0_tx` equal to 4, meaning
each spacer maps to all four isoforms of KRAS. We can retrieve
information about each alignment with the `onTargets` function. Looking
at the on-targets for the first spacer we can see where the target
`pam_site` is relative to the start of the transcript with respect to
each isoform of KRAS.

``` r
onTargets(results["spacer_1"])
```

    ## GRanges object with 4 ranges and 9 metadata columns:
    ##                   seqnames    ranges strand |                  spacer
    ##                      <Rle> <IRanges>  <Rle> |          <DNAStringSet>
    ##   spacer_1 ENST00000256078        24      + | CCGCCGCCGCGGCCGCCGCCTAG
    ##   spacer_1 ENST00000311936        24      + | CCGCCGCCGCGGCCGCCGCCTAG
    ##   spacer_1 ENST00000556131        24      + | CCGCCGCCGCGGCCGCCGCCTAG
    ##   spacer_1 ENST00000557334        31      + | CCGCCGCCGCGGCCGCCGCCTAG
    ##                        protospacer            pam  pam_site n_mismatches
    ##                     <DNAStringSet> <DNAStringSet> <numeric>    <integer>
    ##   spacer_1 CTAGGCGGCGGCCGCGGCGGCGG              A        24            0
    ##   spacer_1 CTAGGCGGCGGCCGCGGCGGCGG              A        24            0
    ##   spacer_1 CTAGGCGGCGGCCGCGGCGGCGG              A        24            0
    ##   spacer_1 CTAGGCGGCGGCCGCGGCGGCGG              A        31            0
    ##            canonical  cut_site         gene_id gene_symbol
    ##            <logical> <numeric>     <character> <character>
    ##   spacer_1      TRUE        NA ENSG00000133703        KRAS
    ##   spacer_1      TRUE        NA ENSG00000133703        KRAS
    ##   spacer_1      TRUE        NA ENSG00000133703        KRAS
    ##   spacer_1      TRUE        NA ENSG00000133703        KRAS
    ##   -------
    ##   seqinfo: 7514 sequences from an unspecified genome; no seqlengths

Note that each annotated alignment is specific to the transcript ID
given under `seqnames`.

Below is a spacer that targets (with no mismatches) multiple genes.

``` r
results["spacer_244"]
```

    ## GuideSet object with 1 range and 10 metadata columns:
    ##              seqnames    ranges strand |             protospacer            pam
    ##                 <Rle> <IRanges>  <Rle> |          <DNAStringSet> <DNAStringSet>
    ##   spacer_244 region_1       267      + | CTTGACGATACAGCTAATTCAGA              A
    ##               pam_site  cut_site      region     n0_tx     n1_tx   n0_gene
    ##              <numeric> <numeric> <character> <numeric> <numeric> <numeric>
    ##   spacer_244       267        NA    region_1         5         0         2
    ##                n1_gene
    ##              <numeric>
    ##   spacer_244         0
    ##                                                                        alignments
    ##                                                                     <GRangesList>
    ##   spacer_244 ENST00000256078:267:+,ENST00000311936:267:+,ENST00000407852:77:+,...
    ##   -------
    ##   seqinfo: 1 sequence from custom genome
    ##   crisprNuclease: CasRx

Upon further inspection of this spacer’s alignments, however, we can see
that the off-target occurs in the pseudogene KRASP1, and should be
harmless.

``` r
onTargets(results["spacer_244"])
```

    ## GRanges object with 5 ranges and 9 metadata columns:
    ##                     seqnames    ranges strand |                  spacer
    ##                        <Rle> <IRanges>  <Rle> |          <DNAStringSet>
    ##   spacer_244 ENST00000256078       267      + | TCTGAATTAGCTGTATCGTCAAG
    ##   spacer_244 ENST00000311936       267      + | TCTGAATTAGCTGTATCGTCAAG
    ##   spacer_244 ENST00000407852        77      + | TCTGAATTAGCTGTATCGTCAAG
    ##   spacer_244 ENST00000556131       254      + | TCTGAATTAGCTGTATCGTCAAG
    ##   spacer_244 ENST00000557334       274      + | TCTGAATTAGCTGTATCGTCAAG
    ##                          protospacer            pam  pam_site n_mismatches
    ##                       <DNAStringSet> <DNAStringSet> <numeric>    <integer>
    ##   spacer_244 CTTGACGATACAGCTAATTCAGA              A       267            0
    ##   spacer_244 CTTGACGATACAGCTAATTCAGA              A       267            0
    ##   spacer_244 CTTGACGATACAGCTAATTCAGA              A        77            0
    ##   spacer_244 CTTGACGATACAGCTAATTCAGA              A       254            0
    ##   spacer_244 CTTGACGATACAGCTAATTCAGA              A       274            0
    ##              canonical  cut_site         gene_id gene_symbol
    ##              <logical> <numeric>     <character> <character>
    ##   spacer_244      TRUE        NA ENSG00000133703        KRAS
    ##   spacer_244      TRUE        NA ENSG00000133703        KRAS
    ##   spacer_244      TRUE        NA ENSG00000220635      KRASP1
    ##   spacer_244      TRUE        NA ENSG00000133703        KRAS
    ##   spacer_244      TRUE        NA ENSG00000133703        KRAS
    ##   -------
    ##   seqinfo: 7514 sequences from an unspecified genome; no seqlengths

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
    ## [11] crisprDesignData_0.99.13          crisprDesign_0.99.112            
    ## [13] crisprBase_1.1.2                 
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] bitops_1.0-7                  matrixStats_0.62.0           
    ##   [3] bit64_4.0.5                   filelock_1.0.2               
    ##   [5] progress_1.2.2                httr_1.4.3                   
    ##   [7] tools_4.2.0                   utf8_1.2.2                   
    ##   [9] R6_2.5.1                      DBI_1.1.3                    
    ##  [11] tidyselect_1.1.2              prettyunits_1.1.1            
    ##  [13] bit_4.0.4                     curl_4.3.2                   
    ##  [15] compiler_4.2.0                crisprBowtie_1.0.0           
    ##  [17] cli_3.3.0                     Biobase_2.56.0               
    ##  [19] basilisk.utils_1.9.1          crisprScoreData_1.1.3        
    ##  [21] xml2_1.3.3                    DelayedArray_0.22.0          
    ##  [23] randomForest_4.7-1.1          readr_2.1.2                  
    ##  [25] rappdirs_0.3.3                stringr_1.4.0                
    ##  [27] digest_0.6.29                 Rsamtools_2.12.0             
    ##  [29] rmarkdown_2.14                crisprScore_1.1.14           
    ##  [31] basilisk_1.9.2                pkgconfig_2.0.3              
    ##  [33] htmltools_0.5.3               MatrixGenerics_1.8.1         
    ##  [35] dbplyr_2.2.1                  fastmap_1.1.0                
    ##  [37] rlang_1.0.4                   rstudioapi_0.13              
    ##  [39] RSQLite_2.2.15                shiny_1.7.2                  
    ##  [41] BiocIO_1.6.0                  generics_0.1.3               
    ##  [43] jsonlite_1.8.0                vroom_1.5.7                  
    ##  [45] BiocParallel_1.30.3           dplyr_1.0.9                  
    ##  [47] VariantAnnotation_1.42.1      RCurl_1.98-1.8               
    ##  [49] magrittr_2.0.3                GenomeInfoDbData_1.2.8       
    ##  [51] Matrix_1.4-1                  Rcpp_1.0.9                   
    ##  [53] fansi_1.0.3                   reticulate_1.25              
    ##  [55] Rbowtie_1.36.0                lifecycle_1.0.1              
    ##  [57] stringi_1.7.8                 yaml_2.3.5                   
    ##  [59] SummarizedExperiment_1.26.1   zlibbioc_1.42.0              
    ##  [61] BiocFileCache_2.4.0           AnnotationHub_3.4.0          
    ##  [63] grid_4.2.0                    blob_1.2.3                   
    ##  [65] promises_1.2.0.1              parallel_4.2.0               
    ##  [67] ExperimentHub_2.4.0           crayon_1.5.1                 
    ##  [69] crisprBwa_1.0.0               dir.expiry_1.4.0             
    ##  [71] lattice_0.20-45               GenomicFeatures_1.48.3       
    ##  [73] hms_1.1.1                     KEGGREST_1.36.3              
    ##  [75] knitr_1.39                    pillar_1.8.0                 
    ##  [77] rjson_0.2.21                  codetools_0.2-18             
    ##  [79] biomaRt_2.52.0                BiocVersion_3.15.2           
    ##  [81] XML_3.99-0.10                 glue_1.6.2                   
    ##  [83] evaluate_0.15                 BiocManager_1.30.18          
    ##  [85] httpuv_1.6.5                  png_0.1-7                    
    ##  [87] vctrs_0.4.1                   tzdb_0.3.0                   
    ##  [89] purrr_0.3.4                   assertthat_0.2.1             
    ##  [91] cachem_1.0.6                  xfun_0.31                    
    ##  [93] mime_0.12                     Rbwa_1.0.0                   
    ##  [95] xtable_1.8-4                  restfulr_0.0.15              
    ##  [97] later_1.3.0                   tibble_3.1.8                 
    ##  [99] GenomicAlignments_1.32.1      AnnotationDbi_1.58.0         
    ## [101] memoise_2.0.1                 interactiveDisplayBase_1.34.0
    ## [103] ellipsis_0.3.2

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
