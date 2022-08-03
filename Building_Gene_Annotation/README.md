Building a gene annotation object
================

-   [Introduction](#introduction)
-   [Installation](#installation)
    -   [Software requirements](#software-requirements)
        -   [OS Requirements](#os-requirements)
    -   [Installation](#installation-1)
        -   [Getting started](#getting-started)
-   [Building gene annotation
    objects](#building-gene-annotation-objects)
    -   [Building a txObject](#building-a-txobject)
    -   [Building a tssObject](#building-a-tssobject)
-   [Using gene annotation objects](#using-gene-annotation-objects)
    -   [Using txObjects](#using-txobjects)
    -   [Using tssObjects](#using-tssobjects)
-   [License](#license)
-   [Reproducibility](#reproducibility)

Authors: Jean-Philippe Fortin, Luke Hoberecht

Date: 02 August, 2022

# Introduction

In this tutorial we describe the process for making and using a gene
annotation object, which is required input for certain `crisprDesign`
functions such as `queryTxObject` and `addGeneAnnotation` (it is also
referred to as a `txObject`, after the argument name the object is
passed to for applicable functions). This object is a `GRangesList` that
is composed from a `TxDb` object to provide convenience and speed when
retrieving gene annotation with `crisprDesign` functions. While a `TxDb`
can be passed to `txObject` arguments, we recommend constructing and
using a gene annotation object described herein. We will also describe
the process for constructing and using a transcription start site (TSS)
annotation object (similarly also referred to as a `tssObject`), which
contains TSS coordinates and is a `GRanges` rather than a `GRangesList`.

# Installation

## Software requirements

### OS Requirements

This package is supported for macOS, Linux and Windows machines. It was
developed and tested on R version 4.2.

## Installation

This tutorial will use the `crisprDesign` package to demonstrate how to
construct and retrieve information from gene annotaion objects. Examples
using such objects will take advantages of the premade gene annotation
objects in the `crisprDesignData` package. Both packages can be
installed by typing the following commands in an R session:

``` r
install.packages("devtools")
devtools::install_github("Jfortin1/crisprDesign")
devtools::install_github("Jfortin1/crisprDesignData")
```

### Getting started

The packages can be loaded into an R session in the usual way:

``` r
library(crisprDesign)
library(crisprDesignData)
```

# Building gene annotation objects

The following demonstration will first create a `TxDb` object from the
Ensembl database with `getTxDb` from `crisprDesign`, which uses the
`GenomicFeatures` function `makeTxDbFromEnsembl`, before converting it
to a gene annotation object. Note that this method requires an internet
connection. If you have a General Feature Format (GFF) file from which
you want to construct the gene annotation object, you can pass this to
the `file` argument of the `crisprDesign` function `getTxDb`; this will
create the `TxDb` object using the `GenomicFeatures` function
`makeTxDbFromGFF`.

## Building a txObject

In this example, we construct a gene annotation object for the human
genome, Ensembl release 104 (hg38). We then show how to subset this
object to only contain annotation for our gene of interest, IQSEC3,
although this method also works for an arbitrary list of genes.
Subsetting is useful if memory constraints are a concern and only
annotation for the gene(s) of interest is(are) desired.

``` r
library(GenomicRanges)

## genome-wide annotation
txdb <- getTxDb(organism="Homo sapiens",
                release=104)
grList <- TxDb2GRangesList(txdb)

## subset for IQSEC3
meta <- metadata(grList)
gene_id <- "ENSG00000120645" #IQSEC3
grListExample <- lapply(grList, function(gr){
    gr[gr$gene_id %in% gene_id]
})
grListExample <- GRangesList(grListExample)
metadata(grListExample) <- meta
GenomeInfoDb::genome(grListExample) <- "hg38"
```

Let’s look at our gene annotation object:

``` r
names(grListExample)
```

    ## [1] "transcripts" "exons"       "cds"         "fiveUTRs"    "threeUTRs"  
    ## [6] "introns"     "tss"

``` r
grListExample$transcripts
```

    ## GRanges object with 2 ranges and 14 metadata columns:
    ##    seqnames       ranges strand |           tx_id         gene_id
    ##       <Rle>    <IRanges>  <Rle> |     <character>     <character>
    ##          12 66767-178455      + | ENST00000538872 ENSG00000120645
    ##          12 77376-171330      + | ENST00000382841 ENSG00000120645
    ##         protein_id        tx_type gene_symbol     exon_id exon_rank cds_start
    ##        <character>    <character> <character> <character> <integer> <integer>
    ##    ENSP00000437554 protein_coding      IQSEC3        <NA>      <NA>      <NA>
    ##    ENSP00000372292 protein_coding      IQSEC3        <NA>      <NA>      <NA>
    ##      cds_end  tx_start    tx_end   cds_len exon_start  exon_end
    ##    <integer> <integer> <integer> <integer>  <integer> <integer>
    ##         <NA>      <NA>      <NA>      <NA>       <NA>      <NA>
    ##         <NA>      <NA>      <NA>      <NA>       <NA>      <NA>
    ##   -------
    ##   seqinfo: 25 sequences (1 circular) from hg38 genome

## Building a tssObject

Building a `tssObject` requires only one additional step after
constructing the `txObject` described above, with the convenience
function `getTssObjectFromTxObject` in `crisprDesign`. Let’s apply the
function to our subset example:

``` r
tssObjectExample <- getTssObjectFromTxObject(grListExample)
tssObjectExample
```

    ## GRanges object with 2 ranges and 5 metadata columns:
    ##         seqnames    ranges strand |           tx_id         gene_id gene_symbol
    ##            <Rle> <IRanges>  <Rle> |     <character>     <character> <character>
    ##   76436       12     66767      + | ENST00000538872 ENSG00000120645      IQSEC3
    ##   76437       12     77376      + | ENST00000382841 ENSG00000120645      IQSEC3
    ##                promoter                     ID
    ##             <character>            <character>
    ##   76436 ENST00000538872 ENSG00000120645_ENST..
    ##   76437 ENST00000382841 ENSG00000120645_ENST..
    ##   -------
    ##   seqinfo: 25 sequences (1 circular) from hg38 genome

# Using gene annotation objects

Gene (or TSS) annotation objects are necessary for the full
characterization of guide-RNAs (gRNAs) and serve as input for several
`crisprDesign` functions, including `queryTxObject`, `queryTssObject`,
`addGeneAnnotation`, `addTssAnnotation`, and `addSpacerAlignments`. The
following sections demonstrate the information it contains and how to
use the object. We will take advantage of precomputed gene annotation
objects available in the `crisprDesignData` package:

| Object name  | Object class  | Version     | Description                                           |
|--------------|---------------|-------------|-------------------------------------------------------|
| `txdb_human` | `GRangesList` | Release 104 | Ensembl gene model for human (hg38/GRCh38)            |
| `txdb_mouse` | `GRangesList` | Release 102 | Ensembl gene model for mouse (mm10/GRCm38)            |
| `tss_human`  | `GRanges`     | Release 104 | Ensembl-based TSS coordinates for human (hg38/GRCh38) |
| `tss_mouse`  | `GRanges`     | Release 102 | Ensembl-based TSS coordinates for human (mm10/GRCm38) |

## Using txObjects

Let’s first look at the `txdb_human` object.

``` r
data(txdb_human, package="crisprDesignData")
```

We can look at metadata information about the gene model by using the
`metadata` function from the `S4Vectors` package:

``` r
head(S4Vectors::metadata(txdb_human))
##                 name                    value
## 1            Db type                     TxDb
## 2 Supporting package          GenomicFeatures
## 3        Data source                  Ensembl
## 4           Organism             Homo sapiens
## 5    Ensembl release                      104
## 6   Ensembl database homo_sapiens_core_104_38
```

The object is a `GRangesList` with 7 elements that contain genomic
coordinates for different levels of the gene model:

``` r
names(txdb_human)
## [1] "transcripts" "exons"       "cds"         "fiveUTRs"    "threeUTRs"  
## [6] "introns"     "tss"
```

As an example, let’s look at the `GRanges` containing genomic
coordinates for all exons represented in the gene model:

``` r
txdb_human$exons
## GRanges object with 796644 ranges and 14 metadata columns:
##     seqnames      ranges strand |           tx_id         gene_id
##        <Rle>   <IRanges>  <Rle> |     <character>     <character>
##         chr1 11869-12227      + | ENST00000456328 ENSG00000223972
##         chr1 12613-12721      + | ENST00000456328 ENSG00000223972
##         chr1 13221-14409      + | ENST00000456328 ENSG00000223972
##         chr1 12010-12057      + | ENST00000450305 ENSG00000223972
##         chr1 12179-12227      + | ENST00000450305 ENSG00000223972
##   .      ...         ...    ... .             ...             ...
##         chrM   5826-5891      - | ENST00000387409 ENSG00000210144
##         chrM   7446-7514      - | ENST00000387416 ENSG00000210151
##         chrM 14149-14673      - | ENST00000361681 ENSG00000198695
##         chrM 14674-14742      - | ENST00000387459 ENSG00000210194
##         chrM 15956-16023      - | ENST00000387461 ENSG00000210196
##          protein_id                tx_type gene_symbol         exon_id
##         <character>            <character> <character>     <character>
##                <NA>   processed_transcript     DDX11L1 ENSE00002234944
##                <NA>   processed_transcript     DDX11L1 ENSE00003582793
##                <NA>   processed_transcript     DDX11L1 ENSE00002312635
##                <NA> transcribed_unproces..     DDX11L1 ENSE00001948541
##                <NA> transcribed_unproces..     DDX11L1 ENSE00001671638
##   .             ...                    ...         ...             ...
##                <NA>                Mt_tRNA       MT-TY ENSE00001544488
##                <NA>                Mt_tRNA      MT-TS1 ENSE00001544487
##     ENSP00000354665         protein_coding      MT-ND6 ENSE00001434974
##                <NA>                Mt_tRNA       MT-TE ENSE00001544476
##                <NA>                Mt_tRNA       MT-TP ENSE00001544473
##     exon_rank cds_start   cds_end  tx_start    tx_end   cds_len exon_start
##     <integer> <integer> <integer> <integer> <integer> <integer>  <integer>
##             1      <NA>      <NA>     11869     14409         0       <NA>
##             2      <NA>      <NA>     11869     14409         0       <NA>
##             3      <NA>      <NA>     11869     14409         0       <NA>
##             1      <NA>      <NA>     12010     13670         0       <NA>
##             2      <NA>      <NA>     12010     13670         0       <NA>
##   .       ...       ...       ...       ...       ...       ...        ...
##             1      <NA>      <NA>      5826      5891         0       <NA>
##             1      <NA>      <NA>      7446      7514         0       <NA>
##             1     14149     14673     14149     14673       525       <NA>
##             1      <NA>      <NA>     14674     14742         0       <NA>
##             1      <NA>      <NA>     15956     16023         0       <NA>
##      exon_end
##     <integer>
##          <NA>
##          <NA>
##          <NA>
##          <NA>
##          <NA>
##   .       ...
##          <NA>
##          <NA>
##          <NA>
##          <NA>
##          <NA>
##   -------
##   seqinfo: 25 sequences (1 circular) from hg38 genome
```

The function `queryTxObject` in `crisprDesign` is a user-friendly
function to work with such objects. For instance we can find the CDS
coordinates for the KRAS transcripts using the following lines of code:

``` r
cds <- queryTxObject(txdb_human,
                     featureType="cds",
                     queryColumn="gene_symbol",
                     queryValue="KRAS")
head(cds)
## GRanges object with 6 ranges and 14 metadata columns:
##            seqnames            ranges strand |           tx_id         gene_id
##               <Rle>         <IRanges>  <Rle> |     <character>     <character>
##   region_1    chr12 25245274-25245384      - | ENST00000256078 ENSG00000133703
##   region_2    chr12 25227234-25227412      - | ENST00000256078 ENSG00000133703
##   region_3    chr12 25225614-25225773      - | ENST00000256078 ENSG00000133703
##   region_4    chr12 25215441-25215560      - | ENST00000256078 ENSG00000133703
##   region_5    chr12 25245274-25245384      - | ENST00000311936 ENSG00000133703
##   region_6    chr12 25227234-25227412      - | ENST00000311936 ENSG00000133703
##                 protein_id        tx_type gene_symbol         exon_id exon_rank
##                <character>    <character> <character>     <character> <integer>
##   region_1 ENSP00000256078 protein_coding        KRAS ENSE00000936617         2
##   region_2 ENSP00000256078 protein_coding        KRAS ENSE00001719809         3
##   region_3 ENSP00000256078 protein_coding        KRAS ENSE00001644818         4
##   region_4 ENSP00000256078 protein_coding        KRAS ENSE00001189807         5
##   region_5 ENSP00000256078 protein_coding        KRAS ENSE00000936617         2
##   region_6 ENSP00000256078 protein_coding        KRAS ENSE00001719809         3
##            cds_start   cds_end  tx_start    tx_end   cds_len exon_start
##            <integer> <integer> <integer> <integer> <integer>  <integer>
##   region_1      <NA>      <NA>  25205246  25250929       570   25245274
##   region_2      <NA>      <NA>  25205246  25250929       570   25227234
##   region_3      <NA>      <NA>  25205246  25250929       570   25225614
##   region_4      <NA>      <NA>  25205246  25250929       570   25215437
##   region_5      <NA>      <NA>  25205246  25250929       567   25245274
##   region_6      <NA>      <NA>  25205246  25250929       567   25227234
##             exon_end
##            <integer>
##   region_1  25245395
##   region_2  25227412
##   region_3  25225773
##   region_4  25215560
##   region_5  25245395
##   region_6  25227412
##   -------
##   seqinfo: 25 sequences (1 circular) from hg38 genome
```

We can also pass `txdb_human` to the `txObject` parameter of
`addGeneAnnotation` (required) and `addSpacerAlignments` (optional) to
add gene annotation information to a `GuideSet`, provided the `genome`
of the `GuideSet` is also `hg38`. Below is a simple example using
`txdb_human` in `addGeneAnnotation`:

``` r
library(BSgenome.Hsapiens.UCSC.hg38)
gs <- findSpacers(head(cds),
                  crisprNuclease=SpCas9,
                  bsgenome=BSgenome.Hsapiens.UCSC.hg38)
gs <- head(gs)
```

``` r
## add gene annnotation
gs1 <- addGeneAnnotation(gs,
                         txObject=txdb_human)
geneAnnotation(gs1)
## DataFrame with 6 rows and 23 columns
##               chr anchor_site   strand gene_symbol         gene_id
##          <factor>   <integer> <factor> <character>     <character>
## spacer_1    chr12    25215441        -        KRAS ENSG00000133703
## spacer_2    chr12    25215480        -        KRAS ENSG00000133703
## spacer_3    chr12    25215474        +        KRAS ENSG00000133703
## spacer_4    chr12    25215517        +        KRAS ENSG00000133703
## spacer_5    chr12    25215538        -        KRAS ENSG00000133703
## spacer_6    chr12    25215556        -        KRAS ENSG00000133703
##                    tx_id      protein_id   cut_cds cut_fiveUTRs cut_threeUTRs
##              <character>     <character> <logical>    <logical>     <logical>
## spacer_1 ENST00000256078 ENSP00000256078      TRUE        FALSE         FALSE
## spacer_2 ENST00000256078 ENSP00000256078      TRUE        FALSE         FALSE
## spacer_3 ENST00000256078 ENSP00000256078      TRUE        FALSE         FALSE
## spacer_4 ENST00000256078 ENSP00000256078      TRUE        FALSE         FALSE
## spacer_5 ENST00000256078 ENSP00000256078      TRUE        FALSE         FALSE
## spacer_6 ENST00000256078 ENSP00000256078      TRUE        FALSE         FALSE
##          cut_introns percentCDS aminoAcidIndex downtreamATG percentTx nIsoforms
##            <logical>  <numeric>      <numeric>    <numeric> <numeric> <integer>
## spacer_1       FALSE      100.0            190            0      14.0         1
## spacer_2       FALSE       93.2            177            1      13.3         1
## spacer_3       FALSE       94.2            179            1      13.4         1
## spacer_4       FALSE       86.7            165            1      12.6         1
## spacer_5       FALSE       83.0            158            1      12.2         1
## spacer_6       FALSE       79.8            152            1      11.9         1
##          totalIsoforms percentIsoforms isCommonExon nCodingIsoforms
##              <numeric>       <numeric>    <logical>       <integer>
## spacer_1             4              25        FALSE               1
## spacer_2             4              25        FALSE               1
## spacer_3             4              25        FALSE               1
## spacer_4             4              25        FALSE               1
## spacer_5             4              25        FALSE               1
## spacer_6             4              25        FALSE               1
##          totalCodingIsoforms percentCodingIsoforms isCommonCodingExon
##                    <numeric>             <numeric>          <logical>
## spacer_1                   4                    25              FALSE
## spacer_2                   4                    25              FALSE
## spacer_3                   4                    25              FALSE
## spacer_4                   4                    25              FALSE
## spacer_5                   4                    25              FALSE
## spacer_6                   4                    25              FALSE
```

For examples using a gene annotation object in `addSpacerAlignments`,
see the tutorials for [CRISPRko design with
Cas9](https://github.com/crisprVerse/Tutorials/tree/master/Design_CRISPRko_Cas9)
or [CRISPRko design with
Cas12a](https://github.com/crisprVerse/Tutorials/tree/master/Design_CRISPRko_Cas12a).

## Using tssObjects

The `tss_human` and `tss_mouse` objects are `GRanges` representing the
TSS coordinates for human and mouse, respectively. The coordinates were
extracted from the transcripts stored in the Ensembl-based models
`txdb_human` and `txdb_mouse` using the function
`getTssObjectFromTxObject` from `crisprDesign`.

Let’s take a look at `tss_human`:

``` r
data(tss_human, package="crisprDesignData")
head(tss_human)
## GRanges object with 6 ranges and 9 metadata columns:
##                      seqnames    ranges strand |     score peak_start  peak_end
##                         <Rle> <IRanges>  <Rle> | <numeric>  <integer> <integer>
##   ENSG00000000003_P1     chrX 100636805      - |   4.35417  100636805 100636805
##   ENSG00000000005_P1     chrX 100584935      + |   3.29137  100584935 100584935
##   ENSG00000000419_P1    chr20  50958531      - |   5.74747   50958531  50958531
##   ENSG00000000457_P1     chr1 169893895      - |   4.75432  169893895 169893895
##   ENSG00000000460_P1     chr1 169795044      + |   4.92777  169795044 169795044
##   ENSG00000000938_P1     chr1  27635184      - |   4.61214   27635184  27635184
##                                tx_id         gene_id      source    promoter
##                          <character>     <character> <character> <character>
##   ENSG00000000003_P1 ENST00000373020 ENSG00000000003     fantom5          P1
##   ENSG00000000005_P1 ENST00000373031 ENSG00000000005     fantom5          P1
##   ENSG00000000419_P1 ENST00000371588 ENSG00000000419     fantom5          P1
##   ENSG00000000457_P1 ENST00000367771 ENSG00000000457     fantom5          P1
##   ENSG00000000460_P1 ENST00000359326 ENSG00000000460     fantom5          P1
##   ENSG00000000938_P1 ENST00000374005 ENSG00000000938     fantom5          P1
##                                      ID gene_symbol
##                             <character> <character>
##   ENSG00000000003_P1 ENSG00000000003_P1      TSPAN6
##   ENSG00000000005_P1 ENSG00000000005_P1        TNMD
##   ENSG00000000419_P1 ENSG00000000419_P1        DPM1
##   ENSG00000000457_P1 ENSG00000000457_P1       SCYL3
##   ENSG00000000460_P1 ENSG00000000460_P1    C1orf112
##   ENSG00000000938_P1 ENSG00000000938_P1         FGR
##   -------
##   seqinfo: 25 sequences from an unspecified genome; no seqlengths
```

The function `queryTss` in `crisprDesign` is a user-friendly function to
work with such objects, accepting an argument called `tss_window` to
specify a number of nucleotides upstream and downstream of the TSS. This
is particularly useful to return genomic regions to target for CRISPRa
and CRISPRi.

For instance, if we want to target the region 500 nucleotides upstream
of any of the KRAS TSSs, one can use the following lines of code:

``` r
library(crisprDesign)
tss <- queryTss(tss_human,
                queryColumn="gene_symbol",
                queryValue="KRAS",
                tss_window=c(-500,0))
tss
## GRanges object with 1 range and 9 metadata columns:
##            seqnames            ranges strand |     score peak_start  peak_end
##               <Rle>         <IRanges>  <Rle> | <numeric>  <integer> <integer>
##   region_1    chr12 25250929-25251428      - |   5.20187   25250928  25250928
##                      tx_id         gene_id      source    promoter
##                <character>     <character> <character> <character>
##   region_1 ENST00000256078 ENSG00000133703     fantom5          P1
##                            ID gene_symbol
##                   <character> <character>
##   region_1 ENSG00000133703_P1        KRAS
##   -------
##   seqinfo: 25 sequences from an unspecified genome; no seqlengths
```

Similar to `txObject`s and the `addGeneAnnotation` function,
`tssObject`s can be passed to the `crisprDesign` function
`addTssAnnotation`, and is also an optional argument for
`addSpacerAlignments`. The functions append TSS annotation to the
`GuideSet` with respect to gRNAs and their alignments, respectively; see
the tutorials [CRISPRa
design](https://github.com/crisprVerse/Tutorials/tree/master/Design_CRISPRa)
or [CRISPRi
design](https://github.com/crisprVerse/Tutorials/tree/master/Design_CRISPRi)
for more information.

# License

The package is licensed under the MIT license.

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
    ##   [7] bit64_4.0.5                   lubridate_1.8.0              
    ##   [9] interactiveDisplayBase_1.34.0 AnnotationDbi_1.58.0         
    ##  [11] fansi_1.0.3                   xml2_1.3.3                   
    ##  [13] codetools_0.2-18              cachem_1.0.6                 
    ##  [15] knitr_1.39                    pkgload_1.3.0                
    ##  [17] jsonlite_1.8.0                Rsamtools_2.12.0             
    ##  [19] dbplyr_2.2.1                  png_0.1-7                    
    ##  [21] shiny_1.7.2                   BiocManager_1.30.18          
    ##  [23] readr_2.1.2                   compiler_4.2.0               
    ##  [25] httr_1.4.3                    basilisk_1.8.0               
    ##  [27] assertthat_0.2.1              Matrix_1.4-1                 
    ##  [29] fastmap_1.1.0                 cli_3.3.0                    
    ##  [31] later_1.3.0                   htmltools_0.5.3              
    ##  [33] prettyunits_1.1.1             tools_4.2.0                  
    ##  [35] glue_1.6.2                    GenomeInfoDbData_1.2.8       
    ##  [37] crisprBowtie_1.0.0            dplyr_1.0.9                  
    ##  [39] rappdirs_0.3.3                Rcpp_1.0.9                   
    ##  [41] Biobase_2.56.0                vctrs_0.4.1                  
    ##  [43] ExperimentHub_2.4.0           crisprScore_1.1.6            
    ##  [45] xfun_0.31                     stringr_1.4.0                
    ##  [47] ps_1.7.1                      mime_0.12                    
    ##  [49] miniUI_0.1.1.1                lifecycle_1.0.1              
    ##  [51] restfulr_0.0.15               devtools_2.4.4               
    ##  [53] XML_3.99-0.10                 AnnotationHub_3.4.0          
    ##  [55] basilisk.utils_1.8.0          zlibbioc_1.42.0              
    ##  [57] VariantAnnotation_1.42.1      hms_1.1.1                    
    ##  [59] promises_1.2.0.1              MatrixGenerics_1.8.1         
    ##  [61] parallel_4.2.0                SummarizedExperiment_1.26.1  
    ##  [63] RMariaDB_1.2.2                yaml_2.3.5                   
    ##  [65] curl_4.3.2                    reticulate_1.25              
    ##  [67] memoise_2.0.1                 biomaRt_2.52.0               
    ##  [69] stringi_1.7.8                 RSQLite_2.2.15               
    ##  [71] BiocVersion_3.15.2            BiocIO_1.6.0                 
    ##  [73] randomForest_4.7-1.1          crisprScoreData_1.1.3        
    ##  [75] GenomicFeatures_1.48.3        filelock_1.0.2               
    ##  [77] pkgbuild_1.3.1                BiocParallel_1.30.3          
    ##  [79] rlang_1.0.4                   pkgconfig_2.0.3              
    ##  [81] matrixStats_0.62.0            bitops_1.0-7                 
    ##  [83] evaluate_0.15                 lattice_0.20-45              
    ##  [85] purrr_0.3.4                   GenomicAlignments_1.32.1     
    ##  [87] htmlwidgets_1.5.4             bit_4.0.4                    
    ##  [89] processx_3.7.0                tidyselect_1.1.2             
    ##  [91] magrittr_2.0.3                R6_2.5.1                     
    ##  [93] generics_0.1.3                profvis_0.3.7                
    ##  [95] DelayedArray_0.22.0           DBI_1.1.3                    
    ##  [97] pillar_1.8.0                  KEGGREST_1.36.3              
    ##  [99] RCurl_1.98-1.8                dir.expiry_1.4.0             
    ## [101] tibble_3.1.8                  crayon_1.5.1                 
    ## [103] utf8_1.2.2                    BiocFileCache_2.4.0          
    ## [105] tzdb_0.3.0                    rmarkdown_2.14               
    ## [107] urlchecker_1.0.1              progress_1.2.2               
    ## [109] usethis_2.1.6                 grid_4.2.0                   
    ## [111] blob_1.2.3                    callr_3.7.1                  
    ## [113] digest_0.6.29                 xtable_1.8-4                 
    ## [115] httpuv_1.6.5                  sessioninfo_1.2.2
