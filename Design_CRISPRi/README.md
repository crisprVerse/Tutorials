Using crisprDesign to design gRNAs for CRISPRi
================

-   [Introduction](#introduction)
-   [Installation](#installation)
-   [Terminology](#terminology)
-   [CRISPRi design](#crispri-design)
    -   [Creating the GuideSet](#creating-the-guideset)
    -   [Annotating the GuideSet](#annotating-the-guideset)
    -   [Adding TSS annotation](#adding-tss-annotation)
    -   [Adding spacer alignments with TSS
        annotation](#adding-spacer-alignments-with-tss-annotation)
    -   [Adding CRISPRai scores](#adding-crisprai-scores)
-   [Session Info](#session-info)
-   [References](#references)

Authors: Jean-Philippe Fortin, Luke Hoberecht

Date: 10 August, 2022

# Introduction

This tutorial will demonstrate how to use `crisprDesign` to design gRNAs
for CRISPR interference (CRISPRi). Specifically, it will target the
human KRAS gene and use the SpCas9 nuclease, however, much of the same
process can be applied to any genomic target(s) and with any CRISPR
nuclease.

# Installation

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

# Terminology

See the [CRISPRko design
vignette](https://github.com/crisprVerse/Tutorials/tree/master/Design_CRISPRko_Cas9)
to get familiar with the terminology used throughout this tutorial.

# CRISPRi design

For CRISPR activation (CRISPRa) and interference (CRISPRi) applications,
the CRISPR nuclease is engineered to lose its endonuclease activity, and
should therefore not introduce double-stranded breaks (DSBs). We will
use the dead SpCas9 (dSpCas9) nuclease as an example here. Note that
users don’t have to distinguish between dSpCas9 and SpCas9 when
specifying the nuclease in `crisprDesign` and `crisprBase` as they do
not differ in terms of the characteristics stored in the
`CrisprNuclease` object.

In CRISPRi, fusing dSpCas9 with a Krüppel-associated box (KRAB) domain
has been shown to be effective at repressing transcription in mammalian
cells (Gilbert et al. et al. 2013). The dSpCas9-KRAB fused protein is a
commonly-used construct to conduct CRISPR inhibition (CRISPRi)
experiments. To achieve optimal inhibition, gRNAs are usually designed
targeting the region directly downstream of the gene transcription
starting site (TSS).

`crisprDesign` provides functionalities to be able to take into account
design rules that are specific to CRISPRi applications. The `queryTss`
function allows for specifying genomic coordinates of promoter regions.
The `addTssAnnotation` function annotates gRNAs for known TSSs, and
includes a column `dist_to_tss` that gives the distance in nucleotides
between the TSS position and the PAM site of the gRNA. For CRISPRi, we
recommend targeting the region 25-75bp region downstream of the TSS for
optimal inhibition; see Sanson et al. et al. (2018) for more
information. Finally, the function `addCrispraiScores` adds
CRISPRi-specific on-target activity scores based on the work of
(Horlbeck et al. et al. 2016).

## Creating the GuideSet

We first start by loading the required packages:

``` r
library(crisprBase)
library(crisprDesign)
library(crisprDesignData)
library(BSgenome.Hsapiens.UCSC.hg38)
```

To demonstrate CRISPRi design, we will design gRNAs to inhibit
expression of the human KRAS gene using the SpCas9 nuclease. To
accomplish this, we want our gRNAs to target the region downstream of
the KRAS TSS; let’s consider the window containing 500bp immediately
downstream of the TSS to explore candidate gRNAs.

We first need to retrieve the TSS coordinates for KRAS. These data are
conveniently stored in the `crisprDesignData` package as the dataset
`tss_human`. For more information on `tss_human` and how to create
similar TSS annotation objects, see the [Building a gene annotation
object](https://github.com/crisprVerse/Tutorials/tree/master/Building_Gene_Annotation)
tutorial.

We load the TSS coordinates stored in the `tss_human` object

``` r
data("tss_human", package="crisprDesignData")
```

and query for KRAS using the `queryTss` function from `crisprDesign`:

``` r
target_window <- c(0, 500)
target_region <- queryTss(tss_human,
                          queryColumn="gene_symbol",
                          queryValue="KRAS",
                          tss_window=target_window)
```

``` r
target_region
## GRanges object with 1 range and 9 metadata columns:
##            seqnames            ranges strand |     score peak_start  peak_end
##               <Rle>         <IRanges>  <Rle> | <numeric>  <integer> <integer>
##   region_1    chr12 25250429-25250928      - |   5.20187   25250928  25250928
##                      tx_id         gene_id      source    promoter
##                <character>     <character> <character> <character>
##   region_1 ENST00000256078 ENSG00000133703     fantom5          P1
##                            ID gene_symbol
##                   <character> <character>
##   region_1 ENSG00000133703_P1        KRAS
##   -------
##   seqinfo: 25 sequences from an unspecified genome; no seqlengths
```

We load the `crisprNuclease` object storing information about the SpCas9
nuclease from the `crisprBase` package:

``` r
data(SpCas9, package="crisprBase")
```

We then find all candidate protospacer sequences in our target region
with `findSpacers`:

``` r
gs <- findSpacers(target_region,
                  crisprNuclease=SpCas9,
                  bsgenome=BSgenome.Hsapiens.UCSC.hg38)
```

``` r
gs
## GuideSet object with 160 ranges and 5 metadata columns:
##              seqnames    ranges strand |          protospacer            pam
##                 <Rle> <IRanges>  <Rle> |       <DNAStringSet> <DNAStringSet>
##     spacer_1    chr12  25250434      - | GCCGCGGCTGGAGGCTTCTG            GGG
##     spacer_2    chr12  25250435      - | AGCCGCGGCTGGAGGCTTCT            GGG
##     spacer_3    chr12  25250436      - | GAGCCGCGGCTGGAGGCTTC            TGG
##     spacer_4    chr12  25250443      - | TCCCCGAGAGCCGCGGCTGG            AGG
##     spacer_5    chr12  25250446      - | TCCTCCCCGAGAGCCGCGGC            TGG
##          ...      ...       ...    ... .                  ...            ...
##   spacer_156    chr12  25250915      - | ATTTTCCTAGGCGGCGGCCG            CGG
##   spacer_157    chr12  25250916      + | CGCTGCTGCCTCCGCCGCCG            CGG
##   spacer_158    chr12  25250921      - | AGCTCGATTTTCCTAGGCGG            CGG
##   spacer_159    chr12  25250924      - | CGGAGCTCGATTTTCCTAGG            CGG
##   spacer_160    chr12  25250928      + | CGCCGCCGCGGCCGCCGCCT            AGG
##               pam_site  cut_site      region
##              <numeric> <numeric> <character>
##     spacer_1  25250434  25250437    region_1
##     spacer_2  25250435  25250438    region_1
##     spacer_3  25250436  25250439    region_1
##     spacer_4  25250443  25250446    region_1
##     spacer_5  25250446  25250449    region_1
##          ...       ...       ...         ...
##   spacer_156  25250915  25250918    region_1
##   spacer_157  25250916  25250913    region_1
##   spacer_158  25250921  25250924    region_1
##   spacer_159  25250924  25250927    region_1
##   spacer_160  25250928  25250925    region_1
##   -------
##   seqinfo: 640 sequences (1 circular) from hg38 genome
##   crisprNuclease: SpCas9
```

## Annotating the GuideSet

Next, we annotate our candidate gRNAs to assess quality. There are
several functions in `crisprDesign` that provide annotation for features
that are nonspecific to CRISPRi, for which we refer the reader to the
[CRISPRko design with
Cas9](https://github.com/crisprVerse/Tutorials/tree/master/Design_CRISPRko_Cas9)
tutorial for more information. The sections below will cover annotation
functions that are of particular interest to CRISPRi applications.

## Adding TSS annotation

As the name implies, the `addTssAnnotation` function annotates gRNAs
with TSS context such as the distance between the gRNA and the TSS, as
well as which TSS is targeted (many genes contain different TSSs
corresponding to different isoforms).

The function requires a `tssObject` object, and the `tss_window` values
that we used earlier to define the target region. We can then retrieve
the appended annotation with the accessor function `tssAnnotation`:

``` r
gs <- addTssAnnotation(gs,
                       tssObject=tss_human,
                       tss_window=target_window)
tssAnnotation(gs)
## DataFrame with 160 rows and 15 columns
##                 chr anchor_site   strand     score peak_start  peak_end
##            <factor>   <integer> <factor> <numeric>  <integer> <integer>
## spacer_1      chr12    25250437        -   5.20187   25250928  25250928
## spacer_2      chr12    25250438        -   5.20187   25250928  25250928
## spacer_3      chr12    25250439        -   5.20187   25250928  25250928
## spacer_4      chr12    25250446        -   5.20187   25250928  25250928
## spacer_5      chr12    25250449        -   5.20187   25250928  25250928
## ...             ...         ...      ...       ...        ...       ...
## spacer_156    chr12    25250918        -   5.20187   25250928  25250928
## spacer_157    chr12    25250913        +   5.20187   25250928  25250928
## spacer_158    chr12    25250924        -   5.20187   25250928  25250928
## spacer_159    chr12    25250927        -   5.20187   25250928  25250928
## spacer_160    chr12    25250925        +   5.20187   25250928  25250928
##                      tx_id         gene_id      source    promoter
##                <character>     <character> <character> <character>
## spacer_1   ENST00000256078 ENSG00000133703     fantom5          P1
## spacer_2   ENST00000256078 ENSG00000133703     fantom5          P1
## spacer_3   ENST00000256078 ENSG00000133703     fantom5          P1
## spacer_4   ENST00000256078 ENSG00000133703     fantom5          P1
## spacer_5   ENST00000256078 ENSG00000133703     fantom5          P1
## ...                    ...             ...         ...         ...
## spacer_156 ENST00000256078 ENSG00000133703     fantom5          P1
## spacer_157 ENST00000256078 ENSG00000133703     fantom5          P1
## spacer_158 ENST00000256078 ENSG00000133703     fantom5          P1
## spacer_159 ENST00000256078 ENSG00000133703     fantom5          P1
## spacer_160 ENST00000256078 ENSG00000133703     fantom5          P1
##                        tss_id gene_symbol  tss_strand   tss_pos dist_to_tss
##                   <character> <character> <character> <integer>   <numeric>
## spacer_1   ENSG00000133703_P1        KRAS           -  25250928         491
## spacer_2   ENSG00000133703_P1        KRAS           -  25250928         490
## spacer_3   ENSG00000133703_P1        KRAS           -  25250928         489
## spacer_4   ENSG00000133703_P1        KRAS           -  25250928         482
## spacer_5   ENSG00000133703_P1        KRAS           -  25250928         479
## ...                       ...         ...         ...       ...         ...
## spacer_156 ENSG00000133703_P1        KRAS           -  25250928          10
## spacer_157 ENSG00000133703_P1        KRAS           -  25250928          15
## spacer_158 ENSG00000133703_P1        KRAS           -  25250928           4
## spacer_159 ENSG00000133703_P1        KRAS           -  25250928           1
## spacer_160 ENSG00000133703_P1        KRAS           -  25250928           3
```

## Adding spacer alignments with TSS annotation

As with all CRISPR applications, potential off-targets effects are an
important concern in assessing gRNA quality. While this concern is
somewhat moderated for CRISPRi, since the dead CRISPR nuclease does not
make DSBs, we should be aware of off-targetsoccuring in the promoter
regions of other genes. This can be handled by passing our `tssObject`
to the `addSpacerAlignments` function. We will search for up to 2
mismatches and increase the size of our `tss_window` (which defines the
promoter region when searching for off-targets) to err on the safe side.

Similar to the CRISPRko design tutorial, we need to specify a Bowtie
index of the human referenge genome; see the [Building genome indices
for short read
aligners](https://github.com/crisprVerse/Tutorials/tree/master/Building_Genome_Indices)
tutorial to learn how to create such an index.

Here, we specify the index that was available to us when generating this
tutorial:

``` r
# Users need to specify the path of their bowtie index
index_path <- "/Users/fortinj2/crisprIndices/bowtie/hg38/hg38"
```

We are ready to add on- and off-target alignments:

``` r
gs <- addSpacerAlignments(gs,
                          aligner="bowtie",
                          aligner_index=index_path,
                          bsgenome=BSgenome.Hsapiens.UCSC.hg38,
                          n_mismatches=2,
                          tssObject=tss_human,
                          tss_window=c(-500, 2000))
```

``` r
gs
## GuideSet object with 160 ranges and 13 metadata columns:
##              seqnames    ranges strand |          protospacer            pam
##                 <Rle> <IRanges>  <Rle> |       <DNAStringSet> <DNAStringSet>
##     spacer_1    chr12  25250434      - | GCCGCGGCTGGAGGCTTCTG            GGG
##     spacer_2    chr12  25250435      - | AGCCGCGGCTGGAGGCTTCT            GGG
##     spacer_3    chr12  25250436      - | GAGCCGCGGCTGGAGGCTTC            TGG
##     spacer_4    chr12  25250443      - | TCCCCGAGAGCCGCGGCTGG            AGG
##     spacer_5    chr12  25250446      - | TCCTCCCCGAGAGCCGCGGC            TGG
##          ...      ...       ...    ... .                  ...            ...
##   spacer_156    chr12  25250915      - | ATTTTCCTAGGCGGCGGCCG            CGG
##   spacer_157    chr12  25250916      + | CGCTGCTGCCTCCGCCGCCG            CGG
##   spacer_158    chr12  25250921      - | AGCTCGATTTTCCTAGGCGG            CGG
##   spacer_159    chr12  25250924      - | CGGAGCTCGATTTTCCTAGG            CGG
##   spacer_160    chr12  25250928      + | CGCCGCCGCGGCCGCCGCCT            AGG
##               pam_site  cut_site      region        tssAnnotation        n0
##              <numeric> <numeric> <character> <SplitDataFrameList> <numeric>
##     spacer_1  25250434  25250437    region_1 chr12:25250437:-:...         1
##     spacer_2  25250435  25250438    region_1 chr12:25250438:-:...         1
##     spacer_3  25250436  25250439    region_1 chr12:25250439:-:...         1
##     spacer_4  25250443  25250446    region_1 chr12:25250446:-:...         1
##     spacer_5  25250446  25250449    region_1 chr12:25250449:-:...         1
##          ...       ...       ...         ...                  ...       ...
##   spacer_156  25250915  25250918    region_1 chr12:25250918:-:...         1
##   spacer_157  25250916  25250913    region_1 chr12:25250913:+:...         1
##   spacer_158  25250921  25250924    region_1 chr12:25250924:-:...         1
##   spacer_159  25250924  25250927    region_1 chr12:25250927:-:...         1
##   spacer_160  25250928  25250925    region_1 chr12:25250925:+:...         1
##                     n1        n2      n0_p      n1_p      n2_p
##              <numeric> <numeric> <numeric> <numeric> <numeric>
##     spacer_1         0         1         1         0         0
##     spacer_2         0         1         1         0         0
##     spacer_3         0         1         1         0         0
##     spacer_4         0         0         1         0         0
##     spacer_5         0         0         1         0         0
##          ...       ...       ...       ...       ...       ...
##   spacer_156         0         0         1         0         0
##   spacer_157         0        27         1         0        18
##   spacer_158         0         0         1         0         0
##   spacer_159         0         0         1         0         0
##   spacer_160        18       160         1        17       121
##                                                         alignments
##                                                      <GRangesList>
##     spacer_1                     chr12:25250434:-,chr7:155971918:-
##     spacer_2                     chr12:25250435:-,chr3:184602035:+
##     spacer_3                     chr12:25250436:-,chr3:184602034:+
##     spacer_4                                      chr12:25250443:-
##     spacer_5                                      chr12:25250446:-
##          ...                                                   ...
##   spacer_156                                      chr12:25250915:-
##   spacer_157 chr12:25250916:+,chr2:55050346:+,chr20:21397165:+,...
##   spacer_158                                      chr12:25250921:-
##   spacer_159                                      chr12:25250924:-
##   spacer_160  chr12:25250928:+,chr17:49361951:+,chr1:3069055:-,...
##   -------
##   seqinfo: 640 sequences (1 circular) from hg38 genome
##   crisprNuclease: SpCas9
```

Including a `tssObject` parameter in the `addSpacerAlignments` function
appends columns to the `GuideSet` that tallies the alignments restricted
to the defined (via `tss_window`) promoter regions: `n0_p`, `n1_p`, and
`n2_p` (the `_p` suffix denotes “promoter”).

## Adding CRISPRai scores

The CRISPRai algorithm was developed by the Weissman lab to score SpCas9
gRNAs for CRISPRa and CRISPRi applications for the human genome
(Horlbeck et al. et al. 2016). The function `addCrispraiScores`
implements this algorithm to add scores to the `GuideSet`. Compared to
other on-target scoring algorithms, it requires several additional
inputs:

-   The `gr` argument is the `GRanges` object derived from the
    `queryTss` function and used to create the `GuideSet` object. In our
    example, this is the object named `target_region`.
-   The `tssObject` argument is a `GRanges` object that contains TSS
    coordinates and annotation. It must also contain the following
    columns: `ID`, `promoter`, `tx_id`, and `gene_symbol`. Our
    `tssObject` in this instance is `tss_human`.
-   `geneCol` indicates which column of `tssObject` should be used as
    the unique gene identifier.
-   `modality` is the modality of the CRISPR application, in our case,
    `CRISPRi`.
-   `fastaFile` is the path of a FASTA file containing the sequence of
    the human reference genome in hg38 coordinates. This file is
    available
    [here](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz).
-   `chromatinFiles` is a vector of length 3 specifying the path of
    files containing the chromatin accessibility data needed for the
    algorithm in hg38 coordinates. The chromatin files can be downloaded
    from Zenodo [here](https://zenodo.org/record/6716721#.YrzCfS-cY4d).

We first prepare all needed inputs for `addCrispraiScores`. We start by
specifying the location of the FASTA file on our local machine:

``` r
fastaPath <- "/Users/fortinj2/crisprIndices/genomes/hg38/hg38.fa"
```

This corresponds to the path where the downloaded file from
[here](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)
is stored. Next, we specify the location of the chromatin files:

``` r
mnasePath <- "/Users/fortinj2/crisprIndices/chromatin/hg38/crispria_mnase_human_K562_hg38.bigWig"
dnasePath <- "/Users/fortinj2/crisprIndices/chromatin/hg38/crispria_dnase_human_K562_hg38.bigWig"
fairePath <- "/Users/fortinj2/crisprIndices/chromatin/hg38/crispria_faire_human_K562_hg38.bigWig"
chromatinFiles <- c(mnase=mnasePath,
                    dnase=dnasePath,
                    faire=fairePath)
```

This should correspond to the files that were downloaded from
[here](https://zenodo.org/record/6716721#.YrzCfS-cY4d).

We are now ready to add the scores:

``` r
results <- addCrispraiScores(gs,
                             gr=target_region,
                             tssObject=tss_human,
                             geneCol="gene_id",
                             modality="CRISPRi",
                             fastaFile=fastaPath,
                             chromatinFiles=chromatinFiles)
```

Let’s look at the results:

``` r
results
## GuideSet object with 160 ranges and 14 metadata columns:
##              seqnames    ranges strand |          protospacer            pam
##                 <Rle> <IRanges>  <Rle> |       <DNAStringSet> <DNAStringSet>
##     spacer_1    chr12  25250434      - | GCCGCGGCTGGAGGCTTCTG            GGG
##     spacer_2    chr12  25250435      - | AGCCGCGGCTGGAGGCTTCT            GGG
##     spacer_3    chr12  25250436      - | GAGCCGCGGCTGGAGGCTTC            TGG
##     spacer_4    chr12  25250443      - | TCCCCGAGAGCCGCGGCTGG            AGG
##     spacer_5    chr12  25250446      - | TCCTCCCCGAGAGCCGCGGC            TGG
##          ...      ...       ...    ... .                  ...            ...
##   spacer_156    chr12  25250915      - | ATTTTCCTAGGCGGCGGCCG            CGG
##   spacer_157    chr12  25250916      + | CGCTGCTGCCTCCGCCGCCG            CGG
##   spacer_158    chr12  25250921      - | AGCTCGATTTTCCTAGGCGG            CGG
##   spacer_159    chr12  25250924      - | CGGAGCTCGATTTTCCTAGG            CGG
##   spacer_160    chr12  25250928      + | CGCCGCCGCGGCCGCCGCCT            AGG
##               pam_site  cut_site      region        tssAnnotation        n0
##              <numeric> <numeric> <character> <SplitDataFrameList> <numeric>
##     spacer_1  25250434  25250437    region_1 chr12:25250437:-:...         1
##     spacer_2  25250435  25250438    region_1 chr12:25250438:-:...         1
##     spacer_3  25250436  25250439    region_1 chr12:25250439:-:...         1
##     spacer_4  25250443  25250446    region_1 chr12:25250446:-:...         1
##     spacer_5  25250446  25250449    region_1 chr12:25250449:-:...         1
##          ...       ...       ...         ...                  ...       ...
##   spacer_156  25250915  25250918    region_1 chr12:25250918:-:...         1
##   spacer_157  25250916  25250913    region_1 chr12:25250913:+:...         1
##   spacer_158  25250921  25250924    region_1 chr12:25250924:-:...         1
##   spacer_159  25250924  25250927    region_1 chr12:25250927:-:...         1
##   spacer_160  25250928  25250925    region_1 chr12:25250925:+:...         1
##                     n1        n2      n0_p      n1_p      n2_p
##              <numeric> <numeric> <numeric> <numeric> <numeric>
##     spacer_1         0         1         1         0         0
##     spacer_2         0         1         1         0         0
##     spacer_3         0         1         1         0         0
##     spacer_4         0         0         1         0         0
##     spacer_5         0         0         1         0         0
##          ...       ...       ...       ...       ...       ...
##   spacer_156         0         0         1         0         0
##   spacer_157         0        27         1         0        18
##   spacer_158         0         0         1         0         0
##   spacer_159         0         0         1         0         0
##   spacer_160        18       160         1        17       121
##                                                         alignments
##                                                      <GRangesList>
##     spacer_1                     chr12:25250434:-,chr7:155971918:-
##     spacer_2                     chr12:25250435:-,chr3:184602035:+
##     spacer_3                     chr12:25250436:-,chr3:184602034:+
##     spacer_4                                      chr12:25250443:-
##     spacer_5                                      chr12:25250446:-
##          ...                                                   ...
##   spacer_156                                      chr12:25250915:-
##   spacer_157 chr12:25250916:+,chr2:55050346:+,chr20:21397165:+,...
##   spacer_158                                      chr12:25250921:-
##   spacer_159                                      chr12:25250924:-
##   spacer_160  chr12:25250928:+,chr17:49361951:+,chr1:3069055:-,...
##              score_crispri
##                  <numeric>
##     spacer_1      0.372821
##     spacer_2      0.356982
##     spacer_3      0.390816
##     spacer_4      0.421704
##     spacer_5      0.408481
##          ...           ...
##   spacer_156      0.599555
##   spacer_157      0.666575
##   spacer_158      0.599259
##   spacer_159      0.565775
##   spacer_160      0.636552
##   -------
##   seqinfo: 640 sequences (1 circular) from hg38 genome
##   crisprNuclease: SpCas9
```

You can see that the column `score_crispri` was added to the `GuideSet`.
Note that this function works identically for CRISPRa applications, with
the `modality` argument replaced by `CRISPRa`.

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
    ## [15] BiocGenerics_0.42.0               crisprDesignData_0.99.13         
    ## [17] crisprDesign_0.99.115             crisprBase_1.1.3                 
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

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-crispri" class="csl-entry">

Gilbert, Luke A, Matthew H Larson, Leonardo Morsut, Zairan Liu, Gloria A
Brar, Sandra E Torres, Noam Stern-Ginossar, et al., et al. 2013.
“CRISPR-Mediated Modular RNA-Guided Regulation of Transcription in
Eukaryotes.” *Cell* 154 (2): 442–51.

</div>

<div id="ref-crisprai" class="csl-entry">

Horlbeck, Max A, Luke A Gilbert, Jacqueline E Villalta, Britt Adamson,
Ryan A Pak, Yuwen Chen, Alexander P Fields, et al., et al. 2016.
“Compact and Highly Active Next-Generation Libraries for CRISPR-Mediated
Gene Repression and Activation.” *Elife* 5.

</div>

<div id="ref-sanson2018optimized" class="csl-entry">

Sanson, Kendall R, Ruth E Hanna, Mudra Hegde, Katherine F Donovan,
Christine Strand, Meagan E Sullender, Emma W Vaimberg, et al., et al.
2018. “Optimized Libraries for CRISPR-Cas9 Genetic Screens with Multiple
Modalities.” *Nature Communications* 9 (1): 1–15.

</div>

</div>
