Using crisprDesign to design gRNAs for CRISPRa
================

-   <a href="#introduction" id="toc-introduction">Introduction</a>
-   <a href="#installation" id="toc-installation">Installation</a>
-   <a href="#terminology" id="toc-terminology">Terminology</a>
-   <a href="#crispra-design" id="toc-crispra-design">CRISPRa design</a>
    -   <a href="#creating-the-guideset" id="toc-creating-the-guideset">Creating
        the GuideSet</a>
    -   <a href="#annotating-the-guideset"
        id="toc-annotating-the-guideset">Annotating the GuideSet</a>
    -   <a href="#adding-tss-annotation" id="toc-adding-tss-annotation">Adding
        TSS annotation</a>
    -   <a href="#adding-spacer-alignments-with-tss-annotation"
        id="toc-adding-spacer-alignments-with-tss-annotation">Adding spacer
        alignments with TSS annotation</a>
    -   <a href="#adding-crisprai-scores" id="toc-adding-crisprai-scores">Adding
        CRISPRai scores</a>
-   <a href="#session-info" id="toc-session-info">Session Info</a>

Authors: Jean-Philippe Fortin, Luke Hoberecht

Date: 03 August, 2022

# Introduction

`crisprDesign` is a comprehensive software package for designing and
annotating CRISPR guide RNA (gRNA) sequences, including the
characterization of on-targets and off-targets, gene context annotation,
and SNP annotation (human only). The software was developed to be as
applicable and generalizable as possible.

This tutorial will demonstrate how to use `crisprDesign` to design gRNAs
for CRISPR activation (CRISPRa). Specifically, it will target the human
KRAS gene and use the SpCas9 nuclease, however, much of the same process
can be applied to any genomic target(s) and with any CRISPR nuclease.

# Installation

`crisprDesign` can be installed from Bioconductor using the following
commands in an R session:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("crisprDesign")
```

# Terminology

CRISPR nucleases are examples of RNA-guided endonucleases. They require
two binding components for cleavage. First, the nuclease needs to
recognize a constant nucleotide motif in the target DNA called the
protospacer adjacent motif (PAM) sequence. Second, the gRNA, which
guides the nuclease to the target sequence, needs to bind to a
complementary sequence adjacent to the PAM sequence, called the
**protospacer** sequence. The latter can be thought of as a variable
binding motif that can be specified by designing corresponding gRNA
sequences.

The **spacer** sequence is used in the gRNA construct to guide the
CRISPR nuclease to the target **protospacer** sequence in the host
genome.

For DNA-targeting nucleases, the nucleotide sequence of the spacer and
protospacer are identical. For RNA-targeting nucleases, they are the
reverse complement of each other.

While a gRNA spacer sequence may not always uniquely target the host
genome (i.e. it may map to multiple protospacers in the host genome), we
can, for a given reference genome, uniquely identify a protospacer
sequence with a combination of 3 attributes:

-   `chr`: chromosome name
-   `strand`: forward (+) or reverse (-)
-   `pam_site`: genomic coordinate of the first nucleotide of the
    nuclease-specific PAM sequence; for SpCas9 this is the “N” in the
    NGG PAM sequence

# CRISPRa design

For CRISPR activation (CRISPRa) and interference (CRISPRi) applications
the CRISPR nuclease is engineered to lose its endonuclease activity, and
should therefore not introduce double-stranded breaks (DSBs). We will
use the dead SpCas9 (dSpCas9) nuclease as an example here. Note that
users don’t have to distinguish between dSpCas9 and SpCas9 when
specifying the nuclease in `crisprDesign` and `crisprBase` as they do
not differ in terms of the characteristics stored in the
`CrisprNuclease` object.

In CRISPRa, dSpCas9 is used to activate gene expression by coupling the
dead nuclease with activation factors. Several CRISPRa systems have been
developed (see Kampmann (2018) for a review). For optimal activation,
gRNAs are usually designed to target the region directly upstream of the
gene transcription start site (TSS).

`crisprDesign` provides functionalities to be able to take into account
design rules that are specific to CRISPRa applications. The `queryTss`
function allows for specifying genomic coordinates of promoter regions.
The `addTssAnnotation` function annotates gRNAs for known TSSs, and
includes a column `dist_to_tss` that gives the distance in nucleotides
between the TSS position and the PAM site of the gRNA. For CRISPRa, we
recommend targeting the region 75-150bp upstream of the TSS for optimal
activation; see Sanson et al. (2018) for more information.

## Creating the GuideSet

As an example to demonstrate selecting gRNAs for CRISPRa, suppose we
want to activate the human KRAS gene using the SpCas9 nuclease. To
accomplish this we will want our gRNAs to target the region upstream of
the KRAS TSS; let’s consider the window containing 500bp immediately
upstream of the TSS. We first need to retrieve the TSS coordinates for
KRAS. These data are conveniently stored in the `crisprDesignData`
package as `tss_human` (for more information on `tss_human` and how to
create similar TSS annotation objects, see the [Building a gene
annotation
object](https://github.com/crisprVerse/Tutorials/tree/master/Building_Gene_Annotation)
tutorial).

Install `crisprDesignData` with

``` r
install.packages("devtools")
devtools::install_github("Jfortin1/crisprDesignData")
```

Then load the TSS coordinates stored in `tss_human` and query for KRAS
using the `queryTss` function from `crisprDesign`:

``` r
library(crisprDesignData)
data("tss_human", package="crisprDesignData")
library(crisprDesign)
target_window <- c(-500, 0)
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

We can find all protospacers in our target region with `findSpacers`:

``` r
library(BSgenome.Hsapiens.UCSC.hg38)
gs <- findSpacers(target_region,
                  crisprNuclease=SpCas9,
                  bsgenome=BSgenome.Hsapiens.UCSC.hg38)
```

``` r
gs
## GuideSet object with 146 ranges and 5 metadata columns:
##              seqnames    ranges strand |          protospacer            pam
##                 <Rle> <IRanges>  <Rle> |       <DNAStringSet> <DNAStringSet>
##     spacer_1    chr12  25250927      - | GCTCGGAGCTCGATTTTCCT            AGG
##     spacer_2    chr12  25250944      - | CCCGAACTCATCGGTGTGCT            CGG
##     spacer_3    chr12  25250953      - | CCGCCCGGCCCCGAACTCAT            CGG
##     spacer_4    chr12  25250961      + | TCCGAGCACACCGATGAGTT            CGG
##     spacer_5    chr12  25250962      + | CCGAGCACACCGATGAGTTC            GGG
##          ...      ...       ...    ... .                  ...            ...
##   spacer_142    chr12  25251419      - | AGGCCGACCCTGAGGGTGGC            GGG
##   spacer_143    chr12  25251420      - | TAGGCCGACCCTGAGGGTGG            CGG
##   spacer_144    chr12  25251423      - | GTATAGGCCGACCCTGAGGG            TGG
##   spacer_145    chr12  25251429      + | AAGAGCACCCCGCCACCCTC            AGG
##   spacer_146    chr12  25251430      + | AGAGCACCCCGCCACCCTCA            GGG
##               pam_site  cut_site      region
##              <numeric> <numeric> <character>
##     spacer_1  25250927  25250930    region_1
##     spacer_2  25250944  25250947    region_1
##     spacer_3  25250953  25250956    region_1
##     spacer_4  25250961  25250958    region_1
##     spacer_5  25250962  25250959    region_1
##          ...       ...       ...         ...
##   spacer_142  25251419  25251422    region_1
##   spacer_143  25251420  25251423    region_1
##   spacer_144  25251423  25251426    region_1
##   spacer_145  25251429  25251426    region_1
##   spacer_146  25251430  25251427    region_1
##   -------
##   seqinfo: 640 sequences (1 circular) from hg38 genome
##   crisprNuclease: SpCas9
```

## Annotating the GuideSet

The next step is to annotate our candidate gRNAs to assess their
quality. There are several functions in `crisprDesign` that provide
annotation for features that are nonspecific to CRISPRa, for which we
refer the reader to the [CRISPRko design with
Cas9](https://github.com/crisprVerse/Tutorials/tree/master/Design_CRISPRko_Cas9)
tutorial for more information. The sections below will cover annotation
functions that are of particular interest to CRISPRa applications.

## Adding TSS annotation

While KRAS has a single TSS, for gene target(s) having multiple TSSs it
is valuable knowing which promoter region each gRNA targets, and
consequently which isoforms will be affected. It is also helpful to know
the exact location our gRNAs target with respect to each TSS. All this
information is appended to the `GuideSet` object with the
`addTssAnnotation` function. In addition to the `GuideSet`, we simply
need to pass the `tssObject` object and `tss_window` values we used
earlier to define the target region. We can then retrieve the appended
annotation with the accessor function `tssAnnotation`:

``` r
gs <- addTssAnnotation(gs,
                       tssObject=tss_human,
                       tss_window=target_window)
tssAnnotation(gs)
## DataFrame with 146 rows and 15 columns
##                 chr anchor_site   strand     score peak_start  peak_end
##            <factor>   <integer> <factor> <numeric>  <integer> <integer>
## spacer_1      chr12    25250930        -   5.20187   25250928  25250928
## spacer_2      chr12    25250947        -   5.20187   25250928  25250928
## spacer_3      chr12    25250956        -   5.20187   25250928  25250928
## spacer_4      chr12    25250958        +   5.20187   25250928  25250928
## spacer_5      chr12    25250959        +   5.20187   25250928  25250928
## ...             ...         ...      ...       ...        ...       ...
## spacer_142    chr12    25251422        -   5.20187   25250928  25250928
## spacer_143    chr12    25251423        -   5.20187   25250928  25250928
## spacer_144    chr12    25251426        -   5.20187   25250928  25250928
## spacer_145    chr12    25251426        +   5.20187   25250928  25250928
## spacer_146    chr12    25251427        +   5.20187   25250928  25250928
##                      tx_id         gene_id      source    promoter
##                <character>     <character> <character> <character>
## spacer_1   ENST00000256078 ENSG00000133703     fantom5          P1
## spacer_2   ENST00000256078 ENSG00000133703     fantom5          P1
## spacer_3   ENST00000256078 ENSG00000133703     fantom5          P1
## spacer_4   ENST00000256078 ENSG00000133703     fantom5          P1
## spacer_5   ENST00000256078 ENSG00000133703     fantom5          P1
## ...                    ...             ...         ...         ...
## spacer_142 ENST00000256078 ENSG00000133703     fantom5          P1
## spacer_143 ENST00000256078 ENSG00000133703     fantom5          P1
## spacer_144 ENST00000256078 ENSG00000133703     fantom5          P1
## spacer_145 ENST00000256078 ENSG00000133703     fantom5          P1
## spacer_146 ENST00000256078 ENSG00000133703     fantom5          P1
##                        tss_id gene_symbol  tss_strand   tss_pos dist_to_tss
##                   <character> <character> <character> <integer>   <numeric>
## spacer_1   ENSG00000133703_P1        KRAS           -  25250928          -2
## spacer_2   ENSG00000133703_P1        KRAS           -  25250928         -19
## spacer_3   ENSG00000133703_P1        KRAS           -  25250928         -28
## spacer_4   ENSG00000133703_P1        KRAS           -  25250928         -30
## spacer_5   ENSG00000133703_P1        KRAS           -  25250928         -31
## ...                       ...         ...         ...       ...         ...
## spacer_142 ENSG00000133703_P1        KRAS           -  25250928        -494
## spacer_143 ENSG00000133703_P1        KRAS           -  25250928        -495
## spacer_144 ENSG00000133703_P1        KRAS           -  25250928        -498
## spacer_145 ENSG00000133703_P1        KRAS           -  25250928        -498
## spacer_146 ENSG00000133703_P1        KRAS           -  25250928        -499
```

## Adding spacer alignments with TSS annotation

As with all CRISPR applications, off-targets is an important concern in
assessing gRNA quality. While this concern is somewhat moderated for
CRISPRa, since the dead CRISPR nuclease does not make DSBs, we should be
aware of off-targets occuring in the promoter regions of other genes.
This can be handled by passing our `tssObject` to the
`addSpacerAlignments`; we will search for up to 2 mismatches and
increase the size of our `tss_window` to err on the safe side. (Note:
this alignment example uses a local bowtie index file; for information
on how to create index files for available aligners, see the [Building
genome indices for short read
aligners](https://github.com/crisprVerse/Tutorials/tree/master/Building_Genome_Indices)
tutorial.)

``` r
index_path <- "/Users/hoberecl/crisprIndices/bowtie/hg38/hg38"
gs <- addSpacerAlignments(gs,
                          aligner="bowtie",
                          aligner_index=index_path,
                          bsgenome=BSgenome.Hsapiens.UCSC.hg38,
                          n_mismatches=2,
                          tssObject=tss_human,
                          tss_window=c(-2000, 500))
```

``` r
gs
## GuideSet object with 146 ranges and 13 metadata columns:
##              seqnames    ranges strand |          protospacer            pam
##                 <Rle> <IRanges>  <Rle> |       <DNAStringSet> <DNAStringSet>
##     spacer_1    chr12  25250927      - | GCTCGGAGCTCGATTTTCCT            AGG
##     spacer_2    chr12  25250944      - | CCCGAACTCATCGGTGTGCT            CGG
##     spacer_3    chr12  25250953      - | CCGCCCGGCCCCGAACTCAT            CGG
##     spacer_4    chr12  25250961      + | TCCGAGCACACCGATGAGTT            CGG
##     spacer_5    chr12  25250962      + | CCGAGCACACCGATGAGTTC            GGG
##          ...      ...       ...    ... .                  ...            ...
##   spacer_142    chr12  25251419      - | AGGCCGACCCTGAGGGTGGC            GGG
##   spacer_143    chr12  25251420      - | TAGGCCGACCCTGAGGGTGG            CGG
##   spacer_144    chr12  25251423      - | GTATAGGCCGACCCTGAGGG            TGG
##   spacer_145    chr12  25251429      + | AAGAGCACCCCGCCACCCTC            AGG
##   spacer_146    chr12  25251430      + | AGAGCACCCCGCCACCCTCA            GGG
##               pam_site  cut_site      region        tssAnnotation        n0
##              <numeric> <numeric> <character> <SplitDataFrameList> <numeric>
##     spacer_1  25250927  25250930    region_1 chr12:25250930:-:...         1
##     spacer_2  25250944  25250947    region_1 chr12:25250947:-:...         1
##     spacer_3  25250953  25250956    region_1 chr12:25250956:-:...         1
##     spacer_4  25250961  25250958    region_1 chr12:25250958:+:...         1
##     spacer_5  25250962  25250959    region_1 chr12:25250959:+:...         1
##          ...       ...       ...         ...                  ...       ...
##   spacer_142  25251419  25251422    region_1 chr12:25251422:-:...         1
##   spacer_143  25251420  25251423    region_1 chr12:25251423:-:...         1
##   spacer_144  25251423  25251426    region_1 chr12:25251426:-:...         1
##   spacer_145  25251429  25251426    region_1 chr12:25251426:+:...         1
##   spacer_146  25251430  25251427    region_1 chr12:25251427:+:...         1
##                     n1        n2      n0_p      n1_p      n2_p
##              <numeric> <numeric> <numeric> <numeric> <numeric>
##     spacer_1         0         0         1         0         0
##     spacer_2         0         0         1         0         0
##     spacer_3         0         0         1         0         0
##     spacer_4         0         0         1         0         0
##     spacer_5         0         0         1         0         0
##          ...       ...       ...       ...       ...       ...
##   spacer_142         0         1         1         0         0
##   spacer_143         0         0         1         0         0
##   spacer_144         0         0         1         0         0
##   spacer_145         0         0         1         0         0
##   spacer_146         0         0         1         0         0
##                                                          alignments
##                                                       <GRangesList>
##     spacer_1 chr12:25250927:-,chr18:21741459:-,chr7:131738185:-,...
##     spacer_2 chr9:130886259:+,chr12:25251105:+,chr9:111897205:+,...
##     spacer_3  chr3:9072318:-,chr12:117052836:-,chr12:31324170:-,...
##     spacer_4   chr5:136560488:+,chr1:1627633:-,chr9:131739852:+,...
##     spacer_5  chr3:101742769:+,chr16:22191335:-,chr1:19485593:-,...
##          ...                                                    ...
##   spacer_142 chr19:7990814:+,chr12:125795462:-,chr15:58933192:-,...
##   spacer_143   chr14:99531164:+,chr9:81862478:-,chr1:22563612:-,...
##   spacer_144   chr19:18597158:-,chr4:5826778:-,chr13:80341047:+,...
##   spacer_145 chr14:105483246:-,chr12:118087564:+,chr10:135759:+,...
##   spacer_146    chr19:8003020:+,chr16:88227567:-,chr11:556165:+,...
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
(Horlbeck et al. 2016). The function `addCrispraiScores` implements this
algorithm to add scores to the `GuideSet`; it requires several inputs in
addition to the `GuideSet` object:

-   `gr` is the `GRanges` object derived from the `queryTss` function
    and used to create the `GuideSet` object. In our example, this is
    `target_region`.
-   `tssObject` is a `GRanges` object that contains TSS coordinates and
    annotation. It must also contain the following columns: `ID`,
    `promoter`, `tx_id`, and `gene_symbol`. Our `tssObject` in this
    instance is `tss_human`.
-   `geneCol` indicates which column of `tssObject` should be used as
    the unique gene identifier
-   `modality` is the modality of the CRISPR application; in this
    example it will be `"CRISPRa"`
-   `fastaFile` is the path of the FASTA file of the human reference
    genome; the file for the human genome (hg38) can be downloaded
    directly from here:
    <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz>
-   `chromatinFiles` is a vector of length 3 specifying the path of
    files containing the chromatin accessibility data needed for the
    algorithm in hg38 coordinates. The chromatin files can be downloaded
    from Zenodo [here](https://zenodo.org/record/6716721#.YrzCfS-cY4d).

Let’s prepare the inputs for `addCrispraiScores` then annotated our
gRNAs with those scores (the required files have been downloaded locally
in advance):

``` r
fastaPath <- "/Users/hoberecl/crispraiScores/hg38.fa.gz"
mnasePath <- "/Users/hoberecl/crispraiScores/crispria_mnase_human_K562_hg38.bigWig"
dnasePath <- "/Users/hoberecl/crispraiScores/crispria_dnase_human_K562_hg38.bigWig"
fairePath <- "/Users/hoberecl/crispraiScores/crispria_faire_human_K562_hg38.bigWig"
results <- addCrispraiScores(gs,
                             gr=target_region,
                             tssObject=tss_human,
                             geneCol="gene_id",
                             modality="CRISPRa",
                             fastaFile=fastaPath,
                             chromatinFiles=c(mnase=mnasePath,
                                              dnase=dnasePath,
                                              faire=fairePath)
                             )
```

``` r
results
```

This function works identically for CRISPRi applications, with modality
replaced by `CRISPRi`.

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
    ## [11] crisprDesign_0.99.110             crisprBase_1.1.2                 
    ## [13] crisprDesignData_0.99.12         
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
    ##  [55] basilisk.utils_1.9.1          zlibbioc_1.42.0              
    ##  [57] vroom_1.5.7                   VariantAnnotation_1.42.1     
    ##  [59] hms_1.1.1                     promises_1.2.0.1             
    ##  [61] MatrixGenerics_1.8.1          parallel_4.2.0               
    ##  [63] SummarizedExperiment_1.26.1   yaml_2.3.5                   
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
    ## [115] httpuv_1.6.5                  Rbwa_1.0.0                   
    ## [117] sessioninfo_1.2.2

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-crisprai" class="csl-entry">

Horlbeck, Max A, Luke A Gilbert, Jacqueline E Villalta, Britt Adamson,
Ryan A Pak, Yuwen Chen, Alexander P Fields, et al. 2016. “Compact and
Highly Active Next-Generation Libraries for CRISPR-Mediated Gene
Repression and Activation.” *Elife* 5.

</div>

<div id="ref-crispracrisprireview" class="csl-entry">

Kampmann, Martin. 2018. “CRISPRi and CRISPRa Screens in Mammalian Cells
for Precision Biology and Medicine.” *ACS Chemical Biology* 13 (2):
406–16.

</div>

<div id="ref-sanson2018optimized" class="csl-entry">

Sanson, Kendall R, Ruth E Hanna, Mudra Hegde, Katherine F Donovan,
Christine Strand, Meagan E Sullender, Emma W Vaimberg, et al. 2018.
“Optimized Libraries for CRISPR-Cas9 Genetic Screens with Multiple
Modalities.” *Nature Communications* 9 (1): 1–15.

</div>

</div>
