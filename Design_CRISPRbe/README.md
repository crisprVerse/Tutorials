Using crisprDesign to design gRNAs for CRISPRbe
================

-   <a href="#introduction" id="toc-introduction">Introduction</a>
-   <a href="#installation" id="toc-installation">Installation</a>
-   <a href="#terminology" id="toc-terminology">Terminology</a>
-   <a href="#crispr-base-editing-with-be4max"
    id="toc-crispr-base-editing-with-be4max">CRISPR base editing with
    BE4max</a>
    -   <a href="#loading-packages" id="toc-loading-packages">Loading
        packages</a>
    -   <a href="#example-workflow" id="toc-example-workflow">Example
        workflow</a>
-   <a href="#session-info" id="toc-session-info">Session Info</a>
-   <a href="#references" id="toc-references">References</a>

Authors: Jean-Philippe Fortin, Luke Hoberecht

Date: 04 August, 2022

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
    nuclease-specific PAM sequence (e.g. for SpCas9, the “N” in the NGG
    PAM sequence; for AsCas12a, the first “T” of the TTTV PAM sequence)

# CRISPR base editing with BE4max

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

## Example workflow

We illustrate the CRISPR base editing (CRISPRbe) functionalities of
`crisprDesign` by designing and characterizing gRNAs targeting the human
gene KRAS using the cytidine base editor BE4max (Koblan et al. 2018).

We first load the BE4max `BaseEditor` object from the `crisprBase`
package:

``` r
data(BE4max, package="crisprBase")
BE4max
```

    ## Class: BaseEditor
    ##   CRISPR Nuclease name: SpCas9
    ##       Target type: DNA
    ##       Metadata: list of length 2
    ##       PAMs: NGG, NAG, NGA
    ##       Weights: 1, 0.2593, 0.0694
    ##       Spacer length: 20
    ##       PAM side: 3prime
    ##         Distance from PAM: 0
    ##       Prototype protospacers: 5'--SSSSSSSSSSSSSSSSSSSS[NGG]--3', 5'--SSSSSSSSSSSSSSSSSSSS[NAG]--3', 5'--SSSSSSSSSSSSSSSSSSSS[NGA]--3'
    ##   Base editor name: BE4max
    ##       Editing strand: original
    ##       Maximum editing weight: C2T at position -15

The editing probabilities of the base editor BE4max are stored in a
matrix where rows correspond to the different nucleotide substitutions,
and columns correspond to the genomic coordinate relative to the PAM
site. The `editingWeights` function from `crisprBase` retrieves those
probabilities. One can see that C to T editing is optimal around 15
nucleotides upstream of the PAM site for the BE4max base editor:

``` r
crisprBase::editingWeights(BE4max)["C2T",]
```

    ##   -36   -35   -34   -33   -32   -31   -30   -29   -28   -27   -26   -25   -24 
    ## 0.007 0.007 0.008 0.018 0.010 0.020 0.014 0.012 0.023 0.013 0.024 0.022 0.034 
    ##   -23   -22   -21   -20   -19   -18   -17   -16   -15   -14   -13   -12   -11 
    ## 0.022 0.021 0.035 0.058 0.162 0.318 0.632 0.903 1.000 0.870 0.620 0.314 0.163 
    ##   -10    -9    -8    -7    -6    -5    -4    -3    -2    -1 
    ## 0.100 0.056 0.033 0.019 0.018 0.024 0.017 0.005 0.002 0.001

Let’s create a `GuideSet` object targeting KRAS and, for the sake of
brevity, retain only a couple gRNAs. We first need to retrieve the
coordinates for KRAS. These data are conveniently stored in the
`crisprDesignData` package as `txdb_human` (for more information on
`txdb_human` and how to create similar gene annotation objects, see the
[Building a gene annotation
object](https://github.com/crisprVerse/Tutorials/tree/master/Building_Gene_Annotation)
tutorial).

``` r
data("txdb_human", package="crisprDesignData")
```

``` r
bsgenome <- BSgenome.Hsapiens.UCSC.hg38

gr <- queryTxObject(txObject=txdb_human,
                    featureType="cds",
                    queryColumn="gene_symbol",
                    queryValue="KRAS")
gs <- findSpacers(gr,
                  bsgenome=bsgenome,
                  crisprNuclease=BE4max)
gs <- gs[57:58]
```

The function `addEditedAlleles` finds, characterizes, and scores
predicted edited alleles for each gRNA and a chosen transcript. It
requires a transcript-specific annotation that can be obtained with the
`getTxInfoDataFrame` function. Here, we perform the analysis using the
primary isoform of KRAS (Ensembl transcript ID: ENST00000311936).

We first get the transcript table for our transcript,

``` r
txid <- "ENST00000311936"
txTable <- getTxInfoDataFrame(tx_id=txid,
                              txObject=txdb_human,
                              bsgenome=bsgenome)
head(txTable)
```

    ## DataFrame with 6 rows and 10 columns
    ##           chr       pos         nuc          aa aa_number      exon  pos_plot
    ##   <character> <numeric> <character> <character> <integer> <integer> <integer>
    ## 1       chr12  25250929           C          NA        NA         1        31
    ## 2       chr12  25250928           T          NA        NA         1        32
    ## 3       chr12  25250927           A          NA        NA         1        33
    ## 4       chr12  25250926           G          NA        NA         1        34
    ## 5       chr12  25250925           G          NA        NA         1        35
    ## 6       chr12  25250924           C          NA        NA         1        36
    ##    pos_mrna   pos_cds      region
    ##   <integer> <integer> <character>
    ## 1         1        NA        5UTR
    ## 2         2        NA        5UTR
    ## 3         3        NA        5UTR
    ## 4         4        NA        5UTR
    ## 5         5        NA        5UTR
    ## 6         6        NA        5UTR

and then add the edited alleles annotation to the `GuideSet`:

``` r
editingWindow <- c(-20,-8)
gs <- addEditedAlleles(gs,
                       baseEditor=BE4max,
                       txTable=txTable,
                       editingWindow=editingWindow)
```

    ## [addEditedAlleles] Obtaining edited alleles at each gRNA target site.

    ## [addEditedAlleles] Adding functional consequences to alleles.

The `editingWindow` argument specifies the window of editing that we are
interested in. When not provided, it uses the default window provided in
the `BaseEditor` object. Note that providing large windows can
exponentially increase computing time as the number of possible alleles
grows exponentially.

Let’s retrieve the edited alleles for the first gRNA:

``` r
alleles <- editedAlleles(gs)[[1]]
```

We get a `DataFrame` object with useful metadata:

``` r
metadata(alleles)
```

    ## $wildtypeAllele
    ##       spacer_57 
    ## "AGGGACCAGTACA" 
    ## 
    ## $start
    ## [1] 25227310
    ## 
    ## $end
    ## [1] 25227322
    ## 
    ## $chr
    ## [1] "chr12"
    ## 
    ## $strand
    ## [1] "-"
    ## 
    ## $editingWindow
    ## [1] -20  -8
    ## 
    ## $wildtypeAmino
    ## [1] "MYYYQQQDDDRRR"

The `wildtypeAllele` reports the unedited nucleotide sequence of the
region specified by the editing window (with respect to the gRNA PAM
site). It is always reported from the 5’ to 3’ direction on the strand
corresponding to the gRNA strand. The `start` and `end` fields specify
the corresponding coordinates on the transcript.

Let’s look at the edited alleles:

``` r
head(alleles)
```

    ## DataFrame with 6 rows and 4 columns
    ##              seq     score     variant            aa
    ##   <DNAStringSet> <numeric> <character>   <character>
    ## 1  AGGGATTAGTACA 0.6987640    nonsense MYYY***DDDRRR
    ## 2  AGGGATCAGTACA 0.1052162      silent MYYYQQQDDDRRR
    ## 3  AGGGATGAGTACA 0.0449779    missense MYYYEEEDDDRRR
    ## 4  AGGGATTAGTATA 0.0382137    nonsense MYYY***DDDRRR
    ## 5  AGGGAGTAGTACA 0.0377333    nonsense MYYY***EEERRR
    ## 6  AGGGATAAGTACA 0.0216858    missense MYYYKKKDDDRRR

The `DataFrame` is ordered by descending values in the `score` column.
This `score` represents the likelihood of the edited allele to occur
relative to all possible edited alleles, and is calculated using the
editing weights stored in the `BE4max` object. The `seq` column
represents the edited nucleotide sequences. As with the `wildtypeAllele`
in the metadata, they are always reported from the 5’ to 3’ direction on
the strand corresponding to the gRNA strand. The `variant` column
describes the functional consequence of the editing event (silent,
nonsense or missense mutation). If an edited allele results in multiple
editing events, as can happen when multiple bases are edited, the most
pronounced mutation (nonsense over missense, missense over silent) is
reported. Finally, the `aa` column reports the resulting edited amino
acid sequence, with each single letter code mapping to its corresponding
nucleotide (`*` for termination).

Note that `addEditedAlleles` also appended several gRNA-level aggregate
scores to the `GuideSet` object:

``` r
head(gs)
```

    ## GuideSet object with 2 ranges and 11 metadata columns:
    ##             seqnames    ranges strand |          protospacer            pam
    ##                <Rle> <IRanges>  <Rle> |       <DNAStringSet> <DNAStringSet>
    ##   spacer_57    chr12  25227302      - | AGGGACCAGTACATGAGGAC            TGG
    ##   spacer_58    chr12  25227307      - | CAATGAGGGACCAGTACATG            AGG
    ##              pam_site  cut_site      region
    ##             <numeric> <numeric> <character>
    ##   spacer_57  25227302  25227305    region_6
    ##   spacer_58  25227307  25227310    region_2
    ##                                                                                                                editedAlleles
    ##                                                                                                                       <list>
    ##   spacer_57 AGGGATTAGTACA:0.6987640:nonsense:...,AGGGATCAGTACA:0.1052162:silent:...,AGGGATGAGTACA:0.0449779:missense:...,...
    ##   spacer_58 CAATGAGGGATCA:0.0820221:silent:...,TAATGAGGGACCA:0.0452522:missense:...,CAATGAGGGACTA:0.0437344:nonsense:...,...
    ##             score_missense score_nonsense score_silent  maxVariant
    ##                  <numeric>      <numeric>    <numeric> <character>
    ##   spacer_57      0.0864560      0.8025738    0.1109702    nonsense
    ##   spacer_58      0.0665998      0.0516591    0.0820221      silent
    ##             maxVariantScore
    ##                   <numeric>
    ##   spacer_57       0.8025738
    ##   spacer_58       0.0820221
    ##   -------
    ##   seqinfo: 640 sequences (1 circular) from hg38 genome
    ##   crisprNuclease: SpCas9

The `score_missense`, `score_nonsense` and `score_silent` columns report
aggregated scores for each mutation type. They are calculated by summing
all scores of a given mutation type across the set of edited alleles for
a given gRNA. The `maxVariant` column indicates the most probable
mutation type for the given gRNA based on the maximum aggregated score,
which is stored in `maxVariantScore`. In our example, the highest score
for `spacer_57` is `score_nonsense`, and so `maxVariant` is set to
`nonsense`.

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
    ##  [43] crisprScore_1.1.14            xfun_0.31                    
    ##  [45] stringr_1.4.0                 ps_1.7.1                     
    ##  [47] mime_0.12                     miniUI_0.1.1.1               
    ##  [49] lifecycle_1.0.1               restfulr_0.0.15              
    ##  [51] devtools_2.4.4                XML_3.99-0.10                
    ##  [53] AnnotationHub_3.4.0           zlibbioc_1.42.0              
    ##  [55] basilisk.utils_1.9.1          VariantAnnotation_1.42.1     
    ##  [57] hms_1.1.1                     promises_1.2.0.1             
    ##  [59] MatrixGenerics_1.8.1          parallel_4.2.0               
    ##  [61] SummarizedExperiment_1.26.1   yaml_2.3.5                   
    ##  [63] curl_4.3.2                    memoise_2.0.1                
    ##  [65] reticulate_1.25               biomaRt_2.52.0               
    ##  [67] stringi_1.7.8                 RSQLite_2.2.15               
    ##  [69] BiocVersion_3.15.2            BiocIO_1.6.0                 
    ##  [71] randomForest_4.7-1.1          crisprScoreData_1.1.3        
    ##  [73] GenomicFeatures_1.48.3        filelock_1.0.2               
    ##  [75] pkgbuild_1.3.1                BiocParallel_1.30.3          
    ##  [77] rlang_1.0.4                   pkgconfig_2.0.3              
    ##  [79] matrixStats_0.62.0            bitops_1.0-7                 
    ##  [81] evaluate_0.15                 lattice_0.20-45              
    ##  [83] purrr_0.3.4                   htmlwidgets_1.5.4            
    ##  [85] GenomicAlignments_1.32.1      bit_4.0.4                    
    ##  [87] tidyselect_1.1.2              processx_3.7.0               
    ##  [89] magrittr_2.0.3                R6_2.5.1                     
    ##  [91] profvis_0.3.7                 generics_0.1.3               
    ##  [93] DelayedArray_0.22.0           DBI_1.1.3                    
    ##  [95] pillar_1.8.0                  KEGGREST_1.36.3              
    ##  [97] RCurl_1.98-1.8                tibble_3.1.8                 
    ##  [99] dir.expiry_1.4.0              crayon_1.5.1                 
    ## [101] utf8_1.2.2                    BiocFileCache_2.4.0          
    ## [103] urlchecker_1.0.1              tzdb_0.3.0                   
    ## [105] rmarkdown_2.14                usethis_2.1.6                
    ## [107] progress_1.2.2                grid_4.2.0                   
    ## [109] blob_1.2.3                    callr_3.7.1                  
    ## [111] digest_0.6.29                 xtable_1.8-4                 
    ## [113] httpuv_1.6.5                  sessioninfo_1.2.2

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-crispracrisprireview" class="csl-entry">

Kampmann, Martin. 2018. “CRISPRi and CRISPRa Screens in Mammalian Cells
for Precision Biology and Medicine.” *ACS Chemical Biology* 13 (2):
406–16.

</div>

<div id="ref-koblan2018improving" class="csl-entry">

Koblan, Luke W, Jordan L Doman, Christopher Wilson, Jonathan M Levy,
Tristan Tay, Gregory A Newby, Juan Pablo Maianti, Aditya Raguram, and
David R Liu. 2018. “Improving Cytidine and Adenine Base Editors by
Expression Optimization and Ancestral Reconstruction.” *Nature
Biotechnology* 36 (9): 843–46.

</div>

</div>
