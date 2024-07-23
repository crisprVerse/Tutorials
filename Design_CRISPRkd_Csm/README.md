Using crisprDesign to design gRNAs for the CRISPR-Csm complex
================
Jean-Philippe Fortin, Luke Hoberecht

-   <a href="#introduction" id="toc-introduction">Introduction</a>
-   <a href="#installation" id="toc-installation">Installation</a>
-   <a href="#terminology" id="toc-terminology">Terminology</a>
-   <a href="#end-to-end-grna-design-workflow"
    id="toc-end-to-end-grna-design-workflow">End-to-end gRNA design
    workflow</a>
    -   <a href="#creating-the-guideset" id="toc-creating-the-guideset">Creating
        the GuideSet</a>
    -   <a href="#annotating-the-guideset"
        id="toc-annotating-the-guideset">Annotating the GuideSet</a>
        -   <a href="#adding-spacer-alignments"
            id="toc-adding-spacer-alignments">Adding spacer alignments</a>
-   <a href="#session-info" id="toc-session-info">Session Info</a>
-   <a href="#references" id="toc-references">References</a>

# Introduction

The CRISPR-Csm complex is a programmable RNA-targeting system that does
not induce indiscriminate trans-cleavage activity, which is an important
advantage in comparison to the CRISPR-Cas13 family of RNA-targeting
nucleases (Colognori, Trinidad, and Doudna 2023). It has recently been
shown that it can be use to perform effective single-molecule live-cell
RNA imaging (Xia et al. 2024).

In this tutorial, we will design gRNAs for the CRISPR-Csm system for the
primary isoform of the human gene KRAS.

# Installation

See the [Installation
tutorial](https://github.com/crisprVerse/Tutorials/tree/master/Installation)
to learn how to install the packages necessary for this tutorial:
`crisprDesign`, `crisprDesignData`

# Terminology

See the [CRISPRko design
vignette](https://github.com/crisprVerse/Tutorials/tree/master/Design_CRISPRko_Cas9)
to get familiar with the terminology used throughout this tutorial.

# End-to-end gRNA design workflow

We first start by loading the crisprVerse packages needed for this
tutorial:

``` r
library(crisprBase)
library(crisprDesign)
library(crisprDesignData)
```

We will also load the `BSgenome` package containing DNA sequences for
the hg38 genome:

``` r
library(BSgenome.Hsapiens.UCSC.hg38)
```

## Creating the GuideSet

We begin by loading the Csm `CrisprNuclease` object from the
`crisprBase` package:

``` r
data(Csm, package="crisprBase")
Csm
```

    ## Class: CrisprNuclease
    ##   Name: Csm
    ##   Target type: RNA
    ##   Metadata: list of length 2
    ##   PFS: N
    ##   Weights: 1
    ##   Spacer length: 32
    ##   PFS side: 3prime
    ##     Distance from PFS: 0
    ##   Prototype protospacers: 5'--SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS[N]--3'

The PFS sequence (the equivalent of a PAM sequence for RNA-targeting
nucleases) for Csm is `N`, meaning there is no specific PFS sequences
preferred by Csm. The default spacer length of the Csm nuclease is 32nt.
This can be changed using `spacerLength` (for instance,
`spacerLength(Csm) <- 36`).

Next, we will extract the mRNA sequence for the KRAS transcript
ENST00000311936 with the function `getMrnaSequences` from
`crisprDesign`. The function requires a gene annotation object. We will
load the Ensembl model from the `crisprDesignData` package stored in the
`GRangesList` object `txdb_human`:

``` r
data("txdb_human", package="crisprDesignData")
```

For more information on `txdb_human` and how to create similar gene
annotation objects, see the [Building a gene annotation
object](https://github.com/crisprVerse/Tutorials/tree/master/Building_Gene_Annotation)
tutorial).

We also need a `BSgenome` object containing the DNA sequences:

``` r
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
```

We are now ready to obtain our mRNA sequence:

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

Similar to the CRISPRko gRNA design, we use the function `findSpacers`
to design our gRNAs:

``` r
gs <- findSpacers(mrna,
                  crisprNuclease=Csm)
head(gs)
```

    ## GuideSet object with 6 ranges and 5 metadata columns:
    ##                   seqnames    ranges strand |             protospacer
    ##                      <Rle> <IRanges>  <Rle> |          <DNAStringSet>
    ##   spacer_1 ENST00000311936        33      + | CTAGGCGGCG...GAGGCAGCAG
    ##   spacer_2 ENST00000311936        34      + | TAGGCGGCGG...AGGCAGCAGC
    ##   spacer_3 ENST00000311936        35      + | AGGCGGCGGC...GGCAGCAGCG
    ##   spacer_4 ENST00000311936        36      + | GGCGGCGGCC...GCAGCAGCGG
    ##   spacer_5 ENST00000311936        37      + | GCGGCGGCCG...CAGCAGCGGC
    ##   spacer_6 ENST00000311936        38      + | CGGCGGCCGC...AGCAGCGGCG
    ##                       pam  pam_site  cut_site          region
    ##            <DNAStringSet> <numeric> <numeric>     <character>
    ##   spacer_1              C        33        NA ENST00000311936
    ##   spacer_2              G        34        NA ENST00000311936
    ##   spacer_3              G        35        NA ENST00000311936
    ##   spacer_4              C        36        NA ENST00000311936
    ##   spacer_5              G        37        NA ENST00000311936
    ##   spacer_6              G        38        NA ENST00000311936
    ##   -------
    ##   seqinfo: 1 sequence from custom genome
    ##   crisprNuclease: Csm

Note that all protospacer sequences are located on the original strand
of the mRNA sequence. For RNA-targeting nucleases, the spacer and
protospacer sequences are the reverse complement of each other. (Compare
the output of the code below with a `GuideSet` that uses a DNA-targeting
nuclease–for such `GuideSet` pbjects, the output of `spacers` and
`protospacers` are identical.)

``` r
head(spacers(gs))
```

    ## DNAStringSet object of length 6:
    ##     width seq                                               names               
    ## [1]    32 CTGCTGCCTCCGCCGCCGCGGCCGCCGCCTAG                  spacer_1
    ## [2]    32 GCTGCTGCCTCCGCCGCCGCGGCCGCCGCCTA                  spacer_2
    ## [3]    32 CGCTGCTGCCTCCGCCGCCGCGGCCGCCGCCT                  spacer_3
    ## [4]    32 CCGCTGCTGCCTCCGCCGCCGCGGCCGCCGCC                  spacer_4
    ## [5]    32 GCCGCTGCTGCCTCCGCCGCCGCGGCCGCCGC                  spacer_5
    ## [6]    32 CGCCGCTGCTGCCTCCGCCGCCGCGGCCGCCG                  spacer_6

``` r
head(protospacers(gs))
```

    ## DNAStringSet object of length 6:
    ##     width seq                                               names               
    ## [1]    32 CTAGGCGGCGGCCGCGGCGGCGGAGGCAGCAG                  spacer_1
    ## [2]    32 TAGGCGGCGGCCGCGGCGGCGGAGGCAGCAGC                  spacer_2
    ## [3]    32 AGGCGGCGGCCGCGGCGGCGGAGGCAGCAGCG                  spacer_3
    ## [4]    32 GGCGGCGGCCGCGGCGGCGGAGGCAGCAGCGG                  spacer_4
    ## [5]    32 GCGGCGGCCGCGGCGGCGGAGGCAGCAGCGGC                  spacer_5
    ## [6]    32 CGGCGGCCGCGGCGGCGGAGGCAGCAGCGGCG                  spacer_6

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
bowtie_index <- "/Users/fortinj2/crisprIndices/bowtie/ensembl_human_104/ensembl_human_104"
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
    ##                   seqnames    ranges strand |             protospacer
    ##                      <Rle> <IRanges>  <Rle> |          <DNAStringSet>
    ##   spacer_1 ENST00000311936        33      + | CTAGGCGGCG...GAGGCAGCAG
    ##   spacer_2 ENST00000311936        34      + | TAGGCGGCGG...AGGCAGCAGC
    ##   spacer_3 ENST00000311936        35      + | AGGCGGCGGC...GGCAGCAGCG
    ##   spacer_4 ENST00000311936        36      + | GGCGGCGGCC...GCAGCAGCGG
    ##   spacer_5 ENST00000311936        37      + | GCGGCGGCCG...CAGCAGCGGC
    ##   spacer_6 ENST00000311936        38      + | CGGCGGCCGC...AGCAGCGGCG
    ##                       pam  pam_site  cut_site          region     n0_tx
    ##            <DNAStringSet> <numeric> <numeric>     <character> <numeric>
    ##   spacer_1              C        33        NA ENST00000311936         4
    ##   spacer_2              G        34        NA ENST00000311936         4
    ##   spacer_3              G        35        NA ENST00000311936         4
    ##   spacer_4              C        36        NA ENST00000311936         4
    ##   spacer_5              G        37        NA ENST00000311936         4
    ##   spacer_6              G        38        NA ENST00000311936         4
    ##                n1_tx   n0_gene   n1_gene
    ##            <numeric> <numeric> <numeric>
    ##   spacer_1         0         1         0
    ##   spacer_2         0         1         0
    ##   spacer_3         0         1         0
    ##   spacer_4         0         1         0
    ##   spacer_5         0         1         0
    ##   spacer_6         0         1         0
    ##                                                                    alignments
    ##                                                                 <GRangesList>
    ##   spacer_1 ENST00000256078:33:+,ENST00000311936:33:+,ENST00000556131:33:+,...
    ##   spacer_2 ENST00000256078:34:+,ENST00000311936:34:+,ENST00000556131:34:+,...
    ##   spacer_3 ENST00000256078:35:+,ENST00000311936:35:+,ENST00000556131:35:+,...
    ##   spacer_4 ENST00000256078:36:+,ENST00000311936:36:+,ENST00000556131:36:+,...
    ##   spacer_5 ENST00000256078:37:+,ENST00000311936:37:+,ENST00000556131:37:+,...
    ##   spacer_6 ENST00000256078:38:+,ENST00000311936:38:+,ENST00000556131:38:+,...
    ##   -------
    ##   seqinfo: 1 sequence from custom genome
    ##   crisprNuclease: Csm

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
    ##   spacer_1 ENST00000256078        33      + | CTGCTGCCTC...CGCCGCCTAG
    ##   spacer_1 ENST00000311936        33      + | CTGCTGCCTC...CGCCGCCTAG
    ##   spacer_1 ENST00000556131        33      + | CTGCTGCCTC...CGCCGCCTAG
    ##   spacer_1 ENST00000557334        40      + | CTGCTGCCTC...CGCCGCCTAG
    ##                        protospacer            pam  pam_site n_mismatches
    ##                     <DNAStringSet> <DNAStringSet> <numeric>    <integer>
    ##   spacer_1 CTAGGCGGCG...GAGGCAGCAG              C        33            0
    ##   spacer_1 CTAGGCGGCG...GAGGCAGCAG              C        33            0
    ##   spacer_1 CTAGGCGGCG...GAGGCAGCAG              C        33            0
    ##   spacer_1 CTAGGCGGCG...GAGGCAGCAG              C        40            0
    ##            canonical  cut_site         gene_id gene_symbol
    ##            <logical> <numeric>     <character> <character>
    ##   spacer_1      TRUE        NA ENSG00000133703        KRAS
    ##   spacer_1      TRUE        NA ENSG00000133703        KRAS
    ##   spacer_1      TRUE        NA ENSG00000133703        KRAS
    ##   spacer_1      TRUE        NA ENSG00000133703        KRAS
    ##   -------
    ##   seqinfo: 2825 sequences from an unspecified genome; no seqlengths

Note that each annotated alignment is specific to the transcript ID
given under `seqnames`.

Below is a spacer that targets (with no mismatches) multiple genes:

``` r
results["spacer_244"]
```

    ## GuideSet object with 1 range and 10 metadata columns:
    ##                     seqnames    ranges strand |             protospacer
    ##                        <Rle> <IRanges>  <Rle> |          <DNAStringSet>
    ##   spacer_244 ENST00000311936       276      + | CTTGACGATA...AATCATTTTG
    ##                         pam  pam_site  cut_site          region     n0_tx
    ##              <DNAStringSet> <numeric> <numeric>     <character> <numeric>
    ##   spacer_244              T       276        NA ENST00000311936         5
    ##                  n1_tx   n0_gene   n1_gene
    ##              <numeric> <numeric> <numeric>
    ##   spacer_244         0         2         0
    ##                                                                        alignments
    ##                                                                     <GRangesList>
    ##   spacer_244 ENST00000256078:276:+,ENST00000311936:276:+,ENST00000407852:86:+,...
    ##   -------
    ##   seqinfo: 1 sequence from custom genome
    ##   crisprNuclease: Csm

Upon further inspection of this spacer’s alignments, however, we can see
that the off-target occurs in the pseudogene KRASP1, and should be
harmless.

``` r
onTargets(results["spacer_244"])
```

    ## GRanges object with 5 ranges and 9 metadata columns:
    ##                     seqnames    ranges strand |                  spacer
    ##                        <Rle> <IRanges>  <Rle> |          <DNAStringSet>
    ##   spacer_244 ENST00000256078       276      + | CAAAATGATT...TATCGTCAAG
    ##   spacer_244 ENST00000311936       276      + | CAAAATGATT...TATCGTCAAG
    ##   spacer_244 ENST00000407852        86      + | CAAAATGATT...TATCGTCAAG
    ##   spacer_244 ENST00000556131       263      + | CAAAATGATT...TATCGTCAAG
    ##   spacer_244 ENST00000557334       283      + | CAAAATGATT...TATCGTCAAG
    ##                          protospacer            pam  pam_site n_mismatches
    ##                       <DNAStringSet> <DNAStringSet> <numeric>    <integer>
    ##   spacer_244 CTTGACGATA...AATCATTTTG              T       276            0
    ##   spacer_244 CTTGACGATA...AATCATTTTG              T       276            0
    ##   spacer_244 CTTGACGATA...AATCATTTTG              T        86            0
    ##   spacer_244 CTTGACGATA...AATCATTTTG              T       263            0
    ##   spacer_244 CTTGACGATA...AATCATTTTG              T       283            0
    ##              canonical  cut_site         gene_id gene_symbol
    ##              <logical> <numeric>     <character> <character>
    ##   spacer_244      TRUE        NA ENSG00000133703        KRAS
    ##   spacer_244      TRUE        NA ENSG00000133703        KRAS
    ##   spacer_244      TRUE        NA ENSG00000220635      KRASP1
    ##   spacer_244      TRUE        NA ENSG00000133703        KRAS
    ##   spacer_244      TRUE        NA ENSG00000133703        KRAS
    ##   -------
    ##   seqinfo: 2825 sequences from an unspecified genome; no seqlengths

# Session Info

``` r
sessionInfo()
```

    ## R version 4.4.1 (2024-06-14)
    ## Platform: x86_64-apple-darwin20
    ## Running under: macOS Ventura 13.6.7
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-x86_64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] BSgenome.Hsapiens.UCSC.hg38_1.4.5 BSgenome_1.73.0                  
    ##  [3] rtracklayer_1.65.0                BiocIO_1.15.0                    
    ##  [5] Biostrings_2.73.1                 XVector_0.45.0                   
    ##  [7] GenomicRanges_1.57.1              GenomeInfoDb_1.41.1              
    ##  [9] IRanges_2.39.2                    S4Vectors_0.43.2                 
    ## [11] BiocGenerics_0.51.0               crisprDesignData_0.99.30         
    ## [13] crisprDesign_1.7.2                crisprBase_1.9.1                 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] DBI_1.2.3                   bitops_1.0-7               
    ##  [3] httr2_1.0.2                 biomaRt_2.61.2             
    ##  [5] rlang_1.1.4                 magrittr_2.0.3             
    ##  [7] Rbowtie_1.45.0              matrixStats_1.3.0          
    ##  [9] compiler_4.4.1              RSQLite_2.3.7              
    ## [11] GenomicFeatures_1.57.0      dir.expiry_1.13.0          
    ## [13] png_0.1-8                   vctrs_0.6.5                
    ## [15] txdbmaker_1.1.1             stringr_1.5.1              
    ## [17] pkgconfig_2.0.3             crayon_1.5.3               
    ## [19] fastmap_1.2.0               dbplyr_2.5.0               
    ## [21] utf8_1.2.4                  Rsamtools_2.21.0           
    ## [23] rmarkdown_2.27              tzdb_0.4.0                 
    ## [25] UCSC.utils_1.1.0            bit_4.0.5                  
    ## [27] xfun_0.46                   randomForest_4.7-1.1       
    ## [29] zlibbioc_1.51.1             cachem_1.1.0               
    ## [31] jsonlite_1.8.8              progress_1.2.3             
    ## [33] blob_1.2.4                  DelayedArray_0.31.9        
    ## [35] BiocParallel_1.39.0         parallel_4.4.1             
    ## [37] prettyunits_1.2.0           R6_2.5.1                   
    ## [39] VariantAnnotation_1.51.0    stringi_1.8.4              
    ## [41] reticulate_1.38.0           Rcpp_1.0.13                
    ## [43] SummarizedExperiment_1.35.1 knitr_1.48                 
    ## [45] readr_2.1.5                 Matrix_1.7-0               
    ## [47] tidyselect_1.2.1            rstudioapi_0.16.0          
    ## [49] abind_1.4-5                 yaml_2.3.9                 
    ## [51] codetools_0.2-20            curl_5.2.1                 
    ## [53] lattice_0.22-6              tibble_3.2.1               
    ## [55] withr_3.0.0                 Biobase_2.65.0             
    ## [57] basilisk.utils_1.17.0       KEGGREST_1.45.1            
    ## [59] evaluate_0.24.0             crisprScoreData_1.9.0      
    ## [61] BiocFileCache_2.13.0        xml2_1.3.6                 
    ## [63] ExperimentHub_2.13.0        pillar_1.9.0               
    ## [65] BiocManager_1.30.23         filelock_1.0.3             
    ## [67] MatrixGenerics_1.17.0       crisprScore_1.9.1          
    ## [69] generics_0.1.3              vroom_1.6.5                
    ## [71] RCurl_1.98-1.16             BiocVersion_3.20.0         
    ## [73] hms_1.1.3                   glue_1.7.0                 
    ## [75] tools_4.4.1                 AnnotationHub_3.13.0       
    ## [77] GenomicAlignments_1.41.0    XML_3.99-0.17              
    ## [79] grid_4.4.1                  AnnotationDbi_1.67.0       
    ## [81] GenomeInfoDbData_1.2.12     basilisk_1.17.0            
    ## [83] restfulr_0.0.15             cli_3.6.3                  
    ## [85] rappdirs_0.3.3              fansi_1.0.6                
    ## [87] S4Arrays_1.5.5              dplyr_1.1.4                
    ## [89] crisprBowtie_1.9.0          digest_0.6.36              
    ## [91] SparseArray_1.5.25          rjson_0.2.21               
    ## [93] memoise_2.0.1               htmltools_0.5.8.1          
    ## [95] lifecycle_1.0.4             httr_1.4.7                 
    ## [97] bit64_4.0.5

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-csm1" class="csl-entry">

Colognori, David, Marena Trinidad, and Jennifer A Doudna. 2023. “Precise
Transcript Targeting by CRISPR-Csm Complexes.” *Nature Biotechnology* 41
(9): 1256–64.

</div>

<div id="ref-csm2" class="csl-entry">

Xia, Chenglong, David Colognori, Xueyang Jiang, Ke Xu, and Jennifer A
Doudna. 2024. “Single-Molecule Live-Cell RNA Imaging with CRISPR-Csm.”
*bioRxiv*, 2024–07.

</div>

</div>
