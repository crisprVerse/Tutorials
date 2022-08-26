Using crisprDesign to design gRNAs for custom sequences
================
Jean-Philippe Fortin, Luke Hoberecht

-   <a href="#introduction" id="toc-introduction">Introduction</a>
-   <a href="#installation" id="toc-installation">Installation</a>
-   <a href="#use-case-designing-grnas-against-egfp"
    id="toc-use-case-designing-grnas-against-egfp">Use case: designing gRNAs
    against EGFP</a>
    -   <a href="#loading-necessary-packages"
        id="toc-loading-necessary-packages">Loading necessary packages</a>
    -   <a href="#obtaining-the-dna-sequence"
        id="toc-obtaining-the-dna-sequence">Obtaining the DNA sequence</a>
    -   <a href="#constructing-the-guideset-object"
        id="toc-constructing-the-guideset-object">Constructing the
        <code>GuideSet</code> object:</a>
    -   <a
        href="#finding-off-targets-in-the-human-genome-to-find-grnas-specific-to-egfp"
        id="toc-finding-off-targets-in-the-human-genome-to-find-grnas-specific-to-egfp">Finding
        off-targets in the human genome to find gRNAs specific to EGFP</a>
    -   <a href="#predicting-on-target-activity"
        id="toc-predicting-on-target-activity">Predicting on-target activity</a>
    -   <a href="#final-selection" id="toc-final-selection">Final selection</a>
-   <a href="#session-info" id="toc-session-info">Session Info</a>
-   <a href="#references" id="toc-references">References</a>

# Introduction

In this tutorial, we illustrate the main functionalities of
`crisprDesign` for designing gRNAs for custom sequences. To design gRNAs
for targets located in an organism genome, see the [introductory
CRISPRko
tutorial](https://github.com/crisprVerse/Tutorials/tree/master/Design_CRISPRko_Cas9).

# Installation

See the [Installation
tutorial](https://github.com/crisprVerse/Tutorials/tree/master/Installation)
to learn how to install the packages necessary for this tutorial:
`crisprDesign`, `crisprDesignData`

# Use case: designing gRNAs against EGFP

Suppose we are engineering a human cell line to express the enhanced
green fluorescent protein (EGFP) marker, and that we want to design
gRNAs that knockout EGFP as experimental controls. Such control gRNAs
should target EGFP with (1) high efficiency, and (2) should be specific
to EGFP, that is, should not target the cell genome (human genome in
this case). Supposed also that the cell line is also stably expressing
SpCas9.

## Loading necessary packages

We first start by loading the necessary packages:

``` r
library(Biostrings)
library(crisprBase)
library(crisprDesign)
library(crisprDesignData)
library(BSgenome.Hsapiens.UCSC.hg38)
```

## Obtaining the DNA sequence

In the folder `data`, we have included a fasta file containing the DNA
sequence of the EGFP marker. The sequence was obtained from the
[SnapGene
website](https://www.snapgene.com/resources/plasmid-files/?set=fluorescent_protein_genes_and_plasmids&plasmid=EGFP)

We can read in the fasta file using the `readDNAStringSet` function from
the package `Biostrings`:

``` r
dna <- Biostrings::readDNAStringSet("data/egfp.fa")
names(dna) <- "EGFP"
dna
```

    ## DNAStringSet object of length 1:
    ##     width seq                                               names               
    ## [1]   720 ATGGTGAGCAAGGGCGAGGAGCT...GCATGGACGAGCTGTACAAGTAA EGFP

This could also be simply constructed from a regular string:

``` r
dna <- "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA"
dna <- DNAStringSet(dna)
names(dna) <- "EGFP"
```

(Note that the function also accepts a simple string, which would be
internally converted into a `DNAStringSet`). This is the custom sequence
input that we will use to design gRNAs.

## Constructing the `GuideSet` object:

Next, we design all possible SpCas9 gRNAs targeting EGFP. First, we load
the SpCas9 object from the `crisprBase` package:

``` r
data(SpCas9, package="crisprBase")
```

and we design gRNAs using the function `findSpacers` from
`crisprDesign`:

``` r
gs <- findSpacers(dna, 
                  crisprNuclease=SpCas9)
head(gs)
```

    ## GuideSet object with 6 ranges and 5 metadata columns:
    ##            seqnames    ranges strand |          protospacer            pam
    ##               <Rle> <IRanges>  <Rle> |       <DNAStringSet> <DNAStringSet>
    ##   spacer_1     EGFP        30      + | AAGGGCGAGGAGCTGTTCAC            CGG
    ##   spacer_2     EGFP        31      + | AGGGCGAGGAGCTGTTCACC            GGG
    ##   spacer_3     EGFP        31      - | GACCAGGATGGGCACCACCC            CGG
    ##   spacer_4     EGFP        32      + | GGGCGAGGAGCTGTTCACCG            GGG
    ##   spacer_5     EGFP        35      + | CGAGGAGCTGTTCACCGGGG            TGG
    ##   spacer_6     EGFP        42      - | CCGTCCAGCTCGACCAGGAT            GGG
    ##             pam_site  cut_site      region
    ##            <numeric> <numeric> <character>
    ##   spacer_1        30        27        EGFP
    ##   spacer_2        31        28        EGFP
    ##   spacer_3        31        34        EGFP
    ##   spacer_4        32        29        EGFP
    ##   spacer_5        35        32        EGFP
    ##   spacer_6        42        45        EGFP
    ##   -------
    ##   seqinfo: 1 sequence from custom genome
    ##   crisprNuclease: SpCas9

The resulting output is a regular `GuideSet` object, and all
functionalities described in the [introductory CRISPRko
tutorial](https://github.com/crisprVerse/Tutorials/tree/master/Design_CRISPRko_Cas9)
can be applied here as well.

There are a few key differences to note with respect to a `GuideSet`
object constructed using a reference genome. First, the name of the
input DNA sequence (EGFP) is used as the chromosome name stored in the
`seqnames` field. Second, the `pam_site` and `cut_site` coordinates are
all relative to the first nucleotide of the custom DNA sequence.
Finally, the `GuideSet` object stores the input sequence, which can be
accessed using the function `customSequences`:

``` r
customSequences(gs)
```

    ## DNAStringSet object of length 1:
    ##     width seq                                               names               
    ## [1]   720 ATGGTGAGCAAGGGCGAGGAGCT...GCATGGACGAGCTGTACAAGTAA EGFP

## Finding off-targets in the human genome to find gRNAs specific to EGFP

Now that we have designed all possible gRNAs targeting EGFP, we will
filter out gRNAs that have on- and off-targets located in the human
genome. We will use the bowtie aligner to find targets, so we need to
first specify the path of a bowtie index constructed on the human
genome:

``` r
# Path of the hg38 bowtie index on my personal laptop:
bowtie_index <- "/Users/fortinj2/crisprIndices/bowtie/hg38/hg38"
```

For instructions on how to build a Bowtie index from a given reference
genome, see the [genome index
tutorial](https://github.com/crisprVerse/Tutorials/tree/master/Building_Genome_Indices).

To annotate off-targets with genomic context, for instance to know
whether or not they are located in coding regions, we will also need a
gene model object. We will use the gene model object `txdb_human` from
`crisprDesignData`, which contains genomic coordinates of all human
protein-coding genes. See the [crisprDesignData
package](https://github.com/crisprVerse/crisprDesignData) for more
details.

``` r
data(txdb_human, package="crisprDesignData")
```

We are now ready to find all on- and off-targets using the
`addSpacerAlignments` function from `crisprDesign`:

``` r
gs <- addSpacerAlignments(gs,
                          aligner="bowtie",
                          aligner_index=bowtie_index,
                          bsgenome=BSgenome.Hsapiens.UCSC.hg38,
                          n_mismatches=3,
                          txObject=txdb_human)
```

    ## [runCrisprBowtie] Using BSgenome.Hsapiens.UCSC.hg38 
    ## [runCrisprBowtie] Searching for SpCas9 protospacers

``` r
gs
```

    ## GuideSet object with 119 ranges and 14 metadata columns:
    ##              seqnames    ranges strand |          protospacer            pam
    ##                 <Rle> <IRanges>  <Rle> |       <DNAStringSet> <DNAStringSet>
    ##     spacer_1     EGFP        30      + | AAGGGCGAGGAGCTGTTCAC            CGG
    ##     spacer_2     EGFP        31      + | AGGGCGAGGAGCTGTTCACC            GGG
    ##     spacer_3     EGFP        31      - | GACCAGGATGGGCACCACCC            CGG
    ##     spacer_4     EGFP        32      + | GGGCGAGGAGCTGTTCACCG            GGG
    ##     spacer_5     EGFP        35      + | CGAGGAGCTGTTCACCGGGG            TGG
    ##          ...      ...       ...    ... .                  ...            ...
    ##   spacer_115     EGFP       684      + | CTGGAGTTCGTGACCGCCGC            CGG
    ##   spacer_116     EGFP       685      + | TGGAGTTCGTGACCGCCGCC            GGG
    ##   spacer_117     EGFP       685      - | GTCCATGCCGAGAGTGATCC            CGG
    ##   spacer_118     EGFP       696      + | ACCGCCGCCGGGATCACTCT            CGG
    ##   spacer_119     EGFP       701      + | CGCCGGGATCACTCTCGGCA            TGG
    ##               pam_site  cut_site      region        n0        n1        n2
    ##              <numeric> <numeric> <character> <numeric> <numeric> <numeric>
    ##     spacer_1        30        27        EGFP         0         0         1
    ##     spacer_2        31        28        EGFP         0         0         2
    ##     spacer_3        31        34        EGFP         0         0         1
    ##     spacer_4        32        29        EGFP         0         0         0
    ##     spacer_5        35        32        EGFP         0         0         0
    ##          ...       ...       ...         ...       ...       ...       ...
    ##   spacer_115       684       681        EGFP         0         0         1
    ##   spacer_116       685       682        EGFP         0         0         0
    ##   spacer_117       685       688        EGFP         0         0         0
    ##   spacer_118       696       693        EGFP         0         0         0
    ##   spacer_119       701       698        EGFP         0         0         0
    ##                     n3      n0_c      n1_c      n2_c      n3_c
    ##              <numeric> <numeric> <numeric> <numeric> <numeric>
    ##     spacer_1         7         0         0         0         2
    ##     spacer_2         6         0         0         0         0
    ##     spacer_3        23         0         0         0         3
    ##     spacer_4         6         0         0         0         1
    ##     spacer_5         4         0         0         0         2
    ##          ...       ...       ...       ...       ...       ...
    ##   spacer_115         5         0         0         0         0
    ##   spacer_116         2         0         0         0         0
    ##   spacer_117         2         0         0         0         0
    ##   spacer_118         2         0         0         0         1
    ##   spacer_119         0         0         0         0         0
    ##                                                          alignments
    ##                                                       <GRangesList>
    ##     spacer_1   chr4:151052480:-,chr1:29350153:-,chrX:44232571:-,...
    ##     spacer_2  chr6:115711229:+,chr6:52987132:+,chr3:186748384:+,...
    ##     spacer_3 chr6:149167095:+,chr17:37844933:-,chr17:82019964:-,...
    ##     spacer_4 chr17:82484917:+,chr11:35064010:-,chr18:48539310:-,...
    ##     spacer_5    chr4:426282:-,chr19:16847037:+,chr19:14471687:+,...
    ##          ...                                                    ...
    ##   spacer_115  chr2:142494626:+,chr1:54281906:+,chr2:117131280:-,...
    ##   spacer_116                      chr4:139453357:-,chr17:82078529:+
    ##   spacer_117                        chr1:53536150:+,chr7:24540197:+
    ##   spacer_118                       chr18:77345044:+,chr3:51662810:-
    ##   spacer_119                                                       
    ##   -------
    ##   seqinfo: 1 sequence from custom genome
    ##   crisprNuclease: SpCas9

## Predicting on-target activity

We also want to make sure to filter out gRNAs that are predicted to have
poor on-target activity. To do so, we annotate gRNAs with the DeepHF
on-target activity score:

``` r
gs <- addOnTargetScores(gs, methods="deephf")
```

    ## [addOnTargetScores] Adding deephf scores.

    ## snapshotDate(): 2022-08-23

    ## see ?crisprScoreData and browseVignettes('crisprScoreData') for documentation

    ## loading from cache

Finally, we characterize the spacer sequences using the
`addSequenceFeatures` function from `crisprDesign`:

``` r
gs <- addSequenceFeatures(gs)
```

## Final selection

For our use case, we will only retain gRNAs that do not map to the human
genome (`n0=0`), don’t have any 1 or 2-mismatch off-targets (`n1=0` and
`n2=0`), and do not have 3-mismatch off-targets located in coding
regions (`n3_c=0`):

``` r
gs <- gs[gs$n0==0 & gs$n1==0 & gs$n2==0 & gs$n3_c==0]
```

We also remove gRNAs that contain polyT sequences

``` r
gs <- gs[!gs$polyT,]
```

and only keep gRNAs that don’t have extreme GC content:

``` r
gs <- gs[gs$percentGC>=20 & gs$percentGC<=80]
```

Finally, we rank gRNAs from the highest to the lowest on-target activity
score:

``` r
gs <- gs[order(-gs$score_deephf)]
head(gs)
```

    ## GuideSet object with 6 ranges and 21 metadata columns:
    ##              seqnames    ranges strand |          protospacer            pam
    ##                 <Rle> <IRanges>  <Rle> |       <DNAStringSet> <DNAStringSet>
    ##    spacer_64     EGFP       359      + | GAAGTTCGAGGGCGACACCC            TGG
    ##   spacer_114     EGFP       682      - | CATGCCGAGAGTGATCCCGG            CGG
    ##    spacer_77     EGFP       446      - | CGGCCATGATATAGACGTTG            TGG
    ##    spacer_16     EGFP        98      + | CAAGTTCAGCGTGTCCGGCG            AGG
    ##    spacer_45     EGFP       229      - | GTCGTGCTGCTTCATGTGGT            CGG
    ##   spacer_107     EGFP       635      - | TGTGATCGCGCTTCTCGTTG            GGG
    ##               pam_site  cut_site      region        n0        n1        n2
    ##              <numeric> <numeric> <character> <numeric> <numeric> <numeric>
    ##    spacer_64       359       356        EGFP         0         0         0
    ##   spacer_114       682       685        EGFP         0         0         0
    ##    spacer_77       446       449        EGFP         0         0         0
    ##    spacer_16        98        95        EGFP         0         0         0
    ##    spacer_45       229       232        EGFP         0         0         0
    ##   spacer_107       635       638        EGFP         0         0         0
    ##                     n3      n0_c      n1_c      n2_c      n3_c
    ##              <numeric> <numeric> <numeric> <numeric> <numeric>
    ##    spacer_64         0         0         0         0         0
    ##   spacer_114         0         0         0         0         0
    ##    spacer_77         1         0         0         0         0
    ##    spacer_16         3         0         0         0         0
    ##    spacer_45         9         0         0         0         0
    ##   spacer_107         2         0         0         0         0
    ##                                                       alignments score_deephf
    ##                                                    <GRangesList>    <numeric>
    ##    spacer_64                                                         0.716188
    ##   spacer_114                                                         0.700199
    ##    spacer_77                                    chr3:140247159:+     0.686111
    ##    spacer_16    chr5:132092368:+,chr6:42782483:+,chrX:97074134:+     0.670439
    ##    spacer_45 chr7:87066799:+,chr4:89369669:-,chrX:82007103:+,...     0.664180
    ##   spacer_107                    chr2:128980117:+,chr1:53729319:-     0.654151
    ##              percentGC     polyA     polyC     polyG     polyT startingGGGGG
    ##              <numeric> <logical> <logical> <logical> <logical>     <logical>
    ##    spacer_64        65     FALSE     FALSE     FALSE     FALSE         FALSE
    ##   spacer_114        65     FALSE     FALSE     FALSE     FALSE         FALSE
    ##    spacer_77        50     FALSE     FALSE     FALSE     FALSE         FALSE
    ##    spacer_16        65     FALSE     FALSE     FALSE     FALSE         FALSE
    ##    spacer_45        55     FALSE     FALSE     FALSE     FALSE         FALSE
    ##   spacer_107        55     FALSE     FALSE     FALSE     FALSE         FALSE
    ##   -------
    ##   seqinfo: 1 sequence from custom genome
    ##   crisprNuclease: SpCas9

Users can select the top gRNAs as their control gRNAs.

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
    ##  [1] crisprScoreData_1.1.3             ExperimentHub_2.5.0              
    ##  [3] AnnotationHub_3.5.0               BiocFileCache_2.5.0              
    ##  [5] dbplyr_2.2.1                      BSgenome.Hsapiens.UCSC.hg38_1.4.4
    ##  [7] BSgenome_1.65.2                   rtracklayer_1.57.0               
    ##  [9] GenomicRanges_1.49.1              crisprDesignData_0.99.17         
    ## [11] crisprDesign_0.99.133             crisprBase_1.1.5                 
    ## [13] Biostrings_2.65.2                 GenomeInfoDb_1.33.5              
    ## [15] XVector_0.37.0                    IRanges_2.31.2                   
    ## [17] S4Vectors_0.35.1                  BiocGenerics_0.43.1              
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bitops_1.0-7                  matrixStats_0.62.0           
    ##  [3] bit64_4.0.5                   filelock_1.0.2               
    ##  [5] progress_1.2.2                httr_1.4.4                   
    ##  [7] tools_4.2.1                   utf8_1.2.2                   
    ##  [9] R6_2.5.1                      DBI_1.1.3                    
    ## [11] tidyselect_1.1.2              prettyunits_1.1.1            
    ## [13] bit_4.0.4                     curl_4.3.2                   
    ## [15] compiler_4.2.1                crisprBowtie_1.1.1           
    ## [17] cli_3.3.0                     Biobase_2.57.1               
    ## [19] basilisk.utils_1.9.1          xml2_1.3.3                   
    ## [21] DelayedArray_0.23.1           randomForest_4.7-1.1         
    ## [23] readr_2.1.2                   rappdirs_0.3.3               
    ## [25] stringr_1.4.1                 digest_0.6.29                
    ## [27] Rsamtools_2.13.4              rmarkdown_2.15.2             
    ## [29] crisprScore_1.1.14            basilisk_1.9.2               
    ## [31] pkgconfig_2.0.3               htmltools_0.5.3              
    ## [33] MatrixGenerics_1.9.1          fastmap_1.1.0                
    ## [35] rlang_1.0.4                   rstudioapi_0.14              
    ## [37] RSQLite_2.2.16                shiny_1.7.2                  
    ## [39] BiocIO_1.7.1                  generics_0.1.3               
    ## [41] jsonlite_1.8.0                vroom_1.5.7                  
    ## [43] BiocParallel_1.31.12          dplyr_1.0.9                  
    ## [45] VariantAnnotation_1.43.3      RCurl_1.98-1.8               
    ## [47] magrittr_2.0.3                GenomeInfoDbData_1.2.8       
    ## [49] Matrix_1.4-1                  Rcpp_1.0.9                   
    ## [51] fansi_1.0.3                   reticulate_1.25              
    ## [53] Rbowtie_1.37.0                lifecycle_1.0.1              
    ## [55] stringi_1.7.8                 yaml_2.3.5                   
    ## [57] SummarizedExperiment_1.27.1   zlibbioc_1.43.0              
    ## [59] grid_4.2.1                    blob_1.2.3                   
    ## [61] promises_1.2.0.1              parallel_4.2.1               
    ## [63] crayon_1.5.1                  crisprBwa_1.1.3              
    ## [65] dir.expiry_1.5.0              lattice_0.20-45              
    ## [67] GenomicFeatures_1.49.6        hms_1.1.2                    
    ## [69] KEGGREST_1.37.3               knitr_1.40                   
    ## [71] pillar_1.8.1                  rjson_0.2.21                 
    ## [73] codetools_0.2-18              biomaRt_2.53.2               
    ## [75] BiocVersion_3.16.0            XML_3.99-0.10                
    ## [77] glue_1.6.2                    evaluate_0.16                
    ## [79] BiocManager_1.30.18           httpuv_1.6.5                 
    ## [81] png_0.1-7                     vctrs_0.4.1                  
    ## [83] tzdb_0.3.0                    purrr_0.3.4                  
    ## [85] assertthat_0.2.1              cachem_1.0.6                 
    ## [87] xfun_0.32                     mime_0.12                    
    ## [89] Rbwa_1.1.0                    xtable_1.8-4                 
    ## [91] restfulr_0.0.15               later_1.3.0                  
    ## [93] tibble_3.1.8                  GenomicAlignments_1.33.1     
    ## [95] AnnotationDbi_1.59.1          memoise_2.0.1                
    ## [97] interactiveDisplayBase_1.35.0 ellipsis_0.3.2

# References
