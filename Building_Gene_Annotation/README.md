Building a gene annotation object
================
Jean-Philippe Fortin, Luke Hoberecht

-   <a href="#introduction" id="toc-introduction">Introduction</a>
-   <a href="#installation" id="toc-installation">Installation</a>
    -   <a href="#getting-started" id="toc-getting-started">Getting started</a>
-   <a href="#building-gene-annotation-objects"
    id="toc-building-gene-annotation-objects">Building gene annotation
    objects</a>
-   <a href="#building-a-grangeslist-from-ensembl"
    id="toc-building-a-grangeslist-from-ensembl">Building a GRangesList from
    Ensembl</a>
    -   <a href="#building-a-tssobject" id="toc-building-a-tssobject">Building a
        tssObject</a>
-   <a href="#using-gene-annotation-objects"
    id="toc-using-gene-annotation-objects">Using gene annotation objects</a>
-   <a href="#building-a-gene-annotation-object-from-a-gff-file"
    id="toc-building-a-gene-annotation-object-from-a-gff-file">Building a
    gene annotation object from a GFF file</a>
-   <a href="#reproducibility" id="toc-reproducibility">Reproducibility</a>

# Introduction

In this tutorial, we describe the process for making and using  
rich gene annotation objects to be used throughout the crisprVerse
ecosystem. Such objects enable users to retrieve coordinates of
transcripts, exons, etc. Those objects are also used by several
functions in the [crisprDesign
package](https://github.com/crisprVerse/crisprDesign) to add gene
annotations to both gRNA on-targets and off-targets. This is what the
`txObject` argument in many of the functions expect.

We will also describe the process for constructing and using a
transcription start site (TSS) annotation object (`tssObject` argument
in many of the functions).

# Installation

See the [Installation
tutorial](https://github.com/crisprVerse/Tutorials/tree/master/Installation)
to learn how to install the packages `crisprDesign` and
`crisprDesignData` required in this tutorial.

### Getting started

The packages can be loaded into an R session in the usual way:

``` r
library(crisprDesign)
library(crisprDesignData)
```

# Building gene annotation objects

In the crisprVerse, we represent gene annotations using `GRangesList`
object, and this can be easily constructed using the commonly-used
Bioconductor objects `TxDb` (see the [GenomicFeatures
package](https://bioconductor.riken.jp/packages/3.2/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf)
to learn more about `TxDb` objects). We will now show several ways of
constructing such objects.

# Building a GRangesList from Ensembl

We construct a gene annotation object for the human genome using the
Ensembl release 104 (hg38). This can be done using the function
`getTxDb` in `crisprDesign`:

``` r
txdb <- getTxDb(organism="Homo sapiens", release=104)
```

This may take several minutes, and note that this requires an internet
connection. In case it times out, one can increase the timeout option
using the following:

``` r
options(timeout = max(10000000, getOption("timeout")))
```

Once obtained, we can convert the object into a `GRangesList` using the
function `TxDb2GRangesList` from `crisprDesign`:

``` r
grList <- TxDb2GRangesList(txdb)
```

We will specify that the genome is hg38:

``` r
GenomeInfoDb::genome(grList) <- "hg38"
```

And that’s it! The `grList` object contains all of the information about
the Ensembl release 104 gene model, and is ready to be used in the
crisprVerse. Let’s take a quick look at our gene annotation object:

``` r
names(grList)
```

    ## [1] "transcripts" "exons"       "cds"         "fiveUTRs"    "threeUTRs"  
    ## [6] "introns"     "tss"

``` r
grList$transcripts
```

    ## GRanges object with 111751 ranges and 14 metadata columns:
    ##     seqnames      ranges strand |           tx_id         gene_id
    ##        <Rle>   <IRanges>  <Rle> |     <character>     <character>
    ##            1 11869-14409      + | ENST00000456328 ENSG00000223972
    ##            1 12010-13670      + | ENST00000450305 ENSG00000223972
    ##            1 29554-31097      + | ENST00000473358 ENSG00000243485
    ##            1 30267-31109      + | ENST00000469289 ENSG00000243485
    ##            1 30366-30503      + | ENST00000607096 ENSG00000284332
    ##   .      ...         ...    ... .             ...             ...
    ##           MT   5826-5891      - | ENST00000387409 ENSG00000210144
    ##           MT   7446-7514      - | ENST00000387416 ENSG00000210151
    ##           MT 14149-14673      - | ENST00000361681 ENSG00000198695
    ##           MT 14674-14742      - | ENST00000387459 ENSG00000210194
    ##           MT 15956-16023      - | ENST00000387461 ENSG00000210196
    ##          protein_id                tx_type gene_symbol     exon_id exon_rank
    ##         <character>            <character> <character> <character> <integer>
    ##                <NA>   processed_transcript     DDX11L1        <NA>      <NA>
    ##                <NA> transcribed_unproces..     DDX11L1        <NA>      <NA>
    ##                <NA>                 lncRNA MIR1302-2HG        <NA>      <NA>
    ##                <NA>                 lncRNA MIR1302-2HG        <NA>      <NA>
    ##                <NA>                  miRNA   MIR1302-2        <NA>      <NA>
    ##   .             ...                    ...         ...         ...       ...
    ##                <NA>                Mt_tRNA       MT-TY        <NA>      <NA>
    ##                <NA>                Mt_tRNA      MT-TS1        <NA>      <NA>
    ##     ENSP00000354665         protein_coding      MT-ND6        <NA>      <NA>
    ##                <NA>                Mt_tRNA       MT-TE        <NA>      <NA>
    ##                <NA>                Mt_tRNA       MT-TP        <NA>      <NA>
    ##     cds_start   cds_end  tx_start    tx_end   cds_len exon_start  exon_end
    ##     <integer> <integer> <integer> <integer> <integer>  <integer> <integer>
    ##          <NA>      <NA>      <NA>      <NA>      <NA>       <NA>      <NA>
    ##          <NA>      <NA>      <NA>      <NA>      <NA>       <NA>      <NA>
    ##          <NA>      <NA>      <NA>      <NA>      <NA>       <NA>      <NA>
    ##          <NA>      <NA>      <NA>      <NA>      <NA>       <NA>      <NA>
    ##          <NA>      <NA>      <NA>      <NA>      <NA>       <NA>      <NA>
    ##   .       ...       ...       ...       ...       ...        ...       ...
    ##          <NA>      <NA>      <NA>      <NA>      <NA>       <NA>      <NA>
    ##          <NA>      <NA>      <NA>      <NA>      <NA>       <NA>      <NA>
    ##          <NA>      <NA>      <NA>      <NA>      <NA>       <NA>      <NA>
    ##          <NA>      <NA>      <NA>      <NA>      <NA>       <NA>      <NA>
    ##          <NA>      <NA>      <NA>      <NA>      <NA>       <NA>      <NA>
    ##   -------
    ##   seqinfo: 25 sequences (1 circular) from hg38 genome

## Building a tssObject

Building a TSS annotation object requires only one additional step after
constructing the `GRangesList` object described above. This can be
obtained using the function `getTssObjectFromTxObject` in
`crisprDesign`:

``` r
tssObject <- getTssObjectFromTxObject(grList)
tssObject
```

    ## GRanges object with 52547 ranges and 5 metadata columns:
    ##          seqnames    ranges strand |           tx_id         gene_id
    ##             <Rle> <IRanges>  <Rle> |     <character>     <character>
    ##    11402        1     65419      + | ENST00000641515 ENSG00000186092
    ##    11442        1    923923      + | ENST00000616016 ENSG00000187634
    ##    11444        1    925731      + | ENST00000342066 ENSG00000187634
    ##    11445        1    960584      + | ENST00000338591 ENSG00000187961
    ##    11446        1    960639      + | ENST00000622660 ENSG00000187961
    ##      ...      ...       ...    ... .             ...             ...
    ##   123058        Y  24047689      - | ENST00000382407 ENSG00000172352
    ##   123073        Y  24813186      - | ENST00000382365 ENSG00000187191
    ##   123074        Y  24813186      - | ENST00000315357 ENSG00000187191
    ##   123075        Y  24813186      - | ENST00000446723 ENSG00000187191
    ##   123080        Y  25052074      - | ENST00000382287 ENSG00000185894
    ##          gene_symbol        promoter                     ID
    ##          <character>     <character>            <character>
    ##    11402       OR4F5 ENST00000641515 ENSG00000186092_ENST..
    ##    11442      SAMD11 ENST00000616016 ENSG00000187634_ENST..
    ##    11444      SAMD11 ENST00000342066 ENSG00000187634_ENST..
    ##    11445      KLHL17 ENST00000338591 ENSG00000187961_ENST..
    ##    11446        <NA> ENST00000622660 ENSG00000187961_ENST..
    ##      ...         ...             ...                    ...
    ##   123058       CDY1B ENST00000382407 ENSG00000172352_ENST..
    ##   123073        DAZ3 ENST00000382365 ENSG00000187191_ENST..
    ##   123074        DAZ3 ENST00000315357 ENSG00000187191_ENST..
    ##   123075        DAZ3 ENST00000446723 ENSG00000187191_ENST..
    ##   123080       BPY2C ENST00000382287 ENSG00000185894_ENST..
    ##   -------
    ##   seqinfo: 25 sequences (1 circular) from hg38 genome

# Using gene annotation objects

The gene (or TSS) annotation objects described above are often necessary
for the full characterization of CRISPR gRNAs as they as inputs for
several of the `crisprDesign` functions, including `queryTxObject`,
`queryTssObject`, `addGeneAnnotation`, `addTssAnnotation`, and
`addSpacerAlignments`.

For convenience, we provide in the [crisprDesignData
package](https://github.com/crisprVerse/crisprDesignData) precomputed
gene annotation for human and mouse:

| Object name  | Object class  | Version     | Description                                           |
|--------------|---------------|-------------|-------------------------------------------------------|
| `txdb_human` | `GRangesList` | Release 104 | Ensembl gene model for human (hg38/GRCh38)            |
| `txdb_mouse` | `GRangesList` | Release 102 | Ensembl gene model for mouse (mm10/GRCm38)            |
| `tss_human`  | `GRanges`     | Release 104 | Ensembl-based TSS coordinates for human (hg38/GRCh38) |
| `tss_mouse`  | `GRanges`     | Release 102 | Ensembl-based TSS coordinates for human (mm10/GRCm38) |

# Building a gene annotation object from a GFF file

If you have a General Feature Format (GFF) file from which you want to
construct the gene annotation object, you can pass this to the `file`
argument of the `crisprDesign` function `getTxDb`; this will create the
`TxDb` object using the `GenomicFeatures` function `makeTxDbFromGFF`.

# Reproducibility

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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] crisprDesignData_0.99.17 crisprDesign_0.99.133    crisprBase_1.1.5        
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] bitops_1.0-7                  matrixStats_0.62.0           
    ##   [3] lubridate_1.8.0               bit64_4.0.5                  
    ##   [5] filelock_1.0.2                progress_1.2.2               
    ##   [7] httr_1.4.4                    GenomeInfoDb_1.33.5          
    ##   [9] tools_4.2.1                   utf8_1.2.2                   
    ##  [11] R6_2.5.1                      DBI_1.1.3                    
    ##  [13] BiocGenerics_0.43.1           tidyselect_1.1.2             
    ##  [15] prettyunits_1.1.1             bit_4.0.4                    
    ##  [17] curl_4.3.2                    compiler_4.2.1               
    ##  [19] crisprBowtie_1.1.1            cli_3.3.0                    
    ##  [21] Biobase_2.57.1                basilisk.utils_1.9.1         
    ##  [23] crisprScoreData_1.1.3         xml2_1.3.3                   
    ##  [25] DelayedArray_0.23.1           rtracklayer_1.57.0           
    ##  [27] randomForest_4.7-1.1          readr_2.1.2                  
    ##  [29] rappdirs_0.3.3                stringr_1.4.1                
    ##  [31] digest_0.6.29                 Rsamtools_2.13.4             
    ##  [33] rmarkdown_2.15.2              crisprScore_1.1.14           
    ##  [35] basilisk_1.9.2                XVector_0.37.0               
    ##  [37] pkgconfig_2.0.3               htmltools_0.5.3              
    ##  [39] MatrixGenerics_1.9.1          dbplyr_2.2.1                 
    ##  [41] fastmap_1.1.0                 BSgenome_1.65.2              
    ##  [43] RMariaDB_1.2.2                rlang_1.0.4                  
    ##  [45] rstudioapi_0.14               RSQLite_2.2.16               
    ##  [47] shiny_1.7.2                   BiocIO_1.7.1                 
    ##  [49] generics_0.1.3                jsonlite_1.8.0               
    ##  [51] BiocParallel_1.31.12          dplyr_1.0.9                  
    ##  [53] VariantAnnotation_1.43.3      RCurl_1.98-1.8               
    ##  [55] magrittr_2.0.3                GenomeInfoDbData_1.2.8       
    ##  [57] Matrix_1.4-1                  Rcpp_1.0.9                   
    ##  [59] S4Vectors_0.35.1              fansi_1.0.3                  
    ##  [61] reticulate_1.25               Rbowtie_1.37.0               
    ##  [63] lifecycle_1.0.1               stringi_1.7.8                
    ##  [65] yaml_2.3.5                    SummarizedExperiment_1.27.1  
    ##  [67] zlibbioc_1.43.0               BiocFileCache_2.5.0          
    ##  [69] AnnotationHub_3.5.0           grid_4.2.1                   
    ##  [71] blob_1.2.3                    promises_1.2.0.1             
    ##  [73] parallel_4.2.1                ExperimentHub_2.5.0          
    ##  [75] crayon_1.5.1                  dir.expiry_1.5.0             
    ##  [77] lattice_0.20-45               Biostrings_2.65.2            
    ##  [79] GenomicFeatures_1.49.6        hms_1.1.2                    
    ##  [81] KEGGREST_1.37.3               knitr_1.40                   
    ##  [83] pillar_1.8.1                  GenomicRanges_1.49.1         
    ##  [85] rjson_0.2.21                  codetools_0.2-18             
    ##  [87] biomaRt_2.53.2                stats4_4.2.1                 
    ##  [89] BiocVersion_3.16.0            XML_3.99-0.10                
    ##  [91] glue_1.6.2                    evaluate_0.16                
    ##  [93] BiocManager_1.30.18           httpuv_1.6.5                 
    ##  [95] png_0.1-7                     vctrs_0.4.1                  
    ##  [97] tzdb_0.3.0                    purrr_0.3.4                  
    ##  [99] assertthat_0.2.1              cachem_1.0.6                 
    ## [101] xfun_0.32                     mime_0.12                    
    ## [103] xtable_1.8-4                  restfulr_0.0.15              
    ## [105] later_1.3.0                   tibble_3.1.8                 
    ## [107] GenomicAlignments_1.33.1      AnnotationDbi_1.59.1         
    ## [109] memoise_2.0.1                 IRanges_2.31.2               
    ## [111] interactiveDisplayBase_1.35.0 ellipsis_0.3.2
