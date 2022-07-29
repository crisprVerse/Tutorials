Building genome indices for short read aligners
================

-   [Introduction](#introduction)
-   [Installation](#installation)
    -   [OS Requirements](#os-requirements)
    -   [R Dependencies](#r-dependencies)
-   [Building a genome index](#building-a-genome-index)
    -   [Bowtie index](#bowtie-index)
    -   [Bwa index](#bwa-index)
-   [Building a transcriptome index](#building-a-transcriptome-index)
-   [Reproducibility](#reproducibility)
-   [References](#references)

Authors: Jean-Philippe Fortin, Luke Hoberecht

Date: 29 July, 2022

# Introduction

This vignette demonstrates how to build genome indices for use with the
short read aligners bowtie (Langmead et al. 2009), as used by the
`Rbowtie` and `crisprBowtie` packages, and BWA-backtrack (Li and Durbin
2009), as used by the `Rbwa` and `crisprBwa` packages.

# Installation

### OS Requirements

The `crisprBowtie` package is supported for macOS, Linux and Windows
machines. The `crisprBwa` package is supported for macOS and Linux only.
Both packages were developed and tested on R version 4.2.

### R Dependencies

-   RBowtie:
    <https://bioconductor.org/packages/release/bioc/html/Rbowtie.html>
-   RBwa: <https://github.com/Jfortin1/Rbwa>

# Building a genome index

A genome index file is necessary to use the aligner functions in the
`crisprBowtie` package (`runBowtie` and `runCrisprBowtie`) and the
`crisprBwa` package (`runBwa` and `runCrisprBwa`). For a given genome,
this step only needs to be done once.

## Bowtie index

The `Rbowtie` package offers the convenient function `bowtie_build`,
which builds a bowtie index for any custom genome from a FASTA file. In
the following example, we build a bowtie index for a small region of the
human chromosome 1 (provided in the `crisprBowtie` package) and save the
index file to a temporary directory:

``` r
library(Rbowtie)
fasta <- file.path(find.package("crisprBowtie"), "example/chr1.fa")
tempDir <- tempdir()
Rbowtie::bowtie_build(fasta,
                      outdir=tempDir,
                      force=TRUE,
                      prefix="myIndex")
```

See the [crisprBowtie](https://github.com/Jfortin1/crisprBowtie) package
for information on how to obtain alignments using this index file.

## Bwa index

Building a BWA index is made simple with the `bwa_build_index` function
from the `Rbwa` package. This function builds the index file for any
custom genome from a FASTA file. As an example, we build a BWA index for
a small portion of the human chromosome 12 (provided in the `crisprBwa`
package) and save the index file to a temporary directory:

``` r
library(Rbwa)
fasta <- system.file(package="crisprBwa", "example/chr12.fa")
outdir <- tempdir()
index <- file.path(outdir, "chr12")
Rbwa::bwa_build_index(fasta,
                      index_prefix=index)
```

See the [crisprBwa](https://github.com/Jfortin1/crisprBwa) package for
information on how to obtain alignments using this index file.

# Building a transcriptome index

For applications using RNA-targeting nucleases such as CasRx, it is
preferable to search for alignments against the transcriptome rather
than the entire genome. To build an index file for the transcriptome we
must first generate a FASTA file containing the transcriptome. This is
easily accomplished with the `crisprDesign` function `getMrnaSequences`,
as shown in the following example:

``` r
library(BSgenome.Hsapiens.UCSC.hg38)
library(crisprDesign)
library(crisprDesignData)
data("txdb_human", package="crisprDesignData")
exon_ids <- unique(txdb_human$exons$tx_id)
mrnasHuman <- getMrnaSequences(exon_ids,
                               bsgenome=BSgenome.Hsapiens.UCSC.hg38,
                               txObject=txdb_human)
library(Biostrings)
writeXStringSet(mrnasHuman,
                file="ensembl_human_104.fasta",
                format="fasta")
```

Note that the `seqnames` of this FASTA file are Ensembl transcript IDs
instead of chromosomes. Once the FASTA file has been generated, the
process for constructing either a bowtie or BWA index file is the same
as described in the above sections.

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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] Rbwa_1.0.0     Rbowtie_1.36.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] digest_0.6.29   magrittr_2.0.3  evaluate_0.15   rlang_1.0.3    
    ##  [5] stringi_1.7.6   cli_3.3.0       rstudioapi_0.13 rmarkdown_2.14 
    ##  [9] tools_4.2.0     stringr_1.4.0   xfun_0.31       yaml_2.3.5     
    ## [13] fastmap_1.1.0   compiler_4.2.0  htmltools_0.5.2 knitr_1.39

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-langmead2009bowtie" class="csl-entry">

Langmead, Ben, Cole Trapnell, Mihai Pop, and Steven L. Salzberg. 2009.
“Ultrafast and Memory-Efficient Alignment of Short DNA Sequences to the
Human Genome.” *Genome Biology* 10 (3): R25.
<https://doi.org/10.1186/gb-2009-10-3-r25>.

</div>

<div id="ref-bwa" class="csl-entry">

Li, Heng, and Richard Durbin. 2009. “Fast and Accurate Short Read
Alignment with Burrows–Wheeler Transform.” *Bioinformatics* 25 (14):
1754–60.

</div>

</div>
