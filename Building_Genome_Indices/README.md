Building genome indices off-target alignment
================
Jean-Philippe Fortin, Luke Hoberecht

-   <a href="#introduction" id="toc-introduction">Introduction</a>
-   <a href="#installation" id="toc-installation">Installation</a>
-   <a href="#building-a-bowtie-index"
    id="toc-building-a-bowtie-index">Building a bowtie index</a>
-   <a href="#building-a-bwa-index" id="toc-building-a-bwa-index">Building a
    BWA index</a>
-   <a href="#building-a-transcriptome-index"
    id="toc-building-a-transcriptome-index">Building a transcriptome
    index</a>
-   <a href="#reproducibility" id="toc-reproducibility">Reproducibility</a>
-   <a href="#references" id="toc-references">References</a>

# Introduction

This vignette demonstrates how to build genome indices for the purpose
of performing on- and off-target alignment. In particular, we show how
to build such indices for the short read aligners bowtie (Langmead et
al. 2009), as used by the `Rbowtie` and `crisprBowtie` packages, and
BWA-backtrack (Li and Durbin 2009), as used by the `Rbwa` and
`crisprBwa` packages. Note that BWA is not available for Windows users.

Generating a genome index file is time consuming, but only needs to be
done once for a given genome.

# Installation

See the [Installation
tutorial](https://github.com/crisprVerse/Tutorials/tree/master/Installation)
to learn how to install the `crisprBowtie` and `crisprBwa` packages.

# Building a bowtie index

In the following example, we build a bowtie index for the human genome
using the hg38 build. First, users will need to donwload the FASTA file
from the UCSC genome browser. Here’s the link:
<https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz>

Next, assuming the `hg38.fa.gz` is located in the current directory, we
build the bowtie genome index using the function `bowtie_build` from the
`Rbowtie` package (which is installed when `crisprBowtie` is installed):

``` r
library(Rbowtie)
fastaFile <- "./hg38.fa.gz"
bowtie_build(fastaFile,
             outdir="./",
             force=TRUE,
             prefix="hg38")
```

This should take a couple of hours to run, and the resulting bowtie
index files will be located in the folder `./hg38` and can be used to
run bowtie alignment. See the
[crisprBowtie](https://github.com/crisprVerse/crisprBowtie) package to
learn how to perform a bowtie alignment within R.

# Building a BWA index

Building a BWA index is similar to building a bowtie index. Assuming the
`hg38.fa.gz` is located in the current directory, we build the BWA
genome index using the function `bwa_build_index` from the `Rbwa`
package (which is installed when `crisprBwa` is installed):

``` r
library(Rbwa)
fastaFile <- "./hg38.fa.gz"
bwa_build_index(fastaFile,
                index_prefix="hg38")
```

This should take a couple of hours to run, and the resulting BWA index
files will be located in the folder `./hg38` and can be used to run BWA
alignment. See the [crisprBwa](https://github.com/crisprVerse/crisprBwa)
package to learn how to perform a BWA alignment within R.

# Building a transcriptome index

For applications using RNA-targeting nucleases such as CasRx, off-target
search is performed against against transcriptomes rather than genomes.
Building a transcriptome index works similar, except that we first need
to generate a FASTA file containing the transcriptome sequences. This is
easily accomplished with the function `getMrnaSequences` from the
`crisprDesign` package, assuming that a gene model is provided, as well
as a `BSgenome` object containing the DNA sequences for the hg38 genome
(`BSgenome.Hsapiens.UCSC.hg38`).

We first load the necessary packages

``` r
library(BSgenome.Hsapiens.UCSC.hg38)
library(crisprDesign)
```

The `crisprDesignData` package (see Installation) contains a gene model
annotation for the hg38 genome, and can be loaded using the following:

``` r
library(crisprDesignData)
data("txdb_human", package="crisprDesignData")
```

See the [Gene annotation
tutorial](https://github.com/crisprVerse/Tutorials/tree/master/Building_Gene_Annotation)
to learn more about how to build such gene annotation objects.

We will now extract mRNA sequences for all available transcripts:

``` r
txids <- unique(txdb_human$exons$tx_id)
mrnasHuman <- getMrnaSequences(txids,
                               bsgenome=BSgenome.Hsapiens.UCSC.hg38,
                               txObject=txdb_human)
```

This should take less than an hour to run. Once completed, we will write
the extracted mRNA sequences to disk using the FASTA format. This can be
accomplished using the `writeXStringSet` function from the `Biostrings`
package:

``` r
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
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_4.2.1   magrittr_2.0.3   fastmap_1.1.0    cli_3.3.0       
    ##  [5] tools_4.2.1      htmltools_0.5.3  rstudioapi_0.14  yaml_2.3.5      
    ##  [9] stringi_1.7.8    rmarkdown_2.15.2 knitr_1.40       stringr_1.4.1   
    ## [13] xfun_0.32        digest_0.6.29    rlang_1.0.4      evaluate_0.16

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
