Installing the crisprVerse and packages necessary for the tutorials
================
Jean-Philippe Fortin, Luke Hoberecht

-   <a href="#installation" id="toc-installation">Installation</a>
    -   <a href="#requirements" id="toc-requirements">Requirements</a>
    -   <a href="#bioconductor-versions"
        id="toc-bioconductor-versions">Bioconductor versions</a>
    -   <a href="#installing-the-core-crisprverse-packages"
        id="toc-installing-the-core-crisprverse-packages">Installing the core
        crisprVerse packages</a>
    -   <a href="#installing-data-packages"
        id="toc-installing-data-packages">Installing data packages</a>
    -   <a href="#installing-optional-packages"
        id="toc-installing-optional-packages">Installing optional packages</a>
-   <a href="#reproducibility" id="toc-reproducibility">Reproducibility</a>
-   <a href="#references" id="toc-references">References</a>

# Installation

We show in this tutorial how to install the crisprVerse packages, as
well as other packages necessary for some of the [crisprVerse
tutorials](https://github.com/crisprVerse/Tutorials).

## Requirements

The crisprVerse is supported for macOS, Linux and Windows machines. It
requires R version \>=4.4. Some of the third-party functionalities are
not available for Windows machines (BWA alignment, and some of the
scoring functions). To download and install R, follow the instructions
on the [R-project website](https://www.r-project.org/).

## Bioconductor versions

The Bioconductor project has 2 concurrent branches: `release` and
`devel`. Currently (July 2024), the release branch is `3.19`, and the
devel branch is `3.20`. See the [Bioconductor install
page](https://www.bioconductor.org/install/) for more information
regarding Bioconductor versions.

The crisprVerse ecosystem is currently available on both the
Bioconductor release and devel branches. The release branch, updated
twice a year, is frozen in time and contains stable releases of the
packages. The devel branch is built over night and contains the latest
changes pushed to the GitHub repositories. We recommend using the
release branch, unless new features are only available in the devel
branch.

## Installing the core crisprVerse packages

Type in the following commands in an R session to install the core
crisprVerse packages from the Bioconductor devel branch:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version="devel")
BiocManager::install("crisprVerse")
```

To install the packages from the release branch instead:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version="release")
BiocManager::install("crisprVerse")
```

This will install the following packages:

-   [crisprBase](https://github.com/crisprVerse/crisprBase) to specify
    and manipulate CRISPR nucleases.
-   [crisprBowtie](https://github.com/crisprVerse/crisprBowtie) to
    perform gRNA spacer sequence alignment with Bowtie.
-   [crisprScore](https://github.com/crisprVerse/crisprScore) to
    annotate gRNAs with on-target and off-target scores.
-   [crisprDesign](https://github.com/crisprVerse/crisprDesign) to
    design and manipulate gRNAs with `GuideSet` objects.
-   [crisprScoreData](https://github.com/crisprVerse/crisprScoreData) to
    use pre-trained models for the `crisprScore` package.

The following command will load all of those packages in an R session:

``` r
library(crisprVerse)
```

You can check that all crisprVerse packages are up-to-date with
`crisprVerse_update()`:

``` r
crisprVerse_update()
```

## Installing data packages

The following genome data packages from Bioconductor are required for
several of the tutorials:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version="devel")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38.dbSNP151.major")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38.dbSNP151.minor")
```

The [crisrpDesignData](https://github.com/crisprVerse/crisprDesignData)
package is also required for most of the tutorials and can be installed
directly from our GitHub page using the `devtools` package:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install.packages("crisprVerse/crisprDesignData")
```

## Installing optional packages

For maxOS and Linux users, the
[crisprBwa](https://github.com/crisprVerse/crisprBwa) can be installed
from Bioconductor using the following:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version="devel")
BiocManager::install("crisprBwa")
```

The [crisprViz](https://github.com/crisprVerse/crisprViz) package is
currently under review at Bioconductor, but can be installed directly
from GitHub:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install.packages("crisprVerse/crisprViz")
```

# Reproducibility

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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_4.4.1    fastmap_1.2.0     cli_3.6.3         tools_4.4.1      
    ##  [5] htmltools_0.5.8.1 rstudioapi_0.16.0 yaml_2.3.9        rmarkdown_2.27   
    ##  [9] knitr_1.48        xfun_0.46         digest_0.6.36     rlang_1.1.4      
    ## [13] evaluate_0.24.0

# References
