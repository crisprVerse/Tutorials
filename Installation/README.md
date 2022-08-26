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
-   <a href="#reproducibility-1"
    id="toc-reproducibility-1">Reproducibility</a>
-   <a href="#references" id="toc-references">References</a>

# Installation

We show in this tutorial how to install the crisprVerse packages, as
well as other packages necessary for some of the [crisprVerse
tutorials](https://github.com/crisprVerse/Tutorials).

## Requirements

The crisprVerse is supported for macOS, Linux and Windows machines. It
requires R version \>=4.2.1. Some of the third-party functionalities are
not available for Windows machines (BWA alignment, and some of the
scoring functions). To download and install R, see the [R-project
website](https://www.r-project.org/).

## Bioconductor versions

The Bioconductor project has 2 concurrent versions: `release` and
`devel`. Currently (August 2022), the release version is`3.15`, and the
devel version is `3.16`. Release versions are created twice a year.

The current version of the crisprVerse was developed on the devel
version of Bioconductor (`3.16`) to make sure it accesses all of the
latest developments. Earlier versions of some of our packages are
available on the release version, but we do not recommend using the
release version as most of the functionalities described in the
tutorials require devel functionalities.

See the [Bioconductor install
page](https://www.bioconductor.org/install/) for more information re.
Bioconductor.

## Installing the core crisprVerse packages

The [crisprVerse package](https://github.com/crisprVerse/crisprVerse)
installs the core crisprVerse packages in a single command from
Bioconductor. Simply type in the following commands in an R session:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(version="devel")
BiocManager::install("crisprVerse")
```

Note that we specify the devel branch of Bioconductor so that we can use
the latest functionalities.

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
