Building a gRNA database for the human genome using Cas9
================

-   [Introduction](#introduction)
-   [Loading necessary packages](#loading-necessary-packages)
    -   [Specifying the genome](#specifying-the-genome)
    -   [Specifying the genome index](#specifying-the-genome-index)
    -   [Specifying a SNP VCF file](#specifying-a-snp-vcf-file)
    -   [Specifying the nuclease](#specifying-the-nuclease)
    -   [Specifying on-target scoring
        methods](#specifying-on-target-scoring-methods)
    -   [Specifying gene models and TSS
        annotations](#specifying-gene-models-and-tss-annotations)
    -   [Specifying repeat elements](#specifying-repeat-elements)
-   [Building a complete annotation for a given
    gene](#building-a-complete-annotation-for-a-given-gene)
    -   [Converting the `GuideSet` object to a list of
        data.frames](#converting-the-guideset-object-to-a-list-of-dataframes)
-   [Building a complete annotation for all
    genes](#building-a-complete-annotation-for-all-genes)
-   [Considerations for
    CRISPRa/CRISPRi](#considerations-for-crispracrispri)
-   [Considerations for CRISPRkd
    (RfxCas13d)](#considerations-for-crisprkd-rfxcas13d)
-   [Reproducibility](#reproducibility)

Authors: Jean-Philippe Fortin

Date: July 16, 2022

# Introduction

In this tutorial, we provide reproducible code to design and annotates
all gRNAs in a transcriptome, for a given organism, and for a specific
nuclease. As an example, we design CRISPRko gRNAs for the human
transcriptome, in GRCh38 coordinates, for the commonly-used wildtype
nuclease SpCas9.

# Loading necessary packages

We first load the necessary packages:

``` r
library(crisprBase)
library(crisprScore)
library(crisprDesign)
library(crisprDesignData)
library(BSgenome.Hsapiens.UCSC.hg38)
```

### Specifying the genome

We specify a `BSGenome` object that contains the reference genome of the
target organism:

``` r
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
```

### Specifying the genome index

We specify the file path of the genome index needed for Bowtie
alignment:

``` r
bowtie_index <- "/Users/fortinj2/crisprIndices/bowtie/hg38/hg38"
```

See the `crisprBowtie` vignette for instructions about how to create a
Bowtie index from a given reference genome.

### Specifying a SNP VCF file

To add a SNP annotation, we specify a path to a VCF file obtained from
the dbSNP website representing common SNPs for a given dbSNP release:

``` r
vcf <- "/Users/fortinj2/crisprIndices/snps/dbsnp151.grch38/00-common_all.vcf.gz"
```

The VCF file can be obtained from [this
website](https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf).

### Specifying the nuclease

The `crisprBase` provides functionalities for creating custom
`crisprNuclease` objects, and also provides already-available
`crisprNuclease` objects such as the commonly-used `SpCas9`,
`enAsCas12a`, `CasRx` nucleases.

Here, we will use the popular wildtype Cas9 nuclease `SpCas9`:

``` r
data(SpCas9, package="crisprBase")
crisprNuclease=SpCas9
```

### Specifying on-target scoring methods

We specify which on-target scoring methods we want to use:

``` r
scoring_methods=c("crisprscan", "ruleset1")
```

More information about which scoring methods are available for a given
nuclease can be obtained using the following:

``` r
crisprScore::scoringMethodsInfo
```

    ##        method   nuclease left right       type      label len
    ## 1    ruleset1     SpCas9  -24     5  On-target   RuleSet1  30
    ## 2     azimuth     SpCas9  -24     5  On-target    Azimuth  30
    ## 3      deephf     SpCas9  -20     2  On-target     DeepHF  23
    ## 4      lindel     SpCas9  -33    31  On-target     Lindel  65
    ## 5         mit     SpCas9  -20     2 Off-target        MIT  23
    ## 6         cfd     SpCas9  -20     2 Off-target        CFD  23
    ## 7    deepcpf1   AsCas12a   -4    29  On-target   DeepCpf1  34
    ## 8     enpamgb enAsCas12a   -4    29  On-target    EnPAMGB  34
    ## 9  crisprscan     SpCas9  -26     8  On-target CRISPRscan  35
    ## 10    casrxrf      CasRx   NA    NA  On-target   CasRx-RF  NA
    ## 11   crisprai     SpCas9  -19     2  On-target   CRISPRai  22
    ## 12 crisprater     SpCas9  -20    -1  On-target CRISPRater  20
    ## 13 deepspcas9     SpCas9  -24     5  On-target DeepSpCas9  30
    ## 14   ruleset3     SpCas9  -24     5  On-target   RuleSet3  30

### Specifying gene models and TSS annotations

To annotate gRNAs with gene context, we need to specify a gene model
formatted as a `TxDbObject`. To annotate gRNAs with TSS information, we
also need to specify a `GRanges` object containg TSS coordinates. The
`crisprDesignData` contains such objects for both the human and mouse
genomes, in GRCh38 (hg38) and GRCm38 (mm10) coordinates, respectively.
Ensembl gene models were used to generate such objects.

``` r
data(txdb_human, package="crisprDesignData")
data(tss_human, package="crisprDesignData")
txObject <- txdb_human
tssObject <- tss_human
```

The `crisprDesignData` vignette shows how to create such objects from
other transcriptomes.

### Specifying repeat elements

To avoid designing gRNAs in repeat elements, we can specify a `GRanges`
object containing repeats coordinates for a given annotation. Here, we
use the object `gr.repeats.hg38` in `crisprDesignData` that contains
genomic coordinates of the RepeatMasker UCSC track, for the hg38
reference genome:

``` r
data(gr.repeats.hg38, package="crisprDesignData")
grRepeats <- gr.repeats.hg38
```

# Building a complete annotation for a given gene

The `designCompleteAnnotation` function in `crisprDesign` provides a
one-step workflow to design and annotate gRNAs targeting a coding gene
for a user-specific combination of parameters.

Here, we design all gRNAs targeting the human KRAS gene
(ENSG00000133703) for CRISPRko applications using SpCas9:

``` r
gs <- designCompleteAnnotation(queryValue="ENSG00000133703",
                               queryColumn="gene_id",
                               modality="CRISPRko",
                               bsgenome=bsgenome,
                               bowtie_index=bowtie_index,
                               crisprNuclease=SpCas9,
                               txObject=txObject,
                               tssObject=tssObject,
                               grRepeats=grRepeats,
                               vcf=vcf,
                               n_mismatches=1,
                               scoring_methods=scoring_methods)
```

    ## [designCompleteAnnotation] Adding sequence statistics 
    ## [designCompleteAnnotation] Adding spacer alignments

    ## Loading required namespace: crisprBwa

    ## [runCrisprBowtie] Using BSgenome.Hsapiens.UCSC.hg38 
    ## [runCrisprBowtie] Searching for SpCas9 protospacers 
    ## [runCrisprBowtie] Using BSgenome.Hsapiens.UCSC.hg38 
    ## [runCrisprBowtie] Searching for SpCas9 protospacers 
    ## [designCompleteAnnotation] Adding gene annotation 
    ## [designCompleteAnnotation] Adding on-target scores

    ## [addOnTargetScores] Adding ruleset1 scores.

    ## [addOnTargetScores] Adding crisprscan scores.

    ## [designCompleteAnnotation] Adding CFD scores annotation 
    ## [designCompleteAnnotation] Adding SNP annotation

The resulting object is a `GuideSet` object containing the fully
annotated gRNAs, as described in the `crisprDesign` vignette. See the
Tutorial `CRISPRko` to learn more about `GuideSet` objects and how to
access their rich annotations.

``` r
gs
```

    ## GuideSet object with 56 ranges and 27 metadata columns:
    ##                      seqnames    ranges strand |          protospacer
    ##                         <Rle> <IRanges>  <Rle> |       <DNAStringSet>
    ##    ENSG00000133703_1    chr12  25209843      - | AAAGAAAAGATGAGCAAAGA
    ##    ENSG00000133703_2    chr12  25209896      + | TTCTCGAACTAATGTATAGA
    ##    ENSG00000133703_3    chr12  25215438      - | AAATGCATTATAATGTAATC
    ##    ENSG00000133703_4    chr12  25215477      - | AGCAAAGAAGAAAAGACTCC
    ##    ENSG00000133703_5    chr12  25215477      + | TTTTTAATTTTCACACAGCC
    ##                  ...      ...       ...    ... .                  ...
    ##   ENSG00000133703_52    chr12  25245349      - | CTTGTGGTAGTTGGAGCTGG
    ##   ENSG00000133703_53    chr12  25245352      - | AAACTTGTGGTAGTTGGAGC
    ##   ENSG00000133703_54    chr12  25245358      - | GAATATAAACTTGTGGTAGT
    ##   ENSG00000133703_55    chr12  25245365      - | AATGACTGAATATAAACTTG
    ##   ENSG00000133703_56    chr12  25245392      + | TATATTCAGTCATTTTCAGC
    ##                                 pam  pam_site  cut_site      region inRepeats
    ##                      <DNAStringSet> <numeric> <numeric> <character> <logical>
    ##    ENSG00000133703_1            TGG  25209843  25209846    region_8     FALSE
    ##    ENSG00000133703_2            AGG  25209896  25209893    region_8     FALSE
    ##    ENSG00000133703_3            TGG  25215438  25215441    region_4     FALSE
    ##    ENSG00000133703_4            TGG  25215477  25215480    region_4     FALSE
    ##    ENSG00000133703_5            AGG  25215477  25215474    region_4     FALSE
    ##                  ...            ...       ...       ...         ...       ...
    ##   ENSG00000133703_52            TGG  25245349  25245352    region_1     FALSE
    ##   ENSG00000133703_53            TGG  25245352  25245355    region_1     FALSE
    ##   ENSG00000133703_54            TGG  25245358  25245361    region_1     FALSE
    ##   ENSG00000133703_55            TGG  25245365  25245368    region_1     FALSE
    ##   ENSG00000133703_56            AGG  25245392  25245389    region_1     FALSE
    ##                      percentGC     polyA     polyC     polyG     polyT
    ##                      <numeric> <logical> <logical> <logical> <logical>
    ##    ENSG00000133703_1        30      TRUE     FALSE     FALSE     FALSE
    ##    ENSG00000133703_2        30     FALSE     FALSE     FALSE     FALSE
    ##    ENSG00000133703_3        20     FALSE     FALSE     FALSE     FALSE
    ##    ENSG00000133703_4        40      TRUE     FALSE     FALSE     FALSE
    ##    ENSG00000133703_5        30     FALSE     FALSE     FALSE      TRUE
    ##                  ...       ...       ...       ...       ...       ...
    ##   ENSG00000133703_52        55     FALSE     FALSE     FALSE     FALSE
    ##   ENSG00000133703_53        45     FALSE     FALSE     FALSE     FALSE
    ##   ENSG00000133703_54        30     FALSE     FALSE     FALSE     FALSE
    ##   ENSG00000133703_55        25     FALSE     FALSE     FALSE     FALSE
    ##   ENSG00000133703_56        30     FALSE     FALSE     FALSE      TRUE
    ##                      startingGGGGG        n0      n0_c      n0_p        n1
    ##                          <logical> <numeric> <numeric> <numeric> <numeric>
    ##    ENSG00000133703_1         FALSE         1         1         0         4
    ##    ENSG00000133703_2         FALSE         1         1         0         1
    ##    ENSG00000133703_3         FALSE         1         1         0         0
    ##    ENSG00000133703_4         FALSE         1         1         0         0
    ##    ENSG00000133703_5         FALSE         1         1         0         0
    ##                  ...           ...       ...       ...       ...       ...
    ##   ENSG00000133703_52         FALSE         1         1         0         1
    ##   ENSG00000133703_53         FALSE         1         1         0         1
    ##   ENSG00000133703_54         FALSE         1         1         0         1
    ##   ENSG00000133703_55         FALSE         2         1         0         2
    ##   ENSG00000133703_56         FALSE         2         0         0         1
    ##                           n1_c      n1_p
    ##                      <numeric> <numeric>
    ##    ENSG00000133703_1         0         0
    ##    ENSG00000133703_2         0         0
    ##    ENSG00000133703_3         0         0
    ##    ENSG00000133703_4         0         0
    ##    ENSG00000133703_5         0         0
    ##                  ...       ...       ...
    ##   ENSG00000133703_52         0         0
    ##   ENSG00000133703_53         0         0
    ##   ENSG00000133703_54         0         0
    ##   ENSG00000133703_55         0         0
    ##   ENSG00000133703_56         0         0
    ##                                                                 alignments
    ##                                                              <GRangesList>
    ##    ENSG00000133703_1  chr12:25209843:-,chr8:68551391:-,chr6:54771089:+,...
    ##    ENSG00000133703_2                      chr12:25209896:+,chr6:54771050:-
    ##    ENSG00000133703_3                                      chr12:25215438:-
    ##    ENSG00000133703_4                                      chr12:25215477:-
    ##    ENSG00000133703_5                                      chr12:25215477:+
    ##                  ...                                                   ...
    ##   ENSG00000133703_52                      chr12:25245349:-,chr6:54770618:+
    ##   ENSG00000133703_53                      chr12:25245352:-,chr6:54770615:+
    ##   ENSG00000133703_54                      chr12:25245358:-,chr6:54770609:+
    ##   ENSG00000133703_55 chr12:25245365:-,chr6:54770602:+,chr13:60822020:-,...
    ##   ENSG00000133703_56     chr12:25245392:+,chr6:54770575:-,chr1:210618123:-
    ##                                                                          geneAnnotation
    ##                                                                    <SplitDataFrameList>
    ##    ENSG00000133703_1 chr12:25209846:-:...,chr12:25209846:-:...,chr12:25209846:-:...,...
    ##    ENSG00000133703_2 chr12:25209893:+:...,chr12:25209893:+:...,chr12:25209893:+:...,...
    ##    ENSG00000133703_3                                               chr12:25215441:-:...
    ##    ENSG00000133703_4                                               chr12:25215480:-:...
    ##    ENSG00000133703_5                                               chr12:25215474:+:...
    ##                  ...                                                                ...
    ##   ENSG00000133703_52 chr12:25245352:-:...,chr12:25245352:-:...,chr12:25245352:-:...,...
    ##   ENSG00000133703_53 chr12:25245355:-:...,chr12:25245355:-:...,chr12:25245355:-:...,...
    ##   ENSG00000133703_54 chr12:25245361:-:...,chr12:25245361:-:...,chr12:25245361:-:...,...
    ##   ENSG00000133703_55 chr12:25245368:-:...,chr12:25245368:-:...,chr12:25245368:-:...,...
    ##   ENSG00000133703_56 chr12:25245389:+:...,chr12:25245389:+:...,chr12:25245389:+:...,...
    ##                           enzymeAnnotation score_ruleset1 score_crisprscan
    ##                       <SplitDataFrameList>      <numeric>        <numeric>
    ##    ENSG00000133703_1 FALSE:FALSE:FALSE:...      0.0432227      0.513389669
    ##    ENSG00000133703_2 FALSE:FALSE:FALSE:...      0.0244329      0.235156658
    ##    ENSG00000133703_3 FALSE:FALSE:FALSE:...      0.0184826      0.110354708
    ##    ENSG00000133703_4 FALSE:FALSE:FALSE:...      0.0839624     -0.000986411
    ##    ENSG00000133703_5 FALSE:FALSE:FALSE:...      0.0913933      0.369056539
    ##                  ...                   ...            ...              ...
    ##   ENSG00000133703_52 FALSE:FALSE:FALSE:...      0.1958908         0.792829
    ##   ENSG00000133703_53 FALSE:FALSE:FALSE:...      0.0379117         0.234116
    ##   ENSG00000133703_54 FALSE:FALSE:FALSE:...      0.0281395         0.311893
    ##   ENSG00000133703_55 FALSE:FALSE:FALSE:...      0.1930296         0.325055
    ##   ENSG00000133703_56 FALSE:FALSE:FALSE:...      0.0737499         0.130166
    ##                      score_cfd score_mit    hasSNP                     snps
    ##                      <numeric> <numeric> <logical>     <SplitDataFrameList>
    ##    ENSG00000133703_1  0.425027  0.426600      TRUE rs1137282:25209843:0:...
    ##    ENSG00000133703_2  0.500000  0.577367     FALSE                 :...,...
    ##    ENSG00000133703_3  1.000000  1.000000     FALSE                 :...,...
    ##    ENSG00000133703_4  1.000000  1.000000     FALSE                 :...,...
    ##    ENSG00000133703_5  1.000000  1.000000     FALSE                 :...,...
    ##                  ...       ...       ...       ...                      ...
    ##   ENSG00000133703_52  0.500000  0.547046     FALSE                 :...,...
    ##   ENSG00000133703_53  0.500000  0.619963     FALSE                 :...,...
    ##   ENSG00000133703_54  0.777778  0.759301     FALSE                 :...,...
    ##   ENSG00000133703_55  0.458599  0.489579     FALSE                 :...,...
    ##   ENSG00000133703_56  0.442623  0.464868     FALSE                 :...,...
    ##   -------
    ##   seqinfo: 640 sequences (1 circular) from hg38 genome
    ##   crisprNuclease: SpCas9

### Converting the `GuideSet` object to a list of data.frames

The `flattenGuideSet` function in `crisprDesign` is a convenience
function to convert the complex `GuideSet` object into a set of 5
data.frames that can be saved as plain text files:

``` r
dfs <- flattenGuideSet(gs)
```

We can look at the names of the data.frames:

``` r
names(dfs)
```

    ## [1] "primary"          "alignments"       "geneAnnotation"   "enzymeAnnotation"
    ## [5] "snps"

As an example, let’s look at the first rows of the primary data.frame:

``` r
head(dfs$primary)
```

    ##                  ID               spacer   chr    start      end strand pam
    ## 1 ENSG00000133703_1 AAAGAAAAGATGAGCAAAGA chr12 25209844 25209863      - TGG
    ## 2 ENSG00000133703_2 TTCTCGAACTAATGTATAGA chr12 25209876 25209895      + AGG
    ## 3 ENSG00000133703_3 AAATGCATTATAATGTAATC chr12 25215439 25215458      - TGG
    ## 4 ENSG00000133703_4 AGCAAAGAAGAAAAGACTCC chr12 25215478 25215497      - TGG
    ## 5 ENSG00000133703_5 TTTTTAATTTTCACACAGCC chr12 25215457 25215476      + AGG
    ## 6 ENSG00000133703_6 TTTTTTTCAATCTGTATTGT chr12 25215500 25215519      + CGG
    ##   pam_site cut_site   region inRepeats percentGC polyA polyC polyG polyT
    ## 1 25209843 25209846 region_8     FALSE        30  TRUE FALSE FALSE FALSE
    ## 2 25209896 25209893 region_8     FALSE        30 FALSE FALSE FALSE FALSE
    ## 3 25215438 25215441 region_4     FALSE        20 FALSE FALSE FALSE FALSE
    ## 4 25215477 25215480 region_4     FALSE        40  TRUE FALSE FALSE FALSE
    ## 5 25215477 25215474 region_4     FALSE        30 FALSE FALSE FALSE  TRUE
    ## 6 25215520 25215517 region_4     FALSE        20 FALSE FALSE FALSE  TRUE
    ##   startingGGGGG n0 n0_c n0_p n1 n1_c n1_p score_ruleset1 score_crisprscan
    ## 1         FALSE  1    1    0  4    0    0     0.04322268     0.5133896687
    ## 2         FALSE  1    1    0  1    0    0     0.02443294     0.2351566583
    ## 3         FALSE  1    1    0  0    0    0     0.01848258     0.1103547080
    ## 4         FALSE  1    1    0  0    0    0     0.08396242    -0.0009864112
    ## 5         FALSE  1    1    0  0    0    0     0.09139327     0.3690565395
    ## 6         FALSE  1    1    0  4    0    0     0.06525962     0.2724537448
    ##   score_cfd score_mit hasSNP
    ## 1 0.4250273 0.4266001   TRUE
    ## 2 0.5000000 0.5773672  FALSE
    ## 3 1.0000000 1.0000000  FALSE
    ## 4 1.0000000 1.0000000  FALSE
    ## 5 1.0000000 1.0000000  FALSE
    ## 6 0.5212645 0.8835838  FALSE

# Building a complete annotation for all genes

We first get all possibles genes from our gene model:

``` r
gene_ids <- unique(txObject$cds$gene_id)
head(gene_ids)
```

    ## [1] "ENSG00000186092" "ENSG00000187634" "ENSG00000187961" "ENSG00000187583"
    ## [5] "ENSG00000187608" "ENSG00000188157"

and build an index to loop over for generating the GuideSet objects:

``` r
gene_index <- seq_along(gene_ids)
head(gene_index)
```

    ## [1] 1 2 3 4 5 6

We specify where to save the GuideSets:

``` r
dir <- "./crisprko_cas9_hg38"
if (!dir.exists(dir)){
  dir.create(dir, recursive=TRUE)
}
```

We are now ready to generate and save all data with the function
`designCompleteAnnotation` from `crisprDesign`. The function was
designed to be as comprehensive as possible to design and annotate gRNAs
in one step. It does the following:

-   Extract the DNA/RNA sequences with `queryTss`/`queryTxDB`
-   Design gRNAs with `findSpacers`
-   Remove gRNAs targeting repeat elements with `removeRepeats`
-   Characterize spacer sequences with `addSequenceFeatures`
-   Find on- and off-targets with `addSpacerAlignmentsIterative`
-   Add gene annotation with `addGeneAnnotation`
-   Add TSS annotation with `addTssAnnotation`
-   Add on-target efficiency scores with `addOnTargetScores`
-   Add off-target specificity scores with `addOffTargetScores`
-   Add SNP annotation with `addSNPAnnotation`
-   Add restriction enzymes information with `addRestrictionEnzymes`

We are now looping over all genes to generate the data:

``` r
lapply(gene_index, function(gene){
    gs <- designCompleteAnnotation(queryValue=gene,
                                   queryColumn="gene_id",
                                   modality="CRISPRko",
                                   bsgenome=bsgenome,
                                   bowtie_index=bowtie_index,
                                   crisprNuclease=SpCas9,
                                   txObject=txObject,
                                   tssObject=tssObject,
                                   grRepeats=grRepeats,
                                   vcf=vcf,
                                   n_mismatches=3,
                                   scoring_methods=scoring_methods)
    write.rds(gs, file=file.path(dir, paste0(gene, ".rds")))
})
```

This loop can be modified by the user to use an embarrassingly-parallel
approach to save time using the `BiocParallel` package, for instance.

# Considerations for CRISPRa/CRISPRi

# Considerations for CRISPRkd (RfxCas13d)

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
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] BSgenome.Hsapiens.UCSC.hg38_1.4.4 BSgenome_1.63.5                  
    ##  [3] rtracklayer_1.55.4                Biostrings_2.63.2                
    ##  [5] XVector_0.35.0                    GenomicRanges_1.47.6             
    ##  [7] GenomeInfoDb_1.31.6               IRanges_2.29.1                   
    ##  [9] S4Vectors_0.33.11                 crisprDesignData_0.99.8          
    ## [11] crisprDesign_0.99.102             crisprScore_1.1.9                
    ## [13] crisprScoreData_1.1.3             ExperimentHub_2.3.5              
    ## [15] AnnotationHub_3.3.9               BiocFileCache_2.3.4              
    ## [17] dbplyr_2.1.1                      BiocGenerics_0.41.2              
    ## [19] crisprBase_1.1.2                 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] rjson_0.2.21                  ellipsis_0.3.2               
    ##  [3] Rbowtie_1.35.0                rstudioapi_0.13              
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
    ## [37] crisprBwa_1.1.2               xfun_0.30                    
    ## [39] stringr_1.4.0                 mime_0.12                    
    ## [41] lifecycle_1.0.1               restfulr_0.0.13              
    ## [43] XML_3.99-0.9                  zlibbioc_1.41.0              
    ## [45] basilisk.utils_1.9.1          vroom_1.5.7                  
    ## [47] VariantAnnotation_1.41.3      hms_1.1.1                    
    ## [49] promises_1.2.0.1              MatrixGenerics_1.7.0         
    ## [51] parallel_4.2.0                SummarizedExperiment_1.25.3  
    ## [53] yaml_2.3.5                    curl_4.3.2                   
    ## [55] memoise_2.0.1                 reticulate_1.24              
    ## [57] biomaRt_2.51.3                stringi_1.7.6                
    ## [59] RSQLite_2.2.12                BiocVersion_3.15.0           
    ## [61] BiocIO_1.5.0                  randomForest_4.7-1           
    ## [63] GenomicFeatures_1.47.13       filelock_1.0.2               
    ## [65] BiocParallel_1.29.18          rlang_1.0.2                  
    ## [67] pkgconfig_2.0.3               matrixStats_0.61.0           
    ## [69] bitops_1.0-7                  evaluate_0.15                
    ## [71] lattice_0.20-45               purrr_0.3.4                  
    ## [73] GenomicAlignments_1.31.2      bit_4.0.4                    
    ## [75] tidyselect_1.1.2              magrittr_2.0.2               
    ## [77] R6_2.5.1                      generics_0.1.2               
    ## [79] DelayedArray_0.21.2           DBI_1.1.2                    
    ## [81] pillar_1.7.0                  KEGGREST_1.35.0              
    ## [83] RCurl_1.98-1.6                tibble_3.1.6                 
    ## [85] dir.expiry_1.3.0              crayon_1.5.0                 
    ## [87] utf8_1.2.2                    tzdb_0.2.0                   
    ## [89] rmarkdown_2.13                progress_1.2.2               
    ## [91] grid_4.2.0                    blob_1.2.2                   
    ## [93] digest_0.6.29                 xtable_1.8-4                 
    ## [95] httpuv_1.6.5                  Rbwa_1.1.0
