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
-   [Building a complete annotation for a given
    gene](#building-a-complete-annotation-for-a-given-gene)
    -   [Flattening](#flattening)
-   [Reproducibility](#reproducibility)

Authors: Jean-Philippe Fortin

Date: July 16, 2022

# Introduction

# Loading necessary packages

``` r
library(crisprDesign)
library(crisprDesignData)
library(crisprBase)
library(BSgenome.Hsapiens.UCSC.hg38)
```

### Specifying the genome

``` r
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
```

### Specifying the genome index

``` r
bowtie_index <- "/Users/fortinj2/crisprIndices/bowtie/hg38/hg38"
```

### Specifying a SNP VCF file

``` r
vcf <- "/Users/fortinj2/crisprIndices/snps/dbsnp151.grch38/00-common_all.vcf.gz"
```

### Specifying the nuclease

``` r
data(SpCas9, package="crisprBase")
crisprNuclease=SpCas9
```

### Specifying on-target scoring methods

``` r
scoring_methods=c("crisprscan", "ruleset1")
```

### Specifying gene models and TSS annotations

The package `crisprDesignData` contains already-processed Ensembl
objects for gene annotation of human and mouse gRNAs that are needed:

``` r
data(txdb_human, package="crisprDesignData")
data(tss_human, package="crisprDesignData")
data(gr.repeats.hg38, package="crisprDesignData")
txObject <- txdb_human
tssObject <- tss_human
grRepeats <- gr.repeats.hg38
```

# Building a complete annotation for a given gene

``` r
gs <- precomputeGuides("ENSG00000133703",
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

    ## [precomputeGuides] Adding sequence statistics 
    ## [precomputeGuides] Adding spacer alignments

    ## Loading required namespace: crisprBwa

    ## [runCrisprBowtie] Using BSgenome.Hsapiens.UCSC.hg38 
    ## [runCrisprBowtie] Searching for SpCas9 protospacers 
    ## [runCrisprBowtie] Using BSgenome.Hsapiens.UCSC.hg38 
    ## [runCrisprBowtie] Searching for SpCas9 protospacers 
    ## [precomputeGuides] Adding gene annotation 
    ## [precomputeGuides] Adding on-target scores

    ## [addOnTargetScores] Adding ruleset1 scores.

    ## [addOnTargetScores] Adding crisprscan scores.

    ## [precomputeGuides] Adding CFD scores annotation 
    ## [precomputeGuides] Adding SNP annotation

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

### Flattening

``` r
tables <- flattenGuideSet(gs)
head(tables)
```

    ## $primary
    ##                    ID   chr    start      end strand          protospacer pam
    ## 1   ENSG00000133703_1 chr12 25209843 25209843      - AAAGAAAAGATGAGCAAAGA TGG
    ## 2   ENSG00000133703_2 chr12 25209896 25209896      + TTCTCGAACTAATGTATAGA AGG
    ## 3   ENSG00000133703_3 chr12 25215438 25215438      - AAATGCATTATAATGTAATC TGG
    ## 4   ENSG00000133703_4 chr12 25215477 25215477      - AGCAAAGAAGAAAAGACTCC TGG
    ## 5   ENSG00000133703_5 chr12 25215477 25215477      + TTTTTAATTTTCACACAGCC AGG
    ## 6   ENSG00000133703_6 chr12 25215520 25215520      + TTTTTTTCAATCTGTATTGT CGG
    ## 7   ENSG00000133703_7 chr12 25215535 25215535      - GGAGGATGCTTTTTATACAT TGG
    ## 8   ENSG00000133703_8 chr12 25215553 25215553      - TTTTACAATGCAGAGAGTGG AGG
    ## 9   ENSG00000133703_9 chr12 25215556 25215556      - GTGTTTTACAATGCAGAGAG TGG
    ## 10 ENSG00000133703_10 chr12 25225615 25225615      - AACATCAGCAAAGACAAGAC AGG
    ## 11 ENSG00000133703_11 chr12 25225644 25225644      + TTTGCTGATGTTTCAATAAA AGG
    ## 12 ENSG00000133703_12 chr12 25225653 25225653      - CAGGACTTAGCAAGAAGTTA TGG
    ## 13 ENSG00000133703_13 chr12 25225672 25225672      - AGTAGACACAAAACAGGCTC AGG
    ## 14 ENSG00000133703_14 chr12 25225678 25225678      - TAGAACAGTAGACACAAAAC AGG
    ## 15 ENSG00000133703_15 chr12 25225701 25225701      + TTTGTGTCTACTGTTCTAGA AGG
    ## 16 ENSG00000133703_16 chr12 25225722 25225722      - GATGTACCTATGGTCCTAGT AGG
    ## 17 ENSG00000133703_17 chr12 25225726 25225726      + AATCACATTTATTTCCTACT AGG
    ## 18 ENSG00000133703_18 chr12 25225732 25225732      - GGACTCTGAAGATGTACCTA TGG
    ## 19 ENSG00000133703_19 chr12 25225734 25225734      + TTATTTCCTACTAGGACCAT AGG
    ## 20 ENSG00000133703_20 chr12 25225753 25225753      - AGAACAAATTAAAAGAGTTA AGG
    ## 21 ENSG00000133703_21 chr12 25225775 25225775      + AACTCTTTTAATTTGTTCTC TGG
    ## 22 ENSG00000133703_22 chr12 25225776 25225776      + ACTCTTTTAATTTGTTCTCT GGG
    ## 23 ENSG00000133703_23 chr12 25227231 25227231      - AGATATTCACCATTATAGGT GGG
    ## 24 ENSG00000133703_24 chr12 25227232 25227232      - AAGATATTCACCATTATAGG TGG
    ## 25 ENSG00000133703_25 chr12 25227235 25227235      - TTGAAGATATTCACCATTAT AGG
    ## 26 ENSG00000133703_26 chr12 25227240 25227240      + CAATTTAAACCCACCTATAA TGG
    ## 27 ENSG00000133703_27 chr12 25227274 25227274      + AAATGATTTAGTATTATTTA TGG
    ## 28 ENSG00000133703_28 chr12 25227296 25227296      - CAGTACATGAGGACTGGGGA GGG
    ## 29 ENSG00000133703_29 chr12 25227297 25227297      - CCAGTACATGAGGACTGGGG AGG
    ## 30 ENSG00000133703_30 chr12 25227300 25227300      - GGACCAGTACATGAGGACTG GGG
    ## 31 ENSG00000133703_31 chr12 25227301 25227301      - GGGACCAGTACATGAGGACT GGG
    ## 32 ENSG00000133703_32 chr12 25227302 25227302      - AGGGACCAGTACATGAGGAC TGG
    ## 33 ENSG00000133703_33 chr12 25227307 25227307      - CAATGAGGGACCAGTACATG AGG
    ## 34 ENSG00000133703_34 chr12 25227315 25227315      + CCTCCCCAGTCCTCATGTAC TGG
    ## 35 ENSG00000133703_35 chr12 25227321 25227321      - AGAGGAGTACAGTGCAATGA GGG
    ## 36 ENSG00000133703_36 chr12 25227322 25227322      - AAGAGGAGTACAGTGCAATG AGG
    ## 37 ENSG00000133703_37 chr12 25227339 25227339      - TCTCGACACAGCAGGTCAAG AGG
    ## 38 ENSG00000133703_38 chr12 25227347 25227347      - TTGGATATTCTCGACACAGC AGG
    ## 39 ENSG00000133703_39 chr12 25227366 25227366      - TGATGGAGAAACCTGTCTCT TGG
    ## 40 ENSG00000133703_40 chr12 25227373 25227373      + GTCGAGAATATCCAAGAGAC AGG
    ## 41 ENSG00000133703_41 chr12 25227383 25227383      - AGGAAGCAAGTAGTAATTGA TGG
    ## 42 ENSG00000133703_42 chr12 25227403 25227403      - TCCCTTCTCAGGATTCCTAC AGG
    ## 43 ENSG00000133703_43 chr12 25227406 25227406      + AATTACTACTTGCTTCCTGT AGG
    ## 44 ENSG00000133703_44 chr12 25227419 25227419      + TTCCTGTAGGAATCCTGAGA AGG
    ## 45 ENSG00000133703_45 chr12 25227420 25227420      + TCCTGTAGGAATCCTGAGAA GGG
    ## 46 ENSG00000133703_46 chr12 25235231 25235231      + TACTGCTTAATAACACCTGT AGG
    ## 47 ENSG00000133703_47 chr12 25245275 25245275      - CGAATATGATCCAACAATAG AGG
    ## 48 ENSG00000133703_48 chr12 25245283 25245283      + ACAAGATTTACCTCTATTGT TGG
    ## 49 ENSG00000133703_49 chr12 25245299 25245299      - GCTAATTCAGAATCATTTTG TGG
    ## 50 ENSG00000133703_50 chr12 25245330 25245330      + CTGAATTAGCTGTATCGTCA AGG
    ## 51 ENSG00000133703_51 chr12 25245343 25245343      - GTAGTTGGAGCTGGTGGCGT AGG
    ## 52 ENSG00000133703_52 chr12 25245349 25245349      - CTTGTGGTAGTTGGAGCTGG TGG
    ## 53 ENSG00000133703_53 chr12 25245352 25245352      - AAACTTGTGGTAGTTGGAGC TGG
    ## 54 ENSG00000133703_54 chr12 25245358 25245358      - GAATATAAACTTGTGGTAGT TGG
    ## 55 ENSG00000133703_55 chr12 25245365 25245365      - AATGACTGAATATAAACTTG TGG
    ## 56 ENSG00000133703_56 chr12 25245392 25245392      + TATATTCAGTCATTTTCAGC AGG
    ##    pam_site cut_site    region inRepeats percentGC polyA polyC polyG polyT
    ## 1  25209843 25209846  region_8     FALSE        30  TRUE FALSE FALSE FALSE
    ## 2  25209896 25209893  region_8     FALSE        30 FALSE FALSE FALSE FALSE
    ## 3  25215438 25215441  region_4     FALSE        20 FALSE FALSE FALSE FALSE
    ## 4  25215477 25215480  region_4     FALSE        40  TRUE FALSE FALSE FALSE
    ## 5  25215477 25215474  region_4     FALSE        30 FALSE FALSE FALSE  TRUE
    ## 6  25215520 25215517  region_4     FALSE        20 FALSE FALSE FALSE  TRUE
    ## 7  25215535 25215538  region_4     FALSE        35 FALSE FALSE FALSE  TRUE
    ## 8  25215553 25215556  region_4     FALSE        40 FALSE FALSE FALSE  TRUE
    ## 9  25215556 25215559  region_4     FALSE        40 FALSE FALSE FALSE  TRUE
    ## 10 25225615 25225618  region_3     FALSE        40 FALSE FALSE FALSE FALSE
    ## 11 25225644 25225641  region_3     FALSE        25 FALSE FALSE FALSE FALSE
    ## 12 25225653 25225656  region_3     FALSE        40 FALSE FALSE FALSE FALSE
    ## 13 25225672 25225675  region_3     FALSE        45  TRUE FALSE FALSE FALSE
    ## 14 25225678 25225681  region_3     FALSE        35  TRUE FALSE FALSE FALSE
    ## 15 25225701 25225698  region_3     FALSE        35 FALSE FALSE FALSE FALSE
    ## 16 25225722 25225725  region_3     FALSE        45 FALSE FALSE FALSE FALSE
    ## 17 25225726 25225723  region_3     FALSE        25 FALSE FALSE FALSE FALSE
    ## 18 25225732 25225735  region_3     FALSE        45 FALSE FALSE FALSE FALSE
    ## 19 25225734 25225731  region_3     FALSE        35 FALSE FALSE FALSE FALSE
    ## 20 25225753 25225756  region_3     FALSE        20  TRUE FALSE FALSE FALSE
    ## 21 25225775 25225772  region_3     FALSE        25 FALSE FALSE FALSE  TRUE
    ## 22 25225776 25225773  region_3     FALSE        25 FALSE FALSE FALSE  TRUE
    ## 23 25227231 25227234  region_2     FALSE        30 FALSE FALSE FALSE FALSE
    ## 24 25227232 25227235  region_2     FALSE        30 FALSE FALSE FALSE FALSE
    ## 25 25227235 25227238  region_2     FALSE        25 FALSE FALSE FALSE FALSE
    ## 26 25227240 25227237  region_2     FALSE        30 FALSE FALSE FALSE FALSE
    ## 27 25227274 25227271  region_2     FALSE        10 FALSE FALSE FALSE FALSE
    ## 28 25227296 25227299  region_2     FALSE        55 FALSE FALSE  TRUE FALSE
    ## 29 25227297 25227300  region_2     FALSE        60 FALSE FALSE  TRUE FALSE
    ## 30 25227300 25227303  region_2     FALSE        55 FALSE FALSE FALSE FALSE
    ## 31 25227301 25227304  region_2     FALSE        55 FALSE FALSE FALSE FALSE
    ## 32 25227302 25227305  region_2     FALSE        55 FALSE FALSE FALSE FALSE
    ## 33 25227307 25227310  region_2     FALSE        50 FALSE FALSE FALSE FALSE
    ## 34 25227315 25227312  region_2     FALSE        60 FALSE  TRUE FALSE FALSE
    ## 35 25227321 25227324  region_2     FALSE        45 FALSE FALSE FALSE FALSE
    ## 36 25227322 25227325  region_2     FALSE        45 FALSE FALSE FALSE FALSE
    ## 37 25227339 25227342  region_2     FALSE        55 FALSE FALSE FALSE FALSE
    ## 38 25227347 25227350  region_2     FALSE        45 FALSE FALSE FALSE FALSE
    ## 39 25227366 25227369  region_2     FALSE        45 FALSE FALSE FALSE FALSE
    ## 40 25227373 25227370  region_2     FALSE        45 FALSE FALSE FALSE FALSE
    ## 41 25227383 25227386  region_2     FALSE        35 FALSE FALSE FALSE FALSE
    ## 42 25227403 25227406  region_2     FALSE        50 FALSE FALSE FALSE FALSE
    ## 43 25227406 25227403  region_2     FALSE        35 FALSE FALSE FALSE FALSE
    ## 44 25227419 25227416  region_2     FALSE        45 FALSE FALSE FALSE FALSE
    ## 45 25227420 25227417  region_2     FALSE        45 FALSE FALSE FALSE FALSE
    ## 46 25235231 25235228 region_12     FALSE        35 FALSE FALSE FALSE FALSE
    ## 47 25245275 25245278  region_1     FALSE        35 FALSE FALSE FALSE FALSE
    ## 48 25245283 25245280  region_1     FALSE        30 FALSE FALSE FALSE FALSE
    ## 49 25245299 25245302  region_1     FALSE        30 FALSE FALSE FALSE  TRUE
    ## 50 25245330 25245327  region_1     FALSE        40 FALSE FALSE FALSE FALSE
    ## 51 25245343 25245346  region_1     FALSE        60 FALSE FALSE FALSE FALSE
    ## 52 25245349 25245352  region_1     FALSE        55 FALSE FALSE FALSE FALSE
    ## 53 25245352 25245355  region_1     FALSE        45 FALSE FALSE FALSE FALSE
    ## 54 25245358 25245361  region_1     FALSE        30 FALSE FALSE FALSE FALSE
    ## 55 25245365 25245368  region_1     FALSE        25 FALSE FALSE FALSE FALSE
    ## 56 25245392 25245389  region_1     FALSE        30 FALSE FALSE FALSE  TRUE
    ##    startingGGGGG n0 n0_c n0_p n1 n1_c n1_p score_ruleset1 score_crisprscan
    ## 1          FALSE  1    1    0  4    0    0    0.043222682     0.5133896687
    ## 2          FALSE  1    1    0  1    0    0    0.024432939     0.2351566583
    ## 3          FALSE  1    1    0  0    0    0    0.018482579     0.1103547080
    ## 4          FALSE  1    1    0  0    0    0    0.083962416    -0.0009864112
    ## 5          FALSE  1    1    0  0    0    0    0.091393271     0.3690565395
    ## 6          FALSE  1    1    0  4    0    0    0.065259621     0.2724537448
    ## 7          FALSE  1    1    0  0    0    0    0.533938710     0.3539000935
    ## 8          FALSE  1    1    0  0    0    0    0.219489163     0.5827207543
    ## 9          FALSE  1    1    0  0    0    0    0.062261532     0.5705158651
    ## 10         FALSE  1    1    0  0    0    0    0.163449646     0.4215891467
    ## 11         FALSE  1    1    0  0    0    0    0.042497973     0.4581810540
    ## 12         FALSE  1    1    0  1    0    0    0.086608395     0.4398528428
    ## 13         FALSE  1    1    0  0    0    0    0.098959664     0.3545034324
    ## 14         FALSE  1    1    0  1    0    0    0.150109449     0.3636388955
    ## 15         FALSE  2    1    0  0    0    0    0.277445350     0.4138651398
    ## 16         FALSE  1    1    0  2    1    0    0.222075141     0.6632815210
    ## 17         FALSE  2    1    0  0    0    0    0.045234486     0.2679694809
    ## 18         FALSE  1    1    0  1    0    0    0.111680307     0.3051037540
    ## 19         FALSE  2    1    0  0    0    0    0.327827733     0.5448407158
    ## 20         FALSE  1    1    0  0    0    0    0.053409998     0.1453487123
    ## 21         FALSE  1    1    0  1    0    0    0.007536694    -0.0442802663
    ## 22         FALSE  1    1    0  0    0    0    0.020817600     0.1456761556
    ## 23         FALSE  1    1    0  0    0    0    0.512341024     0.3504188624
    ## 24         FALSE  1    1    0  0    0    0    0.258675763     0.3923183291
    ## 25         FALSE  1    1    0  0    0    0    0.072213986     0.4102132935
    ## 26         FALSE  1    1    0  0    0    0    0.033552012     0.2022816416
    ## 27         FALSE  1    1    0  1    0    0    0.009571396     0.2047197767
    ## 28         FALSE  1    1    0  0    0    0    0.223314082     0.5848892426
    ## 29         FALSE  1    1    0  1    0    0    0.123079697     0.6300169444
    ## 30         FALSE  1    1    0  0    0    0    0.277529315     0.6121450231
    ## 31         FALSE  1    1    0  1    0    0    0.060000789     0.5150873975
    ## 32         FALSE  1    1    0  0    0    0    0.036750591     0.3638034597
    ## 33         FALSE  1    1    0  0    0    0    0.011575258     0.4524457160
    ## 34         FALSE  1    1    0  0    0    0    0.077578224     0.2317496922
    ## 35         FALSE  1    1    0  0    0    0    0.339797580     0.4874341560
    ## 36         FALSE  1    1    0  0    0    0    0.089310331     0.4758715650
    ## 37         FALSE  1    1    0  0    0    0    0.291064878     0.5096967766
    ## 38         FALSE  1    1    0  0    0    0    0.050474158     0.4048914907
    ## 39         FALSE  2    1    0  0    0    0    0.052816238     0.1865341473
    ## 40         FALSE  1    1    0  1    0    0    0.158064168     0.4437388578
    ## 41         FALSE  2    1    0  0    0    0    0.115418540     0.3717978184
    ## 42         FALSE  1    1    0  0    0    0    0.592452224     0.3388411782
    ## 43         FALSE  2    1    0  0    0    0    0.371080957     0.1340496768
    ## 44         FALSE  1    0    0  0    0    0    0.052443950     0.4399618205
    ## 45         FALSE  1    0    0  0    0    0    0.101322961     0.5668301886
    ## 46         FALSE  1    0    0  0    0    0    0.518885298     0.1249291062
    ## 47         FALSE  1    1    0  1    0    0    0.231280485     0.2836962755
    ## 48         FALSE  1    1    0  0    0    0    0.048728261     0.1828534831
    ## 49         FALSE  2    1    0  0    0    0    0.528400392     0.1566492632
    ## 50         FALSE  2    1    0  0    0    0    0.080909495     0.3837528502
    ## 51         FALSE  2    1    0  0    0    0    0.134949031     0.6763424113
    ## 52         FALSE  1    1    0  1    0    0    0.195890811     0.7928287023
    ## 53         FALSE  1    1    0  1    0    0    0.037911704     0.2341157806
    ## 54         FALSE  1    1    0  1    0    0    0.028139468     0.3118932129
    ## 55         FALSE  2    1    0  2    0    0    0.193029595     0.3250546585
    ## 56         FALSE  2    0    0  1    0    0    0.073749862     0.1301660348
    ##    score_cfd score_mit hasSNP
    ## 1  0.4250273 0.4266001   TRUE
    ## 2  0.5000000 0.5773672  FALSE
    ## 3  1.0000000 1.0000000  FALSE
    ## 4  1.0000000 1.0000000  FALSE
    ## 5  1.0000000 1.0000000  FALSE
    ## 6  0.5212645 0.8835838  FALSE
    ## 7  1.0000000 1.0000000  FALSE
    ## 8  1.0000000 1.0000000  FALSE
    ## 9  1.0000000 1.0000000  FALSE
    ## 10 1.0000000 1.0000000  FALSE
    ## 11 1.0000000 1.0000000  FALSE
    ## 12 0.6060606 0.7168459  FALSE
    ## 13 1.0000000 1.0000000  FALSE
    ## 14 0.9455253 0.9071057  FALSE
    ## 15 0.5000000 0.5000000  FALSE
    ## 16 0.5668016 0.4701457  FALSE
    ## 17 0.5000000 0.5000000  FALSE
    ## 18 0.8387097 1.0000000  FALSE
    ## 19 0.5000000 0.5000000  FALSE
    ## 20 1.0000000 1.0000000  FALSE
    ## 21 0.9422633 0.9700236  FALSE
    ## 22 1.0000000 1.0000000  FALSE
    ## 23 1.0000000 1.0000000  FALSE
    ## 24 1.0000000 1.0000000  FALSE
    ## 25 1.0000000 1.0000000  FALSE
    ## 26 1.0000000 1.0000000  FALSE
    ## 27 0.5000000 0.9267841  FALSE
    ## 28 1.0000000 1.0000000  FALSE
    ## 29 0.7000000 1.0000000  FALSE
    ## 30 1.0000000 1.0000000  FALSE
    ## 31 0.5263158 0.6317119  FALSE
    ## 32 1.0000000 1.0000000  FALSE
    ## 33 1.0000000 1.0000000  FALSE
    ## 34 1.0000000 1.0000000  FALSE
    ## 35 1.0000000 1.0000000  FALSE
    ## 36 1.0000000 1.0000000  FALSE
    ## 37 1.0000000 1.0000000  FALSE
    ## 38 1.0000000 1.0000000  FALSE
    ## 39 0.5000000 0.5000000  FALSE
    ## 40 0.5263158 0.6191950  FALSE
    ## 41 0.5000000 0.5000000  FALSE
    ## 42 1.0000000 1.0000000  FALSE
    ## 43 0.5000000 0.5000000  FALSE
    ## 44 1.0000000 1.0000000  FALSE
    ## 45 1.0000000 1.0000000  FALSE
    ## 46 1.0000000 1.0000000  FALSE
    ## 47 0.9541284 0.9545907  FALSE
    ## 48 1.0000000 1.0000000  FALSE
    ## 49 0.5000000 0.5000000  FALSE
    ## 50 0.9350649 0.9350649  FALSE
    ## 51 0.7941176 0.7941176  FALSE
    ## 52 0.5000000 0.5470460  FALSE
    ## 53 0.5000000 0.6199628  FALSE
    ## 54 0.7777778 0.7593014  FALSE
    ## 55 0.4585987 0.4895794  FALSE
    ## 56 0.4426230 0.4648680  FALSE
    ## 
    ## $alignments
    ##                    ID   chr     start       end strand               spacer
    ## 1   ENSG00000133703_1 chr12  25209843  25209843      - AAAGAAAAGATGAGCAAAGA
    ## 2   ENSG00000133703_1  chr8  68551391  68551391      - AAAGAAAAGATGAGCAAAGA
    ## 3   ENSG00000133703_1  chr6  54771089  54771089      + AAAGAAAAGATGAGCAAAGA
    ## 4   ENSG00000133703_1  chr5   4348033   4348033      + AAAGAAAAGATGAGCAAAGA
    ## 5   ENSG00000133703_1  chr1  48362810  48362810      + AAAGAAAAGATGAGCAAAGA
    ## 6   ENSG00000133703_2 chr12  25209896  25209896      + TTCTCGAACTAATGTATAGA
    ## 7   ENSG00000133703_2  chr6  54771050  54771050      - TTCTCGAACTAATGTATAGA
    ## 8   ENSG00000133703_3 chr12  25215438  25215438      - AAATGCATTATAATGTAATC
    ## 9   ENSG00000133703_4 chr12  25215477  25215477      - AGCAAAGAAGAAAAGACTCC
    ## 10  ENSG00000133703_5 chr12  25215477  25215477      + TTTTTAATTTTCACACAGCC
    ## 11  ENSG00000133703_6 chr12  25215520  25215520      + TTTTTTTCAATCTGTATTGT
    ## 12  ENSG00000133703_6  chrX  34968543  34968543      - TTTTTTTCAATCTGTATTGT
    ## 13  ENSG00000133703_6  chr8  16614003  16614003      - TTTTTTTCAATCTGTATTGT
    ## 14  ENSG00000133703_6  chr2  34098139  34098139      + TTTTTTTCAATCTGTATTGT
    ## 15  ENSG00000133703_6 chr17  32419546  32419546      - TTTTTTTCAATCTGTATTGT
    ## 16  ENSG00000133703_7 chr12  25215535  25215535      - GGAGGATGCTTTTTATACAT
    ## 17  ENSG00000133703_8 chr12  25215553  25215553      - TTTTACAATGCAGAGAGTGG
    ## 18  ENSG00000133703_9 chr12  25215556  25215556      - GTGTTTTACAATGCAGAGAG
    ## 19 ENSG00000133703_10 chr12  25225615  25225615      - AACATCAGCAAAGACAAGAC
    ## 20 ENSG00000133703_11 chr12  25225644  25225644      + TTTGCTGATGTTTCAATAAA
    ## 21 ENSG00000133703_12 chr12  25225653  25225653      - CAGGACTTAGCAAGAAGTTA
    ## 22 ENSG00000133703_12  chr6  54771002  54771002      + CAGGACTTAGCAAGAAGTTA
    ## 23 ENSG00000133703_13 chr12  25225672  25225672      - AGTAGACACAAAACAGGCTC
    ## 24 ENSG00000133703_14 chr12  25225678  25225678      - TAGAACAGTAGACACAAAAC
    ## 25 ENSG00000133703_14  chr8 111613834 111613834      - TAGAACAGTAGACACAAAAC
    ## 26 ENSG00000133703_15 chr12  25225701  25225701      + TTTGTGTCTACTGTTCTAGA
    ## 27 ENSG00000133703_15  chr6  54770956  54770956      - TTTGTGTCTACTGTTCTAGA
    ## 28 ENSG00000133703_16 chr12  25225722  25225722      - GATGTACCTATGGTCCTAGT
    ## 29 ENSG00000133703_16  chr1 114709677 114709677      - GATGTACCTATGGTCCTAGT
    ## 30 ENSG00000133703_16  chr6  54770935  54770935      + GATGTACCTATGGTCCTAGT
    ## 31 ENSG00000133703_17 chr12  25225726  25225726      + AATCACATTTATTTCCTACT
    ## 32 ENSG00000133703_17  chr6  54770931  54770931      - AATCACATTTATTTCCTACT
    ## 33 ENSG00000133703_18 chr12  25225732  25225732      - GGACTCTGAAGATGTACCTA
    ## 34 ENSG00000133703_18  chr6  54770925  54770925      + GGACTCTGAAGATGTACCTA
    ## 35 ENSG00000133703_19 chr12  25225734  25225734      + TTATTTCCTACTAGGACCAT
    ## 36 ENSG00000133703_19  chr6  54770923  54770923      - TTATTTCCTACTAGGACCAT
    ## 37 ENSG00000133703_20 chr12  25225753  25225753      - AGAACAAATTAAAAGAGTTA
    ## 38 ENSG00000133703_21 chr12  25225775  25225775      + AACTCTTTTAATTTGTTCTC
    ## 39 ENSG00000133703_21 chr11 117609656 117609656      + AACTCTTTTAATTTGTTCTC
    ## 40 ENSG00000133703_22 chr12  25225776  25225776      + ACTCTTTTAATTTGTTCTCT
    ## 41 ENSG00000133703_23 chr12  25227231  25227231      - AGATATTCACCATTATAGGT
    ## 42 ENSG00000133703_24 chr12  25227232  25227232      - AAGATATTCACCATTATAGG
    ## 43 ENSG00000133703_25 chr12  25227235  25227235      - TTGAAGATATTCACCATTAT
    ## 44 ENSG00000133703_26 chr12  25227240  25227240      + CAATTTAAACCCACCTATAA
    ## 45 ENSG00000133703_27 chr12  25227274  25227274      + AAATGATTTAGTATTATTTA
    ## 46 ENSG00000133703_27  chr6  54770843  54770843      - AAATGATTTAGTATTATTTA
    ## 47 ENSG00000133703_28 chr12  25227296  25227296      - CAGTACATGAGGACTGGGGA
    ## 48 ENSG00000133703_29 chr12  25227297  25227297      - CCAGTACATGAGGACTGGGG
    ## 49 ENSG00000133703_29  chr6  54770808  54770808      + CCAGTACATGAGGACTGGGG
    ## 50 ENSG00000133703_30 chr12  25227300  25227300      - GGACCAGTACATGAGGACTG
    ## 51 ENSG00000133703_31 chr12  25227301  25227301      - GGGACCAGTACATGAGGACT
    ## 52 ENSG00000133703_31  chr6  54770804  54770804      + GGGACCAGTACATGAGGACT
    ## 53 ENSG00000133703_32 chr12  25227302  25227302      - AGGGACCAGTACATGAGGAC
    ## 54 ENSG00000133703_33 chr12  25227307  25227307      - CAATGAGGGACCAGTACATG
    ## 55 ENSG00000133703_34 chr12  25227315  25227315      + CCTCCCCAGTCCTCATGTAC
    ## 56 ENSG00000133703_35 chr12  25227321  25227321      - AGAGGAGTACAGTGCAATGA
    ## 57 ENSG00000133703_36 chr12  25227322  25227322      - AAGAGGAGTACAGTGCAATG
    ## 58 ENSG00000133703_37 chr12  25227339  25227339      - TCTCGACACAGCAGGTCAAG
    ## 59 ENSG00000133703_38 chr12  25227347  25227347      - TTGGATATTCTCGACACAGC
    ## 60 ENSG00000133703_39 chr12  25227366  25227366      - TGATGGAGAAACCTGTCTCT
    ## 61 ENSG00000133703_39  chr6  54770740  54770740      + TGATGGAGAAACCTGTCTCT
    ## 62 ENSG00000133703_40 chr12  25227373  25227373      + GTCGAGAATATCCAAGAGAC
    ## 63 ENSG00000133703_40  chr6  54770733  54770733      - GTCGAGAATATCCAAGAGAC
    ## 64 ENSG00000133703_41 chr12  25227383  25227383      - AGGAAGCAAGTAGTAATTGA
    ## 65 ENSG00000133703_41  chr6  54770723  54770723      + AGGAAGCAAGTAGTAATTGA
    ## 66 ENSG00000133703_42 chr12  25227403  25227403      - TCCCTTCTCAGGATTCCTAC
    ## 67 ENSG00000133703_43 chr12  25227406  25227406      + AATTACTACTTGCTTCCTGT
    ## 68 ENSG00000133703_43  chr6  54770700  54770700      - AATTACTACTTGCTTCCTGT
    ## 69 ENSG00000133703_44 chr12  25227419  25227419      + TTCCTGTAGGAATCCTGAGA
    ## 70 ENSG00000133703_45 chr12  25227420  25227420      + TCCTGTAGGAATCCTGAGAA
    ## 71 ENSG00000133703_46 chr12  25235231  25235231      + TACTGCTTAATAACACCTGT
    ## 72 ENSG00000133703_47 chr12  25245275  25245275      - CGAATATGATCCAACAATAG
    ## 73 ENSG00000133703_47  chr6  54770692  54770692      + CGAATATGATCCAACAATAG
    ## 74 ENSG00000133703_48 chr12  25245283  25245283      + ACAAGATTTACCTCTATTGT
    ## 75 ENSG00000133703_49 chr12  25245299  25245299      - GCTAATTCAGAATCATTTTG
    ## 76 ENSG00000133703_49  chr6  54770668  54770668      + GCTAATTCAGAATCATTTTG
    ## 77 ENSG00000133703_50 chr12  25245330  25245330      + CTGAATTAGCTGTATCGTCA
    ## 78 ENSG00000133703_50  chr6  54770637  54770637      - CTGAATTAGCTGTATCGTCA
    ## 79 ENSG00000133703_51 chr12  25245343  25245343      - GTAGTTGGAGCTGGTGGCGT
    ## 80 ENSG00000133703_51  chr6  54770624  54770624      + GTAGTTGGAGCTGGTGGCGT
    ## 81 ENSG00000133703_52 chr12  25245349  25245349      - CTTGTGGTAGTTGGAGCTGG
    ## 82 ENSG00000133703_52  chr6  54770618  54770618      + CTTGTGGTAGTTGGAGCTGG
    ## 83 ENSG00000133703_53 chr12  25245352  25245352      - AAACTTGTGGTAGTTGGAGC
    ## 84 ENSG00000133703_53  chr6  54770615  54770615      + AAACTTGTGGTAGTTGGAGC
    ## 85 ENSG00000133703_54 chr12  25245358  25245358      - GAATATAAACTTGTGGTAGT
    ## 86 ENSG00000133703_54  chr6  54770609  54770609      + GAATATAAACTTGTGGTAGT
    ## 87 ENSG00000133703_55 chr12  25245365  25245365      - AATGACTGAATATAAACTTG
    ## 88 ENSG00000133703_55  chr6  54770602  54770602      + AATGACTGAATATAAACTTG
    ## 89 ENSG00000133703_55 chr13  60822020  60822020      - AATGACTGAATATAAACTTG
    ## 90 ENSG00000133703_55 chr14  57999156  57999156      - AATGACTGAATATAAACTTG
    ## 91 ENSG00000133703_56 chr12  25245392  25245392      + TATATTCAGTCATTTTCAGC
    ## 92 ENSG00000133703_56  chr6  54770575  54770575      - TATATTCAGTCATTTTCAGC
    ## 93 ENSG00000133703_56  chr1 210618123 210618123      - TATATTCAGTCATTTTCAGC
    ##             protospacer pam  pam_site n_mismatches canonical  cut_site  cds
    ## 1  AAAGAAAAGATGAGCAAAGA TGG  25209843            0      TRUE  25209846 KRAS
    ## 2  AAAGAAAAGATGAGCAATGA AAG  68551391            1     FALSE  68551394 <NA>
    ## 3  AAATAAAAGATGAGCAAAGA TGG  54771089            1      TRUE  54771086 <NA>
    ## 4  AGAGAAAAGATGAGCAAAGA TGG   4348033            1      TRUE   4348030 <NA>
    ## 5  CAAGAAAAGATGAGCAAAGA TGA  48362810            1     FALSE  48362807 <NA>
    ## 6  TTCTCGAACTAATGTATAGA AGG  25209896            0      TRUE  25209893 KRAS
    ## 7  TTCTCAAACTAATGTATAGA AGG  54771050            1      TRUE  54771053 <NA>
    ## 8  AAATGCATTATAATGTAATC TGG  25215438            0      TRUE  25215441 KRAS
    ## 9  AGCAAAGAAGAAAAGACTCC TGG  25215477            0      TRUE  25215480 KRAS
    ## 10 TTTTTAATTTTCACACAGCC AGG  25215477            0      TRUE  25215474 KRAS
    ## 11 TTTTTTTCAATCTGTATTGT CGG  25215520            0      TRUE  25215517 KRAS
    ## 12 TTTTTTTCAATCAGTATTGT TGG  34968543            1      TRUE  34968546 <NA>
    ## 13 TTTTTTTCAATCTGTATTAT AGA  16614003            1     FALSE  16614006 <NA>
    ## 14 TTTTTTTCAATGTGTATTGT GAG  34098139            1     FALSE  34098136 <NA>
    ## 15 TTTTTTTCATTCTGTATTGT TGA  32419546            1     FALSE  32419549 <NA>
    ## 16 GGAGGATGCTTTTTATACAT TGG  25215535            0      TRUE  25215538 KRAS
    ## 17 TTTTACAATGCAGAGAGTGG AGG  25215553            0      TRUE  25215556 KRAS
    ## 18 GTGTTTTACAATGCAGAGAG TGG  25215556            0      TRUE  25215559 KRAS
    ## 19 AACATCAGCAAAGACAAGAC AGG  25225615            0      TRUE  25225618 KRAS
    ## 20 TTTGCTGATGTTTCAATAAA AGG  25225644            0      TRUE  25225641 KRAS
    ## 21 CAGGACTTAGCAAGAAGTTA TGG  25225653            0      TRUE  25225656 KRAS
    ## 22 CAGGACTTAGCAAGGAGTTA GGG  54771002            1      TRUE  54770999 <NA>
    ## 23 AGTAGACACAAAACAGGCTC AGG  25225672            0      TRUE  25225675 KRAS
    ## 24 TAGAACAGTAGACACAAAAC AGG  25225678            0      TRUE  25225681 KRAS
    ## 25 TAGAACAGTAGACAAAAAAC AAG 111613834            1     FALSE 111613837 <NA>
    ## 26 TTTGTGTCTACTGTTCTAGA AGG  25225701            0      TRUE  25225698 KRAS
    ## 27 TTTGTGTCTACTGTTCTAGA AGG  54770956            0      TRUE  54770959 <NA>
    ## 28 GATGTACCTATGGTCCTAGT AGG  25225722            0      TRUE  25225725 KRAS
    ## 29 GATGTACCTATGGTGCTAGT GGG 114709677            1      TRUE 114709680 NRAS
    ## 30 GATGTGCCTATGGTCCTAGT AGG  54770935            1      TRUE  54770932 <NA>
    ## 31 AATCACATTTATTTCCTACT AGG  25225726            0      TRUE  25225723 KRAS
    ## 32 AATCACATTTATTTCCTACT AGG  54770931            0      TRUE  54770934 <NA>
    ## 33 GGACTCTGAAGATGTACCTA TGG  25225732            0      TRUE  25225735 KRAS
    ## 34 GGACTCTGAAGATGTGCCTA TGG  54770925            1      TRUE  54770922 <NA>
    ## 35 TTATTTCCTACTAGGACCAT AGG  25225734            0      TRUE  25225731 KRAS
    ## 36 TTATTTCCTACTAGGACCAT AGG  54770923            0      TRUE  54770926 <NA>
    ## 37 AGAACAAATTAAAAGAGTTA AGG  25225753            0      TRUE  25225756 KRAS
    ## 38 AACTCTTTTAATTTGTTCTC TGG  25225775            0      TRUE  25225772 KRAS
    ## 39 AACTCTTTTTATTTGTTCTC CGA 117609656            1     FALSE 117609653 <NA>
    ## 40 ACTCTTTTAATTTGTTCTCT GGG  25225776            0      TRUE  25225773 KRAS
    ## 41 AGATATTCACCATTATAGGT GGG  25227231            0      TRUE  25227234 KRAS
    ## 42 AAGATATTCACCATTATAGG TGG  25227232            0      TRUE  25227235 KRAS
    ## 43 TTGAAGATATTCACCATTAT AGG  25227235            0      TRUE  25227238 KRAS
    ## 44 CAATTTAAACCCACCTATAA TGG  25227240            0      TRUE  25227237 KRAS
    ## 45 AAATGATTTAGTATTATTTA TGG  25227274            0      TRUE  25227271 KRAS
    ## 46 AAATGATTTAATATTATTTA TGG  54770843            1      TRUE  54770846 <NA>
    ## 47 CAGTACATGAGGACTGGGGA GGG  25227296            0      TRUE  25227299 KRAS
    ## 48 CCAGTACATGAGGACTGGGG AGG  25227297            0      TRUE  25227300 KRAS
    ## 49 CCAGTACATGAGGACTGGGC GGG  54770808            1      TRUE  54770805 <NA>
    ## 50 GGACCAGTACATGAGGACTG GGG  25227300            0      TRUE  25227303 KRAS
    ## 51 GGGACCAGTACATGAGGACT GGG  25227301            0      TRUE  25227304 KRAS
    ## 52 AGGACCAGTACATGAGGACT GGG  54770804            1      TRUE  54770801 <NA>
    ## 53 AGGGACCAGTACATGAGGAC TGG  25227302            0      TRUE  25227305 KRAS
    ## 54 CAATGAGGGACCAGTACATG AGG  25227307            0      TRUE  25227310 KRAS
    ## 55 CCTCCCCAGTCCTCATGTAC TGG  25227315            0      TRUE  25227312 KRAS
    ## 56 AGAGGAGTACAGTGCAATGA GGG  25227321            0      TRUE  25227324 KRAS
    ## 57 AAGAGGAGTACAGTGCAATG AGG  25227322            0      TRUE  25227325 KRAS
    ## 58 TCTCGACACAGCAGGTCAAG AGG  25227339            0      TRUE  25227342 KRAS
    ## 59 TTGGATATTCTCGACACAGC AGG  25227347            0      TRUE  25227350 KRAS
    ## 60 TGATGGAGAAACCTGTCTCT TGG  25227366            0      TRUE  25227369 KRAS
    ## 61 TGATGGAGAAACCTGTCTCT TGG  54770740            0      TRUE  54770737 <NA>
    ## 62 GTCGAGAATATCCAAGAGAC AGG  25227373            0      TRUE  25227370 KRAS
    ## 63 GTCAAGAATATCCAAGAGAC AGG  54770733            1      TRUE  54770736 <NA>
    ## 64 AGGAAGCAAGTAGTAATTGA TGG  25227383            0      TRUE  25227386 KRAS
    ## 65 AGGAAGCAAGTAGTAATTGA TGG  54770723            0      TRUE  54770720 <NA>
    ## 66 TCCCTTCTCAGGATTCCTAC AGG  25227403            0      TRUE  25227406 KRAS
    ## 67 AATTACTACTTGCTTCCTGT AGG  25227406            0      TRUE  25227403 KRAS
    ## 68 AATTACTACTTGCTTCCTGT AGG  54770700            0      TRUE  54770703 <NA>
    ## 69 TTCCTGTAGGAATCCTGAGA AGG  25227419            0      TRUE  25227416 <NA>
    ## 70 TCCTGTAGGAATCCTGAGAA GGG  25227420            0      TRUE  25227417 <NA>
    ## 71 TACTGCTTAATAACACCTGT AGG  25235231            0      TRUE  25235228 <NA>
    ## 72 CGAATATGATCCAACAATAG AGG  25245275            0      TRUE  25245278 KRAS
    ## 73 CCAATATGATCCAACAATAG AGA  54770692            1     FALSE  54770689 <NA>
    ## 74 ACAAGATTTACCTCTATTGT TGG  25245283            0      TRUE  25245280 KRAS
    ## 75 GCTAATTCAGAATCATTTTG TGG  25245299            0      TRUE  25245302 KRAS
    ## 76 GCTAATTCAGAATCATTTTG TGG  54770668            0      TRUE  54770665 <NA>
    ## 77 CTGAATTAGCTGTATCGTCA AGG  25245330            0      TRUE  25245327 KRAS
    ## 78 CTGAATTAGCTGTATCGTCA AGA  54770637            0     FALSE  54770640 <NA>
    ## 79 GTAGTTGGAGCTGGTGGCGT AGG  25245343            0      TRUE  25245346 KRAS
    ## 80 GTAGTTGGAGCTGGTGGCGT AAG  54770624            0     FALSE  54770621 <NA>
    ## 81 CTTGTGGTAGTTGGAGCTGG TGG  25245349            0      TRUE  25245352 KRAS
    ## 82 CTTGCGGTAGTTGGAGCTGG TGG  54770618            1      TRUE  54770615 <NA>
    ## 83 AAACTTGTGGTAGTTGGAGC TGG  25245352            0      TRUE  25245355 KRAS
    ## 84 AAACTTGCGGTAGTTGGAGC TGG  54770615            1      TRUE  54770612 <NA>
    ## 85 GAATATAAACTTGTGGTAGT TGG  25245358            0      TRUE  25245361 KRAS
    ## 86 GAATATAAACTTGCGGTAGT TGG  54770609            1      TRUE  54770606 <NA>
    ## 87 AATGACTGAATATAAACTTG TGG  25245365            0      TRUE  25245368 KRAS
    ## 88 AATGACTGAATATAAACTTG CGG  54770602            0      TRUE  54770599 <NA>
    ## 89 AATGACTAAATATAAACTTG CGA  60822020            1     FALSE  60822023 <NA>
    ## 90 AATGACTGAATATAAACTTC TAG  57999156            1     FALSE  57999159 <NA>
    ## 91 TATATTCAGTCATTTTCAGC AGG  25245392            0      TRUE  25245389 <NA>
    ## 92 TATATTCAGTCATTTTCAGC AGG  54770575            0      TRUE  54770578 <NA>
    ## 93 AATATTCAGTCATTTTCAGC TAG 210618123            1     FALSE 210618126 <NA>
    ##    fiveUTRs threeUTRs  exons   introns intergenic intergenic_distance promoters
    ## 1      <NA>      KRAS   KRAS      <NA>       <NA>                  NA      <NA>
    ## 2      <NA>      <NA>   <NA>   C8orf34       <NA>                  NA      <NA>
    ## 3      <NA>      <NA> KRASP1      <NA>       <NA>                  NA      <NA>
    ## 4      <NA>      <NA>   <NA>      <NA>                          88819      <NA>
    ## 5      <NA>      <NA>   <NA>    SPATA6       <NA>                  NA      <NA>
    ## 6      <NA>      KRAS   KRAS      <NA>       <NA>                  NA      <NA>
    ## 7      <NA>      <NA> KRASP1      <NA>       <NA>                  NA      <NA>
    ## 8      <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 9      <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 10     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 11     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 12     <NA>      <NA>   <NA>      <NA>     FAM47B               23630      <NA>
    ## 13     <NA>      <NA>   <NA>                 <NA>                  NA      <NA>
    ## 14     <NA>      <NA>   <NA>      <NA>  LINC01318               28585      <NA>
    ## 15     <NA>      <NA>   <NA>      <NA>                           4421      <NA>
    ## 16     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 17     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 18     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 19     <NA>      <NA>  ;KRAS      KRAS       <NA>                  NA      <NA>
    ## 20     <NA>      <NA>  ;KRAS      KRAS       <NA>                  NA      <NA>
    ## 21     <NA>      <NA>  ;KRAS      KRAS       <NA>                  NA      <NA>
    ## 22     <NA>      <NA> KRASP1      <NA>       <NA>                  NA      <NA>
    ## 23     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 24     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 25     <NA>      <NA>   <NA> LINC02237       <NA>                  NA      <NA>
    ## 26     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 27     <NA>      <NA> KRASP1      <NA>       <NA>                  NA      <NA>
    ## 28     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 29     <NA>      <NA>   NRAS      <NA>       <NA>                  NA      <NA>
    ## 30     <NA>      <NA> KRASP1      <NA>       <NA>                  NA      <NA>
    ## 31     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 32     <NA>      <NA> KRASP1      <NA>       <NA>                  NA      <NA>
    ## 33     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 34     <NA>      <NA> KRASP1      <NA>       <NA>                  NA      <NA>
    ## 35     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 36     <NA>      <NA> KRASP1      <NA>       <NA>                  NA      <NA>
    ## 37     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 38     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 39     <NA>      <NA>   <NA>   DSCAML1       <NA>                  NA      <NA>
    ## 40     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 41     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 42     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 43     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 44     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 45     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 46     <NA>      <NA> KRASP1      <NA>       <NA>                  NA      <NA>
    ## 47     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 48     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 49     <NA>      <NA> KRASP1      <NA>       <NA>                  NA      <NA>
    ## 50     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 51     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 52     <NA>      <NA> KRASP1      <NA>       <NA>                  NA      <NA>
    ## 53     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 54     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 55     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 56     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 57     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 58     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 59     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 60     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 61     <NA>      <NA> KRASP1      <NA>       <NA>                  NA      <NA>
    ## 62     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 63     <NA>      <NA> KRASP1      <NA>       <NA>                  NA      <NA>
    ## 64     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 65     <NA>      <NA> KRASP1      <NA>       <NA>                  NA      <NA>
    ## 66     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 67     <NA>      <NA>   KRAS      KRAS       <NA>                  NA      <NA>
    ## 68     <NA>      <NA> KRASP1      <NA>       <NA>                  NA      <NA>
    ## 69     <NA>      <NA>   <NA>      KRAS       <NA>                  NA      <NA>
    ## 70     <NA>      <NA>   <NA>      KRAS       <NA>                  NA      <NA>
    ## 71     <NA>      <NA>   <NA>      KRAS       <NA>                  NA      <NA>
    ## 72     <NA>      <NA>   KRAS      <NA>       <NA>                  NA      <NA>
    ## 73     <NA>      <NA> KRASP1      <NA>       <NA>                  NA      <NA>
    ## 74     <NA>      <NA>   KRAS      <NA>       <NA>                  NA      <NA>
    ## 75     <NA>      <NA>   KRAS      <NA>       <NA>                  NA      <NA>
    ## 76     <NA>      <NA> KRASP1      <NA>       <NA>                  NA      <NA>
    ## 77     <NA>      <NA>   KRAS      <NA>       <NA>                  NA      <NA>
    ## 78     <NA>      <NA> KRASP1      <NA>       <NA>                  NA      <NA>
    ## 79     <NA>      <NA>   KRAS      <NA>       <NA>                  NA      <NA>
    ## 80     <NA>      <NA> KRASP1      <NA>       <NA>                  NA      <NA>
    ## 81     <NA>      <NA>   KRAS      <NA>       <NA>                  NA      <NA>
    ## 82     <NA>      <NA> KRASP1      <NA>       <NA>                  NA      <NA>
    ## 83     <NA>      <NA>   KRAS      <NA>       <NA>                  NA      <NA>
    ## 84     <NA>      <NA> KRASP1      <NA>       <NA>                  NA      <NA>
    ## 85     <NA>      <NA>   KRAS      <NA>       <NA>                  NA      <NA>
    ## 86     <NA>      <NA> KRASP1      <NA>       <NA>                  NA      <NA>
    ## 87     <NA>      <NA>   KRAS      <NA>       <NA>                  NA      <NA>
    ## 88     <NA>      <NA> KRASP1      <NA>       <NA>                  NA      <NA>
    ## 89     <NA>      <NA>   <NA> LINC00378       <NA>                  NA      <NA>
    ## 90     <NA>      <NA>   <NA>      <NA>      ARMH4                1601      <NA>
    ## 91     KRAS      <NA>   KRAS      <NA>       <NA>                  NA      <NA>
    ## 92     <NA>      <NA>   <NA>      <NA>     KRASP1                   4      <NA>
    ## 93     <NA>      <NA>   <NA>      HHAT       <NA>                  NA      <NA>
    ##     score_cfd  score_mit
    ## 1  1.00000000 1.00000000
    ## 2  0.12962963 0.00362963
    ## 3  0.36363636 0.61500000
    ## 4  0.80000000 0.68500000
    ## 5  0.05952381 0.04048611
    ## 6  1.00000000 1.00000000
    ## 7  1.00000000 0.73200000
    ## 8  1.00000000 1.00000000
    ## 9  1.00000000 1.00000000
    ## 10 1.00000000 1.00000000
    ## 11 1.00000000 1.00000000
    ## 12 0.69230769 0.00000000
    ## 13 0.04960317 0.00000000
    ## 14 0.11522634 0.10085185
    ## 15 0.06127451 0.03090278
    ## 16 1.00000000 1.00000000
    ## 17 1.00000000 1.00000000
    ## 18 1.00000000 1.00000000
    ## 19 1.00000000 1.00000000
    ## 20 1.00000000 1.00000000
    ## 21 1.00000000 1.00000000
    ## 22 0.65000000 0.39500000
    ## 23 1.00000000 1.00000000
    ## 24 1.00000000 1.00000000
    ## 25 0.05761317 0.10240741
    ## 26 1.00000000 1.00000000
    ## 27 1.00000000 1.00000000
    ## 28 1.00000000 1.00000000
    ## 29 0.05000000 0.39500000
    ## 30 0.71428571 0.73200000
    ## 31 1.00000000 1.00000000
    ## 32 1.00000000 1.00000000
    ## 33 1.00000000 1.00000000
    ## 34 0.19230769 0.00000000
    ## 35 1.00000000 1.00000000
    ## 36 1.00000000 1.00000000
    ## 37 1.00000000 1.00000000
    ## 38 1.00000000 1.00000000
    ## 39 0.06127451 0.03090278
    ## 40 1.00000000 1.00000000
    ## 41 1.00000000 1.00000000
    ## 42 1.00000000 1.00000000
    ## 43 1.00000000 1.00000000
    ## 44 1.00000000 1.00000000
    ## 45 1.00000000 1.00000000
    ## 46 1.00000000 0.07900000
    ## 47 1.00000000 1.00000000
    ## 48 1.00000000 1.00000000
    ## 49 0.42857143 0.00000000
    ## 50 1.00000000 1.00000000
    ## 51 1.00000000 1.00000000
    ## 52 0.90000000 0.58300000
    ## 53 1.00000000 1.00000000
    ## 54 1.00000000 1.00000000
    ## 55 1.00000000 1.00000000
    ## 56 1.00000000 1.00000000
    ## 57 1.00000000 1.00000000
    ## 58 1.00000000 1.00000000
    ## 59 1.00000000 1.00000000
    ## 60 1.00000000 1.00000000
    ## 61 1.00000000 1.00000000
    ## 62 1.00000000 1.00000000
    ## 63 0.90000000 0.61500000
    ## 64 1.00000000 1.00000000
    ## 65 1.00000000 1.00000000
    ## 66 1.00000000 1.00000000
    ## 67 1.00000000 1.00000000
    ## 68 1.00000000 1.00000000
    ## 69 1.00000000 1.00000000
    ## 70 1.00000000 1.00000000
    ## 71 1.00000000 1.00000000
    ## 72 1.00000000 1.00000000
    ## 73 0.04807692 0.04756944
    ## 74 1.00000000 1.00000000
    ## 75 1.00000000 1.00000000
    ## 76 1.00000000 1.00000000
    ## 77 1.00000000 1.00000000
    ## 78 0.06944444 0.06944444
    ## 79 1.00000000 1.00000000
    ## 80 0.25925926 0.25925926
    ## 81 1.00000000 1.00000000
    ## 82 1.00000000 0.82800000
    ## 83 1.00000000 1.00000000
    ## 84 1.00000000 0.61300000
    ## 85 1.00000000 1.00000000
    ## 86 0.28571429 0.31700000
    ## 87 1.00000000 1.00000000
    ## 88 1.00000000 1.00000000
    ## 89 0.06944444 0.04256944
    ## 90 0.11111111 0.00000000
    ## 91 1.00000000 1.00000000
    ## 92 1.00000000 1.00000000
    ## 93 0.25925926 0.15114815
    ## 
    ## $geneAnnotation
    ##                     ID   chr anchor_site strand gene_symbol         gene_id
    ## 1    ENSG00000133703_1 chr12    25209846      -        KRAS ENSG00000133703
    ## 2    ENSG00000133703_1 chr12    25209846      -        KRAS ENSG00000133703
    ## 3    ENSG00000133703_1 chr12    25209846      -        KRAS ENSG00000133703
    ## 4    ENSG00000133703_2 chr12    25209893      +        KRAS ENSG00000133703
    ## 5    ENSG00000133703_2 chr12    25209893      +        KRAS ENSG00000133703
    ## 6    ENSG00000133703_2 chr12    25209893      +        KRAS ENSG00000133703
    ## 7    ENSG00000133703_3 chr12    25215441      -        KRAS ENSG00000133703
    ## 8    ENSG00000133703_4 chr12    25215480      -        KRAS ENSG00000133703
    ## 9    ENSG00000133703_5 chr12    25215474      +        KRAS ENSG00000133703
    ## 10   ENSG00000133703_6 chr12    25215517      +        KRAS ENSG00000133703
    ## 11   ENSG00000133703_7 chr12    25215538      -        KRAS ENSG00000133703
    ## 12   ENSG00000133703_8 chr12    25215556      -        KRAS ENSG00000133703
    ## 13   ENSG00000133703_9 chr12    25215559      -        KRAS ENSG00000133703
    ## 14  ENSG00000133703_10 chr12    25225618      -             ENSG00000275197
    ## 15  ENSG00000133703_10 chr12    25225618      -        KRAS ENSG00000133703
    ## 16  ENSG00000133703_10 chr12    25225618      -        KRAS ENSG00000133703
    ## 17  ENSG00000133703_11 chr12    25225641      +             ENSG00000275197
    ## 18  ENSG00000133703_11 chr12    25225641      +        KRAS ENSG00000133703
    ## 19  ENSG00000133703_11 chr12    25225641      +        KRAS ENSG00000133703
    ## 20  ENSG00000133703_12 chr12    25225656      -             ENSG00000275197
    ## 21  ENSG00000133703_12 chr12    25225656      -        KRAS ENSG00000133703
    ## 22  ENSG00000133703_12 chr12    25225656      -        KRAS ENSG00000133703
    ## 23  ENSG00000133703_13 chr12    25225675      -        KRAS ENSG00000133703
    ## 24  ENSG00000133703_13 chr12    25225675      -        KRAS ENSG00000133703
    ## 25  ENSG00000133703_14 chr12    25225681      -        KRAS ENSG00000133703
    ## 26  ENSG00000133703_14 chr12    25225681      -        KRAS ENSG00000133703
    ## 27  ENSG00000133703_15 chr12    25225698      +        KRAS ENSG00000133703
    ## 28  ENSG00000133703_15 chr12    25225698      +        KRAS ENSG00000133703
    ## 29  ENSG00000133703_16 chr12    25225725      -        KRAS ENSG00000133703
    ## 30  ENSG00000133703_16 chr12    25225725      -        KRAS ENSG00000133703
    ## 31  ENSG00000133703_17 chr12    25225723      +        KRAS ENSG00000133703
    ## 32  ENSG00000133703_17 chr12    25225723      +        KRAS ENSG00000133703
    ## 33  ENSG00000133703_18 chr12    25225735      -        KRAS ENSG00000133703
    ## 34  ENSG00000133703_18 chr12    25225735      -        KRAS ENSG00000133703
    ## 35  ENSG00000133703_19 chr12    25225731      +        KRAS ENSG00000133703
    ## 36  ENSG00000133703_19 chr12    25225731      +        KRAS ENSG00000133703
    ## 37  ENSG00000133703_20 chr12    25225756      -        KRAS ENSG00000133703
    ## 38  ENSG00000133703_20 chr12    25225756      -        KRAS ENSG00000133703
    ## 39  ENSG00000133703_21 chr12    25225772      +        KRAS ENSG00000133703
    ## 40  ENSG00000133703_21 chr12    25225772      +        KRAS ENSG00000133703
    ## 41  ENSG00000133703_22 chr12    25225773      +        KRAS ENSG00000133703
    ## 42  ENSG00000133703_22 chr12    25225773      +        KRAS ENSG00000133703
    ## 43  ENSG00000133703_23 chr12    25227234      -        KRAS ENSG00000133703
    ## 44  ENSG00000133703_23 chr12    25227234      -        KRAS ENSG00000133703
    ## 45  ENSG00000133703_24 chr12    25227235      -        KRAS ENSG00000133703
    ## 46  ENSG00000133703_24 chr12    25227235      -        KRAS ENSG00000133703
    ## 47  ENSG00000133703_25 chr12    25227238      -        KRAS ENSG00000133703
    ## 48  ENSG00000133703_25 chr12    25227238      -        KRAS ENSG00000133703
    ## 49  ENSG00000133703_26 chr12    25227237      +        KRAS ENSG00000133703
    ## 50  ENSG00000133703_26 chr12    25227237      +        KRAS ENSG00000133703
    ## 51  ENSG00000133703_27 chr12    25227271      +        KRAS ENSG00000133703
    ## 52  ENSG00000133703_27 chr12    25227271      +        KRAS ENSG00000133703
    ## 53  ENSG00000133703_28 chr12    25227299      -        KRAS ENSG00000133703
    ## 54  ENSG00000133703_28 chr12    25227299      -        KRAS ENSG00000133703
    ## 55  ENSG00000133703_29 chr12    25227300      -        KRAS ENSG00000133703
    ## 56  ENSG00000133703_29 chr12    25227300      -        KRAS ENSG00000133703
    ## 57  ENSG00000133703_30 chr12    25227303      -        KRAS ENSG00000133703
    ## 58  ENSG00000133703_30 chr12    25227303      -        KRAS ENSG00000133703
    ## 59  ENSG00000133703_31 chr12    25227304      -        KRAS ENSG00000133703
    ## 60  ENSG00000133703_31 chr12    25227304      -        KRAS ENSG00000133703
    ## 61  ENSG00000133703_32 chr12    25227305      -        KRAS ENSG00000133703
    ## 62  ENSG00000133703_32 chr12    25227305      -        KRAS ENSG00000133703
    ## 63  ENSG00000133703_33 chr12    25227310      -        KRAS ENSG00000133703
    ## 64  ENSG00000133703_33 chr12    25227310      -        KRAS ENSG00000133703
    ## 65  ENSG00000133703_34 chr12    25227312      +        KRAS ENSG00000133703
    ## 66  ENSG00000133703_34 chr12    25227312      +        KRAS ENSG00000133703
    ## 67  ENSG00000133703_35 chr12    25227324      -        KRAS ENSG00000133703
    ## 68  ENSG00000133703_35 chr12    25227324      -        KRAS ENSG00000133703
    ## 69  ENSG00000133703_36 chr12    25227325      -        KRAS ENSG00000133703
    ## 70  ENSG00000133703_36 chr12    25227325      -        KRAS ENSG00000133703
    ## 71  ENSG00000133703_37 chr12    25227342      -        KRAS ENSG00000133703
    ## 72  ENSG00000133703_37 chr12    25227342      -        KRAS ENSG00000133703
    ## 73  ENSG00000133703_38 chr12    25227350      -        KRAS ENSG00000133703
    ## 74  ENSG00000133703_38 chr12    25227350      -        KRAS ENSG00000133703
    ## 75  ENSG00000133703_39 chr12    25227369      -        KRAS ENSG00000133703
    ## 76  ENSG00000133703_39 chr12    25227369      -        KRAS ENSG00000133703
    ## 77  ENSG00000133703_40 chr12    25227370      +        KRAS ENSG00000133703
    ## 78  ENSG00000133703_40 chr12    25227370      +        KRAS ENSG00000133703
    ## 79  ENSG00000133703_41 chr12    25227386      -        KRAS ENSG00000133703
    ## 80  ENSG00000133703_41 chr12    25227386      -        KRAS ENSG00000133703
    ## 81  ENSG00000133703_42 chr12    25227406      -        KRAS ENSG00000133703
    ## 82  ENSG00000133703_42 chr12    25227406      -        KRAS ENSG00000133703
    ## 83  ENSG00000133703_43 chr12    25227403      +        KRAS ENSG00000133703
    ## 84  ENSG00000133703_43 chr12    25227403      +        KRAS ENSG00000133703
    ## 85  ENSG00000133703_47 chr12    25245278      -        KRAS ENSG00000133703
    ## 86  ENSG00000133703_47 chr12    25245278      -        KRAS ENSG00000133703
    ## 87  ENSG00000133703_47 chr12    25245278      -        KRAS ENSG00000133703
    ## 88  ENSG00000133703_47 chr12    25245278      -        KRAS ENSG00000133703
    ## 89  ENSG00000133703_48 chr12    25245280      +        KRAS ENSG00000133703
    ## 90  ENSG00000133703_48 chr12    25245280      +        KRAS ENSG00000133703
    ## 91  ENSG00000133703_48 chr12    25245280      +        KRAS ENSG00000133703
    ## 92  ENSG00000133703_48 chr12    25245280      +        KRAS ENSG00000133703
    ## 93  ENSG00000133703_49 chr12    25245302      -        KRAS ENSG00000133703
    ## 94  ENSG00000133703_49 chr12    25245302      -        KRAS ENSG00000133703
    ## 95  ENSG00000133703_49 chr12    25245302      -        KRAS ENSG00000133703
    ## 96  ENSG00000133703_49 chr12    25245302      -        KRAS ENSG00000133703
    ## 97  ENSG00000133703_50 chr12    25245327      +        KRAS ENSG00000133703
    ## 98  ENSG00000133703_50 chr12    25245327      +        KRAS ENSG00000133703
    ## 99  ENSG00000133703_50 chr12    25245327      +        KRAS ENSG00000133703
    ## 100 ENSG00000133703_50 chr12    25245327      +        KRAS ENSG00000133703
    ## 101 ENSG00000133703_51 chr12    25245346      -        KRAS ENSG00000133703
    ## 102 ENSG00000133703_51 chr12    25245346      -        KRAS ENSG00000133703
    ## 103 ENSG00000133703_51 chr12    25245346      -        KRAS ENSG00000133703
    ## 104 ENSG00000133703_51 chr12    25245346      -        KRAS ENSG00000133703
    ## 105 ENSG00000133703_52 chr12    25245352      -        KRAS ENSG00000133703
    ## 106 ENSG00000133703_52 chr12    25245352      -        KRAS ENSG00000133703
    ## 107 ENSG00000133703_52 chr12    25245352      -        KRAS ENSG00000133703
    ## 108 ENSG00000133703_52 chr12    25245352      -        KRAS ENSG00000133703
    ## 109 ENSG00000133703_53 chr12    25245355      -        KRAS ENSG00000133703
    ## 110 ENSG00000133703_53 chr12    25245355      -        KRAS ENSG00000133703
    ## 111 ENSG00000133703_53 chr12    25245355      -        KRAS ENSG00000133703
    ## 112 ENSG00000133703_53 chr12    25245355      -        KRAS ENSG00000133703
    ## 113 ENSG00000133703_54 chr12    25245361      -        KRAS ENSG00000133703
    ## 114 ENSG00000133703_54 chr12    25245361      -        KRAS ENSG00000133703
    ## 115 ENSG00000133703_54 chr12    25245361      -        KRAS ENSG00000133703
    ## 116 ENSG00000133703_54 chr12    25245361      -        KRAS ENSG00000133703
    ## 117 ENSG00000133703_55 chr12    25245368      -        KRAS ENSG00000133703
    ## 118 ENSG00000133703_55 chr12    25245368      -        KRAS ENSG00000133703
    ## 119 ENSG00000133703_55 chr12    25245368      -        KRAS ENSG00000133703
    ## 120 ENSG00000133703_55 chr12    25245368      -        KRAS ENSG00000133703
    ## 121 ENSG00000133703_56 chr12    25245389      +        KRAS ENSG00000133703
    ## 122 ENSG00000133703_56 chr12    25245389      +        KRAS ENSG00000133703
    ## 123 ENSG00000133703_56 chr12    25245389      +        KRAS ENSG00000133703
    ## 124 ENSG00000133703_56 chr12    25245389      +        KRAS ENSG00000133703
    ##               tx_id      protein_id cut_cds cut_fiveUTRs cut_threeUTRs
    ## 1   ENST00000256078            <NA>   FALSE        FALSE          TRUE
    ## 2   ENST00000311936 ENSP00000308495    TRUE        FALSE         FALSE
    ## 3   ENST00000557334 ENSP00000452512    TRUE        FALSE         FALSE
    ## 4   ENST00000256078            <NA>   FALSE        FALSE          TRUE
    ## 5   ENST00000311936 ENSP00000308495    TRUE        FALSE         FALSE
    ## 6   ENST00000557334 ENSP00000452512    TRUE        FALSE         FALSE
    ## 7   ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 8   ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 9   ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 10  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 11  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 12  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 13  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 14  ENST00000620933            <NA>   FALSE        FALSE         FALSE
    ## 15  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 16  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 17  ENST00000620933            <NA>   FALSE        FALSE         FALSE
    ## 18  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 19  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 20  ENST00000620933            <NA>   FALSE        FALSE         FALSE
    ## 21  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 22  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 23  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 24  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 25  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 26  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 27  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 28  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 29  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 30  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 31  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 32  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 33  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 34  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 35  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 36  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 37  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 38  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 39  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 40  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 41  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 42  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 43  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 44  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 45  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 46  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 47  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 48  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 49  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 50  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 51  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 52  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 53  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 54  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 55  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 56  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 57  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 58  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 59  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 60  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 61  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 62  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 63  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 64  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 65  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 66  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 67  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 68  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 69  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 70  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 71  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 72  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 73  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 74  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 75  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 76  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 77  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 78  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 79  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 80  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 81  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 82  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 83  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 84  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 85  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 86  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 87  ENST00000557334 ENSP00000256078    TRUE        FALSE         FALSE
    ## 88  ENST00000556131 ENSP00000256078    TRUE        FALSE         FALSE
    ## 89  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 90  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 91  ENST00000557334 ENSP00000256078    TRUE        FALSE         FALSE
    ## 92  ENST00000556131 ENSP00000256078    TRUE        FALSE         FALSE
    ## 93  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 94  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 95  ENST00000557334 ENSP00000256078    TRUE        FALSE         FALSE
    ## 96  ENST00000556131 ENSP00000256078    TRUE        FALSE         FALSE
    ## 97  ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 98  ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 99  ENST00000557334 ENSP00000256078    TRUE        FALSE         FALSE
    ## 100 ENST00000556131 ENSP00000256078    TRUE        FALSE         FALSE
    ## 101 ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 102 ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 103 ENST00000557334 ENSP00000256078    TRUE        FALSE         FALSE
    ## 104 ENST00000556131 ENSP00000256078    TRUE        FALSE         FALSE
    ## 105 ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 106 ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 107 ENST00000557334 ENSP00000256078    TRUE        FALSE         FALSE
    ## 108 ENST00000556131 ENSP00000256078    TRUE        FALSE         FALSE
    ## 109 ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 110 ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 111 ENST00000557334 ENSP00000256078    TRUE        FALSE         FALSE
    ## 112 ENST00000556131 ENSP00000256078    TRUE        FALSE         FALSE
    ## 113 ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 114 ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 115 ENST00000557334 ENSP00000256078    TRUE        FALSE         FALSE
    ## 116 ENST00000556131 ENSP00000256078    TRUE        FALSE         FALSE
    ## 117 ENST00000256078 ENSP00000256078    TRUE        FALSE         FALSE
    ## 118 ENST00000311936 ENSP00000256078    TRUE        FALSE         FALSE
    ## 119 ENST00000557334 ENSP00000256078    TRUE        FALSE         FALSE
    ## 120 ENST00000556131 ENSP00000256078    TRUE        FALSE         FALSE
    ## 121 ENST00000256078 ENSP00000256078   FALSE         TRUE         FALSE
    ## 122 ENST00000311936 ENSP00000256078   FALSE        FALSE         FALSE
    ## 123 ENST00000557334 ENSP00000256078   FALSE        FALSE         FALSE
    ## 124 ENST00000556131 ENSP00000256078   FALSE        FALSE         FALSE
    ##     cut_introns percentCDS aminoAcidIndex downtreamATG percentTx nIsoforms
    ## 1         FALSE         NA             NA           NA      15.3         3
    ## 2         FALSE       91.0            172            1      13.3         3
    ## 3         FALSE       77.6             59            1      35.6         3
    ## 4         FALSE         NA             NA           NA      14.4         3
    ## 5         FALSE       82.7            157            2      12.4         3
    ## 6         FALSE       57.0             44            2      31.1         3
    ## 7         FALSE      100.0            190            0      14.0         1
    ## 8         FALSE       93.2            177            1      13.3         1
    ## 9         FALSE       94.2            179            1      13.4         1
    ## 10        FALSE       86.7            165            1      12.6         1
    ## 11        FALSE       83.0            158            1      12.2         1
    ## 12        FALSE       79.8            152            1      11.9         1
    ## 13        FALSE       79.3            151            1      11.8         1
    ## 14        FALSE         NA             NA           NA      91.7         1
    ## 15        FALSE       78.2            149            1      11.7         2
    ## 16        FALSE       78.7            149            2      12.0         2
    ## 17        FALSE         NA             NA           NA      95.7         1
    ## 18        FALSE       74.2            141            1      11.3         2
    ## 19        FALSE       74.6            141            2      11.6         2
    ## 20        FALSE         NA             NA           NA      98.4         1
    ## 21        FALSE       71.6            136            1      11.0         2
    ## 22        FALSE       72.0            136            2      11.3         2
    ## 23        FALSE       68.2            130            1      10.7         2
    ## 24        FALSE       68.6            130            2      10.9         2
    ## 25        FALSE       67.2            128            1      10.6         2
    ## 26        FALSE       67.5            128            2      10.8         2
    ## 27        FALSE       64.2            122            1      10.2         2
    ## 28        FALSE       64.6            122            2      10.5         2
    ## 29        FALSE       59.5            113            1       9.7         2
    ## 30        FALSE       59.8            113            2      10.0         2
    ## 31        FALSE       59.8            114            1       9.8         2
    ## 32        FALSE       60.1            114            2      10.0         2
    ## 33        FALSE       57.7            110            2       9.6         2
    ## 34        FALSE       58.0            110            3       9.8         2
    ## 35        FALSE       58.4            111            1       9.6         2
    ## 36        FALSE       58.7            111            2       9.9         2
    ## 37        FALSE       54.0            103            1       9.2         2
    ## 38        FALSE       54.3            103            2       9.4         2
    ## 39        FALSE       51.2             98            1       8.9         2
    ## 40        FALSE       51.5             98            2       9.1         2
    ## 41        FALSE       51.1             97            1       8.9         2
    ## 42        FALSE       51.3             97            2       9.1         2
    ## 43        FALSE       50.9             97            1       8.8         2
    ## 44        FALSE       51.1             97            2       9.0         2
    ## 45        FALSE       50.7             97            1       8.8         2
    ## 46        FALSE       51.0             97            2       9.0         2
    ## 47        FALSE       50.2             96            1       8.8         2
    ## 48        FALSE       50.4             96            2       9.0         2
    ## 49        FALSE       50.4             96            1       8.8         2
    ## 50        FALSE       50.6             96            2       9.0         2
    ## 51        FALSE       44.4             85            1       8.2         2
    ## 52        FALSE       44.6             85            1       8.3         2
    ## 53        FALSE       39.5             75            1       7.6         2
    ## 54        FALSE       39.7             75            1       7.8         2
    ## 55        FALSE       39.3             75            1       7.6         2
    ## 56        FALSE       39.5             75            1       7.8         2
    ## 57        FALSE       38.8             74            1       7.6         2
    ## 58        FALSE       39.0             74            1       7.7         2
    ## 59        FALSE       38.6             74            1       7.6         2
    ## 60        FALSE       38.8             74            1       7.7         2
    ## 61        FALSE       38.4             73            1       7.5         2
    ## 62        FALSE       38.6             73            1       7.7         2
    ## 63        FALSE       37.5             72            1       7.4         2
    ## 64        FALSE       37.7             72            1       7.6         2
    ## 65        FALSE       37.2             71            2       7.4         2
    ## 66        FALSE       37.4             71            2       7.6         2
    ## 67        FALSE       35.1             67            2       7.2         2
    ## 68        FALSE       35.3             67            2       7.4         2
    ## 69        FALSE       34.9             67            2       7.2         2
    ## 70        FALSE       35.1             67            2       7.3         2
    ## 71        FALSE       31.9             61            3       6.9         2
    ## 72        FALSE       32.1             61            3       7.0         2
    ## 73        FALSE       30.5             58            3       6.7         2
    ## 74        FALSE       30.7             58            3       6.9         2
    ## 75        FALSE       27.2             52            3       6.4         2
    ## 76        FALSE       27.3             52            3       6.5         2
    ## 77        FALSE       27.0             52            3       6.3         2
    ## 78        FALSE       27.2             52            3       6.5         2
    ## 79        FALSE       24.2             46            3       6.0         2
    ## 80        FALSE       24.3             46            3       6.2         2
    ## 81        FALSE       20.7             40            3       5.7         2
    ## 82        FALSE       20.8             40            3       5.8         2
    ## 83        FALSE       21.2             41            3       5.7         2
    ## 84        FALSE       21.3             41            3       5.9         2
    ## 85        FALSE       18.8             36            3       5.5         4
    ## 86        FALSE       18.9             36            3       5.6         4
    ## 87        FALSE       46.9             36            2      28.9         4
    ## 88        FALSE       81.1             36            1      16.7         4
    ## 89        FALSE       18.4             35            3       5.4         4
    ## 90        FALSE       18.5             35            3       5.6         4
    ## 91        FALSE       46.1             35            2      28.7         4
    ## 92        FALSE       79.5             35            1      16.6         4
    ## 93        FALSE       14.6             28            3       5.0         4
    ## 94        FALSE       14.6             28            3       5.1         4
    ## 95        FALSE       36.4             28            2      26.6         4
    ## 96        FALSE       62.9             28            1      15.3         4
    ## 97        FALSE       10.2             20            2       4.6         4
    ## 98        FALSE       10.2             20            2       4.7         4
    ## 99        FALSE       25.4             20            2      24.2         4
    ## 100       FALSE       43.9             20            1      13.9         4
    ## 101       FALSE        6.8             13            2       4.2         4
    ## 102       FALSE        6.9             13            2       4.3         4
    ## 103       FALSE       17.1             13            2      22.4         4
    ## 104       FALSE       29.5             13            1      12.7         4
    ## 105       FALSE        5.8             11            2       4.1         4
    ## 106       FALSE        5.8             11            2       4.2         4
    ## 107       FALSE       14.5             11            2      21.9         4
    ## 108       FALSE       25.0             11            1      12.4         4
    ## 109       FALSE        5.3             10            2       4.1         4
    ## 110       FALSE        5.3             10            2       4.1         4
    ## 111       FALSE       13.2             10            2      21.6         4
    ## 112       FALSE       22.7             10            1      12.2         4
    ## 113       FALSE        4.2              8            2       3.9         4
    ## 114       FALSE        4.2              8            2       4.0         4
    ## 115       FALSE       10.5              8            2      21.0         4
    ## 116       FALSE       18.2              8            1      11.9         4
    ## 117       FALSE        3.0              6            2       3.8         4
    ## 118       FALSE        3.0              6            2       3.9         4
    ## 119       FALSE        7.5              6            2      20.3         4
    ## 120       FALSE       12.9              6            1      11.4         4
    ## 121       FALSE         NA             NA           NA       3.4         4
    ## 122       FALSE         NA             NA           NA       3.5         4
    ## 123       FALSE         NA             NA           NA      18.3         4
    ## 124       FALSE         NA             NA           NA      10.2         4
    ##     totalIsoforms percentIsoforms isCommonExon nCodingIsoforms
    ## 1               4              75        FALSE               3
    ## 2               4              75        FALSE               3
    ## 3               4              75        FALSE               3
    ## 4               4              75        FALSE               3
    ## 5               4              75        FALSE               3
    ## 6               4              75        FALSE               3
    ## 7               4              25        FALSE               1
    ## 8               4              25        FALSE               1
    ## 9               4              25        FALSE               1
    ## 10              4              25        FALSE               1
    ## 11              4              25        FALSE               1
    ## 12              4              25        FALSE               1
    ## 13              4              25        FALSE               1
    ## 14              1             100         TRUE               1
    ## 15              4              50        FALSE               2
    ## 16              4              50        FALSE               2
    ## 17              1             100         TRUE               1
    ## 18              4              50        FALSE               2
    ## 19              4              50        FALSE               2
    ## 20              1             100         TRUE               1
    ## 21              4              50        FALSE               2
    ## 22              4              50        FALSE               2
    ## 23              4              50        FALSE               2
    ## 24              4              50        FALSE               2
    ## 25              4              50        FALSE               2
    ## 26              4              50        FALSE               2
    ## 27              4              50        FALSE               2
    ## 28              4              50        FALSE               2
    ## 29              4              50        FALSE               2
    ## 30              4              50        FALSE               2
    ## 31              4              50        FALSE               2
    ## 32              4              50        FALSE               2
    ## 33              4              50        FALSE               2
    ## 34              4              50        FALSE               2
    ## 35              4              50        FALSE               2
    ## 36              4              50        FALSE               2
    ## 37              4              50        FALSE               2
    ## 38              4              50        FALSE               2
    ## 39              4              50        FALSE               2
    ## 40              4              50        FALSE               2
    ## 41              4              50        FALSE               2
    ## 42              4              50        FALSE               2
    ## 43              4              50        FALSE               2
    ## 44              4              50        FALSE               2
    ## 45              4              50        FALSE               2
    ## 46              4              50        FALSE               2
    ## 47              4              50        FALSE               2
    ## 48              4              50        FALSE               2
    ## 49              4              50        FALSE               2
    ## 50              4              50        FALSE               2
    ## 51              4              50        FALSE               2
    ## 52              4              50        FALSE               2
    ## 53              4              50        FALSE               2
    ## 54              4              50        FALSE               2
    ## 55              4              50        FALSE               2
    ## 56              4              50        FALSE               2
    ## 57              4              50        FALSE               2
    ## 58              4              50        FALSE               2
    ## 59              4              50        FALSE               2
    ## 60              4              50        FALSE               2
    ## 61              4              50        FALSE               2
    ## 62              4              50        FALSE               2
    ## 63              4              50        FALSE               2
    ## 64              4              50        FALSE               2
    ## 65              4              50        FALSE               2
    ## 66              4              50        FALSE               2
    ## 67              4              50        FALSE               2
    ## 68              4              50        FALSE               2
    ## 69              4              50        FALSE               2
    ## 70              4              50        FALSE               2
    ## 71              4              50        FALSE               2
    ## 72              4              50        FALSE               2
    ## 73              4              50        FALSE               2
    ## 74              4              50        FALSE               2
    ## 75              4              50        FALSE               2
    ## 76              4              50        FALSE               2
    ## 77              4              50        FALSE               2
    ## 78              4              50        FALSE               2
    ## 79              4              50        FALSE               2
    ## 80              4              50        FALSE               2
    ## 81              4              50        FALSE               2
    ## 82              4              50        FALSE               2
    ## 83              4              50        FALSE               2
    ## 84              4              50        FALSE               2
    ## 85              4             100         TRUE               4
    ## 86              4             100         TRUE               4
    ## 87              4             100         TRUE               4
    ## 88              4             100         TRUE               4
    ## 89              4             100         TRUE               4
    ## 90              4             100         TRUE               4
    ## 91              4             100         TRUE               4
    ## 92              4             100         TRUE               4
    ## 93              4             100         TRUE               4
    ## 94              4             100         TRUE               4
    ## 95              4             100         TRUE               4
    ## 96              4             100         TRUE               4
    ## 97              4             100         TRUE               4
    ## 98              4             100         TRUE               4
    ## 99              4             100         TRUE               4
    ## 100             4             100         TRUE               4
    ## 101             4             100         TRUE               4
    ## 102             4             100         TRUE               4
    ## 103             4             100         TRUE               4
    ## 104             4             100         TRUE               4
    ## 105             4             100         TRUE               4
    ## 106             4             100         TRUE               4
    ## 107             4             100         TRUE               4
    ## 108             4             100         TRUE               4
    ## 109             4             100         TRUE               4
    ## 110             4             100         TRUE               4
    ## 111             4             100         TRUE               4
    ## 112             4             100         TRUE               4
    ## 113             4             100         TRUE               4
    ## 114             4             100         TRUE               4
    ## 115             4             100         TRUE               4
    ## 116             4             100         TRUE               4
    ## 117             4             100         TRUE               4
    ## 118             4             100         TRUE               4
    ## 119             4             100         TRUE               4
    ## 120             4             100         TRUE               4
    ## 121             4             100         TRUE               4
    ## 122             4             100         TRUE               4
    ## 123             4             100         TRUE               4
    ## 124             4             100         TRUE               4
    ##     totalCodingIsoforms percentCodingIsoforms isCommonCodingExon
    ## 1                     4                    75              FALSE
    ## 2                     4                    75              FALSE
    ## 3                     4                    75              FALSE
    ## 4                     4                    75              FALSE
    ## 5                     4                    75              FALSE
    ## 6                     4                    75              FALSE
    ## 7                     4                    25              FALSE
    ## 8                     4                    25              FALSE
    ## 9                     4                    25              FALSE
    ## 10                    4                    25              FALSE
    ## 11                    4                    25              FALSE
    ## 12                    4                    25              FALSE
    ## 13                    4                    25              FALSE
    ## 14                   NA                    NA                 NA
    ## 15                    4                    50              FALSE
    ## 16                    4                    50              FALSE
    ## 17                   NA                    NA                 NA
    ## 18                    4                    50              FALSE
    ## 19                    4                    50              FALSE
    ## 20                   NA                    NA                 NA
    ## 21                    4                    50              FALSE
    ## 22                    4                    50              FALSE
    ## 23                    4                    50              FALSE
    ## 24                    4                    50              FALSE
    ## 25                    4                    50              FALSE
    ## 26                    4                    50              FALSE
    ## 27                    4                    50              FALSE
    ## 28                    4                    50              FALSE
    ## 29                    4                    50              FALSE
    ## 30                    4                    50              FALSE
    ## 31                    4                    50              FALSE
    ## 32                    4                    50              FALSE
    ## 33                    4                    50              FALSE
    ## 34                    4                    50              FALSE
    ## 35                    4                    50              FALSE
    ## 36                    4                    50              FALSE
    ## 37                    4                    50              FALSE
    ## 38                    4                    50              FALSE
    ## 39                    4                    50              FALSE
    ## 40                    4                    50              FALSE
    ## 41                    4                    50              FALSE
    ## 42                    4                    50              FALSE
    ## 43                    4                    50              FALSE
    ## 44                    4                    50              FALSE
    ## 45                    4                    50              FALSE
    ## 46                    4                    50              FALSE
    ## 47                    4                    50              FALSE
    ## 48                    4                    50              FALSE
    ## 49                    4                    50              FALSE
    ## 50                    4                    50              FALSE
    ## 51                    4                    50              FALSE
    ## 52                    4                    50              FALSE
    ## 53                    4                    50              FALSE
    ## 54                    4                    50              FALSE
    ## 55                    4                    50              FALSE
    ## 56                    4                    50              FALSE
    ## 57                    4                    50              FALSE
    ## 58                    4                    50              FALSE
    ## 59                    4                    50              FALSE
    ## 60                    4                    50              FALSE
    ## 61                    4                    50              FALSE
    ## 62                    4                    50              FALSE
    ## 63                    4                    50              FALSE
    ## 64                    4                    50              FALSE
    ## 65                    4                    50              FALSE
    ## 66                    4                    50              FALSE
    ## 67                    4                    50              FALSE
    ## 68                    4                    50              FALSE
    ## 69                    4                    50              FALSE
    ## 70                    4                    50              FALSE
    ## 71                    4                    50              FALSE
    ## 72                    4                    50              FALSE
    ## 73                    4                    50              FALSE
    ## 74                    4                    50              FALSE
    ## 75                    4                    50              FALSE
    ## 76                    4                    50              FALSE
    ## 77                    4                    50              FALSE
    ## 78                    4                    50              FALSE
    ## 79                    4                    50              FALSE
    ## 80                    4                    50              FALSE
    ## 81                    4                    50              FALSE
    ## 82                    4                    50              FALSE
    ## 83                    4                    50              FALSE
    ## 84                    4                    50              FALSE
    ## 85                    4                   100               TRUE
    ## 86                    4                   100               TRUE
    ## 87                    4                   100               TRUE
    ## 88                    4                   100               TRUE
    ## 89                    4                   100               TRUE
    ## 90                    4                   100               TRUE
    ## 91                    4                   100               TRUE
    ## 92                    4                   100               TRUE
    ## 93                    4                   100               TRUE
    ## 94                    4                   100               TRUE
    ## 95                    4                   100               TRUE
    ## 96                    4                   100               TRUE
    ## 97                    4                   100               TRUE
    ## 98                    4                   100               TRUE
    ## 99                    4                   100               TRUE
    ## 100                   4                   100               TRUE
    ## 101                   4                   100               TRUE
    ## 102                   4                   100               TRUE
    ## 103                   4                   100               TRUE
    ## 104                   4                   100               TRUE
    ## 105                   4                   100               TRUE
    ## 106                   4                   100               TRUE
    ## 107                   4                   100               TRUE
    ## 108                   4                   100               TRUE
    ## 109                   4                   100               TRUE
    ## 110                   4                   100               TRUE
    ## 111                   4                   100               TRUE
    ## 112                   4                   100               TRUE
    ## 113                   4                   100               TRUE
    ## 114                   4                   100               TRUE
    ## 115                   4                   100               TRUE
    ## 116                   4                   100               TRUE
    ## 117                   4                   100               TRUE
    ## 118                   4                   100               TRUE
    ## 119                   4                   100               TRUE
    ## 120                   4                   100               TRUE
    ## 121                   4                   100               TRUE
    ## 122                   4                   100               TRUE
    ## 123                   4                   100               TRUE
    ## 124                   4                   100               TRUE
    ## 
    ## $enzymeAnnotation
    ##                    ID EcoRI  KpnI BsmBI  BsaI  BbsI  PacI  MluI
    ## 1   ENSG00000133703_1 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 2   ENSG00000133703_2 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 3   ENSG00000133703_3 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 4   ENSG00000133703_4 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 5   ENSG00000133703_5 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 6   ENSG00000133703_6 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 7   ENSG00000133703_7 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 8   ENSG00000133703_8 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 9   ENSG00000133703_9 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 10 ENSG00000133703_10 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 11 ENSG00000133703_11 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 12 ENSG00000133703_12 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 13 ENSG00000133703_13 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 14 ENSG00000133703_14 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 15 ENSG00000133703_15 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 16 ENSG00000133703_16 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 17 ENSG00000133703_17 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 18 ENSG00000133703_18 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 19 ENSG00000133703_19 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 20 ENSG00000133703_20 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 21 ENSG00000133703_21 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 22 ENSG00000133703_22 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 23 ENSG00000133703_23 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 24 ENSG00000133703_24 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 25 ENSG00000133703_25 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 26 ENSG00000133703_26 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 27 ENSG00000133703_27 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 28 ENSG00000133703_28 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 29 ENSG00000133703_29 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 30 ENSG00000133703_30 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 31 ENSG00000133703_31 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 32 ENSG00000133703_32 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 33 ENSG00000133703_33 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 34 ENSG00000133703_34 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 35 ENSG00000133703_35 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 36 ENSG00000133703_36 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 37 ENSG00000133703_37 FALSE FALSE  TRUE FALSE FALSE FALSE FALSE
    ## 38 ENSG00000133703_38 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 39 ENSG00000133703_39 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 40 ENSG00000133703_40 FALSE FALSE  TRUE FALSE FALSE FALSE FALSE
    ## 41 ENSG00000133703_41 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 42 ENSG00000133703_42 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 43 ENSG00000133703_43 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 44 ENSG00000133703_44 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 45 ENSG00000133703_45 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 46 ENSG00000133703_46 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 47 ENSG00000133703_47 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 48 ENSG00000133703_48 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 49 ENSG00000133703_49 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 50 ENSG00000133703_50 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 51 ENSG00000133703_51 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 52 ENSG00000133703_52 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 53 ENSG00000133703_53 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 54 ENSG00000133703_54 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 55 ENSG00000133703_55 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 56 ENSG00000133703_56 FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 
    ## $snps
    ##                  ID        rs  rs_site rs_site_rel allele_ref allele_minor
    ## 1 ENSG00000133703_1 rs1137282 25209843           0          A            G
    ##   MAF_1000G MAF_TOPMED type length
    ## 1    0.1755    0.19671  snp      1

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
    ##  [1] crisprScore_1.1.6                 crisprScoreData_1.1.3            
    ##  [3] ExperimentHub_2.3.5               AnnotationHub_3.3.9              
    ##  [5] BiocFileCache_2.3.4               dbplyr_2.1.1                     
    ##  [7] BSgenome.Hsapiens.UCSC.hg38_1.4.4 BSgenome_1.63.5                  
    ##  [9] rtracklayer_1.55.4                Biostrings_2.63.2                
    ## [11] XVector_0.35.0                    GenomicRanges_1.47.6             
    ## [13] GenomeInfoDb_1.31.6               IRanges_2.29.1                   
    ## [15] S4Vectors_0.33.11                 BiocGenerics_0.41.2              
    ## [17] crisprDesignData_0.99.7           crisprDesign_0.99.97             
    ## [19] crisprBase_1.1.2                 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bitops_1.0-7                  matrixStats_0.61.0           
    ##  [3] bit64_4.0.5                   filelock_1.0.2               
    ##  [5] progress_1.2.2                httr_1.4.2                   
    ##  [7] tools_4.2.0                   utf8_1.2.2                   
    ##  [9] R6_2.5.1                      DBI_1.1.2                    
    ## [11] tidyselect_1.1.2              prettyunits_1.1.1            
    ## [13] bit_4.0.4                     curl_4.3.2                   
    ## [15] compiler_4.2.0                crisprBowtie_1.1.1           
    ## [17] cli_3.3.0                     Biobase_2.55.0               
    ## [19] basilisk.utils_1.5.0          xml2_1.3.3                   
    ## [21] DelayedArray_0.21.2           randomForest_4.7-1           
    ## [23] readr_2.1.2                   rappdirs_0.3.3               
    ## [25] stringr_1.4.0                 digest_0.6.29                
    ## [27] Rsamtools_2.11.0              rmarkdown_2.13               
    ## [29] basilisk_1.9.1                pkgconfig_2.0.3              
    ## [31] htmltools_0.5.2               MatrixGenerics_1.7.0         
    ## [33] fastmap_1.1.0                 rlang_1.0.2                  
    ## [35] rstudioapi_0.13               RSQLite_2.2.12               
    ## [37] shiny_1.7.1                   BiocIO_1.5.0                 
    ## [39] generics_0.1.2                jsonlite_1.8.0               
    ## [41] vroom_1.5.7                   BiocParallel_1.29.18         
    ## [43] dplyr_1.0.8                   VariantAnnotation_1.41.3     
    ## [45] RCurl_1.98-1.6                magrittr_2.0.2               
    ## [47] GenomeInfoDbData_1.2.7        Matrix_1.4-0                 
    ## [49] Rcpp_1.0.8.3                  fansi_1.0.2                  
    ## [51] reticulate_1.24               Rbowtie_1.35.0               
    ## [53] lifecycle_1.0.1               stringi_1.7.6                
    ## [55] yaml_2.3.5                    SummarizedExperiment_1.25.3  
    ## [57] zlibbioc_1.41.0               grid_4.2.0                   
    ## [59] blob_1.2.2                    promises_1.2.0.1             
    ## [61] parallel_4.2.0                crayon_1.5.0                 
    ## [63] crisprBwa_1.1.2               dir.expiry_1.3.0             
    ## [65] lattice_0.20-45               GenomicFeatures_1.47.13      
    ## [67] hms_1.1.1                     KEGGREST_1.35.0              
    ## [69] knitr_1.37                    pillar_1.7.0                 
    ## [71] rjson_0.2.21                  biomaRt_2.51.3               
    ## [73] BiocVersion_3.15.0            XML_3.99-0.9                 
    ## [75] glue_1.6.2                    evaluate_0.15                
    ## [77] BiocManager_1.30.16           httpuv_1.6.5                 
    ## [79] png_0.1-7                     vctrs_0.3.8                  
    ## [81] tzdb_0.2.0                    purrr_0.3.4                  
    ## [83] assertthat_0.2.1              cachem_1.0.6                 
    ## [85] xfun_0.30                     mime_0.12                    
    ## [87] Rbwa_1.1.0                    xtable_1.8-4                 
    ## [89] restfulr_0.0.13               later_1.3.0                  
    ## [91] tibble_3.1.6                  GenomicAlignments_1.31.2     
    ## [93] AnnotationDbi_1.57.1          memoise_2.0.1                
    ## [95] interactiveDisplayBase_1.33.0 ellipsis_0.3.2
