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
-   [Session Info](#session-info)

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
    ##    seqnames    start      end strand          protospacer pam pam_site cut_site
    ## 1     chr12 25209843 25209843      - AAAGAAAAGATGAGCAAAGA TGG 25209843 25209846
    ## 2     chr12 25209896 25209896      + TTCTCGAACTAATGTATAGA AGG 25209896 25209893
    ## 3     chr12 25215438 25215438      - AAATGCATTATAATGTAATC TGG 25215438 25215441
    ## 4     chr12 25215477 25215477      - AGCAAAGAAGAAAAGACTCC TGG 25215477 25215480
    ## 5     chr12 25215477 25215477      + TTTTTAATTTTCACACAGCC AGG 25215477 25215474
    ## 6     chr12 25215520 25215520      + TTTTTTTCAATCTGTATTGT CGG 25215520 25215517
    ## 7     chr12 25215535 25215535      - GGAGGATGCTTTTTATACAT TGG 25215535 25215538
    ## 8     chr12 25215553 25215553      - TTTTACAATGCAGAGAGTGG AGG 25215553 25215556
    ## 9     chr12 25215556 25215556      - GTGTTTTACAATGCAGAGAG TGG 25215556 25215559
    ## 10    chr12 25225615 25225615      - AACATCAGCAAAGACAAGAC AGG 25225615 25225618
    ## 11    chr12 25225644 25225644      + TTTGCTGATGTTTCAATAAA AGG 25225644 25225641
    ## 12    chr12 25225653 25225653      - CAGGACTTAGCAAGAAGTTA TGG 25225653 25225656
    ## 13    chr12 25225672 25225672      - AGTAGACACAAAACAGGCTC AGG 25225672 25225675
    ## 14    chr12 25225678 25225678      - TAGAACAGTAGACACAAAAC AGG 25225678 25225681
    ## 15    chr12 25225701 25225701      + TTTGTGTCTACTGTTCTAGA AGG 25225701 25225698
    ## 16    chr12 25225722 25225722      - GATGTACCTATGGTCCTAGT AGG 25225722 25225725
    ## 17    chr12 25225726 25225726      + AATCACATTTATTTCCTACT AGG 25225726 25225723
    ## 18    chr12 25225732 25225732      - GGACTCTGAAGATGTACCTA TGG 25225732 25225735
    ## 19    chr12 25225734 25225734      + TTATTTCCTACTAGGACCAT AGG 25225734 25225731
    ## 20    chr12 25225753 25225753      - AGAACAAATTAAAAGAGTTA AGG 25225753 25225756
    ## 21    chr12 25225775 25225775      + AACTCTTTTAATTTGTTCTC TGG 25225775 25225772
    ## 22    chr12 25225776 25225776      + ACTCTTTTAATTTGTTCTCT GGG 25225776 25225773
    ## 23    chr12 25227231 25227231      - AGATATTCACCATTATAGGT GGG 25227231 25227234
    ## 24    chr12 25227232 25227232      - AAGATATTCACCATTATAGG TGG 25227232 25227235
    ## 25    chr12 25227235 25227235      - TTGAAGATATTCACCATTAT AGG 25227235 25227238
    ## 26    chr12 25227240 25227240      + CAATTTAAACCCACCTATAA TGG 25227240 25227237
    ## 27    chr12 25227274 25227274      + AAATGATTTAGTATTATTTA TGG 25227274 25227271
    ## 28    chr12 25227296 25227296      - CAGTACATGAGGACTGGGGA GGG 25227296 25227299
    ## 29    chr12 25227297 25227297      - CCAGTACATGAGGACTGGGG AGG 25227297 25227300
    ## 30    chr12 25227300 25227300      - GGACCAGTACATGAGGACTG GGG 25227300 25227303
    ## 31    chr12 25227301 25227301      - GGGACCAGTACATGAGGACT GGG 25227301 25227304
    ## 32    chr12 25227302 25227302      - AGGGACCAGTACATGAGGAC TGG 25227302 25227305
    ## 33    chr12 25227307 25227307      - CAATGAGGGACCAGTACATG AGG 25227307 25227310
    ## 34    chr12 25227315 25227315      + CCTCCCCAGTCCTCATGTAC TGG 25227315 25227312
    ## 35    chr12 25227321 25227321      - AGAGGAGTACAGTGCAATGA GGG 25227321 25227324
    ## 36    chr12 25227322 25227322      - AAGAGGAGTACAGTGCAATG AGG 25227322 25227325
    ## 37    chr12 25227339 25227339      - TCTCGACACAGCAGGTCAAG AGG 25227339 25227342
    ## 38    chr12 25227347 25227347      - TTGGATATTCTCGACACAGC AGG 25227347 25227350
    ## 39    chr12 25227366 25227366      - TGATGGAGAAACCTGTCTCT TGG 25227366 25227369
    ## 40    chr12 25227373 25227373      + GTCGAGAATATCCAAGAGAC AGG 25227373 25227370
    ## 41    chr12 25227383 25227383      - AGGAAGCAAGTAGTAATTGA TGG 25227383 25227386
    ## 42    chr12 25227403 25227403      - TCCCTTCTCAGGATTCCTAC AGG 25227403 25227406
    ## 43    chr12 25227406 25227406      + AATTACTACTTGCTTCCTGT AGG 25227406 25227403
    ## 44    chr12 25227419 25227419      + TTCCTGTAGGAATCCTGAGA AGG 25227419 25227416
    ## 45    chr12 25227420 25227420      + TCCTGTAGGAATCCTGAGAA GGG 25227420 25227417
    ## 46    chr12 25235231 25235231      + TACTGCTTAATAACACCTGT AGG 25235231 25235228
    ## 47    chr12 25245275 25245275      - CGAATATGATCCAACAATAG AGG 25245275 25245278
    ## 48    chr12 25245283 25245283      + ACAAGATTTACCTCTATTGT TGG 25245283 25245280
    ## 49    chr12 25245299 25245299      - GCTAATTCAGAATCATTTTG TGG 25245299 25245302
    ## 50    chr12 25245330 25245330      + CTGAATTAGCTGTATCGTCA AGG 25245330 25245327
    ## 51    chr12 25245343 25245343      - GTAGTTGGAGCTGGTGGCGT AGG 25245343 25245346
    ## 52    chr12 25245349 25245349      - CTTGTGGTAGTTGGAGCTGG TGG 25245349 25245352
    ## 53    chr12 25245352 25245352      - AAACTTGTGGTAGTTGGAGC TGG 25245352 25245355
    ## 54    chr12 25245358 25245358      - GAATATAAACTTGTGGTAGT TGG 25245358 25245361
    ## 55    chr12 25245365 25245365      - AATGACTGAATATAAACTTG TGG 25245365 25245368
    ## 56    chr12 25245392 25245392      + TATATTCAGTCATTTTCAGC AGG 25245392 25245389
    ##       region inRepeats percentGC polyA polyC polyG polyT startingGGGGG n0 n0_c
    ## 1   region_8     FALSE        30  TRUE FALSE FALSE FALSE         FALSE  1    1
    ## 2   region_8     FALSE        30 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 3   region_4     FALSE        20 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 4   region_4     FALSE        40  TRUE FALSE FALSE FALSE         FALSE  1    1
    ## 5   region_4     FALSE        30 FALSE FALSE FALSE  TRUE         FALSE  1    1
    ## 6   region_4     FALSE        20 FALSE FALSE FALSE  TRUE         FALSE  1    1
    ## 7   region_4     FALSE        35 FALSE FALSE FALSE  TRUE         FALSE  1    1
    ## 8   region_4     FALSE        40 FALSE FALSE FALSE  TRUE         FALSE  1    1
    ## 9   region_4     FALSE        40 FALSE FALSE FALSE  TRUE         FALSE  1    1
    ## 10  region_3     FALSE        40 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 11  region_3     FALSE        25 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 12  region_3     FALSE        40 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 13  region_3     FALSE        45  TRUE FALSE FALSE FALSE         FALSE  1    1
    ## 14  region_3     FALSE        35  TRUE FALSE FALSE FALSE         FALSE  1    1
    ## 15  region_3     FALSE        35 FALSE FALSE FALSE FALSE         FALSE  2    1
    ## 16  region_3     FALSE        45 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 17  region_3     FALSE        25 FALSE FALSE FALSE FALSE         FALSE  2    1
    ## 18  region_3     FALSE        45 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 19  region_3     FALSE        35 FALSE FALSE FALSE FALSE         FALSE  2    1
    ## 20  region_3     FALSE        20  TRUE FALSE FALSE FALSE         FALSE  1    1
    ## 21  region_3     FALSE        25 FALSE FALSE FALSE  TRUE         FALSE  1    1
    ## 22  region_3     FALSE        25 FALSE FALSE FALSE  TRUE         FALSE  1    1
    ## 23  region_2     FALSE        30 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 24  region_2     FALSE        30 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 25  region_2     FALSE        25 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 26  region_2     FALSE        30 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 27  region_2     FALSE        10 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 28  region_2     FALSE        55 FALSE FALSE  TRUE FALSE         FALSE  1    1
    ## 29  region_2     FALSE        60 FALSE FALSE  TRUE FALSE         FALSE  1    1
    ## 30  region_2     FALSE        55 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 31  region_2     FALSE        55 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 32  region_2     FALSE        55 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 33  region_2     FALSE        50 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 34  region_2     FALSE        60 FALSE  TRUE FALSE FALSE         FALSE  1    1
    ## 35  region_2     FALSE        45 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 36  region_2     FALSE        45 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 37  region_2     FALSE        55 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 38  region_2     FALSE        45 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 39  region_2     FALSE        45 FALSE FALSE FALSE FALSE         FALSE  2    1
    ## 40  region_2     FALSE        45 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 41  region_2     FALSE        35 FALSE FALSE FALSE FALSE         FALSE  2    1
    ## 42  region_2     FALSE        50 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 43  region_2     FALSE        35 FALSE FALSE FALSE FALSE         FALSE  2    1
    ## 44  region_2     FALSE        45 FALSE FALSE FALSE FALSE         FALSE  1    0
    ## 45  region_2     FALSE        45 FALSE FALSE FALSE FALSE         FALSE  1    0
    ## 46 region_12     FALSE        35 FALSE FALSE FALSE FALSE         FALSE  1    0
    ## 47  region_1     FALSE        35 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 48  region_1     FALSE        30 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 49  region_1     FALSE        30 FALSE FALSE FALSE  TRUE         FALSE  2    1
    ## 50  region_1     FALSE        40 FALSE FALSE FALSE FALSE         FALSE  2    1
    ## 51  region_1     FALSE        60 FALSE FALSE FALSE FALSE         FALSE  2    1
    ## 52  region_1     FALSE        55 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 53  region_1     FALSE        45 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 54  region_1     FALSE        30 FALSE FALSE FALSE FALSE         FALSE  1    1
    ## 55  region_1     FALSE        25 FALSE FALSE FALSE FALSE         FALSE  2    1
    ## 56  region_1     FALSE        30 FALSE FALSE FALSE  TRUE         FALSE  2    0
    ##    n0_p n1 n1_c n1_p score_ruleset1 score_crisprscan score_cfd score_mit hasSNP
    ## 1     0  4    0    0    0.043222682     0.5133896687 0.4250273 0.4266001   TRUE
    ## 2     0  1    0    0    0.024432939     0.2351566583 0.5000000 0.5773672  FALSE
    ## 3     0  0    0    0    0.018482579     0.1103547080 1.0000000 1.0000000  FALSE
    ## 4     0  0    0    0    0.083962416    -0.0009864112 1.0000000 1.0000000  FALSE
    ## 5     0  0    0    0    0.091393271     0.3690565395 1.0000000 1.0000000  FALSE
    ## 6     0  4    0    0    0.065259621     0.2724537448 0.5212645 0.8835838  FALSE
    ## 7     0  0    0    0    0.533938710     0.3539000935 1.0000000 1.0000000  FALSE
    ## 8     0  0    0    0    0.219489163     0.5827207543 1.0000000 1.0000000  FALSE
    ## 9     0  0    0    0    0.062261532     0.5705158651 1.0000000 1.0000000  FALSE
    ## 10    0  0    0    0    0.163449646     0.4215891467 1.0000000 1.0000000  FALSE
    ## 11    0  0    0    0    0.042497973     0.4581810540 1.0000000 1.0000000  FALSE
    ## 12    0  1    0    0    0.086608395     0.4398528428 0.6060606 0.7168459  FALSE
    ## 13    0  0    0    0    0.098959664     0.3545034324 1.0000000 1.0000000  FALSE
    ## 14    0  1    0    0    0.150109449     0.3636388955 0.9455253 0.9071057  FALSE
    ## 15    0  0    0    0    0.277445350     0.4138651398 0.5000000 0.5000000  FALSE
    ## 16    0  2    1    0    0.222075141     0.6632815210 0.5668016 0.4701457  FALSE
    ## 17    0  0    0    0    0.045234486     0.2679694809 0.5000000 0.5000000  FALSE
    ## 18    0  1    0    0    0.111680307     0.3051037540 0.8387097 1.0000000  FALSE
    ## 19    0  0    0    0    0.327827733     0.5448407158 0.5000000 0.5000000  FALSE
    ## 20    0  0    0    0    0.053409998     0.1453487123 1.0000000 1.0000000  FALSE
    ## 21    0  1    0    0    0.007536694    -0.0442802663 0.9422633 0.9700236  FALSE
    ## 22    0  0    0    0    0.020817600     0.1456761556 1.0000000 1.0000000  FALSE
    ## 23    0  0    0    0    0.512341024     0.3504188624 1.0000000 1.0000000  FALSE
    ## 24    0  0    0    0    0.258675763     0.3923183291 1.0000000 1.0000000  FALSE
    ## 25    0  0    0    0    0.072213986     0.4102132935 1.0000000 1.0000000  FALSE
    ## 26    0  0    0    0    0.033552012     0.2022816416 1.0000000 1.0000000  FALSE
    ## 27    0  1    0    0    0.009571396     0.2047197767 0.5000000 0.9267841  FALSE
    ## 28    0  0    0    0    0.223314082     0.5848892426 1.0000000 1.0000000  FALSE
    ## 29    0  1    0    0    0.123079697     0.6300169444 0.7000000 1.0000000  FALSE
    ## 30    0  0    0    0    0.277529315     0.6121450231 1.0000000 1.0000000  FALSE
    ## 31    0  1    0    0    0.060000789     0.5150873975 0.5263158 0.6317119  FALSE
    ## 32    0  0    0    0    0.036750591     0.3638034597 1.0000000 1.0000000  FALSE
    ## 33    0  0    0    0    0.011575258     0.4524457160 1.0000000 1.0000000  FALSE
    ## 34    0  0    0    0    0.077578224     0.2317496922 1.0000000 1.0000000  FALSE
    ## 35    0  0    0    0    0.339797580     0.4874341560 1.0000000 1.0000000  FALSE
    ## 36    0  0    0    0    0.089310331     0.4758715650 1.0000000 1.0000000  FALSE
    ## 37    0  0    0    0    0.291064878     0.5096967766 1.0000000 1.0000000  FALSE
    ## 38    0  0    0    0    0.050474158     0.4048914907 1.0000000 1.0000000  FALSE
    ## 39    0  0    0    0    0.052816238     0.1865341473 0.5000000 0.5000000  FALSE
    ## 40    0  1    0    0    0.158064168     0.4437388578 0.5263158 0.6191950  FALSE
    ## 41    0  0    0    0    0.115418540     0.3717978184 0.5000000 0.5000000  FALSE
    ## 42    0  0    0    0    0.592452224     0.3388411782 1.0000000 1.0000000  FALSE
    ## 43    0  0    0    0    0.371080957     0.1340496768 0.5000000 0.5000000  FALSE
    ## 44    0  0    0    0    0.052443950     0.4399618205 1.0000000 1.0000000  FALSE
    ## 45    0  0    0    0    0.101322961     0.5668301886 1.0000000 1.0000000  FALSE
    ## 46    0  0    0    0    0.518885298     0.1249291062 1.0000000 1.0000000  FALSE
    ## 47    0  1    0    0    0.231280485     0.2836962755 0.9541284 0.9545907  FALSE
    ## 48    0  0    0    0    0.048728261     0.1828534831 1.0000000 1.0000000  FALSE
    ## 49    0  0    0    0    0.528400392     0.1566492632 0.5000000 0.5000000  FALSE
    ## 50    0  0    0    0    0.080909495     0.3837528502 0.9350649 0.9350649  FALSE
    ## 51    0  0    0    0    0.134949031     0.6763424113 0.7941176 0.7941176  FALSE
    ## 52    0  1    0    0    0.195890811     0.7928287023 0.5000000 0.5470460  FALSE
    ## 53    0  1    0    0    0.037911704     0.2341157806 0.5000000 0.6199628  FALSE
    ## 54    0  1    0    0    0.028139468     0.3118932129 0.7777778 0.7593014  FALSE
    ## 55    0  2    0    0    0.193029595     0.3250546585 0.4585987 0.4895794  FALSE
    ## 56    0  1    0    0    0.073749862     0.1301660348 0.4426230 0.4648680  FALSE
    ##              spacerID
    ## 1   ENSG00000133703_1
    ## 2   ENSG00000133703_2
    ## 3   ENSG00000133703_3
    ## 4   ENSG00000133703_4
    ## 5   ENSG00000133703_5
    ## 6   ENSG00000133703_6
    ## 7   ENSG00000133703_7
    ## 8   ENSG00000133703_8
    ## 9   ENSG00000133703_9
    ## 10 ENSG00000133703_10
    ## 11 ENSG00000133703_11
    ## 12 ENSG00000133703_12
    ## 13 ENSG00000133703_13
    ## 14 ENSG00000133703_14
    ## 15 ENSG00000133703_15
    ## 16 ENSG00000133703_16
    ## 17 ENSG00000133703_17
    ## 18 ENSG00000133703_18
    ## 19 ENSG00000133703_19
    ## 20 ENSG00000133703_20
    ## 21 ENSG00000133703_21
    ## 22 ENSG00000133703_22
    ## 23 ENSG00000133703_23
    ## 24 ENSG00000133703_24
    ## 25 ENSG00000133703_25
    ## 26 ENSG00000133703_26
    ## 27 ENSG00000133703_27
    ## 28 ENSG00000133703_28
    ## 29 ENSG00000133703_29
    ## 30 ENSG00000133703_30
    ## 31 ENSG00000133703_31
    ## 32 ENSG00000133703_32
    ## 33 ENSG00000133703_33
    ## 34 ENSG00000133703_34
    ## 35 ENSG00000133703_35
    ## 36 ENSG00000133703_36
    ## 37 ENSG00000133703_37
    ## 38 ENSG00000133703_38
    ## 39 ENSG00000133703_39
    ## 40 ENSG00000133703_40
    ## 41 ENSG00000133703_41
    ## 42 ENSG00000133703_42
    ## 43 ENSG00000133703_43
    ## 44 ENSG00000133703_44
    ## 45 ENSG00000133703_45
    ## 46 ENSG00000133703_46
    ## 47 ENSG00000133703_47
    ## 48 ENSG00000133703_48
    ## 49 ENSG00000133703_49
    ## 50 ENSG00000133703_50
    ## 51 ENSG00000133703_51
    ## 52 ENSG00000133703_52
    ## 53 ENSG00000133703_53
    ## 54 ENSG00000133703_54
    ## 55 ENSG00000133703_55
    ## 56 ENSG00000133703_56
    ## 
    ## $alignments
    ##    seqnames     start       end strand               spacer
    ## 1     chr12  25209843  25209843      - AAAGAAAAGATGAGCAAAGA
    ## 2      chr8  68551391  68551391      - AAAGAAAAGATGAGCAAAGA
    ## 3      chr6  54771089  54771089      + AAAGAAAAGATGAGCAAAGA
    ## 4      chr5   4348033   4348033      + AAAGAAAAGATGAGCAAAGA
    ## 5      chr1  48362810  48362810      + AAAGAAAAGATGAGCAAAGA
    ## 6     chr12  25209896  25209896      + TTCTCGAACTAATGTATAGA
    ## 7      chr6  54771050  54771050      - TTCTCGAACTAATGTATAGA
    ## 8     chr12  25215438  25215438      - AAATGCATTATAATGTAATC
    ## 9     chr12  25215477  25215477      - AGCAAAGAAGAAAAGACTCC
    ## 10    chr12  25215477  25215477      + TTTTTAATTTTCACACAGCC
    ## 11    chr12  25215520  25215520      + TTTTTTTCAATCTGTATTGT
    ## 12     chrX  34968543  34968543      - TTTTTTTCAATCTGTATTGT
    ## 13     chr8  16614003  16614003      - TTTTTTTCAATCTGTATTGT
    ## 14     chr2  34098139  34098139      + TTTTTTTCAATCTGTATTGT
    ## 15    chr17  32419546  32419546      - TTTTTTTCAATCTGTATTGT
    ## 16    chr12  25215535  25215535      - GGAGGATGCTTTTTATACAT
    ## 17    chr12  25215553  25215553      - TTTTACAATGCAGAGAGTGG
    ## 18    chr12  25215556  25215556      - GTGTTTTACAATGCAGAGAG
    ## 19    chr12  25225615  25225615      - AACATCAGCAAAGACAAGAC
    ## 20    chr12  25225644  25225644      + TTTGCTGATGTTTCAATAAA
    ## 21    chr12  25225653  25225653      - CAGGACTTAGCAAGAAGTTA
    ## 22     chr6  54771002  54771002      + CAGGACTTAGCAAGAAGTTA
    ## 23    chr12  25225672  25225672      - AGTAGACACAAAACAGGCTC
    ## 24    chr12  25225678  25225678      - TAGAACAGTAGACACAAAAC
    ## 25     chr8 111613834 111613834      - TAGAACAGTAGACACAAAAC
    ## 26    chr12  25225701  25225701      + TTTGTGTCTACTGTTCTAGA
    ## 27     chr6  54770956  54770956      - TTTGTGTCTACTGTTCTAGA
    ## 28    chr12  25225722  25225722      - GATGTACCTATGGTCCTAGT
    ## 29     chr1 114709677 114709677      - GATGTACCTATGGTCCTAGT
    ## 30     chr6  54770935  54770935      + GATGTACCTATGGTCCTAGT
    ## 31    chr12  25225726  25225726      + AATCACATTTATTTCCTACT
    ## 32     chr6  54770931  54770931      - AATCACATTTATTTCCTACT
    ## 33    chr12  25225732  25225732      - GGACTCTGAAGATGTACCTA
    ## 34     chr6  54770925  54770925      + GGACTCTGAAGATGTACCTA
    ## 35    chr12  25225734  25225734      + TTATTTCCTACTAGGACCAT
    ## 36     chr6  54770923  54770923      - TTATTTCCTACTAGGACCAT
    ## 37    chr12  25225753  25225753      - AGAACAAATTAAAAGAGTTA
    ## 38    chr12  25225775  25225775      + AACTCTTTTAATTTGTTCTC
    ## 39    chr11 117609656 117609656      + AACTCTTTTAATTTGTTCTC
    ## 40    chr12  25225776  25225776      + ACTCTTTTAATTTGTTCTCT
    ## 41    chr12  25227231  25227231      - AGATATTCACCATTATAGGT
    ## 42    chr12  25227232  25227232      - AAGATATTCACCATTATAGG
    ## 43    chr12  25227235  25227235      - TTGAAGATATTCACCATTAT
    ## 44    chr12  25227240  25227240      + CAATTTAAACCCACCTATAA
    ## 45    chr12  25227274  25227274      + AAATGATTTAGTATTATTTA
    ## 46     chr6  54770843  54770843      - AAATGATTTAGTATTATTTA
    ## 47    chr12  25227296  25227296      - CAGTACATGAGGACTGGGGA
    ## 48    chr12  25227297  25227297      - CCAGTACATGAGGACTGGGG
    ## 49     chr6  54770808  54770808      + CCAGTACATGAGGACTGGGG
    ## 50    chr12  25227300  25227300      - GGACCAGTACATGAGGACTG
    ## 51    chr12  25227301  25227301      - GGGACCAGTACATGAGGACT
    ## 52     chr6  54770804  54770804      + GGGACCAGTACATGAGGACT
    ## 53    chr12  25227302  25227302      - AGGGACCAGTACATGAGGAC
    ## 54    chr12  25227307  25227307      - CAATGAGGGACCAGTACATG
    ## 55    chr12  25227315  25227315      + CCTCCCCAGTCCTCATGTAC
    ## 56    chr12  25227321  25227321      - AGAGGAGTACAGTGCAATGA
    ## 57    chr12  25227322  25227322      - AAGAGGAGTACAGTGCAATG
    ## 58    chr12  25227339  25227339      - TCTCGACACAGCAGGTCAAG
    ## 59    chr12  25227347  25227347      - TTGGATATTCTCGACACAGC
    ## 60    chr12  25227366  25227366      - TGATGGAGAAACCTGTCTCT
    ## 61     chr6  54770740  54770740      + TGATGGAGAAACCTGTCTCT
    ## 62    chr12  25227373  25227373      + GTCGAGAATATCCAAGAGAC
    ## 63     chr6  54770733  54770733      - GTCGAGAATATCCAAGAGAC
    ## 64    chr12  25227383  25227383      - AGGAAGCAAGTAGTAATTGA
    ## 65     chr6  54770723  54770723      + AGGAAGCAAGTAGTAATTGA
    ## 66    chr12  25227403  25227403      - TCCCTTCTCAGGATTCCTAC
    ## 67    chr12  25227406  25227406      + AATTACTACTTGCTTCCTGT
    ## 68     chr6  54770700  54770700      - AATTACTACTTGCTTCCTGT
    ## 69    chr12  25227419  25227419      + TTCCTGTAGGAATCCTGAGA
    ## 70    chr12  25227420  25227420      + TCCTGTAGGAATCCTGAGAA
    ## 71    chr12  25235231  25235231      + TACTGCTTAATAACACCTGT
    ## 72    chr12  25245275  25245275      - CGAATATGATCCAACAATAG
    ## 73     chr6  54770692  54770692      + CGAATATGATCCAACAATAG
    ## 74    chr12  25245283  25245283      + ACAAGATTTACCTCTATTGT
    ## 75    chr12  25245299  25245299      - GCTAATTCAGAATCATTTTG
    ## 76     chr6  54770668  54770668      + GCTAATTCAGAATCATTTTG
    ## 77    chr12  25245330  25245330      + CTGAATTAGCTGTATCGTCA
    ## 78     chr6  54770637  54770637      - CTGAATTAGCTGTATCGTCA
    ## 79    chr12  25245343  25245343      - GTAGTTGGAGCTGGTGGCGT
    ## 80     chr6  54770624  54770624      + GTAGTTGGAGCTGGTGGCGT
    ## 81    chr12  25245349  25245349      - CTTGTGGTAGTTGGAGCTGG
    ## 82     chr6  54770618  54770618      + CTTGTGGTAGTTGGAGCTGG
    ## 83    chr12  25245352  25245352      - AAACTTGTGGTAGTTGGAGC
    ## 84     chr6  54770615  54770615      + AAACTTGTGGTAGTTGGAGC
    ## 85    chr12  25245358  25245358      - GAATATAAACTTGTGGTAGT
    ## 86     chr6  54770609  54770609      + GAATATAAACTTGTGGTAGT
    ## 87    chr12  25245365  25245365      - AATGACTGAATATAAACTTG
    ## 88     chr6  54770602  54770602      + AATGACTGAATATAAACTTG
    ## 89    chr13  60822020  60822020      - AATGACTGAATATAAACTTG
    ## 90    chr14  57999156  57999156      - AATGACTGAATATAAACTTG
    ## 91    chr12  25245392  25245392      + TATATTCAGTCATTTTCAGC
    ## 92     chr6  54770575  54770575      - TATATTCAGTCATTTTCAGC
    ## 93     chr1 210618123 210618123      - TATATTCAGTCATTTTCAGC
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
    ##     score_cfd  score_mit   spacerID
    ## 1  1.00000000 1.00000000   spacer_1
    ## 2  0.12962963 0.00362963   spacer_1
    ## 3  0.36363636 0.61500000   spacer_1
    ## 4  0.80000000 0.68500000   spacer_1
    ## 5  0.05952381 0.04048611   spacer_1
    ## 6  1.00000000 1.00000000   spacer_3
    ## 7  1.00000000 0.73200000   spacer_3
    ## 8  1.00000000 1.00000000   spacer_5
    ## 9  1.00000000 1.00000000   spacer_6
    ## 10 1.00000000 1.00000000   spacer_7
    ## 11 1.00000000 1.00000000   spacer_8
    ## 12 0.69230769 0.00000000   spacer_8
    ## 13 0.04960317 0.00000000   spacer_8
    ## 14 0.11522634 0.10085185   spacer_8
    ## 15 0.06127451 0.03090278   spacer_8
    ## 16 1.00000000 1.00000000   spacer_9
    ## 17 1.00000000 1.00000000  spacer_10
    ## 18 1.00000000 1.00000000  spacer_11
    ## 19 1.00000000 1.00000000  spacer_12
    ## 20 1.00000000 1.00000000  spacer_14
    ## 21 1.00000000 1.00000000  spacer_16
    ## 22 0.65000000 0.39500000  spacer_16
    ## 23 1.00000000 1.00000000  spacer_18
    ## 24 1.00000000 1.00000000  spacer_20
    ## 25 0.05761317 0.10240741  spacer_20
    ## 26 1.00000000 1.00000000  spacer_22
    ## 27 1.00000000 1.00000000  spacer_22
    ## 28 1.00000000 1.00000000  spacer_24
    ## 29 0.05000000 0.39500000  spacer_24
    ## 30 0.71428571 0.73200000  spacer_24
    ## 31 1.00000000 1.00000000  spacer_26
    ## 32 1.00000000 1.00000000  spacer_26
    ## 33 1.00000000 1.00000000  spacer_28
    ## 34 0.19230769 0.00000000  spacer_28
    ## 35 1.00000000 1.00000000  spacer_30
    ## 36 1.00000000 1.00000000  spacer_30
    ## 37 1.00000000 1.00000000  spacer_32
    ## 38 1.00000000 1.00000000  spacer_34
    ## 39 0.06127451 0.03090278  spacer_34
    ## 40 1.00000000 1.00000000  spacer_36
    ## 41 1.00000000 1.00000000  spacer_38
    ## 42 1.00000000 1.00000000  spacer_40
    ## 43 1.00000000 1.00000000  spacer_42
    ## 44 1.00000000 1.00000000  spacer_44
    ## 45 1.00000000 1.00000000  spacer_46
    ## 46 1.00000000 0.07900000  spacer_46
    ## 47 1.00000000 1.00000000  spacer_48
    ## 48 1.00000000 1.00000000  spacer_50
    ## 49 0.42857143 0.00000000  spacer_50
    ## 50 1.00000000 1.00000000  spacer_52
    ## 51 1.00000000 1.00000000  spacer_54
    ## 52 0.90000000 0.58300000  spacer_54
    ## 53 1.00000000 1.00000000  spacer_56
    ## 54 1.00000000 1.00000000  spacer_58
    ## 55 1.00000000 1.00000000  spacer_60
    ## 56 1.00000000 1.00000000  spacer_62
    ## 57 1.00000000 1.00000000  spacer_64
    ## 58 1.00000000 1.00000000  spacer_66
    ## 59 1.00000000 1.00000000  spacer_68
    ## 60 1.00000000 1.00000000  spacer_70
    ## 61 1.00000000 1.00000000  spacer_70
    ## 62 1.00000000 1.00000000  spacer_72
    ## 63 0.90000000 0.61500000  spacer_72
    ## 64 1.00000000 1.00000000  spacer_74
    ## 65 1.00000000 1.00000000  spacer_74
    ## 66 1.00000000 1.00000000  spacer_76
    ## 67 1.00000000 1.00000000  spacer_78
    ## 68 1.00000000 1.00000000  spacer_78
    ## 69 1.00000000 1.00000000  spacer_80
    ## 70 1.00000000 1.00000000  spacer_82
    ## 71 1.00000000 1.00000000  spacer_84
    ## 72 1.00000000 1.00000000  spacer_85
    ## 73 0.04807692 0.04756944  spacer_85
    ## 74 1.00000000 1.00000000  spacer_89
    ## 75 1.00000000 1.00000000  spacer_93
    ## 76 1.00000000 1.00000000  spacer_93
    ## 77 1.00000000 1.00000000  spacer_97
    ## 78 0.06944444 0.06944444  spacer_97
    ## 79 1.00000000 1.00000000 spacer_101
    ## 80 0.25925926 0.25925926 spacer_101
    ## 81 1.00000000 1.00000000 spacer_105
    ## 82 1.00000000 0.82800000 spacer_105
    ## 83 1.00000000 1.00000000 spacer_109
    ## 84 1.00000000 0.61300000 spacer_109
    ## 85 1.00000000 1.00000000 spacer_113
    ## 86 0.28571429 0.31700000 spacer_113
    ## 87 1.00000000 1.00000000 spacer_117
    ## 88 1.00000000 1.00000000 spacer_117
    ## 89 0.06944444 0.04256944 spacer_117
    ## 90 0.11111111 0.00000000 spacer_117
    ## 91 1.00000000 1.00000000 spacer_121
    ## 92 1.00000000 1.00000000 spacer_121
    ## 93 0.25925926 0.15114815 spacer_121
    ## 
    ## $geneAnnotation
    ##     anchor_site strand gene_symbol         gene_id           tx_id
    ## 1      25209846      -        KRAS ENSG00000133703 ENST00000256078
    ## 2      25209846      -        KRAS ENSG00000133703 ENST00000311936
    ## 3      25209846      -        KRAS ENSG00000133703 ENST00000557334
    ## 4      25209893      +        KRAS ENSG00000133703 ENST00000256078
    ## 5      25209893      +        KRAS ENSG00000133703 ENST00000311936
    ## 6      25209893      +        KRAS ENSG00000133703 ENST00000557334
    ## 7      25215441      -        KRAS ENSG00000133703 ENST00000256078
    ## 8      25215480      -        KRAS ENSG00000133703 ENST00000256078
    ## 9      25215474      +        KRAS ENSG00000133703 ENST00000256078
    ## 10     25215517      +        KRAS ENSG00000133703 ENST00000256078
    ## 11     25215538      -        KRAS ENSG00000133703 ENST00000256078
    ## 12     25215556      -        KRAS ENSG00000133703 ENST00000256078
    ## 13     25215559      -        KRAS ENSG00000133703 ENST00000256078
    ## 14     25225618      -             ENSG00000275197 ENST00000620933
    ## 15     25225618      -        KRAS ENSG00000133703 ENST00000256078
    ## 16     25225618      -        KRAS ENSG00000133703 ENST00000311936
    ## 17     25225641      +             ENSG00000275197 ENST00000620933
    ## 18     25225641      +        KRAS ENSG00000133703 ENST00000256078
    ## 19     25225641      +        KRAS ENSG00000133703 ENST00000311936
    ## 20     25225656      -             ENSG00000275197 ENST00000620933
    ## 21     25225656      -        KRAS ENSG00000133703 ENST00000256078
    ## 22     25225656      -        KRAS ENSG00000133703 ENST00000311936
    ## 23     25225675      -        KRAS ENSG00000133703 ENST00000256078
    ## 24     25225675      -        KRAS ENSG00000133703 ENST00000311936
    ## 25     25225681      -        KRAS ENSG00000133703 ENST00000256078
    ## 26     25225681      -        KRAS ENSG00000133703 ENST00000311936
    ## 27     25225698      +        KRAS ENSG00000133703 ENST00000256078
    ## 28     25225698      +        KRAS ENSG00000133703 ENST00000311936
    ## 29     25225725      -        KRAS ENSG00000133703 ENST00000256078
    ## 30     25225725      -        KRAS ENSG00000133703 ENST00000311936
    ## 31     25225723      +        KRAS ENSG00000133703 ENST00000256078
    ## 32     25225723      +        KRAS ENSG00000133703 ENST00000311936
    ## 33     25225735      -        KRAS ENSG00000133703 ENST00000256078
    ## 34     25225735      -        KRAS ENSG00000133703 ENST00000311936
    ## 35     25225731      +        KRAS ENSG00000133703 ENST00000256078
    ## 36     25225731      +        KRAS ENSG00000133703 ENST00000311936
    ## 37     25225756      -        KRAS ENSG00000133703 ENST00000256078
    ## 38     25225756      -        KRAS ENSG00000133703 ENST00000311936
    ## 39     25225772      +        KRAS ENSG00000133703 ENST00000256078
    ## 40     25225772      +        KRAS ENSG00000133703 ENST00000311936
    ## 41     25225773      +        KRAS ENSG00000133703 ENST00000256078
    ## 42     25225773      +        KRAS ENSG00000133703 ENST00000311936
    ## 43     25227234      -        KRAS ENSG00000133703 ENST00000256078
    ## 44     25227234      -        KRAS ENSG00000133703 ENST00000311936
    ## 45     25227235      -        KRAS ENSG00000133703 ENST00000256078
    ## 46     25227235      -        KRAS ENSG00000133703 ENST00000311936
    ## 47     25227238      -        KRAS ENSG00000133703 ENST00000256078
    ## 48     25227238      -        KRAS ENSG00000133703 ENST00000311936
    ## 49     25227237      +        KRAS ENSG00000133703 ENST00000256078
    ## 50     25227237      +        KRAS ENSG00000133703 ENST00000311936
    ## 51     25227271      +        KRAS ENSG00000133703 ENST00000256078
    ## 52     25227271      +        KRAS ENSG00000133703 ENST00000311936
    ## 53     25227299      -        KRAS ENSG00000133703 ENST00000256078
    ## 54     25227299      -        KRAS ENSG00000133703 ENST00000311936
    ## 55     25227300      -        KRAS ENSG00000133703 ENST00000256078
    ## 56     25227300      -        KRAS ENSG00000133703 ENST00000311936
    ## 57     25227303      -        KRAS ENSG00000133703 ENST00000256078
    ## 58     25227303      -        KRAS ENSG00000133703 ENST00000311936
    ## 59     25227304      -        KRAS ENSG00000133703 ENST00000256078
    ## 60     25227304      -        KRAS ENSG00000133703 ENST00000311936
    ## 61     25227305      -        KRAS ENSG00000133703 ENST00000256078
    ## 62     25227305      -        KRAS ENSG00000133703 ENST00000311936
    ## 63     25227310      -        KRAS ENSG00000133703 ENST00000256078
    ## 64     25227310      -        KRAS ENSG00000133703 ENST00000311936
    ## 65     25227312      +        KRAS ENSG00000133703 ENST00000256078
    ## 66     25227312      +        KRAS ENSG00000133703 ENST00000311936
    ## 67     25227324      -        KRAS ENSG00000133703 ENST00000256078
    ## 68     25227324      -        KRAS ENSG00000133703 ENST00000311936
    ## 69     25227325      -        KRAS ENSG00000133703 ENST00000256078
    ## 70     25227325      -        KRAS ENSG00000133703 ENST00000311936
    ## 71     25227342      -        KRAS ENSG00000133703 ENST00000256078
    ## 72     25227342      -        KRAS ENSG00000133703 ENST00000311936
    ## 73     25227350      -        KRAS ENSG00000133703 ENST00000256078
    ## 74     25227350      -        KRAS ENSG00000133703 ENST00000311936
    ## 75     25227369      -        KRAS ENSG00000133703 ENST00000256078
    ## 76     25227369      -        KRAS ENSG00000133703 ENST00000311936
    ## 77     25227370      +        KRAS ENSG00000133703 ENST00000256078
    ## 78     25227370      +        KRAS ENSG00000133703 ENST00000311936
    ## 79     25227386      -        KRAS ENSG00000133703 ENST00000256078
    ## 80     25227386      -        KRAS ENSG00000133703 ENST00000311936
    ## 81     25227406      -        KRAS ENSG00000133703 ENST00000256078
    ## 82     25227406      -        KRAS ENSG00000133703 ENST00000311936
    ## 83     25227403      +        KRAS ENSG00000133703 ENST00000256078
    ## 84     25227403      +        KRAS ENSG00000133703 ENST00000311936
    ## 85     25245278      -        KRAS ENSG00000133703 ENST00000256078
    ## 86     25245278      -        KRAS ENSG00000133703 ENST00000311936
    ## 87     25245278      -        KRAS ENSG00000133703 ENST00000557334
    ## 88     25245278      -        KRAS ENSG00000133703 ENST00000556131
    ## 89     25245280      +        KRAS ENSG00000133703 ENST00000256078
    ## 90     25245280      +        KRAS ENSG00000133703 ENST00000311936
    ## 91     25245280      +        KRAS ENSG00000133703 ENST00000557334
    ## 92     25245280      +        KRAS ENSG00000133703 ENST00000556131
    ## 93     25245302      -        KRAS ENSG00000133703 ENST00000256078
    ## 94     25245302      -        KRAS ENSG00000133703 ENST00000311936
    ## 95     25245302      -        KRAS ENSG00000133703 ENST00000557334
    ## 96     25245302      -        KRAS ENSG00000133703 ENST00000556131
    ## 97     25245327      +        KRAS ENSG00000133703 ENST00000256078
    ## 98     25245327      +        KRAS ENSG00000133703 ENST00000311936
    ## 99     25245327      +        KRAS ENSG00000133703 ENST00000557334
    ## 100    25245327      +        KRAS ENSG00000133703 ENST00000556131
    ## 101    25245346      -        KRAS ENSG00000133703 ENST00000256078
    ## 102    25245346      -        KRAS ENSG00000133703 ENST00000311936
    ## 103    25245346      -        KRAS ENSG00000133703 ENST00000557334
    ## 104    25245346      -        KRAS ENSG00000133703 ENST00000556131
    ## 105    25245352      -        KRAS ENSG00000133703 ENST00000256078
    ## 106    25245352      -        KRAS ENSG00000133703 ENST00000311936
    ## 107    25245352      -        KRAS ENSG00000133703 ENST00000557334
    ## 108    25245352      -        KRAS ENSG00000133703 ENST00000556131
    ## 109    25245355      -        KRAS ENSG00000133703 ENST00000256078
    ## 110    25245355      -        KRAS ENSG00000133703 ENST00000311936
    ## 111    25245355      -        KRAS ENSG00000133703 ENST00000557334
    ## 112    25245355      -        KRAS ENSG00000133703 ENST00000556131
    ## 113    25245361      -        KRAS ENSG00000133703 ENST00000256078
    ## 114    25245361      -        KRAS ENSG00000133703 ENST00000311936
    ## 115    25245361      -        KRAS ENSG00000133703 ENST00000557334
    ## 116    25245361      -        KRAS ENSG00000133703 ENST00000556131
    ## 117    25245368      -        KRAS ENSG00000133703 ENST00000256078
    ## 118    25245368      -        KRAS ENSG00000133703 ENST00000311936
    ## 119    25245368      -        KRAS ENSG00000133703 ENST00000557334
    ## 120    25245368      -        KRAS ENSG00000133703 ENST00000556131
    ## 121    25245389      +        KRAS ENSG00000133703 ENST00000256078
    ## 122    25245389      +        KRAS ENSG00000133703 ENST00000311936
    ## 123    25245389      +        KRAS ENSG00000133703 ENST00000557334
    ## 124    25245389      +        KRAS ENSG00000133703 ENST00000556131
    ##          protein_id cut_cds cut_fiveUTRs cut_threeUTRs cut_introns percentCDS
    ## 1              <NA>   FALSE        FALSE          TRUE       FALSE         NA
    ## 2   ENSP00000308495    TRUE        FALSE         FALSE       FALSE       91.0
    ## 3   ENSP00000452512    TRUE        FALSE         FALSE       FALSE       77.6
    ## 4              <NA>   FALSE        FALSE          TRUE       FALSE         NA
    ## 5   ENSP00000308495    TRUE        FALSE         FALSE       FALSE       82.7
    ## 6   ENSP00000452512    TRUE        FALSE         FALSE       FALSE       57.0
    ## 7   ENSP00000256078    TRUE        FALSE         FALSE       FALSE      100.0
    ## 8   ENSP00000256078    TRUE        FALSE         FALSE       FALSE       93.2
    ## 9   ENSP00000256078    TRUE        FALSE         FALSE       FALSE       94.2
    ## 10  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       86.7
    ## 11  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       83.0
    ## 12  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       79.8
    ## 13  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       79.3
    ## 14             <NA>   FALSE        FALSE         FALSE       FALSE         NA
    ## 15  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       78.2
    ## 16  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       78.7
    ## 17             <NA>   FALSE        FALSE         FALSE       FALSE         NA
    ## 18  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       74.2
    ## 19  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       74.6
    ## 20             <NA>   FALSE        FALSE         FALSE       FALSE         NA
    ## 21  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       71.6
    ## 22  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       72.0
    ## 23  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       68.2
    ## 24  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       68.6
    ## 25  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       67.2
    ## 26  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       67.5
    ## 27  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       64.2
    ## 28  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       64.6
    ## 29  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       59.5
    ## 30  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       59.8
    ## 31  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       59.8
    ## 32  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       60.1
    ## 33  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       57.7
    ## 34  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       58.0
    ## 35  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       58.4
    ## 36  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       58.7
    ## 37  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       54.0
    ## 38  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       54.3
    ## 39  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       51.2
    ## 40  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       51.5
    ## 41  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       51.1
    ## 42  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       51.3
    ## 43  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       50.9
    ## 44  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       51.1
    ## 45  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       50.7
    ## 46  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       51.0
    ## 47  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       50.2
    ## 48  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       50.4
    ## 49  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       50.4
    ## 50  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       50.6
    ## 51  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       44.4
    ## 52  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       44.6
    ## 53  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       39.5
    ## 54  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       39.7
    ## 55  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       39.3
    ## 56  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       39.5
    ## 57  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       38.8
    ## 58  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       39.0
    ## 59  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       38.6
    ## 60  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       38.8
    ## 61  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       38.4
    ## 62  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       38.6
    ## 63  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       37.5
    ## 64  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       37.7
    ## 65  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       37.2
    ## 66  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       37.4
    ## 67  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       35.1
    ## 68  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       35.3
    ## 69  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       34.9
    ## 70  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       35.1
    ## 71  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       31.9
    ## 72  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       32.1
    ## 73  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       30.5
    ## 74  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       30.7
    ## 75  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       27.2
    ## 76  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       27.3
    ## 77  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       27.0
    ## 78  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       27.2
    ## 79  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       24.2
    ## 80  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       24.3
    ## 81  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       20.7
    ## 82  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       20.8
    ## 83  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       21.2
    ## 84  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       21.3
    ## 85  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       18.8
    ## 86  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       18.9
    ## 87  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       46.9
    ## 88  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       81.1
    ## 89  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       18.4
    ## 90  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       18.5
    ## 91  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       46.1
    ## 92  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       79.5
    ## 93  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       14.6
    ## 94  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       14.6
    ## 95  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       36.4
    ## 96  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       62.9
    ## 97  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       10.2
    ## 98  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       10.2
    ## 99  ENSP00000256078    TRUE        FALSE         FALSE       FALSE       25.4
    ## 100 ENSP00000256078    TRUE        FALSE         FALSE       FALSE       43.9
    ## 101 ENSP00000256078    TRUE        FALSE         FALSE       FALSE        6.8
    ## 102 ENSP00000256078    TRUE        FALSE         FALSE       FALSE        6.9
    ## 103 ENSP00000256078    TRUE        FALSE         FALSE       FALSE       17.1
    ## 104 ENSP00000256078    TRUE        FALSE         FALSE       FALSE       29.5
    ## 105 ENSP00000256078    TRUE        FALSE         FALSE       FALSE        5.8
    ## 106 ENSP00000256078    TRUE        FALSE         FALSE       FALSE        5.8
    ## 107 ENSP00000256078    TRUE        FALSE         FALSE       FALSE       14.5
    ## 108 ENSP00000256078    TRUE        FALSE         FALSE       FALSE       25.0
    ## 109 ENSP00000256078    TRUE        FALSE         FALSE       FALSE        5.3
    ## 110 ENSP00000256078    TRUE        FALSE         FALSE       FALSE        5.3
    ## 111 ENSP00000256078    TRUE        FALSE         FALSE       FALSE       13.2
    ## 112 ENSP00000256078    TRUE        FALSE         FALSE       FALSE       22.7
    ## 113 ENSP00000256078    TRUE        FALSE         FALSE       FALSE        4.2
    ## 114 ENSP00000256078    TRUE        FALSE         FALSE       FALSE        4.2
    ## 115 ENSP00000256078    TRUE        FALSE         FALSE       FALSE       10.5
    ## 116 ENSP00000256078    TRUE        FALSE         FALSE       FALSE       18.2
    ## 117 ENSP00000256078    TRUE        FALSE         FALSE       FALSE        3.0
    ## 118 ENSP00000256078    TRUE        FALSE         FALSE       FALSE        3.0
    ## 119 ENSP00000256078    TRUE        FALSE         FALSE       FALSE        7.5
    ## 120 ENSP00000256078    TRUE        FALSE         FALSE       FALSE       12.9
    ## 121 ENSP00000256078   FALSE         TRUE         FALSE       FALSE         NA
    ## 122 ENSP00000256078   FALSE        FALSE         FALSE       FALSE         NA
    ## 123 ENSP00000256078   FALSE        FALSE         FALSE       FALSE         NA
    ## 124 ENSP00000256078   FALSE        FALSE         FALSE       FALSE         NA
    ##     aminoAcidIndex downtreamATG percentTx nIsoforms totalIsoforms
    ## 1               NA           NA      15.3         3             4
    ## 2              172            1      13.3         3             4
    ## 3               59            1      35.6         3             4
    ## 4               NA           NA      14.4         3             4
    ## 5              157            2      12.4         3             4
    ## 6               44            2      31.1         3             4
    ## 7              190            0      14.0         1             4
    ## 8              177            1      13.3         1             4
    ## 9              179            1      13.4         1             4
    ## 10             165            1      12.6         1             4
    ## 11             158            1      12.2         1             4
    ## 12             152            1      11.9         1             4
    ## 13             151            1      11.8         1             4
    ## 14              NA           NA      91.7         1             1
    ## 15             149            1      11.7         2             4
    ## 16             149            2      12.0         2             4
    ## 17              NA           NA      95.7         1             1
    ## 18             141            1      11.3         2             4
    ## 19             141            2      11.6         2             4
    ## 20              NA           NA      98.4         1             1
    ## 21             136            1      11.0         2             4
    ## 22             136            2      11.3         2             4
    ## 23             130            1      10.7         2             4
    ## 24             130            2      10.9         2             4
    ## 25             128            1      10.6         2             4
    ## 26             128            2      10.8         2             4
    ## 27             122            1      10.2         2             4
    ## 28             122            2      10.5         2             4
    ## 29             113            1       9.7         2             4
    ## 30             113            2      10.0         2             4
    ## 31             114            1       9.8         2             4
    ## 32             114            2      10.0         2             4
    ## 33             110            2       9.6         2             4
    ## 34             110            3       9.8         2             4
    ## 35             111            1       9.6         2             4
    ## 36             111            2       9.9         2             4
    ## 37             103            1       9.2         2             4
    ## 38             103            2       9.4         2             4
    ## 39              98            1       8.9         2             4
    ## 40              98            2       9.1         2             4
    ## 41              97            1       8.9         2             4
    ## 42              97            2       9.1         2             4
    ## 43              97            1       8.8         2             4
    ## 44              97            2       9.0         2             4
    ## 45              97            1       8.8         2             4
    ## 46              97            2       9.0         2             4
    ## 47              96            1       8.8         2             4
    ## 48              96            2       9.0         2             4
    ## 49              96            1       8.8         2             4
    ## 50              96            2       9.0         2             4
    ## 51              85            1       8.2         2             4
    ## 52              85            1       8.3         2             4
    ## 53              75            1       7.6         2             4
    ## 54              75            1       7.8         2             4
    ## 55              75            1       7.6         2             4
    ## 56              75            1       7.8         2             4
    ## 57              74            1       7.6         2             4
    ## 58              74            1       7.7         2             4
    ## 59              74            1       7.6         2             4
    ## 60              74            1       7.7         2             4
    ## 61              73            1       7.5         2             4
    ## 62              73            1       7.7         2             4
    ## 63              72            1       7.4         2             4
    ## 64              72            1       7.6         2             4
    ## 65              71            2       7.4         2             4
    ## 66              71            2       7.6         2             4
    ## 67              67            2       7.2         2             4
    ## 68              67            2       7.4         2             4
    ## 69              67            2       7.2         2             4
    ## 70              67            2       7.3         2             4
    ## 71              61            3       6.9         2             4
    ## 72              61            3       7.0         2             4
    ## 73              58            3       6.7         2             4
    ## 74              58            3       6.9         2             4
    ## 75              52            3       6.4         2             4
    ## 76              52            3       6.5         2             4
    ## 77              52            3       6.3         2             4
    ## 78              52            3       6.5         2             4
    ## 79              46            3       6.0         2             4
    ## 80              46            3       6.2         2             4
    ## 81              40            3       5.7         2             4
    ## 82              40            3       5.8         2             4
    ## 83              41            3       5.7         2             4
    ## 84              41            3       5.9         2             4
    ## 85              36            3       5.5         4             4
    ## 86              36            3       5.6         4             4
    ## 87              36            2      28.9         4             4
    ## 88              36            1      16.7         4             4
    ## 89              35            3       5.4         4             4
    ## 90              35            3       5.6         4             4
    ## 91              35            2      28.7         4             4
    ## 92              35            1      16.6         4             4
    ## 93              28            3       5.0         4             4
    ## 94              28            3       5.1         4             4
    ## 95              28            2      26.6         4             4
    ## 96              28            1      15.3         4             4
    ## 97              20            2       4.6         4             4
    ## 98              20            2       4.7         4             4
    ## 99              20            2      24.2         4             4
    ## 100             20            1      13.9         4             4
    ## 101             13            2       4.2         4             4
    ## 102             13            2       4.3         4             4
    ## 103             13            2      22.4         4             4
    ## 104             13            1      12.7         4             4
    ## 105             11            2       4.1         4             4
    ## 106             11            2       4.2         4             4
    ## 107             11            2      21.9         4             4
    ## 108             11            1      12.4         4             4
    ## 109             10            2       4.1         4             4
    ## 110             10            2       4.1         4             4
    ## 111             10            2      21.6         4             4
    ## 112             10            1      12.2         4             4
    ## 113              8            2       3.9         4             4
    ## 114              8            2       4.0         4             4
    ## 115              8            2      21.0         4             4
    ## 116              8            1      11.9         4             4
    ## 117              6            2       3.8         4             4
    ## 118              6            2       3.9         4             4
    ## 119              6            2      20.3         4             4
    ## 120              6            1      11.4         4             4
    ## 121             NA           NA       3.4         4             4
    ## 122             NA           NA       3.5         4             4
    ## 123             NA           NA      18.3         4             4
    ## 124             NA           NA      10.2         4             4
    ##     percentIsoforms isCommonExon nCodingIsoforms totalCodingIsoforms
    ## 1                75        FALSE               3                   4
    ## 2                75        FALSE               3                   4
    ## 3                75        FALSE               3                   4
    ## 4                75        FALSE               3                   4
    ## 5                75        FALSE               3                   4
    ## 6                75        FALSE               3                   4
    ## 7                25        FALSE               1                   4
    ## 8                25        FALSE               1                   4
    ## 9                25        FALSE               1                   4
    ## 10               25        FALSE               1                   4
    ## 11               25        FALSE               1                   4
    ## 12               25        FALSE               1                   4
    ## 13               25        FALSE               1                   4
    ## 14              100         TRUE               1                  NA
    ## 15               50        FALSE               2                   4
    ## 16               50        FALSE               2                   4
    ## 17              100         TRUE               1                  NA
    ## 18               50        FALSE               2                   4
    ## 19               50        FALSE               2                   4
    ## 20              100         TRUE               1                  NA
    ## 21               50        FALSE               2                   4
    ## 22               50        FALSE               2                   4
    ## 23               50        FALSE               2                   4
    ## 24               50        FALSE               2                   4
    ## 25               50        FALSE               2                   4
    ## 26               50        FALSE               2                   4
    ## 27               50        FALSE               2                   4
    ## 28               50        FALSE               2                   4
    ## 29               50        FALSE               2                   4
    ## 30               50        FALSE               2                   4
    ## 31               50        FALSE               2                   4
    ## 32               50        FALSE               2                   4
    ## 33               50        FALSE               2                   4
    ## 34               50        FALSE               2                   4
    ## 35               50        FALSE               2                   4
    ## 36               50        FALSE               2                   4
    ## 37               50        FALSE               2                   4
    ## 38               50        FALSE               2                   4
    ## 39               50        FALSE               2                   4
    ## 40               50        FALSE               2                   4
    ## 41               50        FALSE               2                   4
    ## 42               50        FALSE               2                   4
    ## 43               50        FALSE               2                   4
    ## 44               50        FALSE               2                   4
    ## 45               50        FALSE               2                   4
    ## 46               50        FALSE               2                   4
    ## 47               50        FALSE               2                   4
    ## 48               50        FALSE               2                   4
    ## 49               50        FALSE               2                   4
    ## 50               50        FALSE               2                   4
    ## 51               50        FALSE               2                   4
    ## 52               50        FALSE               2                   4
    ## 53               50        FALSE               2                   4
    ## 54               50        FALSE               2                   4
    ## 55               50        FALSE               2                   4
    ## 56               50        FALSE               2                   4
    ## 57               50        FALSE               2                   4
    ## 58               50        FALSE               2                   4
    ## 59               50        FALSE               2                   4
    ## 60               50        FALSE               2                   4
    ## 61               50        FALSE               2                   4
    ## 62               50        FALSE               2                   4
    ## 63               50        FALSE               2                   4
    ## 64               50        FALSE               2                   4
    ## 65               50        FALSE               2                   4
    ## 66               50        FALSE               2                   4
    ## 67               50        FALSE               2                   4
    ## 68               50        FALSE               2                   4
    ## 69               50        FALSE               2                   4
    ## 70               50        FALSE               2                   4
    ## 71               50        FALSE               2                   4
    ## 72               50        FALSE               2                   4
    ## 73               50        FALSE               2                   4
    ## 74               50        FALSE               2                   4
    ## 75               50        FALSE               2                   4
    ## 76               50        FALSE               2                   4
    ## 77               50        FALSE               2                   4
    ## 78               50        FALSE               2                   4
    ## 79               50        FALSE               2                   4
    ## 80               50        FALSE               2                   4
    ## 81               50        FALSE               2                   4
    ## 82               50        FALSE               2                   4
    ## 83               50        FALSE               2                   4
    ## 84               50        FALSE               2                   4
    ## 85              100         TRUE               4                   4
    ## 86              100         TRUE               4                   4
    ## 87              100         TRUE               4                   4
    ## 88              100         TRUE               4                   4
    ## 89              100         TRUE               4                   4
    ## 90              100         TRUE               4                   4
    ## 91              100         TRUE               4                   4
    ## 92              100         TRUE               4                   4
    ## 93              100         TRUE               4                   4
    ## 94              100         TRUE               4                   4
    ## 95              100         TRUE               4                   4
    ## 96              100         TRUE               4                   4
    ## 97              100         TRUE               4                   4
    ## 98              100         TRUE               4                   4
    ## 99              100         TRUE               4                   4
    ## 100             100         TRUE               4                   4
    ## 101             100         TRUE               4                   4
    ## 102             100         TRUE               4                   4
    ## 103             100         TRUE               4                   4
    ## 104             100         TRUE               4                   4
    ## 105             100         TRUE               4                   4
    ## 106             100         TRUE               4                   4
    ## 107             100         TRUE               4                   4
    ## 108             100         TRUE               4                   4
    ## 109             100         TRUE               4                   4
    ## 110             100         TRUE               4                   4
    ## 111             100         TRUE               4                   4
    ## 112             100         TRUE               4                   4
    ## 113             100         TRUE               4                   4
    ## 114             100         TRUE               4                   4
    ## 115             100         TRUE               4                   4
    ## 116             100         TRUE               4                   4
    ## 117             100         TRUE               4                   4
    ## 118             100         TRUE               4                   4
    ## 119             100         TRUE               4                   4
    ## 120             100         TRUE               4                   4
    ## 121             100         TRUE               4                   4
    ## 122             100         TRUE               4                   4
    ## 123             100         TRUE               4                   4
    ## 124             100         TRUE               4                   4
    ##     percentCodingIsoforms isCommonCodingExon   spacerID seqnames
    ## 1                      75              FALSE   spacer_1    chr12
    ## 2                      75              FALSE   spacer_1    chr12
    ## 3                      75              FALSE   spacer_1    chr12
    ## 4                      75              FALSE   spacer_3    chr12
    ## 5                      75              FALSE   spacer_3    chr12
    ## 6                      75              FALSE   spacer_3    chr12
    ## 7                      25              FALSE   spacer_5    chr12
    ## 8                      25              FALSE   spacer_6    chr12
    ## 9                      25              FALSE   spacer_7    chr12
    ## 10                     25              FALSE   spacer_8    chr12
    ## 11                     25              FALSE   spacer_9    chr12
    ## 12                     25              FALSE  spacer_10    chr12
    ## 13                     25              FALSE  spacer_11    chr12
    ## 14                     NA                 NA  spacer_12    chr12
    ## 15                     50              FALSE  spacer_12    chr12
    ## 16                     50              FALSE  spacer_12    chr12
    ## 17                     NA                 NA  spacer_14    chr12
    ## 18                     50              FALSE  spacer_14    chr12
    ## 19                     50              FALSE  spacer_14    chr12
    ## 20                     NA                 NA  spacer_16    chr12
    ## 21                     50              FALSE  spacer_16    chr12
    ## 22                     50              FALSE  spacer_16    chr12
    ## 23                     50              FALSE  spacer_18    chr12
    ## 24                     50              FALSE  spacer_18    chr12
    ## 25                     50              FALSE  spacer_20    chr12
    ## 26                     50              FALSE  spacer_20    chr12
    ## 27                     50              FALSE  spacer_22    chr12
    ## 28                     50              FALSE  spacer_22    chr12
    ## 29                     50              FALSE  spacer_24    chr12
    ## 30                     50              FALSE  spacer_24    chr12
    ## 31                     50              FALSE  spacer_26    chr12
    ## 32                     50              FALSE  spacer_26    chr12
    ## 33                     50              FALSE  spacer_28    chr12
    ## 34                     50              FALSE  spacer_28    chr12
    ## 35                     50              FALSE  spacer_30    chr12
    ## 36                     50              FALSE  spacer_30    chr12
    ## 37                     50              FALSE  spacer_32    chr12
    ## 38                     50              FALSE  spacer_32    chr12
    ## 39                     50              FALSE  spacer_34    chr12
    ## 40                     50              FALSE  spacer_34    chr12
    ## 41                     50              FALSE  spacer_36    chr12
    ## 42                     50              FALSE  spacer_36    chr12
    ## 43                     50              FALSE  spacer_38    chr12
    ## 44                     50              FALSE  spacer_38    chr12
    ## 45                     50              FALSE  spacer_40    chr12
    ## 46                     50              FALSE  spacer_40    chr12
    ## 47                     50              FALSE  spacer_42    chr12
    ## 48                     50              FALSE  spacer_42    chr12
    ## 49                     50              FALSE  spacer_44    chr12
    ## 50                     50              FALSE  spacer_44    chr12
    ## 51                     50              FALSE  spacer_46    chr12
    ## 52                     50              FALSE  spacer_46    chr12
    ## 53                     50              FALSE  spacer_48    chr12
    ## 54                     50              FALSE  spacer_48    chr12
    ## 55                     50              FALSE  spacer_50    chr12
    ## 56                     50              FALSE  spacer_50    chr12
    ## 57                     50              FALSE  spacer_52    chr12
    ## 58                     50              FALSE  spacer_52    chr12
    ## 59                     50              FALSE  spacer_54    chr12
    ## 60                     50              FALSE  spacer_54    chr12
    ## 61                     50              FALSE  spacer_56    chr12
    ## 62                     50              FALSE  spacer_56    chr12
    ## 63                     50              FALSE  spacer_58    chr12
    ## 64                     50              FALSE  spacer_58    chr12
    ## 65                     50              FALSE  spacer_60    chr12
    ## 66                     50              FALSE  spacer_60    chr12
    ## 67                     50              FALSE  spacer_62    chr12
    ## 68                     50              FALSE  spacer_62    chr12
    ## 69                     50              FALSE  spacer_64    chr12
    ## 70                     50              FALSE  spacer_64    chr12
    ## 71                     50              FALSE  spacer_66    chr12
    ## 72                     50              FALSE  spacer_66    chr12
    ## 73                     50              FALSE  spacer_68    chr12
    ## 74                     50              FALSE  spacer_68    chr12
    ## 75                     50              FALSE  spacer_70    chr12
    ## 76                     50              FALSE  spacer_70    chr12
    ## 77                     50              FALSE  spacer_72    chr12
    ## 78                     50              FALSE  spacer_72    chr12
    ## 79                     50              FALSE  spacer_74    chr12
    ## 80                     50              FALSE  spacer_74    chr12
    ## 81                     50              FALSE  spacer_76    chr12
    ## 82                     50              FALSE  spacer_76    chr12
    ## 83                     50              FALSE  spacer_78    chr12
    ## 84                     50              FALSE  spacer_78    chr12
    ## 85                    100               TRUE  spacer_85    chr12
    ## 86                    100               TRUE  spacer_85    chr12
    ## 87                    100               TRUE  spacer_85    chr12
    ## 88                    100               TRUE  spacer_85    chr12
    ## 89                    100               TRUE  spacer_89    chr12
    ## 90                    100               TRUE  spacer_89    chr12
    ## 91                    100               TRUE  spacer_89    chr12
    ## 92                    100               TRUE  spacer_89    chr12
    ## 93                    100               TRUE  spacer_93    chr12
    ## 94                    100               TRUE  spacer_93    chr12
    ## 95                    100               TRUE  spacer_93    chr12
    ## 96                    100               TRUE  spacer_93    chr12
    ## 97                    100               TRUE  spacer_97    chr12
    ## 98                    100               TRUE  spacer_97    chr12
    ## 99                    100               TRUE  spacer_97    chr12
    ## 100                   100               TRUE  spacer_97    chr12
    ## 101                   100               TRUE spacer_101    chr12
    ## 102                   100               TRUE spacer_101    chr12
    ## 103                   100               TRUE spacer_101    chr12
    ## 104                   100               TRUE spacer_101    chr12
    ## 105                   100               TRUE spacer_105    chr12
    ## 106                   100               TRUE spacer_105    chr12
    ## 107                   100               TRUE spacer_105    chr12
    ## 108                   100               TRUE spacer_105    chr12
    ## 109                   100               TRUE spacer_109    chr12
    ## 110                   100               TRUE spacer_109    chr12
    ## 111                   100               TRUE spacer_109    chr12
    ## 112                   100               TRUE spacer_109    chr12
    ## 113                   100               TRUE spacer_113    chr12
    ## 114                   100               TRUE spacer_113    chr12
    ## 115                   100               TRUE spacer_113    chr12
    ## 116                   100               TRUE spacer_113    chr12
    ## 117                   100               TRUE spacer_117    chr12
    ## 118                   100               TRUE spacer_117    chr12
    ## 119                   100               TRUE spacer_117    chr12
    ## 120                   100               TRUE spacer_117    chr12
    ## 121                   100               TRUE spacer_121    chr12
    ## 122                   100               TRUE spacer_121    chr12
    ## 123                   100               TRUE spacer_121    chr12
    ## 124                   100               TRUE spacer_121    chr12
    ## 
    ## $enzymeAnnotation
    ##    EcoRI  KpnI BsmBI  BsaI  BbsI  PacI  MluI   spacerID
    ## 1  FALSE FALSE FALSE FALSE FALSE FALSE FALSE   spacer_1
    ## 2  FALSE FALSE FALSE FALSE FALSE FALSE FALSE   spacer_3
    ## 3  FALSE FALSE FALSE FALSE FALSE FALSE FALSE   spacer_5
    ## 4  FALSE FALSE FALSE FALSE FALSE FALSE FALSE   spacer_6
    ## 5  FALSE FALSE FALSE FALSE FALSE FALSE FALSE   spacer_7
    ## 6  FALSE FALSE FALSE FALSE FALSE FALSE FALSE   spacer_8
    ## 7  FALSE FALSE FALSE FALSE FALSE FALSE FALSE   spacer_9
    ## 8  FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_10
    ## 9  FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_11
    ## 10 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_12
    ## 11 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_14
    ## 12 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_16
    ## 13 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_18
    ## 14 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_20
    ## 15 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_22
    ## 16 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_24
    ## 17 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_26
    ## 18 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_28
    ## 19 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_30
    ## 20 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_32
    ## 21 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_34
    ## 22 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_36
    ## 23 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_38
    ## 24 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_40
    ## 25 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_42
    ## 26 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_44
    ## 27 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_46
    ## 28 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_48
    ## 29 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_50
    ## 30 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_52
    ## 31 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_54
    ## 32 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_56
    ## 33 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_58
    ## 34 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_60
    ## 35 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_62
    ## 36 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_64
    ## 37 FALSE FALSE  TRUE FALSE FALSE FALSE FALSE  spacer_66
    ## 38 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_68
    ## 39 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_70
    ## 40 FALSE FALSE  TRUE FALSE FALSE FALSE FALSE  spacer_72
    ## 41 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_74
    ## 42 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_76
    ## 43 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_78
    ## 44 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_80
    ## 45 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_82
    ## 46 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_84
    ## 47 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_85
    ## 48 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_89
    ## 49 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_93
    ## 50 FALSE FALSE FALSE FALSE FALSE FALSE FALSE  spacer_97
    ## 51 FALSE FALSE FALSE FALSE FALSE FALSE FALSE spacer_101
    ## 52 FALSE FALSE FALSE FALSE FALSE FALSE FALSE spacer_105
    ## 53 FALSE FALSE FALSE FALSE FALSE FALSE FALSE spacer_109
    ## 54 FALSE FALSE FALSE FALSE FALSE FALSE FALSE spacer_113
    ## 55 FALSE FALSE FALSE FALSE FALSE FALSE FALSE spacer_117
    ## 56 FALSE FALSE FALSE FALSE FALSE FALSE FALSE spacer_121
    ## 
    ## $snps
    ##          rs  rs_site rs_site_rel allele_ref allele_minor MAF_1000G MAF_TOPMED
    ## 1 rs1137282 25209843           0          A            G    0.1755    0.19671
    ##   type length spacerID
    ## 1  snp      1 spacer_1

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
    ##  [1] crisprScore_1.1.6                 crisprScoreData_1.1.3            
    ##  [3] ExperimentHub_2.3.5               AnnotationHub_3.3.9              
    ##  [5] BiocFileCache_2.3.4               dbplyr_2.1.1                     
    ##  [7] BSgenome.Hsapiens.UCSC.hg38_1.4.4 BSgenome_1.63.5                  
    ##  [9] rtracklayer_1.55.4                Biostrings_2.63.2                
    ## [11] XVector_0.35.0                    GenomicRanges_1.47.6             
    ## [13] GenomeInfoDb_1.31.6               IRanges_2.29.1                   
    ## [15] S4Vectors_0.33.11                 BiocGenerics_0.41.2              
    ## [17] crisprDesignData_0.99.7           crisprDesign_0.99.96             
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
