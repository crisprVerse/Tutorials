Using crisprDesign to design gRNAs for CRISPRko with the SpCas9 nuclease
================

-   [Introduction](#introduction)
-   [Some terminology before we get
    started](#some-terminology-before-we-get-started)
-   [End-to-end gRNA design workflow](#end-to-end-grna-design-workflow)
    -   [Nuclease specification](#nuclease-specification)
    -   [Specification of the target DNA sequence (KRAS
        CDS)](#specification-of-the-target-dna-sequence-kras-cds)
    -   [Finding spacer sequences targeting
        KRAS](#finding-spacer-sequences-targeting-kras)
    -   [Characterizing gRNA spacer
        sequences](#characterizing-grna-spacer-sequences)
    -   [Off-target search with bowtie](#off-target-search-with-bowtie)
    -   [Removing repeat elements](#removing-repeat-elements)
    -   [Off-target scoring (MIT and CFD specificity
        scores)](#off-target-scoring-mit-and-cfd-specificity-scores)
    -   [On-target scoring (gRNA
        efficiency)](#on-target-scoring-grna-efficiency)
    -   [Restriction enzymes](#restriction-enzymes)
    -   [Gene annotation](#gene-annotation)
    -   [TSS annotation](#tss-annotation)
    -   [SNP annotation](#snp-annotation)
    -   [Filtering and ranking gRNAs](#filtering-and-ranking-grnas)
-   [Session Info](#session-info)

Authors: Jean-Philippe Fortin, Luke Hoberecht

Date: 01 August, 2022

# Introduction

In this tutorial, we illustrate the main functionalities of
`crisprDesign`, the central package of the `crisprVerse` ecosystem, by
designing CRISPR/Cas9 gRNAs targeting the coding sequence of the human
gene KRAS. Most steps described in the tutorial are applicable to any
genomic target.

# Some terminology before we get started

Before we start designing gRNAs, we first introduce some terminology
that will be useful throughout this and subsequent tutorials. CRISPR
nucleases require two binding components for cleavage. First, the
nuclease needs to recognize a constant nucleotide motif in the target
DNA called the protospacer adjacent motif (PAM) sequence. Second, the
gRNA, which guides the nuclease to the target sequence, needs to bind to
a complementary sequence adjacent to the PAM sequence, called the
**protospacer** sequence. The latter can be thought of as a variable
binding motif that can be specified by designing corresponding gRNA
sequences.

The **spacer** sequence is used in the gRNA construct to guide the
CRISPR nuclease to the target **protospacer** sequence in the host
genome. While a gRNA spacer sequence may not always uniquely target the
host genome (i.e. it may map to multiple protospacers in the host
genome), we can, for a given reference genome, uniquely identify a
protospacer sequence with a combination of 3 attributes:

-   `chr`: chromosome name
-   `strand`: forward (+) or reverse (-)
-   `pam_site`: genomic coordinate of the first nucleotide of the
    nuclease-specific PAM sequence; for SpCas9 this is the “N” in the
    NGG PAM sequence

For CRISPRko applications, we use an additional genomic coordinate,
called `cut_site`, to represent where the double-stranded break (DSB)
occurs. For SpCas9, the cut site (blunt-ended dsDNA break) is located
4nt upstream of the pam_site (PAM-proximal editing).

# End-to-end gRNA design workflow

We first start by loading the package as usual:

``` r
library(crisprDesign)
```

## Nuclease specification

We load the `SpCas9` nuclease object from the `crisprBase` package (see
the `crisprBase` [vignette](https://github.com/Jfortin1/crisprBase) for
instructions on how to create or load alternative nucleases):

``` r
library(crisprBase)
data(SpCas9, package="crisprBase")
SpCas9
```

    ## Class: CrisprNuclease
    ##   Name: SpCas9
    ##   Target type: DNA
    ##   Metadata: list of length 1
    ##   PAMs: NGG, NAG, NGA
    ##   Weights: 1, 0.2593, 0.0694
    ##   Spacer length: 20
    ##   PAM side: 3prime
    ##     Distance from PAM: 0
    ##   Prototype protospacers: 5'--SSSSSSSSSSSSSSSSSSSS[NGG]--3', 5'--SSSSSSSSSSSSSSSSSSSS[NAG]--3', 5'--SSSSSSSSSSSSSSSSSSSS[NGA]--3'

The three motifs (NGG, NAG and NGA) represent the recognized PAM
sequences by SpCas9, and the weights indicate a recognition score. The
canonical PAM sequence NGG is fully recognized (weight of 1), while the
two non-canonical PAM sequences NAG and NGA are much less tolerated.

The spacer sequence is located on the 5-prime end with respect to the
PAM sequence, and the default spacer sequence length is 20 nucleotides.
If necessary, one can change the spacer length using the function
`crisprBase::spacerLength`. We can inspect the protospacer construct by
using `prototypeSequence`:

``` r
prototypeSequence(SpCas9)
```

    ## [1] "5'--SSSSSSSSSSSSSSSSSSSS[NGG]--3'"

## Specification of the target DNA sequence (KRAS CDS)

Since we aim to design gRNAs that knock out the human KRAS gene, we
first need to retrieve the DNA sequence of the coding region (CDS) of
KRAS. We show in this
[tutorial](https://github.com/crisprVerse/Tutorials/tree/master/Building_Gene_Annotation)
how to build convenient gene model objects that allows to quickly access
gene-specific sequences. Here, we obtain from `crisprDesignData` a
`GRangesList` object that defines the genomic coordinates (in hg38
coordinates) of coding genes in the human genome:

``` r
data(txdb_human, package="crisprDesignData")
```

The `queryTxObject` function allows us to query this object for a
specific gene and feature. Here, we obtain a `GRanges` object containing
the CDS coordinates of KRAS:

``` r
gr <- queryTxObject(txObject=txdb_human,
                    featureType="cds",
                    queryColumn="gene_symbol",
                    queryValue="KRAS")
```

To simplify our design, we will only consider exons that constitute the
primary transcript of the gene (transcript ID: ENST00000311936).

``` r
gr <- gr[gr$tx_id == "ENST00000311936"]
```

Optionally, we could also adjust the arguments in our call to
`queryTxObject` to retrieve those transcript-specific coordinates:

``` r
gr <- queryTxObject(txObject=txObject,
               featureType="cds",
               queryColumn="tx_id",
               queryValue="ENST00000311936")
```

## Finding spacer sequences targeting KRAS

`findSpacers` is the main function of `crisprDesign` for obtaining all
possible spacer sequences that target protospacers located in our target
DNA sequence(s). If a `GRanges` object is provided as input, a
`BSgenome` object (an object that contains sequences of a reference
genome) must be provided as well:

``` r
library(BSgenome.Hsapiens.UCSC.hg38)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38
guideSet <- findSpacers(gr,
                        bsgenome=bsgenome,
                        crisprNuclease=SpCas9)
guideSet
```

    ## GuideSet object with 45 ranges and 5 metadata columns:
    ##             seqnames    ranges strand |          protospacer            pam
    ##                <Rle> <IRanges>  <Rle> |       <DNAStringSet> <DNAStringSet>
    ##    spacer_1    chr12  25209843      - | AAAGAAAAGATGAGCAAAGA            TGG
    ##    spacer_2    chr12  25209896      + | TTCTCGAACTAATGTATAGA            AGG
    ##    spacer_3    chr12  25225615      - | AACATCAGCAAAGACAAGAC            AGG
    ##    spacer_4    chr12  25225644      + | TTTGCTGATGTTTCAATAAA            AGG
    ##    spacer_5    chr12  25225653      - | CAGGACTTAGCAAGAAGTTA            TGG
    ##         ...      ...       ...    ... .                  ...            ...
    ##   spacer_41    chr12  25245343      - | GTAGTTGGAGCTGGTGGCGT            AGG
    ##   spacer_42    chr12  25245349      - | CTTGTGGTAGTTGGAGCTGG            TGG
    ##   spacer_43    chr12  25245352      - | AAACTTGTGGTAGTTGGAGC            TGG
    ##   spacer_44    chr12  25245358      - | GAATATAAACTTGTGGTAGT            TGG
    ##   spacer_45    chr12  25245365      - | AATGACTGAATATAAACTTG            TGG
    ##              pam_site  cut_site      region
    ##             <numeric> <numeric> <character>
    ##    spacer_1  25209843  25209846    region_8
    ##    spacer_2  25209896  25209893    region_8
    ##    spacer_3  25225615  25225618    region_7
    ##    spacer_4  25225644  25225641    region_7
    ##    spacer_5  25225653  25225656    region_7
    ##         ...       ...       ...         ...
    ##   spacer_41  25245343  25245346    region_5
    ##   spacer_42  25245349  25245352    region_5
    ##   spacer_43  25245352  25245355    region_5
    ##   spacer_44  25245358  25245361    region_5
    ##   spacer_45  25245365  25245368    region_5
    ##   -------
    ##   seqinfo: 640 sequences (1 circular) from hg38 genome
    ##   crisprNuclease: SpCas9

This function returns a `GuideSet` object that stores the genomic
coordinates (PAM sites) for all spacer sequences found in the regions
provided by `gr`. The `GuideSet` object is an extension of a
`GenomicRanges` object that stores additional information about gRNAs.

There are several accessor functions we can use to extract information
about the spacer sequences in `guideSet`, and here are a few examples
with their corresponding outputs:

``` r
spacers(guideSet)
```

    ## DNAStringSet object of length 45:
    ##      width seq                                              names               
    ##  [1]    20 AAAGAAAAGATGAGCAAAGA                             spacer_1
    ##  [2]    20 TTCTCGAACTAATGTATAGA                             spacer_2
    ##  [3]    20 AACATCAGCAAAGACAAGAC                             spacer_3
    ##  [4]    20 TTTGCTGATGTTTCAATAAA                             spacer_4
    ##  [5]    20 CAGGACTTAGCAAGAAGTTA                             spacer_5
    ##  ...   ... ...
    ## [41]    20 GTAGTTGGAGCTGGTGGCGT                             spacer_41
    ## [42]    20 CTTGTGGTAGTTGGAGCTGG                             spacer_42
    ## [43]    20 AAACTTGTGGTAGTTGGAGC                             spacer_43
    ## [44]    20 GAATATAAACTTGTGGTAGT                             spacer_44
    ## [45]    20 AATGACTGAATATAAACTTG                             spacer_45

``` r
protospacers(guideSet)
```

    ## DNAStringSet object of length 45:
    ##      width seq                                              names               
    ##  [1]    20 AAAGAAAAGATGAGCAAAGA                             spacer_1
    ##  [2]    20 TTCTCGAACTAATGTATAGA                             spacer_2
    ##  [3]    20 AACATCAGCAAAGACAAGAC                             spacer_3
    ##  [4]    20 TTTGCTGATGTTTCAATAAA                             spacer_4
    ##  [5]    20 CAGGACTTAGCAAGAAGTTA                             spacer_5
    ##  ...   ... ...
    ## [41]    20 GTAGTTGGAGCTGGTGGCGT                             spacer_41
    ## [42]    20 CTTGTGGTAGTTGGAGCTGG                             spacer_42
    ## [43]    20 AAACTTGTGGTAGTTGGAGC                             spacer_43
    ## [44]    20 GAATATAAACTTGTGGTAGT                             spacer_44
    ## [45]    20 AATGACTGAATATAAACTTG                             spacer_45

``` r
pams(guideSet)
```

    ## DNAStringSet object of length 45:
    ##      width seq                                              names               
    ##  [1]     3 TGG                                              spacer_1
    ##  [2]     3 AGG                                              spacer_2
    ##  [3]     3 AGG                                              spacer_3
    ##  [4]     3 AGG                                              spacer_4
    ##  [5]     3 TGG                                              spacer_5
    ##  ...   ... ...
    ## [41]     3 AGG                                              spacer_41
    ## [42]     3 TGG                                              spacer_42
    ## [43]     3 TGG                                              spacer_43
    ## [44]     3 TGG                                              spacer_44
    ## [45]     3 TGG                                              spacer_45

``` r
head(pamSites(guideSet))
```

    ## spacer_1 spacer_2 spacer_3 spacer_4 spacer_5 spacer_6 
    ## 25209843 25209896 25225615 25225644 25225653 25225672

``` r
head(cutSites(guideSet))
```

    ## spacer_1 spacer_2 spacer_3 spacer_4 spacer_5 spacer_6 
    ## 25209846 25209893 25225618 25225641 25225656 25225675

## Characterizing gRNA spacer sequences

There are specific spacer sequence features, independent of the genomic
context of the protospacer sequence, that can reduce or even eliminate
gRNA activity:

-   **Poly-T stretches**: four or more consecutive T nucleotides in the
    spacer sequence may act as a transcriptional termination signal for
    the U6 promoter.
-   **Self-complementarity**: complementary sites with the gRNA backbone
    can compete with the targeted genomic sequence.
-   **Percent GC**: gRNAs with GC content between 20% and 80% are
    preferred.

Use the function `addSequenceFeatures` to evaluate the spacer sequences
with respect to these characteristics and add the results to the
`GuideSet` object:

``` r
guideSet <- addSequenceFeatures(guideSet)
head(guideSet)
```

    ## GuideSet object with 6 ranges and 11 metadata columns:
    ##            seqnames    ranges strand |          protospacer            pam
    ##               <Rle> <IRanges>  <Rle> |       <DNAStringSet> <DNAStringSet>
    ##   spacer_1    chr12  25209843      - | AAAGAAAAGATGAGCAAAGA            TGG
    ##   spacer_2    chr12  25209896      + | TTCTCGAACTAATGTATAGA            AGG
    ##   spacer_3    chr12  25225615      - | AACATCAGCAAAGACAAGAC            AGG
    ##   spacer_4    chr12  25225644      + | TTTGCTGATGTTTCAATAAA            AGG
    ##   spacer_5    chr12  25225653      - | CAGGACTTAGCAAGAAGTTA            TGG
    ##   spacer_6    chr12  25225672      - | AGTAGACACAAAACAGGCTC            AGG
    ##             pam_site  cut_site      region percentGC     polyA     polyC
    ##            <numeric> <numeric> <character> <numeric> <logical> <logical>
    ##   spacer_1  25209843  25209846    region_8        30      TRUE     FALSE
    ##   spacer_2  25209896  25209893    region_8        30     FALSE     FALSE
    ##   spacer_3  25225615  25225618    region_7        40     FALSE     FALSE
    ##   spacer_4  25225644  25225641    region_7        25     FALSE     FALSE
    ##   spacer_5  25225653  25225656    region_7        40     FALSE     FALSE
    ##   spacer_6  25225672  25225675    region_7        45      TRUE     FALSE
    ##                polyG     polyT startingGGGGG
    ##            <logical> <logical>     <logical>
    ##   spacer_1     FALSE     FALSE         FALSE
    ##   spacer_2     FALSE     FALSE         FALSE
    ##   spacer_3     FALSE     FALSE         FALSE
    ##   spacer_4     FALSE     FALSE         FALSE
    ##   spacer_5     FALSE     FALSE         FALSE
    ##   spacer_6     FALSE     FALSE         FALSE
    ##   -------
    ##   seqinfo: 640 sequences (1 circular) from hg38 genome
    ##   crisprNuclease: SpCas9

## Off-target search with bowtie

In order to select gRNAs that are most specific to our target of
interest, it is important to avoid gRNAs that target additional loci in
the genome with either perfect sequence complementarity (multiple
on-targets), or imperfect complementarity through tolerated mismatches
(off-targets). As the SpCas9 nuclease can tolerate mismatches between
the gRNA spacer sequence (RNA) and the protospacer sequence (DNA), it is
necessary to characterize off-targets to minimize the introduction of
double-stranded breaks (DSBs) beyond our intended target.

The `addSpacerAlignments` function appends a list of putative on- and
off-targets to a `GuideSet` object using one of three methods. The first
method uses the fast aligner
[bowtie](http://bowtie-bio.sourceforge.net/index.shtml)
\[@langmead2009bowtie\] via the `crisprBowtie` package to map spacer
sequences to a specified reference genome. This can be done by
specifying `aligner="bowtie` and providing a path to a bowtie index file
to `aligner_index` in `addSpacerAlignments`.

We can control the alignment parameters and output with several function
arguments.

-   `n_mismatches` sets the maximum number of permitted gRNA:DNA
    mismatches (up to 3 mismatches).
-   `n_max_alignments` specifies the maximum number of alignments for a
    given gRNA spacer sequence (1000 by default).
-   `all_alignments`, when set to `TRUE`, overrules the
    `n_max_alignments` and returns all possible alignments.
-   `canonical` filters out protospacer sequences that do not have a
    canonical PAM sequence when `TRUE`.

Let’s search for on- and off-targets having up to 1 mismatch using
bowtie. To use bowtie, we need to specify a bowtie index for the human
genome:

``` r
# Path of the hg38 bowtie index on my personal laptop:
bowtie_index <- "/Users/fortinj2/crisprIndices/bowtie/hg38/hg38"
```

For instructions on how to build a Bowtie index from a given reference
genome, see the [genome index
tutorial](https://github.com/crisprVerse/Tutorials/tree/master/Building_Genome_Indices)
or the [crisprBowtie page](https://github.com/Jfortin1/crisprBowtie) .

We will also specify the gene model object `txdb_human` from
`crisprDesignData` described above for `txObject` argument, which is
needed for the function to annotate genomic alignments with genic
context. This is useful for identifying potentially more problematic
off-targets, such as those located in the CDS of another gene, for
instance.

``` r
guideSet <- addSpacerAlignments(guideSet,
                                aligner="bowtie",
                                aligner_index=bowtie_index,
                                bsgenome=BSgenome.Hsapiens.UCSC.hg38,
                                n_mismatches=1,
                                txObject=txdb_human)
```

    ## [runCrisprBowtie] Using BSgenome.Hsapiens.UCSC.hg38 
    ## [runCrisprBowtie] Searching for SpCas9 protospacers

``` r
guideSet
```

    ## GuideSet object with 45 ranges and 16 metadata columns:
    ##             seqnames    ranges strand |          protospacer            pam
    ##                <Rle> <IRanges>  <Rle> |       <DNAStringSet> <DNAStringSet>
    ##    spacer_1    chr12  25209843      - | AAAGAAAAGATGAGCAAAGA            TGG
    ##    spacer_2    chr12  25209896      + | TTCTCGAACTAATGTATAGA            AGG
    ##    spacer_3    chr12  25225615      - | AACATCAGCAAAGACAAGAC            AGG
    ##    spacer_4    chr12  25225644      + | TTTGCTGATGTTTCAATAAA            AGG
    ##    spacer_5    chr12  25225653      - | CAGGACTTAGCAAGAAGTTA            TGG
    ##         ...      ...       ...    ... .                  ...            ...
    ##   spacer_41    chr12  25245343      - | GTAGTTGGAGCTGGTGGCGT            AGG
    ##   spacer_42    chr12  25245349      - | CTTGTGGTAGTTGGAGCTGG            TGG
    ##   spacer_43    chr12  25245352      - | AAACTTGTGGTAGTTGGAGC            TGG
    ##   spacer_44    chr12  25245358      - | GAATATAAACTTGTGGTAGT            TGG
    ##   spacer_45    chr12  25245365      - | AATGACTGAATATAAACTTG            TGG
    ##              pam_site  cut_site      region percentGC     polyA     polyC
    ##             <numeric> <numeric> <character> <numeric> <logical> <logical>
    ##    spacer_1  25209843  25209846    region_8        30      TRUE     FALSE
    ##    spacer_2  25209896  25209893    region_8        30     FALSE     FALSE
    ##    spacer_3  25225615  25225618    region_7        40     FALSE     FALSE
    ##    spacer_4  25225644  25225641    region_7        25     FALSE     FALSE
    ##    spacer_5  25225653  25225656    region_7        40     FALSE     FALSE
    ##         ...       ...       ...         ...       ...       ...       ...
    ##   spacer_41  25245343  25245346    region_5        60     FALSE     FALSE
    ##   spacer_42  25245349  25245352    region_5        55     FALSE     FALSE
    ##   spacer_43  25245352  25245355    region_5        45     FALSE     FALSE
    ##   spacer_44  25245358  25245361    region_5        30     FALSE     FALSE
    ##   spacer_45  25245365  25245368    region_5        25     FALSE     FALSE
    ##                 polyG     polyT startingGGGGG        n0        n1      n0_c
    ##             <logical> <logical>     <logical> <numeric> <numeric> <numeric>
    ##    spacer_1     FALSE     FALSE         FALSE         1         2         1
    ##    spacer_2     FALSE     FALSE         FALSE         1         1         1
    ##    spacer_3     FALSE     FALSE         FALSE         1         0         1
    ##    spacer_4     FALSE     FALSE         FALSE         1         0         1
    ##    spacer_5     FALSE     FALSE         FALSE         1         1         1
    ##         ...       ...       ...           ...       ...       ...       ...
    ##   spacer_41     FALSE     FALSE         FALSE         1         0         1
    ##   spacer_42     FALSE     FALSE         FALSE         1         1         1
    ##   spacer_43     FALSE     FALSE         FALSE         1         1         1
    ##   spacer_44     FALSE     FALSE         FALSE         1         1         1
    ##   spacer_45     FALSE     FALSE         FALSE         2         0         1
    ##                  n1_c                        alignments
    ##             <numeric>                     <GRangesList>
    ##    spacer_1         0  chr12:25209843:-,chr6:54770740:+
    ##    spacer_2         0 chr12:25225672:-,chr12:25245283:+
    ##    spacer_3         0 chr12:25225678:-,chr12:25245299:-
    ##    spacer_4         0  chr12:25225701:+,chr6:54770668:+
    ##    spacer_5         0  chr6:54770956:-,chr12:25245330:+
    ##         ...       ...                               ...
    ##   spacer_41         0   chr6:54771050:-,chr6:54770723:+
    ##   spacer_42         0 chr12:25225615:-,chr12:25227403:-
    ##   spacer_43         0 chr12:25225644:+,chr12:25227406:+
    ##   spacer_44         0  chr12:25225653:-,chr6:54770700:-
    ##   spacer_45         0  chr6:54771002:+,chr12:25245275:-
    ##   -------
    ##   seqinfo: 640 sequences (1 circular) from hg38 genome
    ##   crisprNuclease: SpCas9

Several columns were added to the `GuideSet` object that summarize the
number of on- and off-targets for each spacer sequence, and take genomic
context into account:

-   **n0, n1, n2, n3**: specify the number of alignments with 0, 1, 2
    and 3 mismatches, respectively.
-   **n0_c, n1_c, n2_c, n3_c**: specify the number of alignments in a
    coding region, with 0, 1, 2 and 3 mismatches, respectively.
-   **n0_p, n1_p, n2_p, n3_p**: specify the number of alignments in a
    promoter region of a coding gene, with 0, 1, 2 and 3 mismatches,
    respectively.

Our `guideSet` now has columns of the first two categories, up to 1
mismatch (the value passed to `n_mismatches`); had we also supplied a
`GRanges` of TSS coordinates to the `tssObject` argument, our `guideSet`
would include columns in the last category.

To inspect individual on- and off-targets and their context, one can use
the `alignments` function, which returns a table of all genomic
alignments stored in the `GuideSet` object:

``` r
alignments(guideSet)
```

    ## GRanges object with 67 ranges and 14 metadata columns:
    ##             seqnames    ranges strand |               spacer
    ##                <Rle> <IRanges>  <Rle> |       <DNAStringSet>
    ##    spacer_1    chr12  25209843      - | AAAGAAAAGATGAGCAAAGA
    ##   spacer_32     chr6  54770740      + | TGATGGAGAAACCTGTCTCT
    ##    spacer_6    chr12  25225672      - | AGTAGACACAAAACAGGCTC
    ##   spacer_38    chr12  25245283      + | ACAAGATTTACCTCTATTGT
    ##    spacer_7    chr12  25225678      - | TAGAACAGTAGACACAAAAC
    ##         ...      ...       ...    ... .                  ...
    ##   spacer_36    chr12  25227406      + | AATTACTACTTGCTTCCTGT
    ##    spacer_5    chr12  25225653      - | CAGGACTTAGCAAGAAGTTA
    ##   spacer_36     chr6  54770700      - | AATTACTACTTGCTTCCTGT
    ##    spacer_5     chr6  54771002      + | CAGGACTTAGCAAGAAGTTA
    ##   spacer_37    chr12  25245275      - | CGAATATGATCCAACAATAG
    ##                      protospacer            pam  pam_site n_mismatches
    ##                   <DNAStringSet> <DNAStringSet> <numeric>    <integer>
    ##    spacer_1 AAAGAAAAGATGAGCAAAGA            TGG  25209843            0
    ##   spacer_32 TGATGGAGAAACCTGTCTCT            TGG  54770740            0
    ##    spacer_6 AGTAGACACAAAACAGGCTC            AGG  25225672            0
    ##   spacer_38 ACAAGATTTACCTCTATTGT            TGG  25245283            0
    ##    spacer_7 TAGAACAGTAGACACAAAAC            AGG  25225678            0
    ##         ...                  ...            ...       ...          ...
    ##   spacer_36 AATTACTACTTGCTTCCTGT            AGG  25227406            0
    ##    spacer_5 CAGGACTTAGCAAGAAGTTA            TGG  25225653            0
    ##   spacer_36 AATTACTACTTGCTTCCTGT            AGG  54770700            0
    ##    spacer_5 CAGGACTTAGCAAGGAGTTA            GGG  54771002            1
    ##   spacer_37 CGAATATGATCCAACAATAG            AGG  25245275            0
    ##             canonical  cut_site         cds    fiveUTRs   threeUTRs       exons
    ##             <logical> <numeric> <character> <character> <character> <character>
    ##    spacer_1      TRUE  25209846        KRAS        <NA>        KRAS        KRAS
    ##   spacer_32      TRUE  54770737        <NA>        <NA>        <NA>      KRASP1
    ##    spacer_6      TRUE  25225675        KRAS        <NA>        <NA>        KRAS
    ##   spacer_38      TRUE  25245280        KRAS        <NA>        <NA>        KRAS
    ##    spacer_7      TRUE  25225681        KRAS        <NA>        <NA>        KRAS
    ##         ...       ...       ...         ...         ...         ...         ...
    ##   spacer_36      TRUE  25227403        KRAS        <NA>        <NA>        KRAS
    ##    spacer_5      TRUE  25225656        KRAS        <NA>        <NA>       ;KRAS
    ##   spacer_36      TRUE  54770703        <NA>        <NA>        <NA>      KRASP1
    ##    spacer_5      TRUE  54770999        <NA>        <NA>        <NA>      KRASP1
    ##   spacer_37      TRUE  25245278        KRAS        <NA>        <NA>        KRAS
    ##                 introns  intergenic intergenic_distance
    ##             <character> <character>           <integer>
    ##    spacer_1        <NA>        <NA>                <NA>
    ##   spacer_32        <NA>        <NA>                <NA>
    ##    spacer_6        KRAS        <NA>                <NA>
    ##   spacer_38        <NA>        <NA>                <NA>
    ##    spacer_7        KRAS        <NA>                <NA>
    ##         ...         ...         ...                 ...
    ##   spacer_36        KRAS        <NA>                <NA>
    ##    spacer_5        KRAS        <NA>                <NA>
    ##   spacer_36        <NA>        <NA>                <NA>
    ##    spacer_5        <NA>        <NA>                <NA>
    ##   spacer_37        <NA>        <NA>                <NA>
    ##   -------
    ##   seqinfo: 25 sequences (1 circular) from hg38 genome

Similarly, the functions `onTargets` and `offTargets` return on-target
alignments (no mismatches) and off-target alignments (having at least
one mismatch), respectively. See `?addSpacerAlignments` for more details
about the different options.

We note that gRNAs that align to hundreds of different locations are
highly unspecific and undesirable. This can also cause
`addSpacerAlignments` to be slow. The function
`addSpacerAlignmentsIterative` is an iterative version of
`addSpacerAlignments` that curtails alignment searches for gRNAs having
more hits than the user-defined threshold. See
`?addSpacerAlignmentsIterative` for more information.

## Removing repeat elements

Many promiscuous protospacer sequences occur in repeats or
low-complexity DNA sequences (regions identified by RepeatMasker). These
sequences are usually not of interest due to their low specificity, and
can be easily removed with `removeRepeats`:

``` r
data("gr.repeats.hg38", package="crisprDesignData")
guideSet <- removeRepeats(guideSet,
                          gr.repeats=gr.repeats.hg38)
```

## Off-target scoring (MIT and CFD specificity scores)

After retrieving a list of putative off-targets and on-targets for a
given spacer sequence, we can use `addOffTargetScores` to predict the
likelihood of the nuclease to cut at the off-target locations based on
mismatch tolerance

``` r
guideSet <- addOffTargetScores(guideSet)
guideSet
```

Note that this will only work after calling `addSpacerAlignments`, as it
requires a list of off-targets for each gRNA entry.

## On-target scoring (gRNA efficiency)

`addOnTargetScores` adds scores from on-target efficiency algorithms
specified by the `methods` argument (or all available methods if `NULL`)
available in the R package
[crisprScore](https://github.com/Jfortin1/crisprScore) and appends them
to the `GuideSet`:

``` r
guideSet <- addOnTargetScores(guideSet,
                              methods=c("deephf"))
```

See the [crisprScore page](https://github.com/Jfortin1/crisprScore) for
a full description of the different scores.

## Restriction enzymes

Since the gRNA library synthesis process usually involves restriction
enzymes, it is often necessary to remove gRNAs that contain restriction
sites of specific enzymes. The function `addRestrictionEnzymes` allows
the user to flag gRNAs containing restriction sites for a user-defined
set of enzymes.

``` r
guideSet <- addRestrictionEnzymes(guideSet)
```

By default (that is, when `includeDefault` is `TRUE`), the function adds
annotation for the following commonly used enzymes: EcoRI, KpnI, BsmBI,
BsaI, BbsI, PacI, ISceI and MluI. Additional enzymes can be included by
name via `enzymeNames`, and custom restriction sites can be defined
using the `patterns` argument. It also accepts arguments to specify the
nucleotide sequence that flanks the spacer sequence on the 5’ end
(`flanking5`) and on the 3’ end (`flanking3`) in the lentiviral cassette
used for gRNA delivery. The function effectively searches for
restriction sites in the full sequence:
`[flanking5][spacer][flanking3]`.

One can use the `enzymeAnnotation` function to retrieve the added
annotation:

``` r
head(enzymeAnnotation(guideSet))
```

    ## DataFrame with 6 rows and 7 columns
    ##              EcoRI      KpnI     BsmBI      BsaI      BbsI      PacI      MluI
    ##          <logical> <logical> <logical> <logical> <logical> <logical> <logical>
    ## spacer_1     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE
    ## spacer_2     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE
    ## spacer_3     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE
    ## spacer_4     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE
    ## spacer_5     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE
    ## spacer_6     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE     FALSE

## Gene annotation

The function `addGeneAnnotation` adds transcript- and gene-level context
to gRNAs from a `TxDb`-like object:

``` r
guideSet <- addGeneAnnotation(guideSet,
                              txObject=txdb_human)
```

The gene annotation can be retrieved using the function
`geneAnnotation`:

``` r
geneAnnotation(guideSet)
```

    ## DataFrame with 113 rows and 23 columns
    ##                chr anchor_site   strand gene_symbol         gene_id
    ##           <factor>   <integer> <factor> <character>     <character>
    ## spacer_1     chr12    25209846        -        KRAS ENSG00000133703
    ## spacer_1     chr12    25209846        -        KRAS ENSG00000133703
    ## spacer_1     chr12    25209846        -        KRAS ENSG00000133703
    ## spacer_2     chr12    25209893        +        KRAS ENSG00000133703
    ## spacer_2     chr12    25209893        +        KRAS ENSG00000133703
    ## ...            ...         ...      ...         ...             ...
    ## spacer_44    chr12    25245361        -        KRAS ENSG00000133703
    ## spacer_45    chr12    25245368        -        KRAS ENSG00000133703
    ## spacer_45    chr12    25245368        -        KRAS ENSG00000133703
    ## spacer_45    chr12    25245368        -        KRAS ENSG00000133703
    ## spacer_45    chr12    25245368        -        KRAS ENSG00000133703
    ##                     tx_id      protein_id   cut_cds cut_fiveUTRs cut_threeUTRs
    ##               <character>     <character> <logical>    <logical>     <logical>
    ## spacer_1  ENST00000256078              NA     FALSE        FALSE          TRUE
    ## spacer_1  ENST00000311936 ENSP00000308495      TRUE        FALSE         FALSE
    ## spacer_1  ENST00000557334 ENSP00000452512      TRUE        FALSE         FALSE
    ## spacer_2  ENST00000256078              NA     FALSE        FALSE          TRUE
    ## spacer_2  ENST00000311936 ENSP00000308495      TRUE        FALSE         FALSE
    ## ...                   ...             ...       ...          ...           ...
    ## spacer_44 ENST00000556131 ENSP00000256078      TRUE        FALSE         FALSE
    ## spacer_45 ENST00000256078 ENSP00000256078      TRUE        FALSE         FALSE
    ## spacer_45 ENST00000311936 ENSP00000256078      TRUE        FALSE         FALSE
    ## spacer_45 ENST00000557334 ENSP00000256078      TRUE        FALSE         FALSE
    ## spacer_45 ENST00000556131 ENSP00000256078      TRUE        FALSE         FALSE
    ##           cut_introns percentCDS aminoAcidIndex downtreamATG percentTx
    ##             <logical>  <numeric>      <numeric>    <numeric> <numeric>
    ## spacer_1        FALSE         NA             NA           NA      15.3
    ## spacer_1        FALSE       91.0            172            1      13.3
    ## spacer_1        FALSE       77.6             59            1      35.6
    ## spacer_2        FALSE         NA             NA           NA      14.4
    ## spacer_2        FALSE       82.7            157            2      12.4
    ## ...               ...        ...            ...          ...       ...
    ## spacer_44       FALSE       18.2              8            1      11.9
    ## spacer_45       FALSE        3.0              6            2       3.8
    ## spacer_45       FALSE        3.0              6            2       3.9
    ## spacer_45       FALSE        7.5              6            2      20.3
    ## spacer_45       FALSE       12.9              6            1      11.4
    ##           nIsoforms totalIsoforms percentIsoforms isCommonExon nCodingIsoforms
    ##           <integer>     <numeric>       <numeric>    <logical>       <integer>
    ## spacer_1          3             4              75        FALSE               3
    ## spacer_1          3             4              75        FALSE               3
    ## spacer_1          3             4              75        FALSE               3
    ## spacer_2          3             4              75        FALSE               3
    ## spacer_2          3             4              75        FALSE               3
    ## ...             ...           ...             ...          ...             ...
    ## spacer_44         4             4             100         TRUE               4
    ## spacer_45         4             4             100         TRUE               4
    ## spacer_45         4             4             100         TRUE               4
    ## spacer_45         4             4             100         TRUE               4
    ## spacer_45         4             4             100         TRUE               4
    ##           totalCodingIsoforms percentCodingIsoforms isCommonCodingExon
    ##                     <numeric>             <numeric>          <logical>
    ## spacer_1                    4                    75              FALSE
    ## spacer_1                    4                    75              FALSE
    ## spacer_1                    4                    75              FALSE
    ## spacer_2                    4                    75              FALSE
    ## spacer_2                    4                    75              FALSE
    ## ...                       ...                   ...                ...
    ## spacer_44                   4                   100               TRUE
    ## spacer_45                   4                   100               TRUE
    ## spacer_45                   4                   100               TRUE
    ## spacer_45                   4                   100               TRUE
    ## spacer_45                   4                   100               TRUE

It provides a great deal of information in describing the genomic
location of the protospacer sequences.

-   Ensembl ID columns are provided for all applicable levels:
    `gene_id`, `tx_id`, `protein_id`, `exon_id`.
-   `exon_rank` gives the order of the exon for the transcript; for
    example `"2"` indicates it is the second exon (from the 5’ end) in
    the mature transcript.
-   several columns describe for which gene the the guide sequence
    overlaps the indicated transcript segment: `cut_cds`,
    `cut_fiveUTRs`, `cut_threeUTRs`, `cut_introns`.
-   `percentCDS` and `percentTx` give the location of the `cut_site`
    within the CDS of the transcript and the entire transcript,
    respectively, as a percent from the 5’ end to the 3’ end.
-   `aminoAcidIndex` gives the number of the specific amino acid in the
    protein where the cut is predicted to occur.
-   `downstreamATG` shows how many in-frame ATGs are downstream of the
    `cut_site` (and upstream from the defined percent transcript cutoff,
    `met_cutoff`), indicating a potential alternative translation
    initiation site that may preserve protein function.
-   isoform coverage is described by four columns:
    -   `nIsoforms` gives the number of isoforms of the target gene
        (from `gene_id`) that overlap with the protospacer sequence.
    -   `totalIsoforms` is the number of isoforms for the target gene.
    -   `percentIsoforms` calculates the percentage of isoforms for the
        target gene that overlap with the protospacer sequence
        (`100*nIsoforms/totalIsoforms`).
    -   `isCommonExon` identifies protospacer sequences that overlap
        with all isoforms for the target gene.
-   isoform coverage when exclusively considering the CDS of the target
    gene is similarly described by the `nCodingIsoforms`,
    `totalCodingIsoforms`, `percentCodingIsoforms`, and
    `isCommonCodingExon` columns.
-   `pfam` gives the ID of Pfam domain(s) overlapping the protospacer
    sequence.

## TSS annotation

Similarly, one might want to know which protospacer sequences are
located within promoter regions of known genes:

``` r
data(tssObjectExample, package="crisprDesign")
guideSet <- addTssAnnotation(guideSet,
                             tssObject=tssObjectExample)
tssAnnotation(guideSet)
```

    ## DataFrame with 0 rows and 11 columns

Not surprisingly, as our `GuideSet` targets the CDS of KRAS, none of our
guides overlap a gene promoter region.

## SNP annotation

Common single-nucleotide polymorphisms (SNPs) can change the on-target
and off-target properties of gRNAs by altering the binding. The function
`addSNPAnnotation` annotates gRNAs with respect to a reference database
of SNPs (stored in a VCF file), specified by the `vcf` argument.

VCF files for common SNPs (dbSNPs) can be downloaded from NCBI on the
[dbSNP
website](https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/).
We will use one of those files, after having downloaded it to our local
machine.

``` r
vcf <- "/Users/fortinj2/crisprIndices/snps/dbsnp151.grch38/00-common_all_snps_only.vcf.gz"
guideSet <- addSNPAnnotation(guideSet, vcf=vcf)
snps(guideSet)
```

    ## DataFrame with 1 row and 9 columns
    ##                   rs   rs_site rs_site_rel     allele_ref   allele_minor
    ##          <character> <integer>   <numeric> <DNAStringSet> <DNAStringSet>
    ## spacer_1   rs1137282  25209843           0              A              G
    ##          MAF_1000G MAF_TOPMED        type    length
    ##          <numeric>  <numeric> <character> <integer>
    ## spacer_1    0.1755    0.19671         snp         1

The `rs_site_rel` gives the relative position of the SNP with respect to
the `pam_site`. `allele_ref` and `allele_minor` report the nucleotide of
the reference and minor alleles, respectively. `MAF_1000G` and
`MAF_TOPMED` report the minor allele frequency (MAF) in the 1000Genomes
and TOPMED populations, respectively.

## Filtering and ranking gRNAs

Once our gRNAs are fully annotated we can filter out unwantd gRNAs with
the function `filterSpacers` then rank the best remaining gRNAs with
`rankSpacers`. But first, let’s check which criteria we can filter and
rank on based on the existing annotation in `guideSet` with the
`validCriteria` function.

``` r
validCriteria(guideSet)
```

    ##                attribute valueType takesTxId takesGeneId
    ## 1                  polyA   logical     FALSE       FALSE
    ## 2                  polyC   logical     FALSE       FALSE
    ## 3                  polyG   logical     FALSE       FALSE
    ## 4                  polyT   logical     FALSE       FALSE
    ## 5          startingGGGGG   logical     FALSE       FALSE
    ## 8                cut_cds   logical      TRUE       FALSE
    ## 9           cut_fiveUTRs   logical      TRUE       FALSE
    ## 10         cut_threeUTRs   logical      TRUE       FALSE
    ## 11           cut_introns   logical      TRUE       FALSE
    ## 12          isCommonExon   logical     FALSE        TRUE
    ## 13    isCommonCodingExon   logical     FALSE        TRUE
    ## 14                hasSNP   logical     FALSE       FALSE
    ## 15             inRepeats   logical     FALSE       FALSE
    ## 16                 EcoRI   logical     FALSE       FALSE
    ## 17                  KpnI   logical     FALSE       FALSE
    ## 18                 BsmBI   logical     FALSE       FALSE
    ## 19                  BsaI   logical     FALSE       FALSE
    ## 20                  BbsI   logical     FALSE       FALSE
    ## 21                  PacI   logical     FALSE       FALSE
    ## 22                  MluI   logical     FALSE       FALSE
    ## 23        aminoAcidIndex       asc      TRUE        TRUE
    ## 25            percentCDS       asc      TRUE        TRUE
    ## 26             percentTx       asc      TRUE        TRUE
    ## 27                    n0       asc     FALSE       FALSE
    ## 28                    n1       asc     FALSE       FALSE
    ## 29                  n0_c       asc     FALSE       FALSE
    ## 30                  n1_c       asc     FALSE       FALSE
    ## 31             nIsoforms      desc     FALSE        TRUE
    ## 32       percentIsoforms      desc     FALSE        TRUE
    ## 33       nCodingIsoforms      desc     FALSE        TRUE
    ## 34 percentCodingIsoforms      desc     FALSE        TRUE
    ## 35             score_cfd      desc     FALSE       FALSE
    ## 36             score_mit      desc     FALSE       FALSE
    ## 37          score_deephf      desc     FALSE       FALSE
    ## 38             percentGC    ranged     FALSE       FALSE

As an example, suppose that we only want to keep gRNAs that meet the
following criteria:

-   has percent GC between 20% and 80%
-   does not contain a polyT strech
-   does not have EcoRI or KpnI restriction sites

``` r
filter_criteria <- list(percentGC=c(20, 80),
                        polyT=FALSE,
                        EcoRI=FALSE,
                        KpnI=FALSE)
guideSet <- filterSpacers(guideSet,
                          criteria=filter_criteria)
```

The arguments for `rankSpacers` are the same, but with a subtle
difference for `criteria`: the order of elements defines the priority
for ranking. As an example, suppose we have the following ranking
criteria, in order of importance:

-   few off-targets with one mismatch; preferably no such off-targets
-   high on-target score (DeepHF)
-   targets an exon common to all KRAS isoforms

In setting up the list for `criteria`, we can see from
`validCriteria(guideSet)` above that `n1` takes ascending values,
`score_deephf` takes descending values, and `isCommonExon` requires a
gene ID (Ensembl ID). For integer values, such as `n1`, we recommend
using non-integer values so the division between ranks is unambiguous.
Here’s an example `criteria` we will use to rank our guides:

``` r
rank_criteria <- list(n1=c(0.5, 1.5, 2.5, 5.5),
                      score_deephf=c(0.8, 0.7, 0.6, 0.5, 0.4),
                      isCommonExon=TRUE)
guideSet <- rankSpacers(guideSet,
                        criteria=rank_criteria,
                        geneId="ENSG00000133703") # required for isCommonExon!
guideSet
```

    ## GuideSet object with 41 ranges and 26 metadata columns:
    ##             seqnames    ranges strand |          protospacer            pam
    ##                <Rle> <IRanges>  <Rle> |       <DNAStringSet> <DNAStringSet>
    ##   spacer_29    chr12  25227322      - | AAGAGGAGTACAGTGCAATG            AGG
    ##   spacer_37    chr12  25245275      - | CGAATATGATCCAACAATAG            AGG
    ##   spacer_41    chr12  25245343      - | GTAGTTGGAGCTGGTGGCGT            AGG
    ##   spacer_45    chr12  25245365      - | AATGACTGAATATAAACTTG            TGG
    ##    spacer_3    chr12  25225615      - | AACATCAGCAAAGACAAGAC            AGG
    ##         ...      ...       ...    ... .                  ...            ...
    ##   spacer_43    chr12  25245352      - | AAACTTGTGGTAGTTGGAGC            TGG
    ##   spacer_44    chr12  25245358      - | GAATATAAACTTGTGGTAGT            TGG
    ##    spacer_2    chr12  25209896      + | TTCTCGAACTAATGTATAGA            AGG
    ##    spacer_9    chr12  25225722      - | GATGTACCTATGGTCCTAGT            AGG
    ##    spacer_1    chr12  25209843      - | AAAGAAAAGATGAGCAAAGA            TGG
    ##              pam_site  cut_site      region percentGC     polyA     polyC
    ##             <numeric> <numeric> <character> <numeric> <logical> <logical>
    ##   spacer_29  25227322  25227325    region_6        45     FALSE     FALSE
    ##   spacer_37  25245275  25245278    region_5        35     FALSE     FALSE
    ##   spacer_41  25245343  25245346    region_5        60     FALSE     FALSE
    ##   spacer_45  25245365  25245368    region_5        25     FALSE     FALSE
    ##    spacer_3  25225615  25225618    region_7        40     FALSE     FALSE
    ##         ...       ...       ...         ...       ...       ...       ...
    ##   spacer_43  25245352  25245355    region_5        45     FALSE     FALSE
    ##   spacer_44  25245358  25245361    region_5        30     FALSE     FALSE
    ##    spacer_2  25209896  25209893    region_8        30     FALSE     FALSE
    ##    spacer_9  25225722  25225725    region_7        45     FALSE     FALSE
    ##    spacer_1  25209843  25209846    region_8        30      TRUE     FALSE
    ##                 polyG     polyT startingGGGGG        n0        n1      n0_c
    ##             <logical> <logical>     <logical> <numeric> <numeric> <numeric>
    ##   spacer_29     FALSE     FALSE         FALSE         1         0         1
    ##   spacer_37     FALSE     FALSE         FALSE         1         0         1
    ##   spacer_41     FALSE     FALSE         FALSE         1         0         1
    ##   spacer_45     FALSE     FALSE         FALSE         2         0         1
    ##    spacer_3     FALSE     FALSE         FALSE         1         0         1
    ##         ...       ...       ...           ...       ...       ...       ...
    ##   spacer_43     FALSE     FALSE         FALSE         1         1         1
    ##   spacer_44     FALSE     FALSE         FALSE         1         1         1
    ##    spacer_2     FALSE     FALSE         FALSE         1         1         1
    ##    spacer_9     FALSE     FALSE         FALSE         1         2         1
    ##    spacer_1     FALSE     FALSE         FALSE         1         2         1
    ##                  n1_c                                        alignments
    ##             <numeric>                                     <GRangesList>
    ##   spacer_29         0                                  chr12:25227322:-
    ##   spacer_37         0                                  chr12:25245275:-
    ##   spacer_41         0                                  chr12:25245343:-
    ##   spacer_45         0                  chr12:25245365:-,chr6:54770602:+
    ##    spacer_3         0                                  chr12:25225615:-
    ##         ...       ...                                               ...
    ##   spacer_43         0                  chr12:25245352:-,chr6:54770615:+
    ##   spacer_44         0                  chr12:25245358:-,chr6:54770609:+
    ##    spacer_2         0                  chr12:25209896:+,chr6:54771050:-
    ##    spacer_9         1 chr12:25225722:-,chr1:114709677:-,chr6:54770935:+
    ##    spacer_1         0   chr12:25209843:-,chr6:54771089:+,chr5:4348033:+
    ##             inRepeats score_cfd score_mit score_deephf      enzymeAnnotation
    ##             <logical> <numeric> <numeric>    <numeric>  <SplitDataFrameList>
    ##   spacer_29     FALSE       1.0       1.0     0.712546 FALSE:FALSE:FALSE:...
    ##   spacer_37     FALSE       1.0       1.0     0.609294 FALSE:FALSE:FALSE:...
    ##   spacer_41     FALSE       1.0       1.0     0.692967 FALSE:FALSE:FALSE:...
    ##   spacer_45     FALSE       0.5       0.5     0.671397 FALSE:FALSE:FALSE:...
    ##    spacer_3     FALSE       1.0       1.0     0.613590 FALSE:FALSE:FALSE:...
    ##         ...       ...       ...       ...          ...                   ...
    ##   spacer_43     FALSE  0.500000  0.619963     0.439317 FALSE:FALSE:FALSE:...
    ##   spacer_44     FALSE  0.777778  0.759301     0.433265 FALSE:FALSE:FALSE:...
    ##    spacer_2     FALSE  0.500000  0.577367     0.428607 FALSE:FALSE:FALSE:...
    ##    spacer_9     FALSE  0.566802  0.470146     0.567114 FALSE:FALSE:FALSE:...
    ##    spacer_1     FALSE  0.462185  0.434783     0.450868 FALSE:FALSE:FALSE:...
    ##                                                                 geneAnnotation
    ##                                                           <SplitDataFrameList>
    ##   spacer_29                      chr12:25227325:-:...,chr12:25227325:-:...,...
    ##   spacer_37 chr12:25245278:-:...,chr12:25245278:-:...,chr12:25245278:-:...,...
    ##   spacer_41 chr12:25245346:-:...,chr12:25245346:-:...,chr12:25245346:-:...,...
    ##   spacer_45 chr12:25245368:-:...,chr12:25245368:-:...,chr12:25245368:-:...,...
    ##    spacer_3 chr12:25225618:-:...,chr12:25225618:-:...,chr12:25225618:-:...,...
    ##         ...                                                                ...
    ##   spacer_43 chr12:25245355:-:...,chr12:25245355:-:...,chr12:25245355:-:...,...
    ##   spacer_44 chr12:25245361:-:...,chr12:25245361:-:...,chr12:25245361:-:...,...
    ##    spacer_2 chr12:25209893:+:...,chr12:25209893:+:...,chr12:25209893:+:...,...
    ##    spacer_9                      chr12:25225725:-:...,chr12:25225725:-:...,...
    ##    spacer_1 chr12:25209846:-:...,chr12:25209846:-:...,chr12:25209846:-:...,...
    ##                    tssAnnotation    hasSNP                     snps
    ##             <SplitDataFrameList> <logical>     <SplitDataFrameList>
    ##   spacer_29             :...,...     FALSE                 :...,...
    ##   spacer_37             :...,...     FALSE                 :...,...
    ##   spacer_41             :...,...     FALSE                 :...,...
    ##   spacer_45             :...,...     FALSE                 :...,...
    ##    spacer_3             :...,...     FALSE                 :...,...
    ##         ...                  ...       ...                      ...
    ##   spacer_43             :...,...     FALSE                 :...,...
    ##   spacer_44             :...,...     FALSE                 :...,...
    ##    spacer_2             :...,...     FALSE                 :...,...
    ##    spacer_9             :...,...     FALSE                 :...,...
    ##    spacer_1             :...,...      TRUE rs1137282:25209843:0:...
    ##                       rankings
    ##                   <data.frame>
    ##   spacer_29  spacer_29:5:1:...
    ##   spacer_37  spacer_37:7:1:...
    ##   spacer_41  spacer_41:7:1:...
    ##   spacer_45  spacer_45:7:1:...
    ##    spacer_3   spacer_3:8:1:...
    ##         ...                ...
    ##   spacer_43 spacer_43:31:2:...
    ##   spacer_44 spacer_44:31:2:...
    ##    spacer_2  spacer_2:32:2:...
    ##    spacer_9  spacer_9:47:3:...
    ##    spacer_1  spacer_1:50:3:...
    ##   -------
    ##   seqinfo: 640 sequences (1 circular) from hg38 genome
    ##   crisprNuclease: SpCas9

Our `guideSet` is now sorted according to our criteria, with the best
guides first. For a more detailed look at the rankings, we can look at
the appended column:

``` r
guideSet$rankings
```

    ##                  id rank bin1_n1 bin1_score_deephf bin1_isCommonExon
    ## spacer_29 spacer_29    5       1                 2                 2
    ## spacer_37 spacer_37    7       1                 3                 1
    ## spacer_41 spacer_41    7       1                 3                 1
    ## spacer_45 spacer_45    7       1                 3                 1
    ## spacer_3   spacer_3    8       1                 3                 2
    ## spacer_6   spacer_6    8       1                 3                 2
    ## spacer_7   spacer_7    8       1                 3                 2
    ## spacer_16 spacer_16    8       1                 3                 2
    ## spacer_17 spacer_17    8       1                 3                 2
    ## spacer_21 spacer_21    8       1                 3                 2
    ## spacer_23 spacer_23    8       1                 3                 2
    ## spacer_25 spacer_25    8       1                 3                 2
    ## spacer_26 spacer_26    8       1                 3                 2
    ## spacer_28 spacer_28    8       1                 3                 2
    ## spacer_34 spacer_34    8       1                 3                 2
    ## spacer_36 spacer_36    8       1                 3                 2
    ## spacer_40 spacer_40   10       1                 4                 1
    ## spacer_30 spacer_30   11       1                 4                 2
    ## spacer_31 spacer_31   11       1                 4                 2
    ## spacer_32 spacer_32   11       1                 4                 2
    ## spacer_38 spacer_38   13       1                 5                 1
    ## spacer_8   spacer_8   14       1                 5                 2
    ## spacer_12 spacer_12   14       1                 5                 2
    ## spacer_13 spacer_13   14       1                 5                 2
    ## spacer_27 spacer_27   14       1                 5                 2
    ## spacer_35 spacer_35   14       1                 5                 2
    ## spacer_4   spacer_4   17       1                 6                 2
    ## spacer_10 spacer_10   17       1                 6                 2
    ## spacer_18 spacer_18   17       1                 6                 2
    ## spacer_19 spacer_19   17       1                 6                 2
    ## spacer_42 spacer_42   25       2                 3                 1
    ## spacer_11 spacer_11   26       2                 3                 2
    ## spacer_22 spacer_22   26       2                 3                 2
    ## spacer_24 spacer_24   26       2                 3                 2
    ## spacer_33 spacer_33   26       2                 3                 2
    ## spacer_5   spacer_5   29       2                 4                 2
    ## spacer_43 spacer_43   31       2                 5                 1
    ## spacer_44 spacer_44   31       2                 5                 1
    ## spacer_2   spacer_2   32       2                 5                 2
    ## spacer_9   spacer_9   47       3                 4                 2
    ## spacer_1   spacer_1   50       3                 5                 2

The data frame contains a column for the bin value of each criterium,
along with an absolute score in the `rank` column, where a rank of `"1"`
indicates a guide that meets the highest specified level (bin value of
1) for each of our criteria (in this example, that would translate to
`n1<=0.5 && score_deephf>0.8 && isCommonExon==TRUE`).

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
    ##  [1] crisprScoreData_1.1.3             ExperimentHub_2.3.5              
    ##  [3] AnnotationHub_3.3.9               BiocFileCache_2.3.4              
    ##  [5] dbplyr_2.1.1                      BSgenome.Hsapiens.UCSC.hg38_1.4.4
    ##  [7] BSgenome_1.63.5                   rtracklayer_1.55.4               
    ##  [9] Biostrings_2.63.2                 XVector_0.35.0                   
    ## [11] GenomicRanges_1.47.6              GenomeInfoDb_1.31.6              
    ## [13] IRanges_2.29.1                    S4Vectors_0.33.11                
    ## [15] BiocGenerics_0.41.2               crisprDesign_0.99.109            
    ## [17] crisprBase_1.1.2                 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bitops_1.0-7                  matrixStats_0.61.0           
    ##  [3] bit64_4.0.5                   filelock_1.0.2               
    ##  [5] progress_1.2.2                httr_1.4.2                   
    ##  [7] tools_4.2.0                   utf8_1.2.2                   
    ##  [9] R6_2.5.1                      DBI_1.1.2                    
    ## [11] withr_2.5.0                   tidyselect_1.1.2             
    ## [13] prettyunits_1.1.1             bit_4.0.4                    
    ## [15] curl_4.3.2                    compiler_4.2.0               
    ## [17] crisprBowtie_1.1.1            cli_3.3.0                    
    ## [19] Biobase_2.55.0                basilisk.utils_1.9.1         
    ## [21] xml2_1.3.3                    DelayedArray_0.21.2          
    ## [23] randomForest_4.7-1            readr_2.1.2                  
    ## [25] rappdirs_0.3.3                stringr_1.4.0                
    ## [27] digest_0.6.29                 Rsamtools_2.11.0             
    ## [29] rmarkdown_2.13                crisprScore_1.1.13           
    ## [31] basilisk_1.9.2                pkgconfig_2.0.3              
    ## [33] htmltools_0.5.2               MatrixGenerics_1.7.0         
    ## [35] fastmap_1.1.0                 rlang_1.0.2                  
    ## [37] rstudioapi_0.13               RSQLite_2.2.12               
    ## [39] shiny_1.7.1                   BiocIO_1.5.0                 
    ## [41] generics_0.1.2                jsonlite_1.8.0               
    ## [43] vroom_1.5.7                   BiocParallel_1.29.18         
    ## [45] dplyr_1.0.8                   VariantAnnotation_1.41.3     
    ## [47] RCurl_1.98-1.6                magrittr_2.0.2               
    ## [49] GenomeInfoDbData_1.2.7        Matrix_1.4-0                 
    ## [51] Rcpp_1.0.8.3                  fansi_1.0.2                  
    ## [53] reticulate_1.25               Rbowtie_1.35.0               
    ## [55] lifecycle_1.0.1               stringi_1.7.6                
    ## [57] yaml_2.3.5                    SummarizedExperiment_1.25.3  
    ## [59] zlibbioc_1.41.0               grid_4.2.0                   
    ## [61] blob_1.2.2                    promises_1.2.0.1             
    ## [63] parallel_4.2.0                crayon_1.5.0                 
    ## [65] crisprBwa_1.1.2               dir.expiry_1.3.0             
    ## [67] lattice_0.20-45               GenomicFeatures_1.47.13      
    ## [69] hms_1.1.1                     KEGGREST_1.35.0              
    ## [71] knitr_1.37                    pillar_1.7.0                 
    ## [73] rjson_0.2.21                  biomaRt_2.51.3               
    ## [75] BiocVersion_3.15.0            XML_3.99-0.9                 
    ## [77] glue_1.6.2                    evaluate_0.15                
    ## [79] BiocManager_1.30.16           httpuv_1.6.5                 
    ## [81] png_0.1-7                     vctrs_0.3.8                  
    ## [83] tzdb_0.2.0                    purrr_0.3.4                  
    ## [85] assertthat_0.2.1              cachem_1.0.6                 
    ## [87] xfun_0.30                     mime_0.12                    
    ## [89] Rbwa_1.1.0                    xtable_1.8-4                 
    ## [91] restfulr_0.0.13               later_1.3.0                  
    ## [93] tibble_3.1.6                  GenomicAlignments_1.31.2     
    ## [95] AnnotationDbi_1.57.1          memoise_2.0.1                
    ## [97] interactiveDisplayBase_1.33.0 ellipsis_0.3.2
