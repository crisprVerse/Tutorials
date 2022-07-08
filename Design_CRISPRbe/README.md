crisprScore: on-target and off-target scoring for CRISPR gRNAs
================

-   [Overview](#overview)
-   [References](#references)

Authors: Jean-Philippe Fortin, Aaron Lun, Luke Hoberecht

Date: July 1, 2022

# Overview

crisprScore provides R wrappers of several on-target and off-target
scoring methods for CRISPR guide RNAs (gRNAs). The following nucleases
are supported: SpCas9, AsCas12a, enAsCas12a, and RfxCas13d (CasRx). The
available on-target cutting efficiency scoring methods are RuleSet1,
RuleSet3, Azimuth, DeepHF, DeepSpCas9, DeepCpf1, enPAM+GB, CRISPRscan
and CRISPRater. Both the CFD and MIT scoring methods are available for
off-target specificity prediction. The package also provides a
Lindel-derived score to predict the probability of a gRNA to produce
indels inducing a frameshift for the Cas9 nuclease. Note that DeepHF,
DeepCpf1 and enPAM+GB are not available on Windows machines.

Our work is described in a recent bioRxiv preprint: [“A comprehensive
Bioconductor ecosystem for the design of CRISPR guide RNAs across
nucleases and
technologies”](https://www.biorxiv.org/content/10.1101/2022.04.21.488824v2)

# References
