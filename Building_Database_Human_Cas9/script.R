# This is the script used in the tutorial:

#0. required packages
library(crisprDatabase)
library(crisprBase)
library(crisprDesign)
library(crisprDesignData)
library(crisprDesignFacilitator)
overwrite <- FALSE


#i=1
#species <- "human"
#species <- "mouse"
#modality <- "crispra"
#modality <- "crisprkd"
#modality <- "crisprko"
#nuclease <- "SpCas9"
#nuclease <- "enAsCas12a"
#nuclease <- "CasRx"
#version  <- "v6"


#1. Tech options
i <- as.numeric(commandArgs(TRUE)[1])
species  <- as.character(commandArgs(TRUE)[2])
modality <- as.character(commandArgs(TRUE)[3])
nuclease <- as.character(commandArgs(TRUE)[4])
version  <- "v6"



#2. gRNA design options:
scoring_methods_cas9 <- c("ruleset1", "azimuth", 
                          "deephf",  "crisprscan", "crisprai",
                          "crisprater", "deepspcas9", "ruleset3")
scoring_methods_cas12a <- c("deepcpf1", "enpamgb")
scoring_methods_cas13d <- "casrxrf"
scoring_methods <- c(scoring_methods_cas9,
                     scoring_methods_cas12a,
                     scoring_methods_cas13d)
n_mismatches <- 3
max_mm <- 2
canonical_offtarget <- FALSE
if (nuclease=="SpCas9"){
    canonical_ontarget=TRUE
} else if (nuclease=="enAsCas12a"){
    canonical_ontarget=FALSE
} else if (nuclease=="CasRx"){
    canonical_ontarget=FALSE
}

### TSS WINDOW
if (modality=="crispko" | modality=="crisprkd"){
    tss_window <- NULL
} else if (modality=="crispra"){
    tss_window <- c(-500,0)
} else if (modality=="crispri"){
    tss_window <- c(0, 500)
} 


### Getting CRISPR nuclease object:
data(SpCas9, package="crisprBase")
data(enAsCas12a, package="crisprBase")
data(CasRx, package="crisprBase")
if (nuclease=="SpCas9"){
    crisprNuclease <- SpCas9
} else if (nuclease=="enAsCas12a"){
    crisprNuclease <- enAsCas12a
} else if (nuclease=="CasRx"){
    crisprNuclease <- CasRx
}


### Necessary annotation files
snpFile <- getSNPFile() 
if (nuclease=="CasRx"){
    bowtie_index <- getBowtieIndex(species=species, what="rna")
} else {
    bowtie_index <- getBowtieIndex(species=species, what="dna")
}
bsgenome  <- getGenomePackage(species=species)

# Stuff for CRISPRai
if (species=="human"){
    chromatinFiles <- getChromatinFiles()
    fastaFile <- getGenomeFasta()
} else {
    chromatinFiles <- NULL
    fastaFile <- NULL
}
if (species=="human" & nuclease=="SpCas9"){
    useDistanceToTss <- FALSE
} else {
    useDistanceToTss <- TRUE
}

# Getting binaries for CasRx
if (nuclease=="CasRx"){
    binaries <- crisprDesignFacilitator::getCasRxRfBinaries()
} else {
    binaries <- NULL
}



# SNP stuff
if (species=="human"){
    vcf <- getSNPFile()
} else {
    vcf <- NULL
}

# Conservation stuff
if (modality=="crisprko"){
    conservationFile <- crisprDesignFacilitator::getConservationFiles(species)
} else {
    conservationFile <- NULL
}



# Isoform stuff:
if (modality=="crisprko"){
    if (species=="human"){
        data(canonicalHuman, package="crisprDesignData")
        canonicalIsoforms <- canonicalHuman
    } else {
        data(canonicalMouse, package="crisprDesignData")
        canonicalIsoforms <- canonicalMouse
    }
} else {
    canonicalIsoforms <- NULL
}


if (species=="human"){
    data(tss_human, package="crisprDesignFacilitator")
    data(txdb_human, package="crisprDesignData")
    data(gr.repeats.hg38, package="crisprDesignData")
    data(pfamTableHuman, package="crisprDesignData")
    txObject <- txdb_human
    tssObject <- tss_human
    grRepeats <- gr.repeats.hg38
    pfamTable <- pfamTableHuman
} else {
    data(tss_mouse, package="crisprDesignFacilitator")
    data(txdb_mouse, package="crisprDesignData")
    data(gr.repeats.mm10, package="crisprDesignData")
    data(pfamTableMouse, package="crisprDesignData")
    txObject <- txdb_mouse
    tssObject <- tss_mouse
    grRepeats <- gr.repeats.mm10
    pfamTable <- pfamTableMouse
}

### Modality for crisprDesign
modality2 <- gsub("crispr", "CRISPR", modality)



#3. gene specification:
genedir <- crisprDatabase::getGeneModelDir(version=version,
                                           species=species)
if (modality=="crisprko"){
    ids <- readRDS(file.path(genedir, "genesids.500chunks.rds"))
} else if (modality=="crisprkd"){
    ids <- readRDS(file.path(genedir, "txids.500chunks.rds"))
} else {
    ids <- readRDS(file.path(genedir, "tssids.500chunks.rds"))
}
ids <- ids[[i]]



if (modality=="crisprko"){
    queryColumn="gene_id"
} else if (modality=="crispra" | modality=="crispri"){
    queryColumn="ID"
} else if (modality=="crisprkd"){
    queryColumn=NULL
}


#gene <- "ENSG00000133703" #KRAS
#gene <- "ENSG00000130270" #With repeats
extdir <- crisprDatabase::getCrisprDir(version=version,
                                       modality=modality,
                                       nuclease=nuclease,
                                       species=species)
if (!dir.exists(extdir)){
    dir.create(extdir, recursive=TRUE)
}
for (k in seq_along(ids)){
    id <- ids[k]
    filename <- file.path(extdir, paste0(id, '.rds'))
    print(modality)
    print(species)
    print(nuclease)
    print(k)
    print(id)
    if (!overwrite & file.exists(filename)){
        cat("Overwrite mode is off, and data already generated for this gene. Skipping :) \n")
    } else { 
        gs <- crisprDesign::designCompleteAnnotation(queryValue=id,
                                                     queryColumn=queryColumn,
                                                     modality=modality2,
                                                     bsgenome=bsgenome,
                                                     vcf=vcf,
                                                     tssObject=tssObject,
                                                     txObject=txObject,
                                                     bowtie_index=bowtie_index,
                                                     crisprNuclease=crisprNuclease,
                                                     n_mismatches=n_mismatches,
                                                     scoring_methods=scoring_methods,
                                                     max_mm=max_mm,
                                                     tss_window=tss_window,
                                                     canonical_ontarget=canonical_ontarget,
                                                     canonical_offtarget=canonical_offtarget,
                                                     grRepeats=grRepeats,
                                                     fastaFile=fastaFile,
                                                     chromatinFiles=chromatinFiles,
                                                     conservationFile=conservationFile,
                                                     geneCol="gene_symbol",
                                                     canonicalIsoforms=canonicalIsoforms,
                                                     binaries=binaries,
                                                     pfamTable=pfamTable)
        if (modality=="crisprko"){
            txid <- canonicalIsoforms$tx_id[match(id, canonicalIsoforms$gene_id)]
        } else {
            txid <- NULL
        }
        if (class(gs)=="GuideSet"){
            gs <- rankSpacers(gs,
                              modality=modality2,
                              commonExon=TRUE,
                              tx_id=txid,
                              useDistanceToTss=useDistanceToTss)
        } 
        saveRDS(gs, file=filename)
    }
}
q("no")









