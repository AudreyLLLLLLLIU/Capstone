## Assume we are at the top of the project


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("proActiv",quietly=TRUE)) {
      BiocManager::install("proActiv")
}

if (!requireNamespace("BiocFileCache",quietly=TRUE)) {
  BiocManager::install("BiocFileCache")
}

## make sure I have the data. I want to use relative path names so this
## will work on any computer

library(proActiv)
library(BiocFileCache)
library(here)


bfc <- BiocFileCache()
bfcinfo(bfc)

url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.basic.annotation.gtf.gz"

gencode46 <- BiocFileCache::bfcadd(bfc,"Gencode46",url)

## From GTF file path
gtf_file <- here("data","gencode.v46.basic.annotation.gtf.gz")

suppressWarnings(dir.create(here("data")))

if (!file.exists(gtf_file)) {
    download.file(url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.basic.annotation.gtf.gz",
                  destfile = gtf_file)
}


promoterAnnotation_gencode_46 <-
    preparePromoterAnnotation(file = gtf_file,
                              species = 'Homo_sapiens')

txdbmaker::makeTxDbFromGFF(file = gtf_file, format = "gtf", organism = 'Homo sapiens')
#Warning message:In .get_cds_IDX(mcols0$type, mcols0$phase) :
#The "phase" metadata column contains non-NA values for features of type stop_codon.This information was ignored.
#do not need to handle stop_codon's phase information specifically, so maybe can ignore this warning

options(timeout = 6000)

##mapping of transcription factors to genomic sites
tf2loci <- get_tf2loci("https://unibind.uio.no/static/data/20220914/bulk_Robust/Homo_sapiens/damo_hg38_TFBS_per_TF.tar.gz")


aggregated_tf_to_promoter <- aggregate_tf_to_promoter(promoterAnnotation_gencode_46, tf2loci)

head(aggregated_tf_to_promoter)
