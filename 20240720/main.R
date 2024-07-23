if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("proActiv")

library(proActiv)


## From GTF file path
gtf.file <- "C:/Users/Jiaxuan/Desktop/Capstone/R/20240720/data/gencode.v46.basic.annotation.gtf.gz"

promoterAnnotation.gencode.v46.subset <- preparePromoterAnnotation(file = gtf.file,
                                                                   species = 'Homo_sapiens')
#Warning message:In .get_cds_IDX(mcols0$type, mcols0$phase) :
#The "phase" metadata column contains non-NA values for features of type stop_codon.This information was ignored.
#do not need to handle stop_codon's phase information specifically, so maybe can ignore this warning


##mapping of transcription factors to genomic sites
tf2loci <- get_tf2loci()

aggregated_tf_to_promoter <- aggregate_tf_to_promoter(promoterAnnotation.gencode.v46.subset, tf2loci)

head(aggregated_tf_to_promoter)
