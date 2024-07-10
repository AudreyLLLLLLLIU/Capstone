#################################
###tx to tss###
#' Get a mapping of transcripts to transcription start sites.
#'
#' @param upstream Nucleotide distance for how far upstream of a transcript
#'   should be considered as part of its transcription start site.
#' @param downstream Nucleotide distance for how far downstream of a transcript
#'   should be considered as part of its transcription start site.
#' @param annotation_gtf_file File path or URL to a GTF annotation file in the
#'   style of ENSEMBL.
#' @param save_to_cache Should the results be cached? (\code{overwrite_cache}
#'   modifies this behaviour.)
#' @param overwrite_cache If there was a pre-existing cache from the same
#'   \code{annotation_gtf_file} with \code{save_to_cache}, should this new
#'   function call overwrite that cached result? Note that different
#'   \code{annotation_gtf_file} files are cached differently (they are
#'   version-aware).
#' @param read_from_cache If there was a pre-existing cache from the same
#'   \code{annotation_gtf_file}, should this new function call just read from
#'   that, rather than try to download and parse it again?
#'
#' @return tibble with two columns: \code{target_id} which contains ENST
#'   transcript accessions, and \code{tss_id} which is a unique ID string for
#'   each trascription start site (TSS). \code{tss_id} is a comma separated
#'   string with chromosome name, start of TSS window, end of TSS window, and
#'   strand (\code{"+"} or \code{"-"}). Note that \code{start <= end} regardless
#'   of strand.
#'
#'   This tibble can be used as \code{sleuth_prep}'s \code{target_mapping}. Just
#'   remebmer to set \code{aggregation_column = "tss_id"}!
#'
#' @importFrom GenomicRanges GRanges resize promoters reduce findOverlaps mcols
#'   seqnames start end strand
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors DataFrame queryHits subjectHits
#' @importFrom assertthat assert_that is.string is.count is.flag
#' @importFrom dplyr filter select mutate rename
#' @importFrom fs dir_create path file_exists
#' @importFrom magrittr %>%
#' @importFrom readr read_rds write_rds
#' @importFrom rlang .data
#' @importFrom rtracklayer import
#' @importFrom tibble tibble as_tibble
#' @importFrom tidyr unite
#' @importFrom tools R_user_dir
#' @importFrom utils packageName
get_tx2tss <- function(
    upstream=500, downstream=100,
    annotation_gtf_file="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz",
    save_to_cache=TRUE,
    overwrite_cache=FALSE,
    read_from_cache=TRUE
) {
  
  assert_that(is.count(upstream) && is.count(downstream))
  assert_that(is.string(annotation_gtf_file))
  assert_that(is.flag(save_to_cache))
  assert_that(is.flag(overwrite_cache))
  assert_that(is.flag(read_from_cache))
  
  # Does a cache file already exist? Do we read_from_cache?
  cache_dir <- R_user_dir('capstone', "cache") %>% dir_create()
  cache_fn <- fs::path(
    cache_dir,
    paste(
      "tx2tss-u", upstream, "-d", downstream, "-",
      digest(annotation_gtf_file, "md5", serialize=TRUE) %>% substr(1, 8),
      ".rds.gz",
      sep=""))
  if (read_from_cache && file_exists(cache_fn)) {
    message("Using saved tx2tss at ", cache_fn, "...")
    return(read_rds(cache_fn))
  }
  
#annotation_gtf_file <- download_if_url(annotation_gtf_file)
library(rtracklayer)
annotation_gtf_file <- import.gff("gencode.v19.annotation.gtf.gz", format="gtf")
 #gtf <- import.gff("gencode.v19.annotation.gtf.gz", format="gtf")
 #head(gtf)
 #length(gtf)
 #str(gtf)
  # From annotation_gtf_file, which has 2619444 records
  #if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
   # BiocManager::install("GenomicRanges")}
  #if (!requireNamespace("dplyr", quietly = TRUE)) {
   # install.packages("dplyr")}
  
  #library(GenomicRanges)
  #library(dplyr)
  #transcripts <- GRanges(
    #ranges = ranges(transcripts),
    #seqnames = seqnames(transcripts),
    #strand = strand(transcripts),
    #mcols = transcripts_mcols)
  # 怎么filter啊，我弄出来length是0 filter for only transcript records, and only the fields we care about.
 

## take every protein-coding transcript and define
## all the promoter regions
annotation_gtf_file |> 
  plyranges::filter(type=='transcript',gene_type=='protein_coding') |> 
  GenomicFeatures::promoters(upstream=1000, downstream=1000) -> 
  protein_coding_promoters

#length(protein_coding_promoters)
annotation_gtf_file |> 
  plyranges::filter(type=='exon',gene_type=='protein_coding') ->
  protein_coding_exons

## The TxID in this corresponds to the Tx ID in supplementary Table 1 of Goeke.
protein_coding_exons |> 
  filter(exon_number==1) -> 
  first_exons


first_exons$transcript_id

## Let's get the TSSs, meaning the *exact* duplicate starts. We ignore overlap
## To get this, we resize to a width of 1 and then reduce. These only overlap
## if they have the exact same start.

TSSs <- resize(first_exons,width=1) |>
  GenomicRanges::reduce(with.revmap=TRUE) 

TSSs


## Now we get the promoters 
promoter_starts <-
  GenomicRanges::reduce(first_exons,with.revmap=TRUE)

promoter_starts

first_exons |>
  filter(gene_id == "ENSG00000000003.10")

find_overlapping_first_exons <- function(promoters, exons) {
  exons <- exons[order(start(exons))]
  overlap <- findOverlaps(promoters, exons, type="within")
  first_exons <- exons[subjectHits(overlap)]
  first_exons <- first_exons[!duplicated(first_exons$gene_id)]
  return(first_exons)
}

#first_exons <- protein_coding_exons[which.min(start(protein_coding_exons$range)), ]
first_exons <- find_overlapping_first_exons(protein_coding_promoters, protein_coding_exons)
# first_exons GRanges object with 0 ranges and 21 metadata columns

transcript_table <- data.frame(
  transcriptId = mcols(first_exons)$transcript_id,
  tssId = paste0("tss.", start(first_exons)),
  promoterId = mcols(first_exons)$promoter_id,
  geneId = mcols(first_exons)$gene_id
)


print(transcript_table)







