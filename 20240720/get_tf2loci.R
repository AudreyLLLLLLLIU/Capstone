library(GenomicRanges)
library(assertthat)
library(digest)
library(fs)
library(magrittr)
library(readr)
library(rtracklayer)
library(tools)
library(utils)
library(tibble)
library(stringr)

get_tf2loci <- function(
    unibind_bed_dir = "https://unibind.uio.no/static/data/20220914/bulk_Robust/Homo_sapiens/damo_hg38_TFBS_per_TF.tar.gz",
    save_to_cache = TRUE,
    overwrite_cache = FALSE,
    read_from_cache = TRUE,
    cache_dir = "C:/Users/Jiaxuan/Desktop/Capstone/R/20240720/cache"
) {
  
  assert_that(is.string(unibind_bed_dir))
  assert_that(is.flag(save_to_cache))
  assert_that(is.flag(overwrite_cache))
  assert_that(is.flag(read_from_cache))
  
  dir_create(cache_dir)
  
  cache_fn <- fs::path(
    cache_dir,
    paste(
      "tf2loci-",
      digest(unibind_bed_dir, "md5", serialize = TRUE) %>% substr(1, 8),
      ".rds.gz",
      sep = ""))
  
  if (read_from_cache && file.exists(cache_fn)) {
    message("Using saved tf2loci at ", cache_fn, "...")
    return(read_rds(cache_fn))
  }
  
  unibind_bed_dir <- unibind_bed_dir %>%
    download_if_url() %>%
    extract_if_targz()
  
  dir_table <- dir(unibind_bed_dir)
  dir_of_tfs <- fs::path(unibind_bed_dir, dir_table)
  tfs <- dir(dir_of_tfs) %>%
    Filter(function(x) is_dir(fs::path(dir_of_tfs, x)), .)
  
  message(sprintf(
    "Loading. Each dot represents a transcription factor (total: %d).",
    length(tfs)))
  
  list_of_granges <- sapply(tfs, function(tf) {
    message(".", appendLF = FALSE)
    dirp <- fs::path(dir_of_tfs, tf)
    
    if (length(dir(dirp)) == 0) {
      warning(dirp, " is empty. Skipping")
      return(NULL)
    }
    
    granges <- lapply(dir(dirp), function(bedfn) {
      granges_sub <- import(file.path(dirp, bedfn))
    }) %>%
      { suppressWarnings(do.call(c, .)) }
    
    granges
  }, simplify = FALSE)
  message()
  
  not_nulls <- list_of_granges %>%
    sapply(is.null) %>%
    `!`
  
  granges_list <- GRangesList(list_of_granges[not_nulls])
  names(granges_list) <- tfs[not_nulls]
  
  if (!save_to_cache) return(granges_list)
  
  if (save_to_cache && file.exists(cache_fn) && !overwrite_cache) {
    message(
      "You indicated save_to_cache=TRUE, but a cache file already exists. ",
      "To overwrite, set overwrite_cache=TRUE.")
    return(granges_list)
  }
  
  message("Caching tf2loci to ", cache_fn, "...")
  write_rds(granges_list, cache_fn)
}


tf2loci <- get_tf2loci(
  unibind_bed_dir = "https://unibind.uio.no/static/data/20220914/bulk_Robust/Homo_sapiens/damo_hg38_TFBS_per_TF.tar.gz",
  save_to_cache = TRUE,
  overwrite_cache = FALSE,
  read_from_cache = TRUE,
  cache_dir = "C:/Users/Jiaxuan/Desktop/Capstone/R/20240720/cache"
)

