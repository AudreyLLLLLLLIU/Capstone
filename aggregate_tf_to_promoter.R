library(GenomicRanges)
library(dplyr)

aggregate_tf_to_promoter <- function(promoterAnnotation, tf2loci) {
  promoterCoordinates <- promoterCoordinates(promoterAnnotation)
  results <- data.frame(promoterId = character(), tf = character(), stringsAsFactors = FALSE)
  for (tf in names(tf2loci)) {
    tf_binding_sites <- tf2loci[[tf]]
    overlaps <- findOverlaps(promoterCoordinates, tf_binding_sites)
    if (length(overlaps) > 0) {
      overlapping_promoters <- queryHits(overlaps)
      promoter_ids <- mcols(promoterCoordinates)$promoterId[overlapping_promoters]
      result <- data.frame(promoterId = promoter_ids, tf = tf, stringsAsFactors = FALSE)
      results <- rbind(results, result)
    }
  }
  return(results)
}
