library(here)
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(plyranges)

## first read in the gencode v19 annotation
gencode <- rtracklayer::import(here('gencode.v19.annotation.gtf.gz'))

transcripts <- gencode |> 
  filter(type=='transcript') #,transcript_status=='KNOWN')

length(transcripts) # 196520, same as in Goke supplement table 1
## clearly they were working with *all* transcripts, not just known.
## I tried both, but this is to replicate their results.

first_exons <- 
  gencode |> 
  filter(type=='exon',exon_number==1) # transcript_status=='KNOWN')

#promoters <- GenomicRanges::reduce(first_exons)

promoters <- first_exons |>
  GenomicRanges::reduce(with.revmap=TRUE)

## let's give the promoters ids. These won't be the same ids as in Goke's paper, but
## we'll still need them.
promoters$promoter_id <- sprintf("promoter.%05d",1:length(promoters))

## we need the number of first exons (or transcripts) for each set of overlapping exons.
promoters$len <- unlist(as.list(promoters$revmap) |> map(length))

## Now build a data frame with two columns
promoter_df <- 
  tibble(
    # column 1 is going to be the transcript_id, which is the rowfrom the row of transcripts 
    # that comes from the revmap
    transcript_id=transcripts$transcript_id[unlist(promoters$revmap)],
    # Then we need to repeat each promoter id for the number of times it is used
         promoter_id=rep(promoters$promoter_id,promoters$len))

dim(promoter_df) ## 196520 rows.
promoter_df ## looks OK

## let's keep a granges object so we remember where each promoter is!
promoter_gr <- select(promoters,promoter_id)

## now go through the same thing with the TSSs. We can use either transcripts
## or first exons for this
tss <- first_exons |> 
  resize(width=1) |>
  GenomicRanges::reduce(with.revmap=TRUE) 

## same process as with promoters
tss$tss_id <- sprintf("tss.%06d",1:length(tss))
tss$len <- unlist(as.list(tss$revmap) |> map(length))

tss_df <- tibble(transcript_id=first_exons$transcript_id[unlist(tss$revmap)],
                 tss_id=rep(tss$tss_id,tss$len),
                 exon_id = first_exons$exon_id[unlist(tss$revmap)]
)
tss_gr <- select(tss,tss_id)

## these should be the same size, and match the number of transcripts.
dim(promoter_df) 
dim(tss_df)

## Now let's put it together into one big tibble.
master_table <- inner_join(tss_df,promoter_df) |> 
  inner_join(as.data.frame(mcols(transcripts)) |> 
               select(transcript_id,gene_id)) |> 
  arrange(gene_id)

if (!dir.exists(here('results'))) {
  dir.create(here('results'))
}

save(master_table,promoter_gr,tss_gr,file=here('results','master_table_etc.RData'))

view(master_table)
#########################################################
install.packages("data.table")
library(data.table)

exon_metadata <- mcols(first_exons)

# Merge the first exon
merged_exons <- first_exons %>%
  GenomicRanges::reduce(with.revmap = TRUE)

## Then I want to add metadata back to merged exons
merged_exons$exon_id <- sapply(merged_exons$revmap, function(x) paste(exon_metadata$exon_id[x], collapse = ","))
merged_exons$transcript_id <- sapply(merged_exons$revmap, function(x) paste(exon_metadata$transcript_id[x], collapse = ","))
merged_exons$gene_id <- sapply(merged_exons$revmap, function(x) paste(exon_metadata$gene_id[x], collapse = ","))

# Add a fixed upstream area
promoters_with_upstream <- promoters(merged_exons, upstream = 2000, downstream = width(merged_exons))
# Assign a new promoter ID
promoters_with_upstream$promoter_id <- sprintf("promoter.%05d", 1:length(promoters_with_upstream))


#Column lengths do not match because of the way strsplit and unlist are used,
#I try to use the metadata of merged_exons directly to generate promoters_with_upstream_dt
#ensuring that all columns are the same length,but this will take a long time
promoter_revmap <- promoters_with_upstream$revmap


promoters_with_upstream_dt <- map_dfr(seq_along(promoter_revmap), function(i) {
  revmap_indices <- promoter_revmap[[i]]
  tibble(
    transcript_id = rep(merged_exons$transcript_id[i], length(revmap_indices)),
    tss_id = sprintf("tss.%06d", revmap_indices),
    promoter_id = promoters_with_upstream$promoter_id[i],
    seqnames = rep(as.character(seqnames(promoters_with_upstream)[i]), length(revmap_indices)),
    start = rep(start(promoters_with_upstream)[i], length(revmap_indices)),
    end = rep(end(promoters_with_upstream)[i], length(revmap_indices)),
    strand = rep(as.character(strand(promoters_with_upstream)[i]), length(revmap_indices))
  )
})

#promoters_with_upstream_df <- do.call(rbind, lapply(seq_along(promoter_revmap), function(i) {
  #revmap_indices <- promoter_revmap[[i]]
  #tibble(
    #transcript_id = rep(merged_exons$transcript_id[i], length(revmap_indices)),
    #tss_id = sprintf("tss.%06d", revmap_indices),
    #promoter_id = promoters_with_upstream$promoter_id[i],
    #seqnames = rep(as.character(seqnames(promoters_with_upstream)[i]), length(revmap_indices)),
    #start = rep(start(promoters_with_upstream)[i], length(revmap_indices)),
    #end = rep(end(promoters_with_upstream)[i], length(revmap_indices)),
    #strand = rep(as.character(strand(promoters_with_upstream)[i]), length(revmap_indices))
#  )
#}))

save(promoters_with_upstream, promoters_with_upstream_dt, file = here('results', 'promoters_with_upstream.RData'))

view(promoters_with_upstream_dt)
tail(promoters_with_upstream_dt, 1000)


#################################################################################################################
##################Merge the first exons and then incorporate overlapping information from cCREs.tsv##############
#从最前面first exon跑完，跳到这里，promoters这个名字重名了，这里的promoter应该是，好了不用了改好了##

cCREs <- read_tsv('data/cCREs.tsv')
cCREs_gr <- makeGRangesFromDataFrame(cCREs, keep.extra.columns = TRUE)

# Find overlaps between promoters_with_upstream and cCREs
overlaps <- findOverlaps(promoters_with_upstream, cCREs_gr)

# Extract overlapping promoters and cCREs
overlapping_promoters <- promoters_with_upstream[queryHits(overlaps)]
overlapping_cCREs <- cCREs_gr[subjectHits(overlaps)]

# Merge the exons and cCREs information
merged_df <- tibble(
  exon_id = overlapping_promoters$exon_id,
  tss_id = overlapping_promoters$tss_id,
  transcript_id = overlapping_promoters$transcript_id,
  promoter_id = overlapping_promoters$promoter_id,
  cCRE_id = overlapping_cCREs$name,
  seqnames = as.character(seqnames(overlapping_promoters)),
  start = start(overlapping_promoters),
  end = end(overlapping_promoters),
  strand = as.character(strand(overlapping_promoters))
)

# Save the result
save(merged_df, file = here('results', 'overlapping_promoters_with_cCREs.RData'))
view(merged_df)











