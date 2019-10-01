library(data.table)
library(GenomicRanges)

tag.transcripts <- readRDS("data/gen26.sel1x.rds")
g2t <- as.data.table(mcols(tag.transcripts[tag.transcripts$type=="transcript"])[,c("gene_name", "transcript_id")])
setkey(g2t, gene_name)
tag.probes <- fread("tagv3.5/all_tag_probes.txt")
tile.probes <- fread("tagv3.5/all_tile_probes.txt")
tag.onco.genes <- readLines("tagv3.5/tagv3.5_tag_list.txt")
tile.onco.genes <- readLines("tagv3.5/tagv3.5_tile_list.txt")
sel.tag <- tag.probes$V1 %in% g2t[tag.onco.genes]$transcript_id
sel.tile <- tile.probes$V1 %in% g2t[tile.onco.genes]$transcript_id
tt.probes <- rbind(
    tag.probes[sel.tag][,1:5],
    tile.probes[sel.tile][,1:5]
)
fwrite(tt.probes, "tagv3.5/tagtile_onco3.5_probes.txt", col.names=FALSE)
