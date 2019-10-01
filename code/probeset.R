library(data.table)
library(stringr)
library(rtracklayer)
library(GenomicFeatures)

#### previously selected transcripts
gen26.sel1 <- readRDS("data/gen26.sel1x.rds")
gen26.sel1.ex <- gen26.sel1[gen26.sel1$type=="exon"]
gen26.sel1.tx <- split(gen26.sel1.ex, gen26.sel1.ex$transcript_id)

#### input probes from Agilent
tagv3.raw <- fread("tagv2/5x_combined_dupesRemoved.txt")
setnames(tagv3.raw, c("transcript_id", "probe_id", "seq", "beg", "end"))
tagv3.raw.gr <- with(tagv3.raw, GRanges(transcript_id, IRanges(beg+1, end), probe_id=probe_id)) # 2M probes

#### COSMIC
cos81.n.rng.red <- readRDS("data/cos81.n.rng.red.rds")
cos81.n.rng.red$select <- as.data.table(mcols(cos81.n.rng.red))[,(n.total>20 & n.max>6)]
cos81.n.rng.red.sel <- cos81.n.rng.red[(cos81.n.rng.red$select)]
cos81.n.rng.red.sel.med <- copy(cos81.n.rng.red.sel)
start(cos81.n.rng.red.sel.med) <- (end(cos81.n.rng.red.sel) + start(cos81.n.rng.red.sel)) %/% 2
end(cos81.n.rng.red.sel.med) <- (end(cos81.n.rng.red.sel) + start(cos81.n.rng.red.sel)) %/% 2 # 4401
cos81.sel1.rng <- reduce(mapToTranscripts(cos81.n.rng.red.sel, gen26.sel1.tx))

#### project exon boundaries onto transcript sequence
gen26.sel1.ex.B <- sort(gen26.sel1[gen26.sel1$type=="exon"])
end(gen26.sel1.ex.B) <- start(gen26.sel1.ex.B)
gen26.sel1.ex.B$edge <- "B"
gen26.sel1.ex.E <- sort(gen26.sel1[gen26.sel1$type=="exon"])
start(gen26.sel1.ex.E) <- end(gen26.sel1.ex.E)
gen26.sel1.ex.E$edge <- "E"
gen26.sel1.ex.BE <- sort(c(gen26.sel1.ex.B, gen26.sel1.ex.E))
gen26.sel1.ex.BE.onTx.hits <- mapToTranscripts(gen26.sel1.ex.BE, gen26.sel1.tx)
gen26.sel1.ex.BE.onTx.hits$transcript_source <- gen26.sel1.ex.BE[gen26.sel1.ex.BE.onTx.hits$xHits]$transcript_id
gen26.sel1.ex.BE.onTx <- gen26.sel1.ex.BE.onTx.hits[as.character(seqnames(gen26.sel1.ex.BE.onTx.hits))==gen26.sel1.ex.BE.onTx.hits$transcript_source]
gen26.sel1.ex.BE.onTx$edge <- gen26.sel1.ex.BE[gen26.sel1.ex.BE.onTx$xHits]$edge
gen26.sel1.ex.BE.onTx$exon_id <- gen26.sel1.ex.BE[gen26.sel1.ex.BE.onTx$xHits]$exon_id
gen26.sel1.ex.BE.onTx$exon_number <- gen26.sel1.ex.BE[gen26.sel1.ex.BE.onTx$xHits]$exon_number

## 1. split coordinate sorted exon beginnings by transcript (
## gen26.sel1.ex.B.noF <- unname(unlist(endoapply(split(gen26.sel1.ex.B, gen26.sel1.ex.B$transcript_id), tail, -1)))
## gen26.sel1.ex.E.noL <- unname(unlist(endoapply(split(gen26.sel1.ex.E, gen26.sel1.ex.E$transcript_id), head, -1)))

####
gen26.sel1.ex.B.onTx <- gen26.sel1.ex.BE.onTx[gen26.sel1.ex.BE.onTx$edge=="B"]
sel.B <- precede(gen26.sel1.ex.B.onTx, tagv3.raw.gr-1)
sel.B <- sel.B[!is.na(sel.B)]
gen26.sel1.ex.E.onTx <- gen26.sel1.ex.BE.onTx[gen26.sel1.ex.BE.onTx$edge=="E"]
sel.E <- follow(gen26.sel1.ex.E.onTx, tagv3.raw.gr-1)
sel.E <- sel.E[!is.na(sel.E)]
sel.BE <- unique(c(sel.B, sel.E))

####
tagv3.sel <- copy(tagv3.raw)
tagv3.sel$edge <- FALSE
tagv3.sel$grid <- FALSE
tagv3.sel$cosmic <- FALSE
tagv3.sel$cosmic <- FALSE
tagv3.sel$idx <- -1L

#### cover each COSMIC hotspot with at least one probe
cos81.sel.tx.ovr <- cos81.sel1.rng[cos81.sel1.rng %over% (tagv3.raw.gr+60)]
tagv3.sel[nearest(cos81.sel.tx.ovr, (tagv3.raw.gr-60)), cosmic:=TRUE] # 3679
tagv3.sel.BE <- tagv3.sel[sel.BE]
tagv3.sel.BE$edge <- TRUE
tagv3.sel.MD <- tagv3.sel[setdiff(1:nrow(tagv3.sel), sel.BE)]
sel.over.BE <- with(tagv3.sel.MD, GRanges(transcript_id, IRanges(beg+1, end))) %over% with(tagv3.sel.BE, GRanges(transcript_id, IRanges(beg+1, end)))
tagv3.sel.MD <- tagv3.sel.MD[!sel.over.BE | cosmic]
tagv3.sel.MD[order(transcript_id,beg,end), idx:=1:nrow(.SD), by=transcript_id]
tagv3.sel.MD[, grid:=(idx %% 5==1)]
tagv3.sel <- rbind(tagv3.sel.BE, tagv3.sel.MD)

#### output
tagv3.sel.tile <- tagv3.sel[(edge | grid | cosmic)]
tagv3.sel.tile.gr <- with(tagv3.sel.tile, GRanges(transcript_id, IRanges(beg+1, end), probe_id=probe_id))
tagv3.sel.tile.gr.onGm <- mapFromTranscripts(tagv3.sel.tile.gr, gen26.sel1.tx)
export(tagv3.sel.tile.gr.onGm, "tagv3.5/all_tile.bed")
fwrite(tagv3.sel.tile, "tagv3.5/all_tile_probes.txt", col.names =FALSE)
