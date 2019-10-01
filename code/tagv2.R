library(data.table)
library(stringr)
library(rtracklayer)
library(GenomicFeatures)
library(VariantAnnotation)

#### COSMIC
cos81.n.rng.red <- readRDS("data/cos81.n.rng.red.rds")
cos81.n.rng.red$select <- as.data.table(mcols(cos81.n.rng.red))[,(n.total>20 & n.max>6)]
cos81.n.rng.red.sel <- cos81.n.rng.red[(cos81.n.rng.red$select)]
cos81.n.rng.red.sel.med <- copy(cos81.n.rng.red.sel)
start(cos81.n.rng.red.sel.med) <- (end(cos81.n.rng.red.sel) + start(cos81.n.rng.red.sel)) %/% 2
end(cos81.n.rng.red.sel.med) <- (end(cos81.n.rng.red.sel) + start(cos81.n.rng.red.sel)) %/% 2

#### ExAC
## param <- ScanVcfParam(fixed="ALT", geno=NA_character_, info=c("AF"))
## exac <- readVcf("./tagv2/ExAC.r1.sites.vep.vcf.gz", "hg19", param, row.names=FALSE)
## exac.af <- sapply(info(exac)$AF, sum)
## exac.hf <- exac[exac.af>0.05 & exac.af<0.95]
## exac.hf.1 <- exac.hf[width(fixed(exac.hf)$REF)==1 & sapply(fixed(exac.hf)$ALT, function(ss) max(str_length(ss)))==1]
## seqlevelsStyle(exac.hf.1) <- "UCSC"
## tmp <- liftOver(rowRanges(exac.hf.1), import.chain("/mctp/users/mcieslik/proj/study/foxa1/data/raw/hg19ToHg38.over.chain"))
## exac.hf.gr <- unlist(tmp[lengths(tmp)==1])
## exac.sel.tx <- reduce(mapToTranscripts(exac.hf.gr, gen26.sel1.tx))
## tagv2.sel$exac <- tagv2.raw.gr %over% exac.sel.tx

#### TRANSCRIPTS
gen26.sel1 <- readRDS("data/gen26.sel1.rds")
gen26.sel1.ex <- gen26.sel1[gen26.sel1$type=="exon"]
gen26.sel1.cd <- gen26.sel1[gen26.sel1$type=="CDS"]
gen26.sel1.st <- gen26.sel1[gen26.sel1$type=="start_codon"]
gen26.sel1.sp <- gen26.sel1[gen26.sel1$type=="stop_codon"]
gen26.sel1.tx <- split(gen26.sel1.ex, gen26.sel1.ex$transcript_id)

#### map onto transcripts
cos81.sel1.rng <- reduce(mapToTranscripts(cos81.n.rng.red.sel, gen26.sel1.tx))
gen26.sel1.cds <- reduce(mapToTranscripts(gen26.sel1.cd, gen26.sel1.tx))
gen26.sel1.txs <- reduce(mapToTranscripts(gen26.sel1.ex, gen26.sel1.tx))
gen26.sel1.sts <- reduce(mapToTranscripts(gen26.sel1.st, gen26.sel1.tx))
gen26.sel1.sps <- reduce(mapToTranscripts(gen26.sel1.sp, gen26.sel1.tx))
gen26.sel1.mds <- gen26.sel1.txs
start(gen26.sel1.mds) <- (end(gen26.sel1.txs) + start(gen26.sel1.txs)) %/% 2
end(gen26.sel1.mds) <- (end(gen26.sel1.txs) + start(gen26.sel1.txs)) %/% 2
    
#### probes
tagv2.raw <- fread("./tagv2/5x_combined_dupesRemoved.txt")
setnames(tagv2.raw, c("transcript_id", "probe_id", "seq", "beg", "end"))
tagv2.raw.gr <- with(tagv2.raw, GRanges(transcript_id, IRanges(beg+1, end), probe_id=probe_id))

#### 1x coverage of TILE transcripts
tagv2.sel <- copy(tagv2.raw[,-3])
tagv2.sel[order(transcript_id,beg,end), idx:=1:nrow(.SD), by=transcript_id]
tagv2.sel[, grid:=(idx %% 5==1)]

#### cover each COSMIC hotspot with at least one probe
tagv2.sel[,cosmic.best:=FALSE]
cos81.sel.tx.ovr <- cos81.sel1.rng[cos81.sel1.rng %over% (tagv2.raw.gr+60)]
tagv2.sel[nearest(cos81.sel.tx.ovr, (tagv2.raw.gr-60)), cosmic.best:=TRUE]

#### probe overlap indices
tagv2.sel$cds <- tagv2.raw.gr %over% gen26.sel1.cds
tagv2.sel$cosmic <- tagv2.raw.gr %over% cos81.sel1.rng
tagv2.sel$starts <- tagv2.raw.gr %over% gen26.sel1.sts
tagv2.sel$stops <- tagv2.raw.gr %over% gen26.sel1.sps
tagv2.sel$mids <- tagv2.raw.gr %over% gen26.sel1.mds

#### 
setkey(tagv2.sel, transcript_id)
tmp <- as.data.table(mcols(gen26.sel1[gen26.sel1$type=="transcript"])[,c("transcript_id", "gene_name")])
setkey(tmp, transcript_id)
tagv2.sel <- tmp[tagv2.sel]
tile <- readLines("tagv2/tagv2_tile_list.txt")
tagv2.sel[,tile:=gene_name %in% tile]

####
tagv2.sel.stop <- tagv2.sel[order(stops, mids, cds, decreasing=TRUE), .SD[1], by=transcript_id]
tagv2.sel.start <- tagv2.sel[order(starts, mids, cds, decreasing=TRUE), .SD[1], by=transcript_id]
tagv2.sel.tag <- unique(rbind(tagv2.sel.start, tagv2.sel.stop))
tagv2.sel.grid <- tagv2.sel[(tile & grid & cds)]
tagv2.sel.cosmic <- tagv2.sel[(!tile & cosmic.best)]
tagv2.sel.tile <- rbind(tagv2.sel.grid, tagv2.sel.cosmic)
tagv2.sel.tagtile <- unique(rbind(tagv2.sel.tag, tagv2.sel.tile))

fwrite(tagv2.raw[probe_id %in% tagv2.sel.tagtile$probe_id], "out2/tagv2.tagtile.txt")
