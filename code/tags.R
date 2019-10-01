library(data.table)
library(stringr)
library(rtracklayer)
library(GenomicFeatures)
library(matrixStats)

#### TRANSCRIPTS
gen26.sel1 <- readRDS("data/gen26.sel1x.rds")
gen26.sel1.ex <- gen26.sel1[gen26.sel1$type=="exon"]
gen26.sel1.cd <- gen26.sel1[gen26.sel1$type=="CDS"]
gen26.sel1.tx <- split(gen26.sel1.ex, gen26.sel1.ex$transcript_id)
gen26.sel1.cx <- split(gen26.sel1.cd, gen26.sel1.cd$transcript_id)
gen26.sel1.cds <- reduce(mapToTranscripts(gen26.sel1.cd, gen26.sel1.tx))
gen26.sel1.txs <- reduce(mapToTranscripts(gen26.sel1.ex, gen26.sel1.tx))
gen26.sel1.iqs <- gen26.sel1.cds
start(gen26.sel1.iqs) <- start(gen26.sel1.cds) + floor((end(gen26.sel1.cds) - start(gen26.sel1.cds)) * 1/4)
end(gen26.sel1.iqs) <- start(gen26.sel1.cds) + floor((end(gen26.sel1.cds) - start(gen26.sel1.cds)) * 3/4)

#### EXPRESSION
meta <- readRDS("/mctp/users/mcieslik/proj/study/met500/data/prod/expression/M.meta.rds")
meta.pair <- meta[sample_source %in% meta[,.N,by=sample_source][N==2]$sample_source]
meta.pair$assay <- str_match(meta.pair$run.id, "(.*)-(capt|poly)-.*")[,3]
meta.pair <- meta.pair[order(sample_source, assay)]
all.bws <- list.files("/mctp/projects/rnascape/data/samples/crisp/2.0.0/mctp/", pattern="*alig_csort.bw",
                      recursive=TRUE, full.name=TRUE)
names(all.bws) <- basename(dirname(all.bws))
meta.pair.sel <- meta.pair
bws <- all.bws[meta.pair.sel$run.id]
res <- mclapply(bws, function(bw) {
    suppressWarnings(mcols(unlist(summary(BigWigFile(bw), gen26.sel1.cd, type="mean", defaultValue=0))))$score
}, mc.cores=8)
resb <- do.call(cbind, res)

## correlation coefficient
tmp <- lapply(split(t(resb), meta.pair.sel$assay), matrix, nrow=nrow(meta.pair.sel)/2, byrow=TRUE)
tmp.cc <- mclapply(1:nrow(resb), function(i) cor(tmp[[1]][,i], tmp[[2]][,i]), mc.cores = 8)
mcols(gen26.sel1.cd) <- cbind(mcols(gen26.sel1.cd), cor=unlist(tmp.cc), cov=colQuantiles(tmp$poly, probs=0.75))
gen26.sel1.cd$cor[is.na(gen26.sel1.cd$cor)] <- -1 ## "data/rebs_exon.rds"
gen26.sel1.cor <- gen26.sel1.cd[gen26.sel1.cd$cor > 0.5]
gen26.sel1.cors <- reduce(mapToTranscripts(gen26.sel1.cor, gen26.sel1.tx))

#### probes
tagv3.raw <- fread("./tagv2/5x_combined_dupesRemoved.txt")
setnames(tagv3.raw, c("transcript_id", "probe_id", "seq", "beg", "end"))
tagv3.raw.gr <- with(tagv3.raw, GRanges(transcript_id, IRanges(beg+1, end), probe_id=probe_id)) # 2M probes

#### probe overlap indices
tagv3.sel <- copy(tagv3.raw[,-3])
tagv3.sel$cds <- tagv3.raw.gr %over% gen26.sel1.cds
tagv3.sel$iqs <- tagv3.raw.gr %over% gen26.sel1.iqs
tagv3.sel$cors <- tagv3.raw.gr %over% gen26.sel1.cors
tagv3.sel.iq <- tagv3.sel[order( iqs, cds, cors, decreasing=TRUE), .SD[1], by=transcript_id]
tagv3.sel.oq <- tagv3.sel[order(!iqs, cds, cors, decreasing=TRUE), .SD[1], by=transcript_id]
tagv3.sel <- unique(rbind(tagv3.sel.iq, tagv3.sel.oq))

### output
tagv3.sel.gr <- with(tagv3.sel, GRanges(transcript_id, IRanges(beg+1, end), probe_id=probe_id))
tagv3.sel.gr.onGm <- sort(mapFromTranscripts(tagv3.sel.gr, gen26.sel1.tx))
export(tagv3.sel.gr.onGm, "tagv3.5/all_tag.bed")
fwrite(tagv3.raw[probe_id %in% tagv3.sel$probe_id], "tagv3.5/all_tag_probes.txt", col.names =FALSE)
