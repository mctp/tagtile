{
suppressMessages({
    library(data.table)
    library(stringr)
    library(GenomicRanges)
    library(GenomicFeatures)
    library(rtracklayer)
    library(BSgenome.Hsapiens.UCSC.hg38)
})

{
goi <- c(readLines("data/gene_list.txt"), "T", "PTPN14") ## genes of interest
cos81.n.rng.red <- readRDS("data/cos81.n.rng.red.rds") ## cosmic hotspots 1.5M total select 10k
agi_v4.hg38  <- readRDS("data/agi_v4.hg38.rds") ## AGI v4 probes in HG38 coordinates
gen26.sel1 <- readRDS("data/gen26.sel1.rds") ## 1 transcript per protein-coding gene
gen26.sel3 <- readRDS("data/gen26.sel3.rds") ## 3 transcrpts per protein-coding gene
gen26.sel1.txw <- readRDS("data/gen26.sel1.txw.rds") ## valid transcript windows
gen26.sel1.txw.map <- readRDS("data/gen26.sel1.txw.map.rds") ## valid transcript windows in genome space (xHits - windw, transcriptHits, includes introns)
t4.score.avg.nrm <- readRDS("data/t4.score.avg.nrm.rds")

}

{ ## get transcript sequences and lengths
gen26.sel1.ex <- gen26.sel1[gen26.sel1$type=="exon"]
gen26.sel1.exl <- split(gen26.sel1.ex, gen26.sel1.ex$transcript_id)
tx.seqs <- extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg38, gen26.sel1.exl)
tx.lens <- sum(width(gen26.sel1.exl))
## all(tx.lens==nchar(tx.seqs))
}


## COV
{
aaa <- cos81.n.rng.red[cos81.n.rng.red$select]
xxx <- findOverlaps(gen26.sel1.txw.map, aaa)
yyy <- data.table(queryHits(xxx), aaa[subjectHits(xxx)]$n.max)[,max(V2),by=V1]
setkey(yyy, V1)
zzz <- yyy[J(1:length(gen26.sel1.txw.map)),]$V1.1
cos <- ifelse(is.na(zzz), 0, zzz)
tmp <- data.table(
    transcript_id = as.character(seqnames(gen26.sel1.txw)),
    probe_id = 1:length(gen26.sel1.txw.map),
    cov = sqrt(t4.score.avg.nrm.cov$capt.cov*t4.score.avg.nrm.cov$poly.cov),
    cds = gen26.sel1.txw.map %over% gen26.sel1[gen26.sel1$type=="CDS"],
    agi = gen26.sel1.txw.map %over% agi_v4.hg38,
    cos = cos
)
tmp[,probe_n:=1:.N, by=transcript_id]
setkey(tmp, transcript_id)
ttt <- unique(as.data.table(mcols(gen26.sel1)[,c("gene_id", "gene_name", "transcript_id")]))
setkey(ttt, transcript_id)
tmp <- tmp[ttt]
setkey(tmp, transcript_id, probe_n)
cov.probes.1 <- tmp[order(-cov, -agi, -cds), head(.SD, 1), by=transcript_id]
cov.probes.1.filter <- cov.probes.1[,.(transcript_id, probe_n,
                                       probe_nm1=probe_n-1, probe_np1=probe_n+1,
                                       probe_nm2=probe_n-2, probe_np2=probe_n+2
                                       )]
setkey(cov.probes.1.filter, transcript_id, probe_n)
tmp <- tmp[!cov.probes.1.filter]
setkey(cov.probes.1.filter, transcript_id, probe_nm1)
tmp <- tmp[!cov.probes.1.filter]
setkey(cov.probes.1.filter, transcript_id, probe_np1)
tmp <- tmp[!cov.probes.1.filter]
setkey(cov.probes.1.filter, transcript_id, probe_nm2)
tmp <- tmp[!cov.probes.1.filter]
setkey(cov.probes.1.filter, transcript_id, probe_np2)
tmp <- tmp[!cov.probes.1.filter]
cov.probes.2 <- tmp[order(-cov, -agi, -cds), head(.SD, 1), by=transcript_id]
cov.probes <- rbind(cov.probes.1, cov.probes.2)
}

## COS
cos.probes <- tmp[!(probe_id %in% cov.probes$probe_id)][gene_name %in% goi][cos>0]



gen26.sel1.txw.cov180 <- gen26.sel1.txw[cov.probes$probe_id]+120 # was 60
gen26.sel1.txw.cos180 <- gen26.sel1.txw[cos.probes$probe_id]+120
gen26.sel1.txw.cov180.map <- mapFromTranscripts(gen26.sel1.txw.cov180, gen26.sel1.exl)
gen26.sel1.txw.cos180.map <- mapFromTranscripts(gen26.sel1.txw.cos180, gen26.sel1.exl)
gen26.sel1.txw.cov180.map$name <- cov.probes[, paste(gene_name, gene_id, transcript_id, "C", probe_id, sep="_")]
gen26.sel1.txw.cos180.map$name <- cos.probes[, paste(gene_name, gene_id, transcript_id, "M", probe_id, sep="_")]


gen26.sel1.txw.cov.seq <- subseq(tx.seqs[seqnames(gen26.sel1.txw.cov180)],
                                 start(gen26.sel1.txw.cov180),
                                 end(gen26.sel1.txw.cov180))
names(gen26.sel1.txw.cov.seq) <- cov.probes[, paste(gene_name, gene_id, transcript_id, "C", probe_id, sep="_")]
gen26.sel1.txw.cos.seq <- subseq(tx.seqs[seqnames(gen26.sel1.txw.cos180)],
                                 start(gen26.sel1.txw.cos180),
                                 end(gen26.sel1.txw.cos180))
names(gen26.sel1.txw.cos.seq) <- cos.probes[, paste(gene_name, gene_id, transcript_id, "M", probe_id, sep="_")]



fin.map <- sort(c(gen26.sel1.txw.cov180.map, gen26.sel1.txw.cos180.map)[,3])
fin.seq <- c(gen26.sel1.txw.cov.seq, gen26.sel1.txw.cos.seq)[fin.map$name]
writeXStringSet(fin.seq, "out/tag_v1.fa")
export.bed(fin.map, "out/tag_v1.bed")

}

##
gen26.sel1.txw.cov300 <- gen26.sel1.txw[cov.probes$probe_id]+120
start(gen26.sel1.txw.cov300[start(gen26.sel1.txw.cov300)<0]) <- 1
tx.len.cov300 <- tx.lens[as.character(seqnames(gen26.sel1.txw.cov300))]
end(gen26.sel1.txw.cov300[end(gen26.sel1.txw.cov300)>tx.len.cov300]) <- tx.len.cov300[end(gen26.sel1.txw.cov300)>tx.len.cov300]
gen26.sel1.txw.cov300.map <- mapFromTranscripts(gen26.sel1.txw.cov300, gen26.sel1.exl)
gen26.sel1.txw.cov300.map$name <- cov.probes[, paste(gene_name, gene_id, transcript_id, "C", probe_id, sep="_")]
gen26.sel1.txw.cov.seq <- subseq(tx.seqs[seqnames(gen26.sel1.txw.cov300)],
                                 start(gen26.sel1.txw.cov300),
                                 end(gen26.sel1.txw.cov300))
names(gen26.sel1.txw.cov.seq) <- cov.probes[, paste(gene_name, gene_id, transcript_id, "C", probe_id, sep="_")]

##
gen26.sel1.txw.cos300 <- gen26.sel1.txw[cos.probes$probe_id]+120
start(gen26.sel1.txw.cos300[start(gen26.sel1.txw.cos300)<0]) <- 1
tx.len.cos300 <- tx.lens[as.character(seqnames(gen26.sel1.txw.cos300))]
end(gen26.sel1.txw.cos300[end(gen26.sel1.txw.cos300)>tx.len.cos300]) <- tx.len.cos300[end(gen26.sel1.txw.cos300)>tx.len.cos300]
gen26.sel1.txw.cos300.map <- mapFromTranscripts(gen26.sel1.txw.cos300, gen26.sel1.exl)
gen26.sel1.txw.cos300.map$name <- cos.probes[, paste(gene_name, gene_id, transcript_id, "M", probe_id, sep="_")]
gen26.sel1.txw.cos.seq <- subseq(tx.seqs[seqnames(gen26.sel1.txw.cos300)],
                                 start(gen26.sel1.txw.cos300),
                                 end(gen26.sel1.txw.cos300))
names(gen26.sel1.txw.cos.seq) <- cos.probes[, paste(gene_name, gene_id, transcript_id, "M", probe_id, sep="_")]



fin.map <- sort(c(gen26.sel1.txw.cov300.map, gen26.sel1.txw.cos300.map)[,3])
fin.seq <- c(gen26.sel1.txw.cov.seq, gen26.sel1.txw.cos.seq)[fin.map$name]
writeXStringSet(fin.seq, "out/tag_v1.fa")
writeXStringSet(tx.seqs, "out/tx_sequences.fa")
export.bed(fin.map, "out/tag_v1.bed")




