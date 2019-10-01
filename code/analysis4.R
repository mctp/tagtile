suppressMessages({
    library(data.table)
    library(stringr)
    library(GenomicRanges)
    library(GenomicFeatures)
    library(rtracklayer)
    library(BSgenome.Hsapiens.UCSC.hg38)
})

{

{ ## prior data
goi <- c(readLines("data/gene_list.txt"), "T", "PTPN14") ## genes of interest
cos81.n.rng.red <- readRDS("data/cos81.n.rng.red.rds") ## cosmic hotspots 1.5M total select 10k
agi_v4.hg38  <- readRDS("data/agi_v4.hg38.rds") ## AGI v4 probes in HG38 coordinates
gen26.sel1 <- readRDS("data/gen26.sel1.rds") ## 1 transcript per protein-coding gene
gen26.sel1.txw <- readRDS("data/gen26.sel1.txw.rds") ## valid transcript windows
gen26.sel1.txw.map <- readRDS("data/gen26.sel1.txw.map.rds") ## valid transcript windows on genome (xHits-window, w/ introns)
t4.score.avg.nrm <- readRDS("data/t4.score.avg.nrm.rds") ## per window coverage
}

{ ## get transcript sequences and lengths
gen26.sel1.ex <- gen26.sel1[gen26.sel1$type=="exon"]
gen26.sel1.exl <- split(gen26.sel1.ex, gen26.sel1.ex$transcript_id)
tx.seqs <- extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg38, gen26.sel1.exl)
tx.lens <- sum(width(gen26.sel1.exl))
}


{ ## coverage data
tmp.1 <- cos81.n.rng.red[cos81.n.rng.red$select]
tmp.2 <- findOverlaps(gen26.sel1.txw.map, tmp.1)
tmp.3 <- data.table(queryHits(tmp.2), tmp.1[subjectHits(tmp.2)]$n.max)[,max(V2),by=V1]
setkey(tmp.3, V1)
tmp.4 <- tmp.3[J(1:length(gen26.sel1.txw.map)),]$V1.1
tmp.5 <- ifelse(is.na(tmp.4), 0, tmp.4)
tmp.6 <- data.table(
    capt.cov=rowMeans(t4.score.avg.nrm[,seq(1,544,2)]),
    poly.cov=rowMeans(t4.score.avg.nrm[,seq(2,544,2)])
)
tmp.7 <- data.table(
    transcript_id = as.character(seqnames(gen26.sel1.txw)),
    probe_id = 1:length(gen26.sel1.txw.map),
    cov = sqrt(tmp.6$capt.cov*tmp.6$poly.cov),
    cds = gen26.sel1.txw.map %over% gen26.sel1[gen26.sel1$type=="CDS"],
    agi = gen26.sel1.txw.map %over% agi_v4.hg38,
    cos = tmp.5
)
tmp.7[,probe_n:=1:.N, by=transcript_id]
setkey(tmp.7, transcript_id)
tmp.8 <- unique(as.data.table(mcols(gen26.sel1)[,c("gene_id", "gene_name", "transcript_id")]))
setkey(tmp.8, transcript_id)
cov.probes <- tmp.7[tmp.8]
setkey(cov.probes, transcript_id, probe_n)
cov.probes.1 <- cov.probes[order(-cov, -agi, -cds), head(.SD, 1), by=transcript_id]
cov.probes.1.filter <- cov.probes.1[,.(transcript_id, probe_n,
                                       probe_nm1=probe_n-1, probe_np1=probe_n+1,
                                       probe_nm2=probe_n-2, probe_np2=probe_n+2
                                       )]
setkey(cov.probes.1.filter, transcript_id, probe_n)
cov.probes <- cov.probes[!cov.probes.1.filter]
setkey(cov.probes.1.filter, transcript_id, probe_nm1)
cov.probes <- cov.probes[!cov.probes.1.filter]
setkey(cov.probes.1.filter, transcript_id, probe_np1)
cov.probes <- cov.probes[!cov.probes.1.filter]
setkey(cov.probes.1.filter, transcript_id, probe_nm2)
cov.probes <- cov.probes[!cov.probes.1.filter]
setkey(cov.probes.1.filter, transcript_id, probe_np2)
cov.probes <- cov.probes[!cov.probes.1.filter]
cov.probes.2 <- cov.probes[order(-cov, -agi, -cds), head(.SD, 1), by=transcript_id]
cov.probes.x <- rbind(cov.probes.1, cov.probes.2)
cos.probes.x <- cov.probes[!(probe_id %in% cov.probes.x$probe_id)][gene_name %in% goi][cos>0]
}

{
gen26.sel1.txw.cov300 <- gen26.sel1.txw[cov.probes.x$probe_id]+120
start(gen26.sel1.txw.cov300[start(gen26.sel1.txw.cov300)<0]) <- 1
tx.len.cov300 <- tx.lens[as.character(seqnames(gen26.sel1.txw.cov300))]
end(gen26.sel1.txw.cov300[end(gen26.sel1.txw.cov300)>tx.len.cov300]) <- tx.len.cov300[end(gen26.sel1.txw.cov300)>tx.len.cov300]
gen26.sel1.txw.cov300.map <- mapFromTranscripts(gen26.sel1.txw.cov300, gen26.sel1.exl)
gen26.sel1.txw.cov300.map$name <- cov.probes.x[, paste(gene_name, gene_id, transcript_id, "C", probe_id, sep="_")]
gen26.sel1.txw.cov300$name <- cov.probes.x[, paste(gene_name, gene_id, transcript_id, "C", probe_id, sep="_")]
gen26.sel1.txw.cov.seq <- subseq(tx.seqs[seqnames(gen26.sel1.txw.cov300)],
                                 start(gen26.sel1.txw.cov300),
                                 end(gen26.sel1.txw.cov300))
names(gen26.sel1.txw.cov.seq) <- cov.probes.x[, paste(gene_name, gene_id, transcript_id, "C", probe_id, sep="_")]
gen26.sel1.txw.cos300 <- gen26.sel1.txw[cos.probes.x$probe_id]+120
start(gen26.sel1.txw.cos300[start(gen26.sel1.txw.cos300)<0]) <- 1
tx.len.cos300 <- tx.lens[as.character(seqnames(gen26.sel1.txw.cos300))]
end(gen26.sel1.txw.cos300[end(gen26.sel1.txw.cos300)>tx.len.cos300]) <- tx.len.cos300[end(gen26.sel1.txw.cos300)>tx.len.cos300]
gen26.sel1.txw.cos300.map <- mapFromTranscripts(gen26.sel1.txw.cos300, gen26.sel1.exl)
gen26.sel1.txw.cos300.map$name <- cos.probes.x[, paste(gene_name, gene_id, transcript_id, "M", probe_id, sep="_")]
gen26.sel1.txw.cos300$name <- cos.probes.x[, paste(gene_name, gene_id, transcript_id, "M", probe_id, sep="_")]
gen26.sel1.txw.cos.seq <- subseq(tx.seqs[seqnames(gen26.sel1.txw.cos300)],
                                 start(gen26.sel1.txw.cos300),
                                 end(gen26.sel1.txw.cos300))
names(gen26.sel1.txw.cos.seq) <- cos.probes.x[, paste(gene_name, gene_id, transcript_id, "M", probe_id, sep="_")]
}

{
##
tag.v1.map <- sort(c(gen26.sel1.txw.cov300.map, gen26.sel1.txw.cos300.map)[,3])
tag.v1.seq <- c(gen26.sel1.txw.cov.seq, gen26.sel1.txw.cos.seq)[tag.v1.map$name]
tag.v1 <- c(gen26.sel1.txw.cov300, gen26.sel1.txw.cos300)
names(tag.v1) <- tag.v1$name
tag.v1 <- tag.v1[tag.v1.map$name]
names(tag.v1) <- NULL
tx.v1.seq <- extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg38, gen26.sel1.exl)
}

}





tag.v1.agi <- fread("./review/from_agilent/tag_v1_Unique_Probes_1-Probe_per_target_CombineRest_Suredesign.txt")
tmp.1 <- as.data.table(as.data.frame(tag.v1))[,.(transcript_id=seqnames, tx_start=1, tx_end=tx.lens[seqnames],target_id=str_replace_all(name, "\\.", "_"), target_start=start, target_end=end)]
tmp.2 <- tag.v1.agi[, .(TargetID, ProbeID, Start, Stop, Sequence)]
tmp.2 <- tmp.2[,.(target_id=TargetID, probe_id=ProbeID, probe_start=Start+1, probe_end=Stop, probe_seq=Sequence)]
setkey(tmp.1, target_id)
setkey(tmp.2, target_id)
tmp.3 <- tmp.1[tmp.2]
tmp.4 <- data.table(transcript_id=names(tx.lens), tx_start=1, tx_end=tx.lens, target_id=names(tx.lens), target_start=1, target_end=tx.lens)
tx.v1.agi <- fread("./review/from_agilent/TX_Dust_GC_unique_Probe1_Probe2_CombineRest_Suredesign.txt")
tmp.5 <- tx.v1.agi[,.(target_id=str_replace(TargetID, "_", "."), probe_id=ProbeID, probe_start=Start+1, probe_end=Stop, probe_seq=Sequence)]
setkey(tmp.4, target_id)
setkey(tmp.5, target_id)
tmp.6 <- tmp.4[tmp.5]
tmp.7 <- rbind(tmp.3, tmp.6)
tmp.7[,":="(probe_tx_start=probe_start+target_start-1, probe_tx_end=probe_end+target_start-1)]
tmp.7$transcript_seq <- as.character(tx.seqs[tmp.7$transcript_id])
tmp.7$probe_seq_verify <- str_sub(tmp.7$transcript_seq, tmp.7$probe_tx_start, tmp.7$probe_tx_end)
tmp.7$transcript_seq <- NULL
all(tmp.7[,probe_seq==probe_seq_verify])
tmp.7$probe_seq_verify <- NULL
tag.v1.fin <- tmp.7


##
tag.v1.fin.probe <- with(tag.v1.fin, GRanges(transcript_id, IRanges(probe_tx_start, probe_tx_end)))
tag.v1.fin.probe.map <- mapFromTranscripts(tag.v1.fin.probe, gen26.sel1.exl)
tag.v1.fin.probe.exl <- gen26.sel1.exl[as.character(seqnames(tag.v1.fin.probe))]
tag.v1.fin.probe.tgt <- pintersect(tag.v1.fin.probe.map, tag.v1.fin.probe.exl, drop.nohit.ranges=TRUE)
tag.v1.fin.probe.tgt.flat <- unlist(tag.v1.fin.probe.tgt)
tag.v1.fin.probe.tgt.flat$probe_id <- rep(tag.v1.fin$probe_id, times=lengths(tag.v1.fin.probe.tgt))
tag.v1.fin.probe.tgt.flat <- unname(sort(tag.v1.fin.probe.tgt.flat)[,c("probe_id", "gene_id", "gene_name", "transcript_id")])
export(tag.v1.fin.probe.tgt.flat, "review/tag_v1_fin_probe_tgt_flat.gtf")
fwrite(as.data.table(tag.v1.fin.probe.tgt.flat)[,.(GeneID=gene_id, Chr=seqnames, Start=start, End=end, Strand=strand)],
       "visit/data/tagv1.sel.tagtile.genome.tgt.flat.gene.saf", sep="\t")
fwrite(as.data.table(tag.v1.fin.probe.tgt.flat)[,.(GeneID=probe_id, Chr=seqnames, Start=start, End=end, Strand=strand)],
       "visit/data/tagv1.sel.tagtile.genome.tgt.flat.probe.saf", sep="\t")
