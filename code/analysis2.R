suppressMessages({
    library(data.table)
    library(stringr)
    library(GenomicRanges)
    library(GenomicFeatures)
    library(rtracklayer)
})

## Agilent V4
{
agi_v4.hg19 <- import("./data/agi_v4.hg19.bed")
tmp <- liftOver(agi_v4.hg19, import.chain("./data/GRCh37To38.over.chain"))
ok.map <- lengths(width(tmp))==1
tmp <- unlist(tmp[ok.map])
agi_v4.hg19.ok <- agi_v4.hg19[ok.map]
agi_v4.hg38 <- tmp[(width(agi_v4.hg19.ok)==width(tmp)) & (seqnames(agi_v4.hg19.ok)==seqnames(tmp))]
saveRDS(agi_v4.hg38, "./data/agi_v4.hg38.rds")
}

## COSMIC hotspots
{
cos81 <- fread("data/CosmicMutantExport.tsv")
cos81 <- cos81[GRCh==38 & `Mutation somatic status` %in%  c("Confirmed somatic variant",
                                                            "Reported in another cancer sample as somatic")]
cos81[,gene_name:=str_match(`Gene name`, "[^_]*")[,1]]
cos81.n <- cos81[,.N,by=.(gene_name, study_id=ID_STUDY, pos=`Mutation genome position`, site=`Primary site`)]
tmp <- str_match(cos81.n$pos, "(.*):(.*)-(.*)")[,2:4]
storage.mode(tmp) <- "integer"
width.ok <- tmp[,3]-tmp[,2] <= 3
cos81.n <- cos81.n[width.ok]
tmp <- tmp[width.ok,]
chr <- paste0("chr", str_replace_all(tmp[,1], c("23"="X", "24"="Y", "25"="M")))
cos81.n.rng <- GRanges(chr, IRanges(tmp[,2], tmp[,3]))
cos81.n.rng.red <- reduce(cos81.n.rng, min.gapwidth=3)
cos81.hits <- findOverlaps(cos81.n.rng, cos81.n.rng.red)
cos81.n$cluster <- subjectHits(cos81.hits)
cos81.n.tally <- cos81.n[,.(n.total=sum(N),
                            n.max=max(N),
                            n.pos=.N,
                            n.site=length(unique(site)),
                            n.study=length(unique(study_id))), by=.(cluster)]
cos81.n.tally$select <- FALSE
cos81.n.tally[n.study>2 & ( (n.site >2 & (n.max>3) & n.total>10) |
                            (n.site==2 & (n.max>6)) |
                            (n.site==1 & (n.max>9))
),  select:=TRUE]

cos81.n.rng.red <- cos81.n.rng.red[cos81.n.tally$cluster]
mcols(cos81.n.rng.red) <- cos81.n.tally
rm(cos81.n.tally)
rm(cos81.n)
rm(cos81)
saveRDS(cos81.n.rng.red, "data/cos81.n.rng.red.rds")
}

## transcript set
{
gen26 <- import("/mctp/projects/rnascape/build/refs/build/motr.v2/output/gencode.v26.basic.annotation.tags.gtf")
gen26.tx <- gen26[gen26$type=="transcript"]
gen26.tx.prot <- gen26.tx[gen26.tx$gene_type=="protein_coding" & gen26.tx$transcript_type=="protein_coding",
                     c("gene_id", "gene_name", "transcript_id", "level", "transcript_support_level", "tags")]
gen26.prot.tx.tbl <- as.data.table(mcols(gen26.tx.prot))
gen26.prot.tx.tbl[,transcript_support_level:=ifelse(transcript_support_level=="NA", NA_character_, transcript_support_level)]

rem0 <- gen26.prot.tx.tbl
sel1 <- grepl("appris_principal_1", rem0$tags)
set1 <- rem0[sel1]
rem1 <- rem0[!(gene_id %in% unique(set1$gene_id))]
set1$set <- 1

sel2 <- grepl("appris_principal_2", rem1$tags)
set2 <- rem1[sel2]
rem2 <- rem1[!(gene_id %in% unique(set2$gene_id))]
set2$set <- 2

sel3 <- grepl("appris_principal_3", rem2$tags)
set3 <- rem2[sel3]
rem3 <- rem2[!(gene_id %in% unique(set3$gene_id))]
set3$set <- 3

sel4 <- grepl("appris_principal_4", rem3$tags)
set4 <- rem3[sel4]
rem4 <- rem3[!(gene_id %in% unique(set4$gene_id))]
set4$set <- 4

sel5 <- grepl("appris_principal_5", rem4$tags)
set5 <- rem4[sel5]
rem5 <- rem4[!(gene_id %in% unique(set5$gene_id))]
set5$set <- 5

sel6 <- grepl("CCDS", rem5$tags)
set6 <- rem5[sel6]
rem6 <- rem5[!(gene_id %in% unique(set6$gene_id))]
set6$set <- 6

sel7 <- rem6$transcript_support_level %in% c("1", "2")
set7 <- rem6[sel7]
rem7 <- rem6[!(gene_id %in% unique(set7$gene_id))]
set7$set <- 7

gen26.ex <- gen26[gen26$type=="exon"]
gen26.exl <- split(gen26.ex, gen26.ex$transcript_id)
tx.lens <- sum(width(gen26.exl))

setX <- rbind(set1, set2, set3, set4, set5, set6, set7)
sel1 <- setdiff(setX[,head(.SD, 1), by=gene_id]$transcript_id, names(tx.lens)[tx.lens<=180])
sel3 <- setdiff(setX[,head(.SD, 3), by=gene_id]$transcript_id, names(tx.lens)[tx.lens<=180])
gen26.sel1 <- gen26[(gen26$transcript_id %in% sel1)]
gen26.sel3 <- gen26[(gen26$transcript_id %in% sel3)]
saveRDS(gen26.sel1, "data/gen26.sel1.rds")
saveRDS(gen26.sel3, "data/gen26.sel3.rds")

}

##
{
gen26.sel1.ex <- gen26.sel1[gen26.sel1$type=="exon"]
gen26.sel1.exl <- split(gen26.sel1.ex, gen26.sel1.ex$transcript_id)
tx.lens <- sum(width(gen26.sel1.exl))

tmp <- lapply(seq_along(tx.lens), function(i) {
    tx.len <- tx.lens[[i]]
    GRanges(names(tx.lens)[i], IRanges(seq(61, tx.len-120, by=60), width=60))
})
gen26.sel1.txw <- unlist(GRangesList(tmp))
gen26.sel1.txw.map <- mapFromTranscripts(gen26.sel1.txw, gen26.sel1.exl)
saveRDS(gen26.sel1.txw, "data/gen26.sel1.txw.rds")
saveRDS(gen26.sel1.txw.map, "data/gen26.sel1.txw.map.rds")
}


## Expression
meta <- readRDS("/proj/study/met500/data/prod/expression/M.meta.rds")
meta.pair <- meta[sample_source %in% meta[,.N,by=sample_source][N==2]$sample_source]
meta.pair$assay <- str_match(meta.pair$run.id, "(.*)-(capt|poly)-.*")[,3]
meta.pair <- meta.pair[order(sample_source, assay)]
all.bws <- list.files("/mctp/projects/rnascape/data/samples/pipeline/grch38.2/mctp/", pattern="*.bw",
                      recursive=TRUE, full.name=TRUE)
names(all.bws) <- basename(dirname(all.bws))
bws <- all.bws[meta.pair$run.id]

t1 <- split(gen26.sel1.txw.map, seqnames(gen26.sel1.txw))
t2 <- relist(unlist(t1)-1, t1)
t3 <- intersect(t2, gen26.sel1.exl)
t4 <- unlist(t3)
t1t <- unlist(t1)
t1t <- GRanges(names(t1t), ranges(t1t), strand(t1t))
t3t <- unlist(t3)
t3t <- GRanges(names(t3t), ranges(t3t), strand(t3t))
t13h <- findOverlaps(t1t, t3t)
saveRDS(t13h, "data/t13h.rds")
saveRDS(meta.pair, "data/meta.pair.rds")

{
system2("date")
tmp <- mclapply(bws, function(bw) {
    score <- unlist(suppressWarnings(summary(BigWigFile(bw), t4, type="mean", defaultValue=0)))$score
}, mc.cores=8)
system2("date")
}
t4.score <- data.table(data.frame(tmp))
saveRDS(t4.score, "data/t4.score.rds")
rm(tmp)


{ ## coverage
t4.score.sum <- as.matrix(rowsum(t4.score, queryHits(t13h)))
t4.score.n <- rowsum(rep(1, length(t13h)), queryHits(t13h))[,1]
t4.score.avg <- t4.score.sum / t4.score.n
rm(t4.score.sum)
rm(t4.score)
gc()
t4.score.avg.tx <- rowsum(copy(t4.score.avg), as.character(seqnames(gen26.sel1.txw))) # R bug
t4.score.avg.tx.idx <- rep(1:nrow(t4.score.avg.tx), runLength(seqnames(gen26.sel1.txw)))
t4.score.avg.tx.long <- t4.score.avg.tx[t4.score.avg.tx.idx,]
t4.score.avg.nrm <- t4.score.avg / t4.score.avg.tx.long
t4.score.avg.nrm[is.na(t4.score.avg.nrm)] <- 0.0
rm(t4.score.avg.tx.long)
rm(t4.score.avg.tx)
rm(t4.score.avg)
t4.score.avg.nrm.cov <- data.table(
    capt.cov=rowMeans(t4.score.avg.nrm[,seq(1,544,2)]),
    poly.cov=rowMeans(t4.score.avg.nrm[,seq(2,544,2)])
)
saveRDS(t4.score.avg.nrm, "data/t4.score.avg.nrm.rds")
saveRDS(t4.score.avg.nrm, "data/t4.score.avg.nrm.cov.rds")

}
t4.score <- readRDS("data/t4.score.rds") # counts of reads overlapping sub-windows in 744 samples
t13h <- readRDS("data/t13h.rds") # mapping from windows to sub-windows
meta.pair <- readRDS("data/meta.pair.rds") # metadata for the 744 samples

## restart
