library(data.table)
library(stringr)
library(GenomicRanges)
library(rtracklayer)

## Agilent V4
agi_v4.hg19 <- import("./data/agi_v4.hg19.bed")
tmp <- liftOver(agi_v4.hg19, import.chain("/refs/genomes/GRCh37To38.over.chain"))
ok.map <- lengths(width(tmp))==1
tmp <- unlist(tmp[ok.map])
agi_v4.hg19.ok <- agi_v4.hg19[ok.map]
agi_v4.hg38 <- tmp[(width(agi_v4.hg19.ok)==width(tmp)) & (seqnames(agi_v4.hg19.ok)==seqnames(tmp))]

## cosmic recurrence
cos81 <- fread("data/CosmicMutantExport.tsv")
cos81 <- cos81[GRCh==38 & `Mutation somatic status` %in%  c("Confirmed somatic variant",
                                                            "Reported in another cancer sample as somatic")]
cos81[,gene_name:=str_match(`Gene name`, "[^_]*")[,1]]
{
cos81.n <- cos81[,.N,by=.(gene_name, study_id=ID_STUDY, pos=`Mutation genome position`, site=`Primary site`)]
tmp <- str_match(cos81.n$pos, "(.*):(.*)-(.*)")[,2:4]
storage.mode(tmp) <- "integer"
width.ok <- tmp[,3]-tmp[,2] <= 6
cos81.n <- cos81.n[width.ok]
tmp <- tmp[width.ok,]
chr <- paste0("chr", str_replace_all(tmp[,1], c("23"="X", "24"="Y", "25"="M")))
cos81.n.rng <- GRanges(chr, IRanges(tmp[,2], tmp[,3]))
cos81.n.rng.red <- reduce(cos81.n.rng, min.gapwidth=6)
cos81.hits <- findOverlaps(cos81.n.rng, cos81.n.rng.red)
}
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
mcols(cos81.n.rng.red) <- cos81.n.tally
tmp <- data.table(cluster=subjectHits(cos81.hits), gene_name=cos81.n[queryHits(cos81.hits)]$gene_name)
tmp <- tmp[,.(gene_names=paste(unique(gene_name), collapse=":")),keyby=cluster]
cos81.n.rng.red[order(cos81.n.rng.red$cluster)]$gene_names <- tmp$gene_names

## exon coverage
gtf.prot <- import("/proj/code/crisp/0_references/build/motr.v2/output/motr.v2-prot.gtf")
gtf.ex <- gtf.prot[gtf.prot$type=="exon"]
gtf.tx <- gtf.prot[gtf.prot$type=="transcript"]
gtf.ex.red <- reduce(gtf.ex)

gtf.ex.g <- split(gtf.ex, gtf.ex$gene_id)
gtf.ex.g.d <- disjoin(gtf.ex.g)
gtf.ex.g.d <- unlist(disjoin(gtf.ex.g))
gtf.ex.g.d$gene_id <- names(gtf.ex.g.d)
names(gtf.ex.g.d) <- NULL
hits <- findOverlaps(gtf.ex, gtf.ex.g.d)
ok.self <- gtf.ex[queryHits(hits)]$gene_id==gtf.ex.g.d[subjectHits(hits)]$gene_id
hits.self <- hits[ok.self]
hits.self.dt <- data.table(as.data.frame(hits.self))
gtf.ex.g.d$exon.cov <- hits.self.dt[,.N,by=subjectHits
                                    ][order(subjectHits)]$N
exon.targets.hg38 <- gtf.ex.g.d

## expression
mx <- readRDS("/proj/study/met500/data/prod/expression/M.mx.rds")
meta <- readRDS("/proj/study/met500/data/prod/expression/M.meta.rds")
meta.pair <- meta[sample_source %in% meta[,.N,by=sample_source][N==2]$sample_source]
meta.pair$assay <- str_match(meta.pair$run.id, "(.*)-(capt|poly)-.*")[,3]
meta.pair <- meta.pair[order(sample_source, assay)]
all.bws <- list.files("/mctp/projects/rnascape/data/samples/pipeline/grch38.2/mctp/", pattern="*.bw",
                      recursive=TRUE, full.name=TRUE)
names(all.bws) <- basename(dirname(all.bws))
bws <- all.bws[meta.pair$run.id]
res <- mclapply(bws[1:400], function(bw) {
    suppressWarnings(mcols(unlist(summary(BigWigFile(bw), agi_v4.hg38, type="mean", defaultValue=0))))$score
}, mc.cores=8)
resb <- do.call(cbind, res)

## correlation coefficient
tmp <- lapply(split(t(resb), meta.pair$assay[1:400]), matrix, nrow=200, byrow=TRUE)
cc <- mclapply(1:nrow(resb), function(i) cor(tmp[[1]][,i], tmp[[2]][,i]), mc.cores = 8)
fin.cor <- unlist(cc)

## probe exon coverage
hit.ex <- as.data.table(findOverlaps(agi_v4.hg38, gtf.ex.red))
hit.ex$cov <- width(pintersect(agi_v4.hg38[hit.ex$queryHits], gtf.ex.red[hit.ex$subjectHits]))/120
hit.ex.rnk <- hit.ex[order(-cov),.SD[1],by=queryHits]

## cosmic
hit.cos <- as.data.table(findOverlaps(agi_v4.hg38, cos81.n.rng.red))
hit.cos$n.max <- cos81.n.rng.red[hit.cos$subjectHits]$n.max
hit.cos$select <- cos81.n.rng.red[hit.cos$subjectHits]$select
hit.cos.rnk <- hit.cos[order(-select, -n.max),.SD[1],by=queryHits]

## n overlapping exons
hit.ovr <- as.data.table(findOverlaps(agi_v4.hg38, exon.targets.hg38))
hit.ovr$exn <- exon.targets.hg38[hit.ovr$subjectHits]$exon.cov
hit.ovr.rnk <- hit.ovr[order(-exn),.SD[1],by=queryHits]

##
agi_v4.hg38$cor <- fin.cor
agi_v4.hg38$cov <- 0.0
mcols(agi_v4.hg38[hit.ex.rnk$queryHits])$cov <- hit.ex.rnk$cov
agi_v4.hg38$cos <- 0
mcols(agi_v4.hg38[hit.cos.rnk$queryHits])$cos <- hit.cos.rnk$n.max*hit.cos.rnk$select
agi_v4.hg38$exn <- 0
mcols(agi_v4.hg38[hit.ovr.rnk$queryHits])$exn <- hit.ovr.rnk$exn
tmp <- reduce(split(gtf.ex, gtf.ex$gene_id))
hits.gene <- as.data.table(findOverlaps(agi_v4.hg38, tmp))
hits.gene$gene_id <- names(tmp[hits.gene$subjectHits])
hits.gene.rnk <- hits.gene[,paste(gene_id, sep=":"),by=queryHits]
agi_v4.hg38$gene_id <- NA_character_
mcols(agi_v4.hg38[hits.gene.rnk$queryHits])$gene_id <- hits.gene.rnk$V1


tmp <- as.data.table(agi_v4.hg38)
tmp$idx <- 1:nrow(tmp)
agi_v4.hg38.sel <- agi_v4.hg38[tmp$idx %in% tmp[!is.na(gene_id)][order(-cov, -(cor>0.8), -exn, -(cor>0.9)),
                                                                 .SD[1], by=gene_id]$idx]
agi_v4.hg38.sel <- unique(c(agi_v4.hg38.sel, agi_v4.hg38[agi_v4.hg38$cos>10]))

tmp <- as.data.table(mcols(agi_v4.hg38.sel)[,c("gene_id", "exn")])
tmp <- tmp[complete.cases(tmp)]
setkey(tmp, gene_id)
tmp2 <- as.data.table(mcols(gtf.tx)[,c("gene_id", "transcript_id")])[,.N,by=gene_id]
setkey(tmp2, gene_id)

tmp3 <- as.data.table(mcols(agi_v4.hg38.sel)[,c("gene_id")])


tmp4 <- as.data.table(mcols(gtf.gn))[,.(gene_id, gene_name)]
setkey(tmp4, gene_id)
tmp4[tmp3[,.N,by=V1][order(-N)][!is.na(V1)][1:20]$V1]

tmp5 <- agi_v4.hg38.sel[!is.na(agi_v4.hg38.sel$gene_id)]




xxx <- ggplot(data.table(tmp[[2]][,515225], tmp[[1]][,515225]))+ aes(x=log(V1+0.01), y=log(V2+0.01)) +
    geom_point() +
    theme_bw(base_size=16) +
    xlab("polyA") +
    ylab("capture")
ggsave("blah.pdf", xxx, useDingbats=FALSE, width=4.5, height=4.5)
