suppressWarnings({
    library(data.table)
    library(stringr)
    library(GenomicRanges)
    library(GenomicFeatures)
    library(rtracklayer)
})



gen26.sel1.txp.map <- readRDS("points.rds")

## Expression
meta <- readRDS("/proj/study/met500/data/prod/expression/M.meta.rds")
meta.pair <- meta[sample_source %in% meta[,.N,by=sample_source][N==2]$sample_source]
meta.pair$assay <- str_match(meta.pair$run.id, "(.*)-(capt|poly)-.*")[,3]
meta.pair <- meta.pair[order(sample_source, assay)]
all.bws <- list.files("/mctp/projects/rnascape/data/samples/pipeline/grch38.2/mctp/", pattern="*.bw",
                      recursive=TRUE, full.name=TRUE)
names(all.bws) <- basename(dirname(all.bws))
bws <- all.bws[meta.pair$run.id]
{
system2("date")
tmp <- mclapply(bws, function(bw) {
    score <- unlist(suppressWarnings(summary(BigWigFile(bw),
             gen26.sel1.txp.map, type="min", defaultValue=0)))$score
    score <- data.table(score,
           ((1:length(score))-1) %/% 6)[, median(score),keyby=V2]$V1
}, mc.cores=32)
system2("date")
}
