library(biomaRt)
library(data.table)
library(stringr)
library(rtracklayer)

refseq.pref <- readLines("/mctp/users/mcieslik/proj/code/tpo/rlibs/carat/inst/extdata/canonical_refseq.txt")
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
ens.pref <- getBM(attributes=c("refseq_mrna", "ensembl_transcript_id"), filters="refseq_mrna", values=refseq.pref,
                  mart=ensembl)

{
gen26 <- import("/mctp/projects/rnascape/build/refs/build/motr.v2/output/gencode.v26.basic.annotation.tags.gtf")
gen26.tx <- gen26[gen26$type=="transcript"]
gen26.tx.prot <- gen26.tx[gen26.tx$gene_type=="protein_coding" & gen26.tx$transcript_type=="protein_coding",
                     c("gene_id", "gene_name", "transcript_id", "level", "transcript_support_level", "tags")]
gen26.prot.tx.tbl <- as.data.table(mcols(gen26.tx.prot))
gen26.prot.tx.tbl[,transcript_support_level:=ifelse(transcript_support_level=="NA", NA_character_, transcript_support_level)]
gen26.prot.tx.tbl[,txid:=str_match(transcript_id, "[^.]*")]

tagv3.raw <- fread("tagv2/5x_combined_dupesRemoved.txt")
setnames(tagv3.raw, c("transcript_id", "probe_id", "seq", "beg", "end"))

####
rem0 <- gen26.prot.tx.tbl
selx <- rem0$txid %in% ens.pref$ensembl_transcript_id & rem0$transcript_id %in% unique(tagv3.raw$transcript_id)
setx <- rem0[selx]
remx <- rem0[!(gene_id %in% unique(setx$gene_id))]
setx$set <- 0

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

setX <- rbind(setx, set1, set2, set3, set4, set5, set6, set7)
sel1 <- setdiff(setX[,head(.SD, 1), by=gene_id]$transcript_id, names(tx.lens)[tx.lens<=180])
sel3 <- setdiff(setX[,head(.SD, 3), by=gene_id]$transcript_id, names(tx.lens)[tx.lens<=180])
gen26.sel1 <- gen26[(gen26$transcript_id %in% sel1)]
gen26.sel3 <- gen26[(gen26$transcript_id %in% sel3)]
saveRDS(gen26.sel1, "data/gen26.sel1x.rds")
saveRDS(gen26.sel3, "data/gen26.sel3x.rds")
}
