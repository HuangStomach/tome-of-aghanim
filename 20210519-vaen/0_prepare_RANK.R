### transcriptome data
ccle <- read.table("./Data/CCLE/CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct",
    skip = 2, sep = "\t", header = T, as.is = T)
any_tcga_rpkm <- read.delim("./Data/TCGA/ACC/HiSeqV2", as.is = T)
genes <- intersect(any_tcga_rpkm[, 1], ccle[, 2]) # 寻找两文件数据中具有交集的基因名称
# > print(length(genes)) [1] 18217

new_ccle <- ccle[match(genes, ccle[, 2]), ] # 筛选出存在于交集中的行
ccle_gene_mat <- as.matrix(new_ccle[, c(-1, -2)]) # 将数据转换为矩阵
rownames(ccle_gene_mat) <- new_ccle[, 2] # 将第二行的内容设置为矩阵的行名(基因名称)

apply(ccle_gene_mat, 1, mean) -> row_mean # 按行进行平均值计算
which(row_mean > 1) -> expressed # 筛选出平均值（基因表达）大于1的行索引
ccle_gene_mat <- ccle_gene_mat[expressed, ] # 数据过滤
# > dim(ccle_gene_mat) [1] 12406  1156

# choose tissues with >= 20 cell lines
sapply(colnames(ccle_gene_mat), function(x) { # 癌细胞类型进行格式化，保留部分信息
    strsplit(x, split = "\\.")[[1]][1] -> u
    strsplit(u, split = "_")[[1]] -> v
    v <- v[-1]
    paste(v, collapse = "_")
}) -> tt

names(tt) <- NULL
table(tt) -> t1 # 构建计算每个元素出现次数的table
names(which(t1 >= 20)) -> tissues_int # 保留数据量多余20的组织信息
which(tt %in% tissues_int) -> ii # 筛选出符合条件的索引
ccle_gene_mat <- ccle_gene_mat[, ii]
# > dim(ccle_gene_mat) [1] 12406  1100


log2.ccle_gene_mat <- log2(ccle_gene_mat + 1) # 基因表达有较大的数据 使用log2-transform规范化

# choose genes that are most variably expressed
apply(log2.ccle_gene_mat, 1, var) -> rowVar
names(which(rowVar > median(rowVar))) -> mad.genes
length(mad.genes)
mad5000.ccle_gene_mat <- log2.ccle_gene_mat[mad.genes, ] ### gene by sample
dim(mad5000.ccle_gene_mat)

# > length(mad.genes)
# [1] 6203

# > dim(mad5000.ccle_gene_mat)
# [1] 6203 1100

##########################################################################
##########################################################################
##########################################################################

cancer.types <- dir("./TCGA/")
sapply(cancer.types, nchar) -> ii
cancer.types <- cancer.types[which(ii <= 4)]
cancer.types <- setdiff(cancer.types, c("FPPP", "LUNG"))


RPKM.mat <- c()
cancer.type.list <- list()
for (k in 1:length(cancer.types)) {
    cancer <- cancer.types[k]
    original.TCGA.RPKM <- read.delim(paste("./TCGA/", cancer, "/HiSeqV2", sep = ""), as.is = T)

    ### exclude genes with rowSum == 0
    apply(original.TCGA.RPKM[, -1], 1, sum) -> rowCheck
    non0.TCGA.RPKM <- original.TCGA.RPKM[which(rowCheck != 0), ]

    shared.TCGA.RPKM <- non0.TCGA.RPKM[match(mad.genes, non0.TCGA.RPKM[, 1]), -1]
    rownames(shared.TCGA.RPKM) <- mad.genes

    apply(shared.TCGA.RPKM, 1, sum) -> check
    shared.TCGA.RPKM <- shared.TCGA.RPKM[!is.na(check), ]

    if (k == 1) {
        cur.genes <- rownames(shared.TCGA.RPKM)
    } else {
        cur.genes <- intersect(cur.genes, rownames(shared.TCGA.RPKM))
    }

    t.shared.TCGA.RPKM <- t(shared.TCGA.RPKM)
    RPKM.mat <- rbind(RPKM.mat[, cur.genes], t.shared.TCGA.RPKM[, cur.genes]) ### cat by samples, columns are mad.genes
    cancer.type.list[[cancer]] <- colnames(shared.TCGA.RPKM)
    cat(cancer, "\t", ncol(shared.TCGA.RPKM), " ", ncol(RPKM.mat), " ", nrow(RPKM.mat), "\n", sep = "")
}

# ACC     79 6203 79
# BLCA    426 6203 505
# BRCA    1218 6203 1723
# CESC    308 6203 2031
# CHOL    45 6203 2076
# COAD    329 6203 2405
# DLBC    48 6203 2453
# ESCA    196 6203 2649
# GBM     172 6203 2821
# HNSC    566 6203 3387
# KICH    91 6203 3478
# KIRC    606 6203 4084
# KIRP    323 6203 4407
# LAML    173 6203 4580
# LGG     530 6203 5110
# LIHC    423 6203 5533
# LUAD    576 6203 6109
# LUSC    553 6203 6662
# MESO    87 6203 6749
# OV      308 6203 7057
# PAAD    183 6203 7240
# PCPG    187 6203 7427
# PRAD    550 6203 7977
# READ    105 6203 8082
# SARC    265 6203 8347
# SKCM    474 6203 8821
# STAD    450 6203 9271
# TGCT    156 6203 9427
# THCA    572 6203 9999
# THYM    122 6203 10121
# UCEC    201 6203 10322
# UCS     57 6203 10379
# UVM     80 6203 10459

##########################################################################
genes2 <- intersect(mad.genes, colnames(RPKM.mat))

ccle.train.mat <- mad5000.ccle_gene_mat[genes2, ]
print(dim(ccle.train.mat))

### rank, per-sample
rank.ccle_gene_mat <- apply(ccle.train.mat, 2, function(u) rank(u) / length(u))
rank.ccle_gene_mat <- apply(rank.ccle_gene_mat, 1, function(u) {
    u[which(u == 1)] <- 6162.5 / 6163
    u
})

### p to z
scaled.ccle_gene_mat <- apply(rank.ccle_gene_mat, 2, function(u) {
    qnorm(u)
})
print(dim(scaled.ccle_gene_mat))

# > length(genes2)
# [1] 6163
# > print(dim(scaled.ccle_gene_mat))
# [1] 1100 6163

##########################################################################
rank.RPKM.mat <- apply(RPKM.mat, 1, function(u) rank(u) / length(u))
rank.RPKM.mat <- apply(rank.RPKM.mat, 1, function(u) {
    u[which(u == 1)] <- 6162.5 / 6163
    u
})

scaled.RPKM.mat <- apply(rank.RPKM.mat, 2, function(u) qnorm(u))

dimnames(scaled.RPKM.mat) <- dimnames(RPKM.mat)

##########################################################################
### dataset for NOPEER.NO01.Sigmoid
write.table(scaled.ccle_gene_mat, file = paste("V15.CCLE.4VAE.RANK.tsv", sep = ""), row.names = T, quote = F, sep = "\t")
write.table(scaled.RPKM.mat, file = paste("V15.TCGA.4VAE.RANK.tsv", sep = ""), row.names = T, quote = F, sep = "\t")
##########################################################################