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
names(which(t1 >= 20)) -> tissues_int # 保留数据量多余20的癌症组织信息
which(tt %in% tissues_int) -> ii # 筛选出符合条件的癌症组织列索引
ccle_gene_mat <- ccle_gene_mat[, ii]
# > dim(ccle_gene_mat) [1] 12406  1100

# choose genes that are most variably expressed
ccle_gene_mat <- log2(ccle_gene_mat + 1) # 基因表达有较大的数据 使用log2-transform规范化
apply(ccle_gene_mat, 1, var) -> row_var # 对每一行计算方差
names(which(row_var > median(row_var))) -> mad_genes # 方差较大的基因
mad5000_ccle_gene_mat <- ccle_gene_mat[mad_genes, ] # 筛选出方差较大的一部分数据
# > length(mad_genes) [1] 6203
# > dim(mad5000_ccle_gene_mat) [1] 6203 1100

cancer_types <- dir("./Data/TCGA/")
sapply(cancer_types, nchar) -> ii
cancer_types <- cancer_types[which(ii <= 4)]
cancer_types <- setdiff(cancer_types, c("FPPP", "LUNG")) # 筛选特定类别的癌症

rpkm_mat <- c()
for (k in seq_len(length(cancer_types))) {
    print
    cancer <- cancer_types[k]
    tcga_rpkm <- read.delim(
        paste("./Data/TCGA/", cancer, "/HiSeqV2", sep = ""),
        as.is = T
    ) # 根据癌症类别读入tcga的基因表达数据

    apply(tcga_rpkm[, -1], 1, sum) -> row_check # 对每一行的基因表达求和
    # 排除基因表达求和为0或nan的行
    non0_tcga_rpkm <- tcga_rpkm[which(row_check != 0 & !is.na(row_check)), ]

    # 筛选出其中方差较大的基因数据，并且删除其第一列
    shared_tcga_rpkm <- non0_tcga_rpkm[match(
        mad_genes, non0_tcga_rpkm[, 1]
    ), -1]
    rownames(shared_tcga_rpkm) <- mad_genes # 并把mad_genes中的基因名称作为行名

    # 如果是第一次 则把shared_tcga_rpkm的基因名称赋予cur_genes
    # 否则选出当前shared_tcga_rpkm涉及的基因名称和上一轮基因名称的交集
    if (k == 1) cur_genes <- rownames(shared_tcga_rpkm)
    else cur_genes <- intersect(cur_genes, rownames(shared_tcga_rpkm))

    t_shared_tcga_rpkm <- t(shared_tcga_rpkm) # 转置 原本行是基因列是样本
    # cat by samples, columns are mad_genes
    # 进行数据的累加，将列置为基因，行置为样本
    rpkm_mat <- rbind(rpkm_mat[, cur_genes], t_shared_tcga_rpkm[, cur_genes])
    cat(cancer, "\t",
        ncol(shared_tcga_rpkm), " ",
        ncol(rpkm_mat), " ",
        nrow(rpkm_mat), "\n",
        sep = "")
}
# 最后一直筛选有交集的基因（也就是都有数据的基因）然后对所有的癌细胞类型进行数据合并
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

genes2 <- intersect(mad_genes, colnames(rpkm_mat)) # 将ccle中的基因和tcga中的取交集
# > length(genes2) [1] 6163

ccle_train_mat <- mad5000_ccle_gene_mat[genes2, ]
# 对每一种癌症样本做rank标准化
rank_ccle_gene_mat <- apply(ccle_train_mat, 2, function(u) rank(u) / length(u))
rank_ccle_gene_mat <- apply(rank_ccle_gene_mat, 1, function(u) {
    u[which(u == 1)] <- 6162.5 / 6163
    u
}) # 没懂 可能是消除极端值

# 对样本数据做正态分布 转换为Z-scores
scaled_ccle_gene_mat <- apply(rank_ccle_gene_mat, 2, function(u) qnorm(u))
# > print(dim(scaled_ccle_gene_mat)) [1] 1100 6163

rank_rpkm_mat <- apply(rpkm_mat, 1, function(u) rank(u) / length(u))
rank_rpkm_mat <- apply(rank_rpkm_mat, 1, function(u) {
    u[which(u == 1)] <- 6162.5 / 6163
    u
})
scaled_rpkm_mat <- apply(rank_rpkm_mat, 2, function(u) qnorm(u))

dimnames(scaled_rpkm_mat) <- dimnames(rpkm_mat)

# dataset for NOPEER.NO01.Sigmoid
write.table(scaled_ccle_gene_mat,
    file = paste("V15.CCLE.4VAE.RANK.tsv", sep = ""),
    row.names = T, quote = F, sep = "\t"
)
write.table(scaled_rpkm_mat,
    file = paste("V15.TCGA.4VAE.RANK.tsv", sep = ""),
    row.names = T, quote = F, sep = "\t"
)
