## 补充脚本 生成正确的 tcga_ss_mat.RData
cancer_types <- dir("./Data/TCGA/")
sapply(cancer_types, nchar) -> ii
cancer_types <- cancer_types[which(ii <= 4)]
cancer_types <- setdiff(cancer_types, c("FPPP", "LUNG"))

tcga_ss_mat <- matrix(ncol = 2, nrow = 0)
for (cancer in cancer_types) {
    print(cancer)
    tcga_rpkm <- read.delim(
        paste("./Data/TCGA/", cancer, "/HiSeqV2", sep = ""),
    as.is = T)
    for (rpkm in colnames(tcga_rpkm)[-1]) {
        tcga_ss_mat <- rbind(tcga_ss_mat, c(rpkm, cancer))
    }
}

gsub("\\.", "-", tcga_ss_mat[, 1]) -> ss
ss -> tcga_ss_mat[, 1]
save(tcga_ss_mat, file = "./Output/1/tcga_ss_mat.RData")