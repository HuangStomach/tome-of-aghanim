library("tsne")
# 进行tSNE降维可视化分析
print(date())
latent <- read.table("../../Output/1/1.TCGA_latent.tsv", as.is = T, header = T)
tsne(latent) -> tpc

load("../../Output/1/tcga_ss_mat.RData")
cancer_types <- dir("../../Data/TCGA/")
sapply(cancer_types, nchar) -> ii
cancer_types <- cancer_types[which(ii <= 4)]
cancer_types <- setdiff(cancer_types, c("FPPP", "LUNG"))

##### create the color vector
library("RColorBrewer")
c(
    brewer.pal(n = 9, name = "Paired"),
    brewer.pal(n = 8, name = "BrBG"),
    brewer.pal(n = 9, name = "Set3"),
    brewer.pal(n = 8, name = "RdBu"),
    brewer.pal(n = 8, name = "RdGy"),
    brewer.pal(n = 8, name = "PiYG"),
    brewer.pal(n = 8, name = "PuBu")
) -> cc

png("C.TCGA-tSNE.png", width = 1366, height = 768, bg = "white")
plot(
    tpc[, 1], tpc[, 2],
    pch = 20, cex = 1.5,
    xlim = c(-100, 150), ylim = c(-100, 100),
    xlab = "tSNE_1", ylab = "tSNE_2",
    col = "white"
)

samples <- tcga_ss_mat[, 2]
for (k in seq_len(length(cancer_types))) {
    which(samples == cancer_types[k]) -> ii
    points(
        tpc[ii, 1], tpc[ii, 2],
        col = cc[k], pch = 20, cex = 1.5
    )
}

legend("topright", fill = cc, legend = cancer_types)
dev.off()
print(date())
