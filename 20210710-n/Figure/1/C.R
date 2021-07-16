library("tsne")
# 进行tSNE降维可视化分析

#latent <- read.table("../../Output/1/1.CCLE_latent.tsv", as.is = T, header = T)
latent <- read.table("../../Output/0/V15.CCLE.4VAE.RANK.tsv", as.is = T, header = T)
tsne(latent) -> tpc
sapply(rownames(latent), function(x) {
    strsplit(x, split = "\\.")[[1]][1] -> u
    strsplit(u, split = "_")[[1]] -> v
    v <- v[-1]
    paste(v, collapse = "_")
}) -> tt
tissues <- sort(unique(tt))

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

png("C.CCLE-tSNE.png", width = 1366, height = 768, bg = "white")
plot(
    tpc[, 1], tpc[, 2],
    pch = 20, cex = 1.5,
    xlim = c(-80, 140), ylim = c(-80, 80),
    xlab = "tSNE_1", ylab = "tSNE_2",
    col = "white"
)

for (k in seq_len(length(tissues))) {
    which(tt == tissues[k]) -> ii
    points(
        tpc[ii, 1], tpc[ii, 2],
        col = cc[k], pch = 20, cex = 1.5
    )
}

legend("topright", fill = cc, legend = tissues)
dev.off()