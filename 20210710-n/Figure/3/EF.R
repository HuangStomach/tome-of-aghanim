library(ggplot2)
load("EF.data.RData")
sapply(rownames(tcga_sensitive_mat),
    function(u) gsub("\\.", "-", u)) -> new_name
rownames(tcga_sensitive_mat) <- unlist(new_name)
sapply(rownames(tcga_resistant_mat),
    function(u) gsub("\\.", "-", u)) -> new_name
rownames(tcga_resistant_mat) <- unlist(new_name)

two_drugs_match <- read.table(
    "../../Data/Match/drugs_match_2.txt",
    as.is = T, sep = "\t"
)

gdsc_anno <- read.delim("../../Data/GDSC/Screened_Compounds.txt", as.is = T)
which(gdsc_anno$TARGET_PATHWAY == "ERK MAPK signaling") -> ii_1
which(gdsc_anno$TARGET_PATHWAY == "EGFR signaling") -> ii_2
mapk_inhibitor <- sort(unique(gdsc_anno[ii_1, 2]))
egfr_inhibitor <- sort(unique(gdsc_anno[ii_2, 2]))

mek_inhibitor <- c(mapk_inhibitor, egfr_inhibitor)

mek_inhibitor <- c(
    mapk_inhibitor, egfr_inhibitor,
    setdiff(two_drugs_match[, 3], c(mapk_inhibitor, egfr_inhibitor))
)
mek_inhibitor[grep("Nutlin", mek_inhibitor)] <- rownames(
    tcga_sensitive_mat
)[grep("Nutlin", rownames(tcga_sensitive_mat))]

tcga_cutoff <- 0.05

### way 2
match(mek_inhibitor, rownames(tcga_sensitive_mat)) -> ii
### no overlap
new_tcga_sensitive_mat <- tcga_sensitive_mat[ii, ]
new_tcga_resistant_mat <- tcga_resistant_mat[ii, ]

new_tcga_sensitive_mat <- new_tcga_sensitive_mat[
    nrow(new_tcga_sensitive_mat):1,
]
new_tcga_resistant_mat <- new_tcga_resistant_mat[
    nrow(new_tcga_resistant_mat):1,
]

png("F.GDSC.MEKonly.lineage.png", width = 600, height = 600)
par(mar = c(1, 7, 10, 3))

m <- -log10(new_tcga_sensitive_mat)
m[which(abs(new_tcga_sensitive_mat) > tcga_cutoff)] <- 0
m[which(is.na(m))] <- 0
m <- t(m)
m[which(m > 100)] <- 100

new_m <- (m - min(m, na.rm = T)) / (max(m, na.rm = T) - min(m, na.rm = T))

white2red <- colorRampPalette(c("white", "red"))
cc <- c(white2red(500))


image(new_m, xaxt = "n", yaxt = "n", col = cc)

### row labels: drug names
gene_count <- ncol(m)
path_col <- rep("black", ncol(m))
lab <- colnames(m)
lab[grep("Nutlin", lab)] <- "Nutlin-3a (-)"
mtext(
    text = lab, side = 2,
    line = 0.3, at = (0:(gene_count - 1)) / (gene_count - 1),
    las = 1, cex = .8, col = path_col
)

### col labels: cancer type
gene_count <- nrow(m)
text(
    (0:(gene_count - 1)) / (gene_count - 1),
    1.05, srt = 90, adj = 0, cex = .8,
    labels = rownames(m), xpd = TRUE
)

gene_count <- nrow(m)
1 / (2 * gene_count - 2) + (0:(gene_count - 1)) / (gene_count - 1) -> tt
for (x in tt) {
    abline(v = x, col = grey(0.9), lwd = .2)
}
gene_count <- ncol(m)
1 / (2 * gene_count - 2) + (0:(gene_count - 1)) / (gene_count - 1) -> tt
for (x in tt) {
    abline(h = x, col = grey(0.8), lwd = .2)
}
box()

m <- -log10(new_tcga_resistant_mat)
m[which(abs(new_tcga_resistant_mat) > tcga_cutoff)] <- 0
m[which(is.na(m))] <- 0
m <- t(m)
m[which(m > 100)] <- 100

new_m <- (m - min(m, na.rm = T)) / (max(m, na.rm = T) - min(m, na.rm = T))

white2blue <- colorRampPalette(c("white", "blue"))
cc <- c(white2blue(500))

gene_count <- nrow(new_m)
-1 / (2 * gene_count - 2) + (0:(gene_count)) / (gene_count - 1) -> tt_x
gene_count <- ncol(m)
-1 / (2 * gene_count - 2) + (0:(gene_count)) / (gene_count - 1) -> tt_y
for (k1 in 1:(length(tt_x) - 1)) {
    for (k2 in 1:(length(tt_y) - 1)) {
        rownames(m)[k1] -> gene
        colnames(m)[k2] -> pathway
        if (!(is.element(gene, rownames(new_m)) &
            is.element(pathway, colnames(new_m)))) next
        if (new_m[gene, pathway] != 1) {
            cc_idx <- round(new_m[gene, pathway] / 2 * 1e3)
            if (cc_idx < 1) cc_idx <- 1
            polygon(
                c(tt_x[k1], tt_x[k1 + 1], tt_x[k1 + 1]),
                c(tt_y[k2], tt_y[k2], tt_y[k2 + 1]),
                col = cc[cc_idx], border = NA
            )
        }
    }
}
box()

# gene_count <- ncol(m)
# rect(
#     1, 0, 1.01, length(ii_1) / (gene_count - 1),
#     col = "red", border = NA
# )
# rect(
#     0.9, 0, 0.91, length(ii_2) / (gene_count - 1),
#     col = "lightblue", border = NA
# )
# rect(
#     0.8, 0, 0.81, length(c(ii_2, ii_1)) / (gene_count - 1),
#     col = "lightgreen", border = NA
# )
dev.off()

#
rbind(
    cbind(
        Cancer = names(sensitive_prop),
        Prop = sensitive_prop, grp = "Sensitive"
    ),
    cbind(
        Cancer = names(resistant_prop),
        Prop = -resistant_prop, grp = "Insensitive"
    )
) -> dat
dat <- as.data.frame(dat)
dat[, 2] <- as.numeric(as.character(dat[, 2]))
the_order <- names(sensitive_prop)
the_order <- the_order[length(the_order):1]
dat[, 3] <- factor(dat[, 3], levels = c("Sensitive", "Insensitive"))

g <- ggplot(dat, aes(x = Cancer, y = Prop, group = grp, fill = grp)) +
    geom_bar(stat = "identity", width = 0.75) +
    coord_flip() +
    scale_x_discrete(limits = the_order) +
    scale_y_continuous(
        breaks = seq(-0.5, 1, 0.1),
        labels = abs(seq(-0.5, 1, 0.1))
    ) +
    xlab("") +
    ylab("") +
    ggtitle("GDSC model") +
    theme(
        legend.position = c(0.9, 0.9),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        plot.title = element_text(hjust = 0.5, size = 7),
        panel.background = element_rect(fill = "grey90"),
        axis.text.y = element_text(size = 6)
    )

ggsave(g, file = "./E.GDSC.PLX4720.png", width = 4.5, height = 3)
