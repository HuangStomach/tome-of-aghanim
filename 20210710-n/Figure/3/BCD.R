library(ggplot2)
load("./BCD.data.RData")

png("C.CCLE.lineage.png", width = 1400, height = 1200)
tcga_cutoff <- 0.05
par(mar = c(1, 6, 10, 4))

tcga_sensitive_mat <- tcga_sensitive_mat[nrow(tcga_sensitive_mat):1, ]
tcga_resistant_mat <- tcga_resistant_mat[nrow(tcga_resistant_mat):1, ]

m <- -log10(tcga_sensitive_mat)

which(is.infinite(m)) -> ii
if (length(ii) > 0) {
    m[which(is.infinite(m))] <- max(m[!is.infinite(m)])
}

m[which(abs(m) < -log10(tcga_cutoff))] <- 0
m[which(abs(m) > 100)] <- 100
m[which(is.na(m))] <- 0
m <- t(m)

new_m <- (m - min(m, na.rm = T)) / (max(m, na.rm = T) - min(m, na.rm = T))

white2red <- colorRampPalette(c("white", "red"))
cc <- c(white2red(500))

image(new_m, xaxt = "n", yaxt = "n", col = cc)
gene_count <- ncol(m)
path_col <- rep("black", ncol(m))
mtext(
    text = colnames(m), side = 2,
    line = 0.3, at = (0:(gene_count - 1)) / (gene_count - 1),
    las = 1, cex = .8, col = path_col
)
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

m <- -log10(tcga_resistant_mat)

which(is.infinite(m))

m[which(abs(m) < -log10(tcga_cutoff))] <- 0
m[which(abs(m) > 100)] <- 100
m[which(is.na(m))] <- 0
m <- t(m)

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
mtext("(C)", side = 2, at = 1.1, las = 1, line = 2)

dev.off()
rbind(
    cbind(
        Cancer = names(sensitive_prop),
        Prop = sensitive_prop,
        grp = "Sensitive"
    ),
    cbind(
        Cancer = names(resistant_prop),
        Prop = -resistant_prop,
        grp = "Insensitive"
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
    ggtitle("CCLE model") +
    theme(
        legend.position = c(0.9, 0.9), legend.title = element_blank(),
        legend.text = element_text(size = 7),
        plot.title = element_text(hjust = 0.5, size = 7),
        panel.background = element_rect(fill = "grey90"),
        axis.text.y = element_text(size = 6)
    )

ggsave(g, file = "./D.CCLE.PLX4720.1.png", width = 4.5, height = 3)

#####
png("B.CCLE.PLX4720.2.png", width = 800, height = 600)
h <- hist(plx4720_pred_dr, breaks = 100, plot = FALSE)

cuts <- cut(
    h$breaks,
    c(-Inf, quantile(plx4720_pred_dr, probs = .05),
    quantile(plx4720_pred_dr, probs = .95), Inf)
)
cc <- rep(c("skyblue", "grey99", "red"), table(cuts))
plot(
    h, col = cc,
    xlab = "PLX4720 predicted ActArea",
    ylab = "",
    main = "Distribution of predicted ActArea\nPLX4720 (CCLE model)",
    cex.main = .7, cex.axis = .7, cex.lab = .7
)

dev.off()
