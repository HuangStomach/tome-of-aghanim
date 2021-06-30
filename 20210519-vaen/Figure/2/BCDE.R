source("../../Lib/multiplot.R")

library(MASS)
library(glmnet)
library(ggplot2)
library(magrittr)

one_drugs_match <- read.table(
    "../../Data/Match/drugs_match.txt", as.is = T
)
two_drugs_match <- read.table(
    "../../Data/Match/drugs_match_2.txt", as.is = T, sep = "\t"
)

ccle_anno <- read.csv(
    "../../Data/CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv", as.is = T
)
gdsc_anno <- read.delim(
    "../../Data/GDSC/v17.3_fitted_dose_response.txt", as.is = T
)
cell_line_anno <- read.csv(
    "../../Data/CCLE/DepMap-2018q3-celllines.csv", as.is = T
)

pps <- read.table(paste("../../Output/1/1.CCLE_latent.tsv", sep = ""))
original_ss_pp <- rownames(pps)
sapply(original_ss_pp, function(x) {
    new_u <- u <- strsplit(x, split = "\\.")[[1]][1]
    if (grepl("^X", u)) {
        substr(u, 2, nchar(u)) -> new_u
    }
    new_u
}) -> ss_pp
names(ss_pp) <- NULL # 先将CCLE中的名称格式化

original_ss_pp <- rownames(pps)
sapply(original_ss_pp, function(x) {
    new_u <- u <- strsplit(x, split = "\\.")[[1]]
    paste(u[3], u[4], sep = "-") -> new_u
    new_u
}) -> ss_ach
names(ss_ach) <- NULL


ccle_pred <- read.table(
    "../../Output/3/VAEN_CCLE.A.pred_CCLE.txt", header = T, as.is = T
)
gdsc_pred <- read.table(
    "../../Output/3/VAEN_GDSC.A.pred_GDSC.txt", header = T, as.is = T, sep = "\t"
)


data4plot_list <- list()
shared_drugs_ori_cor <- c()
for (k in seq_len(nrow(two_drugs_match))) {
    ccle_anno_1 <- ccle_anno[
        which(ccle_anno$Compound == two_drugs_match[k, 1]),
    ]
    gdsc_anno_1 <- gdsc_anno[
        which(gdsc_anno$DRUG_NAME == two_drugs_match[k, 3]),
    ]
    print(c(two_drugs_match[k, 1], nrow(ccle_anno_1), nrow(gdsc_anno_1)))

    match(gdsc_anno_1$COSMIC_ID, cell_line_anno$COSMIC_ID) -> idx

    gdsc_anno_1.match1 <- gdsc_anno_1[which(gdsc_anno_1$COSMIC_ID %in% cell_line_anno$COSMIC_ID), ]
    gdsc_anno_1.match2 <- gdsc_anno_1[which(!gdsc_anno_1$COSMIC_ID %in% cell_line_anno$COSMIC_ID), ]

    match(gdsc_anno_1.match1$COSMIC_ID, cell_line_anno$COSMIC_ID) -> idx1
    match(cell_line_anno[idx1, 1], ss_ach) -> part1.ach.idx
    gdsc_anno_1.match1 <- gdsc_anno_1.match1[which(!is.na(part1.ach.idx)), ]
    match(gdsc_anno_1.match1$COSMIC_ID, cell_line_anno$COSMIC_ID) -> idx1
    match(cell_line_anno[idx1, 1], ss_ach) -> part1.ach.idx
    gdsc.Y <- -gdsc_anno_1.match1[, "LN_IC50"]
    names(gdsc.Y) <- gdsc_anno_1.match1[, "CELL_LINE_NAME"]


    sapply(gdsc_anno_1.match2$CELL_LINE_NAME, function(u) {
        gsub("-", "", u)
    }) -> gdsc.cell.name
    union(which(gdsc.cell.name %in% cell_line_anno$Aliases), which(gdsc_anno_1.match2$CELL_LINE_NAME %in% cell_line_anno$Aliases)) -> ii
    gdsc_anno_1.match2 <- gdsc_anno_1.match2[ii, ]
    gdsc.cell.name <- gdsc.cell.name[ii]

    idx2 <- c()
    for (kk in 1:nrow(gdsc_anno_1.match2)) {
        if (gdsc_anno_1.match2$CELL_LINE_NAME[kk] %in% cell_line_anno$Aliases) {
            idx2 <- c(idx2, match(gdsc_anno_1.match2$CELL_LINE_NAME[kk], cell_line_anno$Aliases))
        } else {
            idx2 <- c(idx2, match(gdsc.cell.name[kk], cell_line_anno$Aliases))
        }
    }
    match(cell_line_anno[idx2, 1], ss_ach) -> part2.ach.idx
    gdsc_anno_1.match2 <- gdsc_anno_1.match2[which(!is.na(part2.ach.idx)), ]

    idx2 <- c()
    for (kk in 1:nrow(gdsc_anno_1.match2)) {
        if (gdsc_anno_1.match2$CELL_LINE_NAME[kk] %in% cell_line_anno$Aliases) {
            idx2 <- c(idx2, match(gdsc_anno_1.match2$CELL_LINE_NAME[kk], cell_line_anno$Aliases))
        } else {
            idx2 <- c(idx2, match(gsub("-", "", gdsc_anno_1.match2$CELL_LINE_NAME[kk]), cell_line_anno$Aliases))
        }
    }


    match(cell_line_anno[idx2, 1], ss_ach) -> part2.ach.idx
    new.Y <- -gdsc_anno_1.match2[, "LN_IC50"]
    names(new.Y) <- gdsc_anno_1.match2[, "CELL_LINE_NAME"]
    gdsc.Y <- c(gdsc.Y, new.Y)
    gdsc.train.data <- pps[c(part1.ach.idx, part2.ach.idx), ]


    intersect(ccle_anno_1[, 1], ss_pp) -> shared.samples
    match(shared.samples, ccle_anno_1[, 1]) -> ii
    ccle_anno.2 <- ccle_anno_1[ii, ]
    ccle.Y <- ccle_anno.2[, "ActArea"]
    ccle.Y.IC50 <- ccle_anno.2[, "IC50..uM."]
    names(ccle.Y) <- ccle_anno.2[, "Primary.Cell.Line.Name"]
    match(shared.samples, ss_pp) -> ii
    ccle.train.data <- pps[ii, ]

    shared.samples <- intersect(rownames(gdsc.train.data), rownames(ccle.train.data))
    match(shared.samples, rownames(ccle.train.data)) -> ccle.ii
    match(shared.samples, rownames(gdsc.train.data)) -> gdsc.ii

    shared.ccle.Y <- ccle.Y[ccle.ii]
    shared.ccle.Y.IC50 <- ccle.Y.IC50[ccle.ii]
    shared.gdsc.Y <- gdsc.Y[gdsc.ii]

    ccle_pred.Y <- ccle_pred[match(shared.samples, ccle_pred[, 1]), two_drugs_match[k, 2]]
    gdsc_pred.Y <- gdsc_pred[match(shared.samples, gdsc_pred[, 1]), one_drugs_match[k, 3]]

    cur.list <- list()
    cur.list[["shared.ccle.Y"]] <- shared.ccle.Y
    cur.list[["shared.gdsc.Y"]] <- shared.gdsc.Y
    cur.list[["ccle_pred.Y"]] <- ccle_pred.Y
    cur.list[["gdsc_pred.Y"]] <- gdsc_pred.Y
    data4plot_list[[two_drugs_match[k, 1]]] <- cur.list
    shared_drugs_ori_cor <- rbind(shared_drugs_ori_cor, c(two_drugs_match[k, ], cor(shared.ccle.Y, shared.gdsc.Y), cor(shared.ccle.Y.IC50, shared.gdsc.Y), length(shared.ccle.Y), cor(ccle_pred.Y, gdsc_pred.Y)))
}
save(shared_drugs_ori_cor, data4plot_list, file = "A.best.shared.drugs.ori.RData")

######
ccle <- read.table("../../Output/3/VAEN_CCLE.A.pred_TCGA.txt", header = T, as.is = T)
gdsc <- read.table("../../Output/3/VAEN_GDSC.A.pred_TCGA.txt", header = T, as.is = T, sep = "\t")

cancer.types <- sort(unique(ccle[, 2]))
shared.drugs.pred.cor <- c()
for (k1 in 1:nrow(two_drugs_match)) {
    df <- as.data.frame(cbind(x = ccle[, two_drugs_match[k1, 2]], y = gdsc[, one_drugs_match[k1, 3]]))
    df[, 1] <- as.numeric(as.character(df[, 1]))
    df[, 2] <- as.numeric(as.character(df[, 2]))
    shared.drugs.pred.cor <- rbind(shared.drugs.pred.cor, c(two_drugs_match[k1, ], cor(df[, 1], df[, 2]), nrow(df)))
}

pdf("./BC.shared.drugs.pdf", width = 8, height = 4)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

plot(x = shared_drugs_ori_cor[, 4], y = shared_drugs_ori_cor[, 7], cex = unlist(shared_drugs_ori_cor[, 6]) / 100, xlab = "Cell line original (PCC)", ylab = "Cell line prediction (PCC)", xlim = c(0.1, 0.9))
text(x = shared_drugs_ori_cor[, 4], y = shared_drugs_ori_cor[, 7], labels = shared_drugs_ori_cor[, 1], pos = 3)

plot(x = shared_drugs_ori_cor[, 4], y = shared.drugs.pred.cor[, 4], cex = unlist(shared_drugs_ori_cor[, 6]) / 100, xlab = "Cell line original (PCC)", ylab = "TCGA prediction (PCC)", xlim = c(0.1, 0.9))
text(x = shared_drugs_ori_cor[, 4], y = shared.drugs.pred.cor[, 4], labels = shared_drugs_ori_cor[, 1], pos = 4)

dev.off()

write.table(shared.drugs.pred.cor, file = "./B.txt", row.names = F, quote = F, sep = "\t")
write.table(shared_drugs_ori_cor, file = "./C.txt", row.names = F, quote = F, sep = "\t")

#####

pdf("./D.pdf", width = 10, height = 10)
for (k in 1:nrow(two_drugs_match)) {
    ############ 1. observed
    data4plot_list[[two_drugs_match[k, 1]]] -> cur.list
    dat.df <- as.data.frame(cbind(x = cur.list$shared.ccle.Y, y = cur.list$shared.gdsc.Y))
    p1 <- ggplot(data = dat.df, aes(x, y)) +
        stat_density2d(aes(fill = ..level.., alpha = ..level..), geom = "polygon", colour = "black") +
        scale_fill_continuous(low = "green", high = "red") +
        geom_smooth(method = lm, linetype = 2, colour = "red", se = F) +
        guides(alpha = "none") +
        ggtitle(paste(one_drugs_match[k, 1], "\nObserved DR, r = ", format(cor(dat.df[, 1], dat.df[, 2]), digits = 3), ", n = ", nrow(dat.df), sep = "")) +
        geom_point() +
        labs(color = "Density", fill = "Density", x = "Observed ActArea, CCLE", y = "Observed -LN(IC50), GDSC") +
        theme(legend.position = c(0, 1), legend.justification = c(0, 1), plot.title = element_text(hjust = 0.5))

    #################### 2. CCLE predicted
    load(paste("../../Output/2/1.CCLE.model.list.RData", sep = ""))
    ccle.res.list <- model.list[[two_drugs_match[k, 1]]]
    Y <- ccle.res.list$Ys[, 1]
    pred_Y <- ccle_pred[, one_drugs_match[k, 2]]
    which(Y != -9) -> ii
    dat.df <- as.data.frame(cbind(x = Y[ii], y = pred_Y[ii]))
    p2 <- ggplot(data = dat.df, aes(x, y)) +
        stat_density2d(aes(fill = ..level.., alpha = ..level..), geom = "polygon", colour = "black") +
        scale_fill_continuous(low = "green", high = "red") +
        geom_smooth(method = lm, linetype = 2, colour = "red", se = F) +
        guides(alpha = "none") +
        ggtitle(paste(two_drugs_match[k, 1], "\nCCLE, r = ", format(cor(dat.df[, 1], dat.df[, 2]), digits = 3), ", n = ", nrow(dat.df), sep = "")) +
        geom_point() +
        labs(color = "Density", fill = "Density", x = "Observed ActArea, CCLE", y = "Predicted ActArea, CCLE") +
        theme(legend.position = c(0, 1), legend.justification = c(0, 1), plot.title = element_text(hjust = 0.5))

    #################### 3. GDSC predicted
    load(paste("../../Output/2/1.GDSC.model.list.RData", sep = ""))
    gdsc.res.list <- model.list[[two_drugs_match[k, 3]]]
    Y <- gdsc.res.list$Ys[, 1]
    pred_Y <- gdsc_pred[, one_drugs_match[k, 3]]
    which(Y != -9) -> ii
    dat.df <- as.data.frame(cbind(x = Y[ii], y = pred_Y[ii]))
    p3 <- ggplot(data = dat.df, aes(x, y)) +
        stat_density2d(aes(fill = ..level.., alpha = ..level..), geom = "polygon", colour = "black") +
        scale_fill_continuous(low = "green", high = "red") +
        geom_smooth(method = lm, linetype = 2, colour = "red", se = F) +
        guides(alpha = "none") +
        ggtitle(paste(two_drugs_match[k, 3], "\nGDSC, r = ", format(cor(dat.df[, 1], dat.df[, 2]), digits = 3), ", n = ", nrow(dat.df), sep = "")) +
        geom_point() +
        labs(color = "Density", fill = "Density", x = "Observed -LN(IC50), GDSC", y = "Predicted -LN(IC50), GDSC") +
        theme(legend.position = c(0, 1), legend.justification = c(0, 1), plot.title = element_text(hjust = 0.5))

    #################### 4. predicted
    dat.df <- as.data.frame(cbind(x = cur.list$ccle_pred.Y, y = cur.list$gdsc_pred.Y))
    p4 <- ggplot(data = dat.df, aes(x, y)) +
        stat_density2d(aes(fill = ..level.., alpha = ..level..), geom = "polygon", colour = "black") +
        scale_fill_continuous(low = "green", high = "red") +
        geom_smooth(method = lm, linetype = 2, colour = "red", se = F) +
        guides(alpha = "none") +
        ggtitle(paste(two_drugs_match[k, 3], "\nPredicted DR, r = ", format(cor(dat.df[, 1], dat.df[, 2]), digits = 3), ", n = ", nrow(dat.df), sep = "")) +
        geom_point() +
        labs(color = "Density", fill = "Density", x = "Predicted ActArea, CCLE", y = "Predicted -LN(IC50), GDSC") +
        theme(legend.position = c(0, 1), legend.justification = c(0, 1), plot.title = element_text(hjust = 0.5))

    multiplot(plotlist = list(p1, p2, p3, p4), layout = matrix(c(1:4), nrow = 2))
    cat(two_drugs_match[k, 1], ",", sep = "")
}
dev.off()


###########################
library("Hmisc")

ccle <- read.table("../../Output/3/VAEN_CCLE.A.pred_TCGA.txt", header = T, as.is = T)
gdsc <- read.table("../../Output/3/VAEN_GDSC.A.pred_TCGA.txt", header = T, as.is = T, sep = "\t")

cancer.types <- sort(unique(ccle[, 2]))

shared.drugs.mat <- matrix(0, nrow = nrow(one_drugs_match), ncol = length(cancer.types))
pdf("./E.percancer.pdf", width = 21, height = 15)
for (k1 in 1:nrow(one_drugs_match)) {
    pp.list <- list()
    for (k in 1:length(cancer.types)) {
        which(ccle[, 2] == cancer.types[k]) -> ii

        df <- as.data.frame(cbind(x = ccle[ii, one_drugs_match[k1, 2]], y = gdsc[ii, one_drugs_match[k1, 3]]))
        df[, 1] <- as.numeric(as.character(df[, 1]))
        df[, 2] <- as.numeric(as.character(df[, 2]))

        rcorr(as.matrix(df), type = "pearson") -> a
        shared.drugs.mat[k1, k] <- ifelse(a$r[1, 2] > 0, a$P[1, 2], 1)

        p1 <- ggplot(df, aes(x = x, y = y)) +
            geom_point() +
            ggtitle(cancer.types[k]) +
            xlab("CCLE") +
            ylab("GDSC") +
            theme(legend.position = "none", text = element_text(size = 10), axis.text.x = element_text(color = "black", size = 10), axis.text.y = element_text(color = "black", size = 10), plot.title = element_text(size = 10, hjust = 0.5)) +
            geom_smooth(method = "lm", color = "red", data = df, aes(x = x, y = y))

        pp.list[[k]] <- p1
    }

    multiplot(plotlist = pp.list, layout = matrix(c(1:35), ncol = 7, byrow = T))
    cat(k1, ".", sep = "")
}
dev.off()


shared.drugs.pcc.mat <- matrix(0, nrow = nrow(one_drugs_match), ncol = length(cancer.types))
for (k1 in 1:nrow(one_drugs_match)) {
    for (k in 1:length(cancer.types)) {
        which(ccle[, 2] == cancer.types[k]) -> ii

        df <- as.data.frame(cbind(
            x = ccle[ii, one_drugs_match[k1, 2]],
            y = gdsc[ii, one_drugs_match[k1, 3]]
        ))
        df[, 1] <- as.numeric(as.character(df[, 1]))
        df[, 2] <- as.numeric(as.character(df[, 2]))

        rcorr(as.matrix(df), type = "pearson") -> a
        shared.drugs.pcc.mat[k1, k] <- a$r[1, 2]
    }
    cat(k1, ".", sep = "")
}

###############################################################################
### 2E

library(reshape2)
rownames(shared.drugs.mat) <- one_drugs_match[, 1]
colnames(shared.drugs.mat) <- cancer.types
log.shared.drugs.mat <- t(-log(shared.drugs.mat + 1e-16))
new <- log.shared.drugs.mat[, order(apply(log.shared.drugs.mat, 2, mean))]

write.table(shared.drugs.mat, file = "./E.txt", sep = "\t", quote = F)
log.shared.drugs.mat <- t(-log(shared.drugs.mat + 1e-16))
new <- log.shared.drugs.mat[, order(apply(log.shared.drugs.mat, 2, mean))]

melt(new) -> dat
pdf("./E.pdf", width = 4, height = 4)
p1 <- ggplot(dat, aes(x = Var2, y = value)) +
    geom_boxplot(outlier.shape = NA) +
    guides(fill = FALSE) +
    geom_jitter(shape = 21, position = position_jitter(width = 0.4), size = 1) +
    guides(colour = FALSE) +
    xlab("") +
    ylab("-log(p)") +
    theme(
        text = element_text(size = 8),
        axis.text.x = element_text(
            color = "black",
            size = 8,
            angle = 90,
            hjust = 1
        ),
        axis.text.y = element_text(
            color = "black",
            size = 8
        ),
        plot.margin = unit(c(2, 2, 2, 2), "mm"),
        plot.title = element_text(size = 8, hjust = 0.5)
    )
print(p1)
dev.off()