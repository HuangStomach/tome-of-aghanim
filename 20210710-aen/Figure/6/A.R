library("gplots")

load("../../Output/2/1.CCLE.model.list.RData")
ccle_model_list <- model_list

obsd_ccle_mat <- c()
for (k in seq_len(length(ccle_model_list))) {
    res_list <- ccle_model_list[[k]]
    # drug <- names(ccle_model_list)[k]
    ys <- res_list$ys # 现在才想明白这个ys可能是y的复数形式……

    obsd_ccle_mat <- cbind(obsd_ccle_mat, ys[, 1])
}
colnames(obsd_ccle_mat) <- names(ccle_model_list)

apply(obsd_ccle_mat, 1, function(u) sum(u == -9)) -> check
obsd_ccle_mat <- obsd_ccle_mat[check < 1, ] # 筛选数据

### Plot 1
x <- as.matrix(obsd_ccle_mat)
apply(x, 2, scale) -> x1 # 数据规范化
heatmap.2(x1, trace = "none", col = bluered(75)) -> h2_2

### predicted CCLE
self_pred <- read.table(
    paste("../../Output/3/VAEN_CCLE.A.pred_CCLE.txt", sep = ""),
    as.is = T, header = T, sep = "\t"
)
dim(self_pred)
colnames(self_pred) -> dd
dd[dd == "X17.AAG"] <- "17-AAG"
colnames(self_pred) <- gsub("\\.", "-", dd) # 读取ccLe的自我预测 格式化

apply(self_pred[, 2:ncol(self_pred)], 1, function(u) sum(u == -9)) -> check
self_pred <- self_pred[check < 1, ]

x <- as.matrix(self_pred[, 2:25])
apply(x, 2, scale) -> x1
heatmap.2(x1, trace = "none", col = bluered(75)) -> pred_ccle_hm_2

### imputed CCLE
self_pred <- read.table(
    paste("../../Output/3/VAEN_CCLE.A.pred_CCLE.full.txt", sep = ""),
    as.is = T, header = T, sep = "\t"
) # 这不是自我预测不知道为什么叫这个
dim(self_pred)
colnames(self_pred) -> dd
dd[dd == "X17.AAG"] <- "17-AAG"
colnames(self_pred) <- gsub("\\.", "-", dd)

apply(self_pred[, 2:ncol(self_pred)], 1, function(u) sum(u == -9)) -> check
self_pred <- self_pred[check < 1, ]

x <- as.matrix(self_pred[, 2:25])
apply(x, 2, scale) -> x1
heatmap.2(x1, trace = "none", col = bluered(75)) -> imputed_ccle_hm_2


### predicted TCGA
drug_ccle <- read.table(
    file = "../../Output/3/VAEN_CCLE.A.pred_TCGA.txt",
    header = T, as.is = T, sep = "\t"
)
colnames(drug_ccle) -> dd
dd[dd == "X17.AAG"] <- "17-AAG"
colnames(drug_ccle) <- gsub("\\.", "-", dd)

cancer_types <- unique(drug_ccle[, 2])
sample_type <- substr(drug_ccle[, 1], 14, 15)

cancer_drug_ccle <- c()
for (ct in seq_len(length(cancer_types))) {
    cancer <- cancer_types[ct]

    type_code <- "01"
    if (cancer == "LAML") {
        type_code <- "03"
    }
    if (cancer == "SKCM") {
        type_code <- "06"
    }

    blca_ccle <- drug_ccle[
        which(drug_ccle[, 2] == cancer & sample_type == type_code),
    ]
    cancer_drug_ccle <- rbind(cancer_drug_ccle, blca_ccle)
}
ccle <- cancer_drug_ccle


cc <- read.table(
    "../../Data/TCGA.color.txt",
    as.is = T, sep = "\t", comment.char = ""
)
rowbar <- rep("", nrow(ccle))
for (ct in seq_len(length(cancer_types))) {
    cancer <- cancer_types[ct]
    rowbar[which(ccle[, 2] == cancer)] <- cc[which(cc[, 1] == cancer), 2]
}

x <- as.matrix(ccle[, 3:26])
apply(x, 2, scale) -> x1
heatmap.2(x1,
    trace = "none", col = bluered(75),
    RowSideColors = rowbar
) -> ccle_h_2

cc <- rep("black", 24)
cc[which(colnames(obsd_ccle_mat) %in%
    c("Lapatinib", "ZD-6474", "Erlotinib", "AZD0530"))] <- "orange"
cc[which(colnames(obsd_ccle_mat) %in%
    c("PD.0325901", "PD-0325901", "AZD6244"))] <- "red"
cc[which(colnames(obsd_ccle_mat) %in%
    c("RAF265", "PLX4720"))] <- "pink"
cc[which(colnames(obsd_ccle_mat) %in%
    c("PF2341066", "PHA-665752"))] <- "cyan"
cc[which(colnames(obsd_ccle_mat) %in%
    c("PD-0332991", "Nutlin-3"))] <- "green"
cc[which(colnames(obsd_ccle_mat) %in%
    c("Irinotecan", "Topotecan", "Paclitaxel"))] <- "blue"
drug_cc <- cc

library(ape)
png("A.png", width = 800, height = 1200)
par(mfrow = c(4, 1), mar = c(1, 5, 2, 5))
plot(
    as.phylo(as.hclust(h2_2$colDendrogram)),
    tip.color = drug_cc,
    direction = "d",
    font = 1, label.offset = 0,
    cex = 2, srt = -180, adj = 1
)
mtext("Cell lines, observed", cex = 2)

plot(
    as.phylo(as.hclust(pred_ccle_hm_2$colDendrogram)),
    tip.color = drug_cc,
    direction = "d",
    font = 1, label.offset = 0.5,
    cex = 2, srt = -180, adj = 1
)
mtext("Cell lines, predicted", cex = 2)

plot(
    as.phylo(as.hclust(imputed_ccle_hm_2$colDendrogram)),
    tip.color = drug_cc,
    direction = "d",
    font = 1, label.offset = 0.5,
    cex = 2, srt = -180, adj = 1
)
mtext("Cell lines, imputed", cex = 2)

plot(
    as.phylo(as.hclust(ccle_h_2$colDendrogram)),
    tip.color = drug_cc,
    direction = "d",
    font = 1, label.offset = 0.5,
    cex = 2, srt = -180, adj = 1
)
mtext("TCGA, predicted", cex = 2)

segments(1, 1, 3.5, 1, lwd = 2)
text(2.5, 0.5, "EGFRi")

segments(5, 1, 6.5, 1, lwd = 2)
text(6.5, 0.5, "MEKi")

segments(15, 1, 16.5, 1, lwd = 2)
text(16.5, 0.5, "BRAFi")

segments(9, 1, 11.8, 1, lwd = 2)
text(10, 0.5, "Cytotoxic")

segments(12, 1, 14, 1, lwd = 2)
text(13, 0.5, "c-METi")

dev.off()