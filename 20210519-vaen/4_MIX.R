library("MASS")
library("magrittr")
library("glmnet")
library("modEvA")
library("vegan")

anno <- read.csv(
    "./Data/CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv", as.is = T
)
drugs <- sort(unique(anno$Compound))

load("./Output/3/CCLE.A.info.RData")
load("./Output/3/CCLE.S.info.RData")

solid_drugs <- c("Erlotinib", "AZD0530", "PLX4720", "TKI258", "ZD.6474")
pdf("./Output/4/CCLE.MIX-F1-W5-PCC.ROC.pdf", width = 8, height = 10.5)
par(mfrow = c(4, 3), mar = c(4, 4, 2, 1))
for (k in seq_len(length(drugs))) { # 这什么写法？ 确定就24种药物？
    drug <- drugs[k]
    l <- min(c(all_avg_cv_r2.mat[, k], solid_avg_cv_r2_mat[, k]), na.rm = T)
    h <- max(c(all_avg_cv_r2.mat[, k], solid_avg_cv_r2_mat[, k]), na.rm = T)
    lx <- min(c(
        all_in_sample_r2_mat[, k],
        solid_in_sample_r2_mat[, k]), na.rm = T
    )
    hx <- max(c(
        all_in_sample_r2_mat[, k],
        solid_in_sample_r2_mat[, k]), na.rm = T
    )

    ### model: All
    plot(
        x = all_in_sample_r2_mat[, k],
        y = all_avg_cv_r2.mat[, k],
        main = , xlab = "", ylab = "",
        ylim = c(L, H), xlim = c(lx, hx)
    )
    points(
        x = solid_in_sample_r2_mat[, k],
        y = solid_avg_cv_r2_mat[, k],
        pch = 3
    )
    ordiellipse(
        rbind(
            cbind(all_in_sample_r2_mat[, k], all_avg_cv_r2.mat[, k]),
            cbind(solid_in_sample_r2_mat[, k], solid_avg_cv_r2_mat[, k])
        ),
        groups = c(rep(1, 100), rep(2, 100)),
        col = c(1:2), display = "sites", kind = "sd",
        label = F, conf = 0.95, lty = 5, lwd = 0.5
    )

    which.max(all_avg_cv_r2.mat[, k]) -> idx
    points(
        all_in_sample_r2_mat[idx, k],
        all_avg_cv_r2.mat[idx, k],
        pch = 19, col = "cyan"
    )

    which.max(solid_avg_cv_r2_mat[, k]) -> idx
    points(
        solid_in_sample_r2_mat[idx, k],
        solid_avg_cv_r2_mat[idx, k],
        pch = 3, col = "cyan"
    )

    mtext(text = "In-sample PCC", side = 1, line = 2.3, cex = .9)
    mtext(text = "CV-R2", side = 2, line = 2.3, cex = .9)
    mtext(text = drugs[k], side = 3, line = .6)
}
dev.off()
# > solid_drugs [1] "AZD0530"   "Lapatinib" "LBW242"    "PLX4720"

all_tcga_pred_mat <- read.table(
    "./Output/3/VAEN_CCLE.A.pred_TCGA.txt", as.is = T, header = T
)
solid_tcga_pred_mat <- read.table(
    "./Output/3/VAEN_CCLE.S.pred_TCGA.txt", as.is = T, header = T
)

tcga_pred_mat <- all_tcga_pred_mat
immune_cancer <- c("LAML", "DLBC", "THYM")
which(!(tcga_pred_mat[, 2] %in% immune_cancer)) -> ii

for (k in seq_len(length(solid_drugs))) {
    drug <- solid_drugs[k]
    cat("Updated ", drug, "\n", sep = "")
    tcga_pred_mat[ii, drug] <- solid_tcga_pred_mat[ii, drug]
}
write.table(tcga_pred_mat,
    file = paste("./Output/4/VAEN_CCLE.MIX.pred_TCGA.txt", sep = ""),
    quote = F, sep = "\t", row.names = FALSE
) # 将A和S(我猜测是不同组织类型的癌症)的数据对tcga的预测进行覆盖混合

all_ccle_pred_mat <- read.table(
    "./Output/3/VAEN_CCLE.A.pred_CCLE.txt", as.is = T, header = T
)
solid_ccle_pred_mat <- read.table(
    "./Output/3/VAEN_CCLE.S.pred_CCLE.txt", as.is = T, header = T
)

pdf("./Output/4/MIX-F1-W5-PCC.obsd.vs.pred.pdf", width = 8, height = 4)
par(mfrow = c(1, 2), cex = 1, mar = c(4, 4, 3, 1))
for (k in seq_len(length(drugs))) {
    drug <- drugs[k]

    load(paste("./Output/2/1.CCLE.model_list.RData", sep = ""))
    model_list[[drug]] -> res_list
    ys <- res_list$ys
    which(ys[, 1] != -9) -> ii
    r <- cor(ys[ii, 1], all_ccle_pred_mat[ii, k + 1])
    plot(
        ys[ii, 1],
        all_ccle_pred_mat[ii, k + 1],
        main = paste(drug, ", all, PCC = ", format(r, digits = 3), sep = ""),
        xlab = "Observed CCLE", ylab = "Predicted CCLE (all)"
    )

    load(paste("./Output/2/1.CCLE.model_list.S.RData", sep = ""))
    model_list[[drug]] -> res_list
    ys <- res_list$ys
    which(ys[, 1] != -9) -> ii
    r <- cor(ys[ii, 1], solid_ccle_pred_mat[ii, k + 1])
    plot(
        ys[ii, 1],
        solid_ccle_pred_mat[ii, k + 1],
        main = paste(drug, ", solid, PCC = ", format(r, digits = 3), sep = ""),
        xlab = "Observed CCLE", ylab = "Predicted CCLE (solid)"
    )
}
dev.off() # 将A和S的训练模型和预测进行图示

all_ccle_pred_mat <- read.table(
    "./Output/3/VAEN_CCLE.A.pred_CCLE.txt", as.is = T, header = T
)
solid_ccle_pred_mat <- read.table(
    "./Output/3/VAEN_CCLE.S.pred_CCLE.txt", as.is = T, header = T
)

ccle_pred_mat <- all_ccle_pred_mat
for (k in seq_len(length(solid_drugs))) {
    drug <- solid_drugs[k]
    ccle_pred_mat[, drug] <- solid_ccle_pred_mat[, drug]
}
write.table(ccle_pred_mat,
    file = paste("./Output/4/VAEN_CCLE.MIX.pred_CCLE.txt", sep = ""),
    quote = F, sep = "\t", row.names = FALSE
)
# 将A和S(我猜测是不同组织类型的癌症)的数据对ccle的预测进行覆盖混合

all_ccle_pred_mat <- read.table(
    "./Output/3/VAEN_CCLE.A.pred_CCLE.full.txt", as.is = T, header = T
)
solid_ccle_pred_mat <- read.table(
    "./Output/3/VAEN_CCLE.S.pred_CCLE.full.txt", as.is = T, header = T
)

ccle_pred_mat <- all_ccle_pred_mat
for (k in seq_len(length(solid_drugs))) {
    drug <- solid_drugs[k]
    ccle_pred_mat[, drug] <- solid_ccle_pred_mat[, drug]
}
write.table(ccle_pred_mat,
    file = paste("./Output/4/VAEN_CCLE.MIX.pred_CCLE.full.txt", sep = ""),
    quote = F, sep = "\t", row.names = FALSE
)
# 将A和S(我猜测是不同组织类型的癌症)的数据对ccle_latent的预测进行覆盖混合
