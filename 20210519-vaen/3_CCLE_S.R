library("MASS")
library("magrittr")
library("glmnet")
library("modEvA")
library("vegan")

load("./Output/1/tcga_ss_mat.RData")
anno <- read.csv(
    "./DATA/CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv", as.is = T
)
drugs <- sort(unique(anno$Compound))
matrix(0, nrow = 100, ncol = length(drugs)) ->
solid_f1_r2_mat -> solid_avg_cv_r2_mat ->
solid_in_sample_r2_mat -> solid_sample_size

drugs -> colnames(solid_f1_r2_mat) -> colnames(solid_in_sample_r2_mat)
colnames(solid_in_sample_r2_mat) -> colnames(solid_sample_size)

solid_mat <- c()
for (ksigmoid in 1:100) {
    load(paste("./Output/2/", ksigmoid, ".CCLE.model.list.S.RData", sep = ""))
    for (k in seq_len(length(drugs))) {
        drug <- drugs[k]
        model_list[[drug]] -> res_list
        if (length(res_list) == 0) next
        fit <- res_list$model
        ys <- res_list$ys
        which(ys[, 1] != -9) -> ii
        ys <- ys[ii, ]
        if (sd(ys[, 2]) == 0) next
        solid_mat <- rbind(
            solid_mat,
            c(ksigmoid, res_list$model_summary, cor(ys[, 1], ys[, 2]))
        )

        #### way 4, PCC
        recall <- cor(ys[, 1], ys[, 2])
        precision <- as.numeric(res_list$model_summary[5])

        solid_sample_size[ksigmoid, k] <- nrow(ys)
        solid_in_sample_r2_mat[ksigmoid, k] <- recall
        solid_avg_cv_r2_mat[ksigmoid, k] <- precision
        solid_f1_r2_mat[ksigmoid, k] <- as.numeric(res_list$model_summary[7])
    }
    cat(ksigmoid, ".", sep = "")
}
save(solid_mat, solid_sample_size,
    solid_in_sample_r2_mat, solid_avg_cv_r2_mat,
    solid_f1_r2_mat,
    file = "./Output/3/CCLE.S.info.RData"
)

pdf("./Output/3/CCLE.S.ROC.pdf", width = 5, height = 5)
for (k in seq_len(length(drugs))) {
    drug <- drugs[k]
    plot(
        x = solid_in_sample_r2_mat[, k],
        y = solid_avg_cv_r2_mat[, k],
        main = drugs[k],
        xlab = "Self in_sample PCC",
        ylab = "avg PCC (in_sample)",
        col = rep("blue", 200), pch = 20, cex = .6
    )
    tmp <- cbind(
        idx = c(1:100),
        solid_f1_r2_mat[, drug],
        solid_in_sample_r2_mat[, drug],
        solid_avg_cv_r2_mat[, drug]
    )
    tmp <- tmp[order(tmp[, 4], decreasing = T), ]
    idx <- tmp[1:10, 1]
    points(solid_in_sample_r2_mat[idx, k], solid_avg_cv_r2_mat[idx, k],
        pch = 4, col = "red"
    )
}
dev.off()

tcga_pred_mat <- c()
solid_model_summary <- c()
for (k in seq_len(length(drugs))) {
    drug <- drugs[k]

    tmp <- cbind(
        idx = c(1:100),
        solid_f1_r2_mat[, drug],
        solid_in_sample_r2_mat[, drug],
        solid_avg_cv_r2_mat[, drug]
    )
    tmp <- tmp[order(tmp[, 4], decreasing = T), ] ### avg_cv_r2
    best_index <- tmp[1, 1]

    load(paste("./Output/2/", best_index, ".CCLE.model.list.S.RData", sep = ""))
    model_list[[drug]] -> res_list
    fit <- res_list$model

    tcga_pred <- read.table(
        paste("./Output/1/", best_index, ".TCGA_latent.tsv", sep = ""), 
        header = T, sep = "\t", as.is = T
    )
    tcga_test_data <- tcga_pred[, -1]
    tcga_probabilities <- predict(
        fit, as.matrix(tcga_test_data), s = "lambda.min"
    )

    tcga_pred_mat <- cbind(tcga_pred_mat, tcga_probabilities)
    cat("...", drug, ".", sep = "")
}

tcga_pred_mat <- cbind(tcga_pred[, 1], "A", tcga_pred_mat)
gsub("\\.", "-", tcga_pred_mat[, 1]) -> ss
tcga_pred_mat[, 1] <- ss
match(tcga_pred_mat[, 1], tcga_ss_mat[, 1]) -> ii
tcga_pred_mat[, 2] <- tcga_ss_mat[ii, 2]
colnames(tcga_pred_mat) <- c("Sample", "Cancer", drugs)
write.table(tcga_pred_mat,
    file = "./Output/3/VAEN_CCLE.S.pred_TCGA.txt",
    quote = F, sep = "\t", row.names = FALSE
)

for (kdrug in seq_len(length(drugs))) {
    drug <- drugs[kdrug]
    if (drug == "X17.AAG") drug <- "17-AAG"
    gsub("\\.", "-", drug) -> drug

    tmp <- cbind(
        idx = c(1:100),
        solid_f1_r2_mat[, drug],
        solid_in_sample_r2_mat[, drug],
        solid_avg_cv_r2_mat[, drug]
    )
    tmp <- tmp[order(tmp[, 4], decreasing = T), ]

    pred_mat <- c()
    best_index <- tmp[1, 1]

    load(paste("./Output/2/", best_index, ".CCLE.model.list.S.RData", sep = ""))
    model_list[[drug]] -> res_list
    ys <- res_list$ys

    if (kdrug == 1) {
        self_prediction_mat <- ys
    } else {
        self_prediction_mat <- cbind(self_prediction_mat, ys[, 2])
    }
}

self_prediction_mat[, 1] <- rownames(ys)
colnames(self_prediction_mat) <- c("CELLLINE", drugs)
write.table(self_prediction_mat,
    file = "./Output/3/VAEN_CCLE.S.pred_CCLE.txt",
    quote = F, sep = "\t", row.names = FALSE
)

ccle_pred_full_mat <- c()
for (k in seq_len(length(drugs))) {
    drug <- drugs[k]

    tmp <- cbind(
        idx = c(1:100),
        solid_f1_r2_mat[, drug],
        solid_in_sample_r2_mat[, drug],
        solid_avg_cv_r2_mat[, drug]
    )
    tmp <- tmp[order(tmp[, 4], decreasing = T), ]
    best_index <- tmp[1, 1]

    load(paste("./Output/2/", best_index, ".CCLE.model.list.S.RData", sep = ""))
    model_list[[drug]] -> res_list

    ccle_latent <- read.table(
        paste("./Output/1/", best_index, ".CCLE_latent.tsv", sep = ""),
        header = T, sep = "\t", as.is = T
    )
    ccle_latent_data <- ccle_latent[, -1]
    fit <- res_list$model
    ccle_probabilities <- predict(
        fit, as.matrix(ccle_latent_data), s = "lambda.min"
    )

    ccle_pred_full_mat <- cbind(ccle_pred_full_mat, ccle_probabilities)
    cat(drug, ".", sep = "")
}

self_pred_mat <- cbind(CELLINE = ccle_latent[, 1], ccle_pred_full_mat)
colnames(self_pred_mat) <- c("CELLINE", drugs)
write.table(self_pred_mat,
    file = paste("./Output/3/VAEN_CCLE.S.pred_CCLE.full.txt", sep = ""),
    quote = F, sep = "\t", row.names = FALSE
)