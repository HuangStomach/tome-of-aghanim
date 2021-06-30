library("MASS")
library("magrittr")
library("glmnet")
library("modEvA")
library("vegan")

load("./Output/1/tcga_ss_mat.RData")
anno <- read.delim("./Data/GDSC/v17.3_fitted_dose_response.txt", as.is = T)
drugs <- sort(unique(anno$DRUG_NAME))

# 初始化矩阵
matrix(0, nrow = 100, ncol = length(drugs)) -> all_f1_r2_mat
all_f1_r2_mat -> all_avg_cv_r2_mat -> all_in_sample_r2_mat -> all_sample_size

drugs -> colnames(all_f1_r2_mat) -> colnames(all_sample_size) ->
colnames(all_in_sample_r2_mat) -> colnames(all_avg_cv_r2_mat)

all_mat <- c()
for (ksigmoid in 1:100) {
    load(paste("./Output/2/", ksigmoid, ".GDSC.model.list.RData", sep = ""))
    for (k in seq_len(length(drugs))) {
        drug <- drugs[k]
        model_list[[drug]] -> res_list
        if (length(res_list) == 0) next
        fit <- res_list$model
        ys <- res_list$ys
        which(ys[, 1] != -9) -> ii # 筛选出合法数据（不合法的被设置为-9
        ys <- ys[ii, ]
        if (sd(ys[, 2]) == 0) next
        all_mat <- rbind(
            all_mat,
            c(ksigmoid, res_list$model_summary, cor(ys[, 1], ys[, 2]))
        ) # 将符合条件的模型追加到矩阵 ys: ['真实反映', '预测反应']

        #### way 4, PCC
        recall <- cor(ys[, 1], ys[, 2])
        precision <- as.numeric(res_list$model_summary[5]) # R2_avg

        all_sample_size[ksigmoid, k] <- nrow(ys)
        all_in_sample_r2_mat[ksigmoid, k] <- recall
        all_avg_cv_r2_mat[ksigmoid, k] <- precision
        all_f1_r2_mat[ksigmoid, k] <- as.numeric(res_list$model_summary[7])
        # cv_R2_avg
    }
    cat(ksigmoid, ".", sep = "")
}

save(all_mat, all_sample_size,
    all_in_sample_r2_mat, all_avg_cv_r2_mat,
    all_f1_r2_mat,
    file = "./Output/3/GDSC.A.info.RData"
)

pdf("./Output/3/GDSC.A.ROC.pdf", width = 5, height = 5)
for (k in seq_len(length(drugs))) {
    drug <- drugs[k]
    plot(
        x = all_in_sample_r2_mat[, k],
        y = all_avg_cv_r2_mat[, k],
        main = drug,
        xlab = "Self in_sample PCC",
        ylab = "avg PCC (in_sample)",
        col = rep("blue", 200), pch = 20, cex = .6
    )
    tmp <- cbind(
        idx = c(1:100),
        all_f1_r2_mat[, drug],
        all_in_sample_r2_mat[, drug],
        all_avg_cv_r2_mat[, drug]
    )
    tmp <- tmp[order(tmp[, 4], decreasing = T), ]
    idx <- tmp[1:10, 1]
    points(all_in_sample_r2_mat[idx, k], all_avg_cv_r2_mat[idx, k],
        pch = 4, col = "red"
    )
}
dev.off()

# 选择误差最小的模型对TCGA进行预测记录
tcga_pred_mat <- c()
all_model_summary <- c()
for (k in seq_len(length(drugs))) {
    drug <- drugs[k]

    tmp <- cbind(
        idx = c(1:100),
        all_f1_r2_mat[, drug], # cv_R2_avg
        all_in_sample_r2_mat[, drug], # 预测误差（协方差
        all_avg_cv_r2_mat[, drug] # R2_avg
    )
    tmp <- tmp[order(tmp[, 4], decreasing = T), ] # 排序 cv_R2_avg
    best_index <- tmp[1, 1] # 记录误差最小的索引

    load(paste("./Output/2/", best_index, ".GDSC.model_list.RData", sep = ""))
     # 根据表现最佳的索引获取对应药物的训练模型
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

    # 将tcga_latent的数据预测后的结果进行追加
    tcga_pred_mat <- cbind(tcga_pred_mat, tcga_probabilities)
    cat("...", drug, ".", sep = "")
}

tcga_pred_mat <- cbind(tcga_pred[, 1], "A", tcga_pred_mat) # 先用A填充
gsub("\\.", "-", tcga_pred_mat[, 1]) -> ss
tcga_pred_mat[, 1] <- ss
match(tcga_pred_mat[, 1], tcga_ss_mat[, 1]) -> ii
tcga_pred_mat[, 2] <- tcga_ss_mat[ii, 2] # 用TCGA癌症名代替之前的占位符A
colnames(tcga_pred_mat) <- c("Sample", "Cancer", drugs)
write.table(tcga_pred_mat,
    file = "./Output/3/VAEN_GDSC.A.pred_TCGA.txt",
    quote = F, sep = "\t", row.names = FALSE
) # 输出每种癌症和药物下的最好预测

for (kdrug in seq_len(length(drugs))) {
    drug <- drugs[kdrug]
    if (drug == "X17.AAG") drug <- "17-AAG"
    gsub("\\.", "-", drug) -> drug

    tmp <- cbind(
        idx = c(1:100),
        all_f1_r2_mat[, drug],
        all_in_sample_r2_mat[, drug],
        all_avg_cv_r2_mat[, drug]
    ) # 各种损失 排序
    tmp <- tmp[order(tmp[, 4], decreasing = T), ]

    pred_mat <- c()
    best_index <- tmp[1, 1]

    load(paste("./Output/2/", best_index, ".GDSC.model_list.RData", sep = ""))
    model_list[[drug]] -> res_list
    ys <- res_list$ys

    if (kdrug == 1) {
        self_prediction_mat <- ys
    } else {
        self_prediction_mat <- cbind(self_prediction_mat, ys[, 2])
    }
    # 同上把每种药物最好的模型自我预测追加至self_prediction_mat
}

self_prediction_mat[, 1] <- rownames(ys)
colnames(self_prediction_mat) <- c("CELLLINE", drugs)
write.table(self_prediction_mat,
    file = "./Output/3/VAEN_GDSC.A.pred_GDSC.txt",
    quote = F, sep = "\t", row.names = FALSE
)

ccle_pred_full_mat <- c()
for (k in seq_len(length(drugs))) {
    drug <- drugs[k]
    tmp <- cbind(
        idx = c(1:100),
        all_f1_r2_mat[, drug],
        all_in_sample_r2_mat[, drug],
        all_avg_cv_r2_mat[, drug]
    )
    tmp <- tmp[order(tmp[, 4], decreasing = T), ]

    pred_mat <- c()
    best_index <- tmp[1, 1]

    load(paste("./Output/2/", best_index, ".GDSC.model_list.RData", sep = ""))
    model_list[[drug]] -> res_list # 同样找到表现最好的模型

    ccle_latent <- read.table(
        paste("./Output/1/", best_index, ".CCLE_latent.tsv", sep = ""),
    header = T, sep = "\t", as.is = T)
    ccle_latent_data <- ccle_latent[, -1]
    fit <- res_list$model
    ccle_probabilities <- predict(fit,
        as.matrix(ccle_latent_data), s = "lambda.min"
    )

    # 将ccle_latent的数据预测后的结果进行追加
    ccle_pred_full_mat <- cbind(ccle_pred_full_mat, ccle_probabilities)
    cat(drug, ".", sep = "")
}

self_pred_mat <- cbind(CELLINE = ccle_latent[, 1], ccle_pred_full_mat)
colnames(self_pred_mat) <- c("CELLINE", drugs)
write.table(self_pred_mat,
    file = paste("./Output/3/VAEN_GDSC.A.pred_CCLE.full.txt", sep = ""),
    quote = F, sep = "\t", row.names = FALSE
)
