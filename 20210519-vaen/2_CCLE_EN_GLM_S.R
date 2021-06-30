#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
    stop("Must provide start and end\n", call. = FALSE)
} else if (length(args) == 2) {
    start <- args[1]
    end <- args[2]
    scale.factor <- 1
} else if (length(args) == 3) {
    start <- args[1]
    end <- args[2]
    scale.factor <- as.numeric(args[3])
}
print(c(start, end))

library("MASS")
library("magrittr")
library("glmnet")
library("modEvA")
library(parallel)

source("./Lib/nested_EN.R")
parallel_main <- function(kk, train_data, y,
                          n_folds = 10, n_train_test_folds = 5,
                          seed = NA, alpha = 0.5, null_testing = FALSE, drug) {
    main(train_data, y,
        n_folds = n_folds,
        n_train_test_folds = n_train_test_folds,
        seed = seed, alpha = alpha,
        null_testing = null_testing, drug = drug
    ) -> res_list
    res_list
}

#####
load("./Output/1/tcga_ss_mat.RData")

anno <- read.csv(
    "./Data/CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv",
as.is = T)
drugs <- sort(unique(anno$Compound)) # 提取其中所有的抗癌药物

for (ksigmoid in start:end) {
    cat("ksigmoid = ", ksigmoid, "\n", sep = "")

    model_summary_file <- paste(
        "./Output/2/", ksigmoid,
        ".model_summary.txt", sep = ""
    )
    model_summary_cols <- c(
        "Drug", "alpha", "n_snps_in_model", "lambda_min_mse",
        "test_R2_avg", "test_R2_sd", "cv_r2_avg", "cv_R2_sd", "in_sample_R2",
        "nested_cv_fisher_pval", "rho_avg", "rho_se", "rho_zscore",
        "rho_avg_squared", "zscore_pval", "cv_rho_avg", "cv_rho_se",
        "cv_rho_avg_squared", "cv_zscore_est", "cv_zscore_pval", "cv_pval_est"
    ) # 声明一些统计信息

    tcga_pred <- read.table(
        paste("./Output/1/", ksigmoid, ".TCGA_latent.tsv", sep = ""),
    header = T, sep = "\t", as.is = T) # 读取经vae编码后的tcga数据
    tcga_test_data <- tcga_pred[, -1] # 去除第一列（样本名称）

    pps <- read.table(
        paste("./Output/1/", ksigmoid, ".CCLE_latent.tsv", sep = "")
    ) # 读取vae编码后的ccle数据
    original_ss_pp <- rownames(pps) # 获取其行名称（癌症类型）
    sapply(original_ss_pp, function(x) {
        new_u <- u <- strsplit(x, split = "\\.")[[1]][1]
        if (grepl("^X", u)) {
            substr(u, 2, nchar(u)) -> new_u
        }
        new_u
    }) -> ss_pp
    names(ss_pp) <- NULL # 对癌症类型简单处理留下主要信息

    self_prediction_mat <- matrix(-9,
        nrow = length(unique(anno[, 1])),
        ncol = length(drugs) + 2
    )
    self_prediction_mat[, 1] <- unique(anno[, "Primary.Cell.Line.Name"])
    self_prediction_mat[, 2] <- unique(anno[, "CCLE.Cell.Line.Name"])
    colnames(self_prediction_mat) <- c("CELLLINE", "Type", drugs)
    # 初始化预测矩阵 细胞系 类别 药物1 药物2...

    tcga_drug_response_mat <- c()
    model_list <- list()
    for (k in seq_len(length(drugs))) {
        drug <- drugs[k]
        cat(drugs[k], " ======== start\n", sep = "")

        anno_1 <- anno[
            which(anno$Compound == drugs[k] & anno[, "ActArea"] != -9),
        ] # 筛选出该化合物有效数据
        intersect(anno_1[, 1], ss_pp) -> shared_samples # 筛选出和CCLE有交集的数据

        grep("HAEMATOPOIETIC", shared_samples) -> haema_ii
        shared_samples <- shared_samples[-haema_ii] # 排除 海马特佩蒂克

        match(shared_samples, anno_1[, 1]) -> ii_y
        anno_2 <- anno_1[ii_y, ] # 提取有交集的数据
        y <- anno_2[, "ActArea"] # 提取ActArea(药物反应)作为以观测数据

        match(shared_samples, ss_pp) -> ii
        train_data <- pps[ii, ] # 提取CCLE中相对应数据为训练数据

        tmp_list <- list()
        mclapply(1:10, parallel_main,
            train_data, y,
            n_folds = 10, n_train_test_folds = 5,
            seed = NA, alpha = 0.5, null_testing = FALSE,
            drug = drugs[k], mc.cores = 10
        ) -> test # 对训练数据做 nested elasticnet

        for (kk in 1:10) {
            res_list <- test[[kk]]
            model <- res_list$model
            beta <- coef(model, model$lambda.min)
            n <- sum(beta != 0)
            cat("PCC = ", cor(y, res_list$self_pred[, 1]),
                "; n = ", n, "\n", sep = "") # 做参数输出
            tmp_list[[kk]] <- res_list
        }

        unlist(lapply(tmp_list,
            function(u) as.numeric(u$model_summary[5]))) -> cv_r2_avg
        which.max(cv_r2_avg) -> idx

        unlist(lapply(tmp_list,
            function(u) as.numeric(u$model_summary[4]))) -> lambda
        unlist(lapply(tmp_list,
            function(u) as.numeric(u$model_summary[9]))) -> self_r2
        unlist(lapply(tmp_list,
            function(u) as.numeric(u$model_summary[3]))) -> n
        cbind(n, lambda, cv_r2_avg, self_r2) -> x
        x[order(x[, 3]), ]

        if (sum(is.na(cv_r2_avg)) == length(cv_r2_avg)) { # 如果全无效
            model_summary <- c(drug, 0.5, 0, NA, NA, NA,
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
            tcga_probabilities <- rep(0, nrow(tcga_test_data))
        } else {
            res_list <- tmp_list[[idx]] # 取出 cv_r2_avg 最小的训练模型进行预测
            cat("selected ", idx, ", ", sep = "")
            fit <- res_list$model
            tcga_probabilities <- predict(
                fit, as.matrix(tcga_test_data), s = "lambda.min"
            )
        }
        tcga_drug_response_mat <- cbind(
            tcga_drug_response_mat, tcga_probabilities
        ) # 将不同药物的预测进行合并

        ys <- matrix(-9, nrow = nrow(pps), ncol = 2)
        rownames(ys) <- rownames(pps)
        ys[match(shared_samples, ss_pp), 1] <- y # 将药物反应放置在第一列
        # 将模型的自我预测放在第二列
        ys[match(shared_samples, ss_pp), 2] <- res_list$self_pred[, 1]
        res_list$ys <- ys
        model_list[[drugs[k]]] <- res_list # 将对应的模型和其预测按药物归纳

        # match(anno_2[, 1], self_prediction_mat[, 2]) -> pii_1
        # self_prediction_mat[pii_1, k + 2] <- res_list$self_pred

        cat(drugs[k], " end \n", sep = "")
    }

    # 整理tcga_latent的数据 使其靠拢tcga_ss_mat match
    # gsub("\\.", "-", tcga_pred[, 1]) -> ss
    # tcga_pred[, 1] <- ss
    # match(tcga_pred[, 1], tcga_ss_mat[, 1]) -> ii
    # tcga_drug_response_mat <- cbind(
    #     tcga_pred[, 1], tcga_ss_mat[ii, 2], tcga_drug_response_mat
    # )
    # colnames(tcga_drug_response_mat) <- c("TCGA", "Cancer", drugs)
    # write.table(tcga_drug_response_mat, file=paste(ksigmoid,".pred_TCGA.txt", sep=""), quote=F, sep="\t", row.names=FALSE)
    # write.table(self_prediction_mat, file=paste(ksigmoid,".pred_CCLE.txt", sep=""), quote=F, sep="\t", row.names=FALSE)
    save(model_list, file = paste(
        "./Output/2/", ksigmoid, ".CCLE.model.list.S.RData",
        sep = ""
    ))
}