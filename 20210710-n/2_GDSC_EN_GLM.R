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

# load("./Output/1/tcga_ss_mat.RData")
anno <- read.delim("./Data/GDSC/v17.3_fitted_dose_response.txt", as.is = T)
drugs <- sort(unique(anno$DRUG_NAME))
cell_line_anno <- read.csv("./Data/CCLE/DepMap-2018q3-celllines.csv", as.is = T)

for (ksigmoid in start:end) {
    cat("ksigmoid = ", ksigmoid, "\n", sep = "")

    # model_summary_file <- paste(
    #     "./Output/2/", ksigmoid,
    #     ".model_summary.txt",
    #     sep = ""
    # )
    # model_summary_cols <- c(
    #     "Drug", "alpha", "n_snps_in_model", "lambda_min_mse",
    #     "test_R2_avg", "test_R2_sd", "cv_r2_avg", "cv_R2_sd", "in_sample_R2",
    #     "nested_cv_fisher_pval", "rho_avg", "rho_se", "rho_zscore",
    #     "rho_avg_squared", "zscore_pval", "cv_rho_avg", "cv_rho_se",
    #     "cv_rho_avg_squared", "cv_zscore_est", "cv_zscore_pval", "cv_pval_est"
    # ) # 声明一些统计信息

    # write(model_summary_cols, file = model_summary_file, ncol = 21, sep = "\t")

    # tcga_pred <- read.table(
    #     paste("./Output/1/", ksigmoid, ".TCGA_latent.tsv", sep = ""),
    #     header = T, sep = "\t", as.is = T
    # )
    # tcga_test_data <- tcga_pred[, -1]

    # Prediction
    pps <- read.table(
        paste("./Output/1/", ksigmoid, ".CCLE_latent.tsv", sep = "")
    )
    original_ss_pp <- rownames(pps)
    sapply(original_ss_pp, function(x) {
        new_u <- u <- strsplit(x, split = "\\.")[[1]][1]
        if (grepl("^X", u)) {
            substr(u, 2, nchar(u)) -> new_u
        }
        new_u
    }) -> ss_pp
    names(ss_pp) <- NULL

    original_ss_pp <- rownames(pps)
    sapply(original_ss_pp, function(x) {
        new_u <- u <- strsplit(x, split = "\\.")[[1]]
        paste(u[3], u[4], sep = "-") -> new_u
        new_u
    }) -> ss_ach
    names(ss_ach) <- NULL

    # original drug data
    # self_prediction_mat <- matrix(-9,
    #     nrow = length(unique(anno[, 4])),
    #     ncol = length(drugs) + 2
    # )
    # self_prediction_mat[, 1] <- unique(anno[, "COSMIC_ID"])
    # self_prediction_mat[, 2] <- unique(anno[, "CELL_LINE_NAME"])
    # colnames(self_prediction_mat) <- c("CELLLINE", "Type", drugs)
    # 初始化预测矩阵 细胞系 类别 药物1 药物2...

    tcga_drug_response_mat <- c()
    model_list <- list()
    for (k in seq_len(length(drugs))) {
        drug <- drugs[k]
        cat(drug, " ======== start\n", sep = "")

        anno_1 <- anno[which(anno$DRUG_NAME == drug), ] # 找到GDSC中符合drug的数据
        # 找到和CCLE的COSMIC_ID相符的数据索引
        # match(anno_1$COSMIC_ID, cell_line_anno$COSMIC_ID) -> idx
        anno_1_match1 <- anno_1[
            which(anno_1$COSMIC_ID %in% cell_line_anno$COSMIC_ID),
        ]

        gdsc_anno_1_match2 <- anno_1[
            which(!anno_1$COSMIC_ID %in% cell_line_anno$COSMIC_ID),
        ]

        match(anno_1_match1$COSMIC_ID, cell_line_anno$COSMIC_ID) -> idx1
        match(cell_line_anno[idx1, 1], ss_ach) -> part1_ach_idx
        part1_ach_idx <- part1_ach_idx[which(!is.na(part1_ach_idx))]
        anno_1_match1 <- anno_1_match1[part1_ach_idx, ]

        y <- -anno_1_match1[, "LN_IC50"]

        sapply(gdsc_anno_1_match2$CELL_LINE_NAME, function(u) {
            gsub("-", "", u)
        }) -> gdsc_cell_name
        union(
            which(gdsc_cell_name %in% cell_line_anno$Aliases),
            which(gdsc_anno_1_match2$CELL_LINE_NAME %in% cell_line_anno$Aliases)
        ) -> ii
        gdsc_anno_1_match2 <- gdsc_anno_1_match2[ii, ]
        gdsc_cell_name <- gdsc_cell_name[ii]

        idx2 <- c()
        for (kk in seq_len(nrow(gdsc_anno_1_match2))) {
            if (gdsc_anno_1_match2$CELL_LINE_NAME[kk]
            %in% cell_line_anno$Aliases) {
                idx2 <- c(idx2, match(
                    gdsc_anno_1_match2$CELL_LINE_NAME[kk],
                    cell_line_anno$Aliases
                ))
            } else {
                idx2 <- c(idx2, match(
                    gdsc_cell_name[kk],
                    cell_line_anno$Aliases
                ))
            }
        }

        match(cell_line_anno[idx2, 1], ss_ach) -> part2_ach_idx
        gdsc_anno_1_match2 <- gdsc_anno_1_match2[which(!is.na(part2_ach_idx)), ]

        idx2 <- c()
        for (kk in seq_len(nrow(gdsc_anno_1_match2))) {
            if (gdsc_anno_1_match2$CELL_LINE_NAME[kk]
            %in% cell_line_anno$Aliases) {
                idx2 <- c(idx2, match(
                    gdsc_anno_1_match2$CELL_LINE_NAME[kk],
                    cell_line_anno$Aliases
                ))
            } else {
                idx2 <- c(idx2, match(
                    gsub("-", "", gdsc_anno_1_match2$CELL_LINE_NAME[kk]),
                    cell_line_anno$Aliases
                ))
            }
        }

        match(cell_line_anno[idx2, 1], ss_ach) -> part2_ach_idx
        new_y <- -gdsc_anno_1_match2[, "LN_IC50"]
        names(new_y) <- gdsc_anno_1_match2[, "CELL_LINE_NAME"]
        # 将gdsc中ID可以对上的数据的反应(LN_IC50)值
        # 和不能对上，只能使用别名对应的数据进行合并
        y <- c(y, new_y)

        train_data <- pps[c(part1_ach_idx, part2_ach_idx), ]

        tmp_list <- list()
        if (drug == "VNLG/124") drug <- "VNLG.124"

        mclapply(1:10, parallel_main,
            train_data, y,
            n_folds = 10, n_train_test_folds = 5,
            seed = NA, alpha = 0.5, null_testing = FALSE,
            drug = drug, mc.cores = 10
        ) -> test # 对训练数据做 nested elasticnet 回归

        for (kk in 1:10) {
            res_list <- test[[kk]]
            model <- res_list$model
            beta <- coef(model, model$lambda.min)
            n <- sum(beta != 0)
            cat("PCC = ", cor(y, res_list$self_pred[, 1]),
                "; n = ", n, "\n", sep = "")
            tmp_list[[kk]] <- res_list
        }

        unlist(lapply(tmp_list,
            function(u) as.numeric(u$model_summary[5]))) -> cv_r2_avg
        which.max(cv_r2_avg) -> idx

        if (sum(is.na(cv_r2_avg)) == length(cv_r2_avg)) {
            model_summary <- c(drug, 0.5, 0, NA, NA, NA, NA,
                NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
            # write(model_summary, file = model_summary_file,
            #     append = TRUE, ncol = 24, sep = "\t")
            # tcga_probabilities <- rep(0, nrow(tcga_test_data))
        } else {
            res_list <- tmp_list[[idx]]
            cat("selected ", idx, ", ", sep = "")
            # write(res_list$model_summary, file = model_summary_file,
            #     append = TRUE, ncol = 24, sep = "\t")
            # fit <- res_list$model
            # tcga_probabilities <- predict(
            #     fit, as.matrix(tcga_test_data), s = "lambda.min"
            # )
        }

        # tcga_drug_response_mat <- cbind(
        #     tcga_drug_response_mat, tcga_probabilities
        # ) # 将不同药物的预测进行合并

        ys <- matrix(-9, nrow = nrow(pps), ncol = 2)
        rownames(ys) <- rownames(pps)
        ys[c(part1_ach_idx, part2_ach_idx), 1] <- y
        ys[c(part1_ach_idx, part2_ach_idx), 2] <- res_list$self_pred
        res_list$ys <- ys

        model_list[[drugs[k]]] <- res_list

        # match(anno_1_match1$COSMIC_ID, self_prediction_mat[, 1]) -> pii_1
        # match(
        #     gdsc_anno_1_match2$CELL_LINE_NAME,
        #     self_prediction_mat[, 2]
        # ) -> pii_2
        # self_prediction_mat[c(pii_1, pii_2), k + 2] <- res_list$self_pred

        cat(drugs[k], " end \n", sep = "")
    }

    # gsub("\\.", "-", tcga_pred[, 1]) -> ss
    # tcga_pred[, 1] <- ss
    # match(tcga_pred[, 1], tcga_ss_mat[, 1]) -> ii
    # tcga_drug_response_mat <- cbind(
    #     tcga_pred[, 1],
    #     tcga_ss_mat[ii, 2],
    #     tcga_drug_response_mat
    # )
    # colnames(tcga_drug_response_mat) <- c("TCGA", "Cancer", drugs)
    # write.table(tcga_drug_response_mat,
    #     file = paste("./Output/2/", ksigmoid, ".GDSC.pred_TCGA.txt", sep = ""),
    #     quote = F, sep = "\t", row.names = FALSE
    # )
    # write.table(self_prediction_mat,
    #     file = paste("./Output/2/", ksigmoid, ".pred_GDSC.txt", sep = ""),
    #     quote = F, sep = "\t", row.names = FALSE
    # )
    save(model_list, file = paste(
        "./Output/2/", ksigmoid, ".GDSC.model.list.RData", sep = ""
    ))
}
