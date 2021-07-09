end <- 2
load("./Output/3/GDSC.A.info.RData")

anno <- read.delim("./DATA/GDSC/v17.3_fitted_dose_response.txt", as.is = T)
drugs <- sort(unique(anno$DRUG_NAME)) # 只为了取药品名称

dr_gdsc_models <- list()
for (kdrug in seq_len(length(drugs))) {
    drug <- drugs[kdrug]
    tmp <- cbind(idx = c(1:end), all_f1_r2_mat[, drug],
        all_in_sample_r2_mat[, drug], all_avg_cv_r2_mat[, drug])
    tmp <- tmp[order(tmp[, 4], decreasing = T), ]
    best_index <- tmp[1, 1]

    load(paste("./Output/2/", best_index, ".GDSC.model.list.RData", sep = ""))
    model_list[[drug]] -> res_list
    res_list[["best_index"]] <- best_index
    dr_gdsc_models[[drug]] <- res_list
}

# 将每种药品下表现最好的模型进行存储
save(dr_gdsc_models, file = "./Output/4/dr.GDSC.A.models.RData")
