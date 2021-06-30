end <- 1
load("./Output/3/CCLE.A.info.RData")
drugs <- colnames(all_f1_r2_mat)
dr_ccle_models <- list()
for (kdrug in seq_len(length(drugs))) {
    drug <- drugs[kdrug]
    if (drug == "X17.AAG") drug <- "17-AAG"
    gsub("\\.", "-", drug) -> drug

    tmp <- cbind(
        idx = c(1:end), all_f1_r2_mat[, drug],
        all_in_sample_r2_mat[, drug], all_avg_cv_r2_mat[, drug]
    )
    tmp <- tmp[order(tmp[, 4], decreasing = T), ]
    best_index <- tmp[1, 1] # 再次选出每个药品表现最好的模型

    load(paste("./Output/2/", best_index, ".CCLE.model.list.RData", sep = ""))
    model_list[[drug]] -> res_list
    res_list[["best_index"]] <- best_index
    dr_ccle_models[[drug]] <- res_list
}
# 将每种药品下表现最好的模型进行存储
save(dr_ccle_models, file = "./Output/4/dr.CCLE.A.models.RData")

load("./Output/3/CCLE.S.info.RData")
dr_ccle_models <- list()
for (kdrug in seq_len(length(drugs))) {
    drug <- drugs[kdrug]
    if (drug == "X17.AAG") drug <- "17-AAG"
    gsub("\\.", "-", drug) -> drug
    tmp <- cbind(
        idx = c(1:end), solid_f1_r2_mat[, drug],
        solid_in_sample_r2_mat[, drug], solid_avg_cv_r2_mat[, drug]
    )
    tmp <- tmp[order(tmp[, 4], decreasing = T), ]
    best_index <- tmp[1, 1]

    load(paste("./Output/2/", best_index, ".CCLE.model.list.S.RData", sep = ""))
    model_list[[drug]] -> res_list
    res_list[["best_index"]] <- best_index
    dr_ccle_models[[drug]] <- res_list
}

save(dr_ccle_models, file = "./Output/4/dr.CCLE.S.models.RData") # 同上