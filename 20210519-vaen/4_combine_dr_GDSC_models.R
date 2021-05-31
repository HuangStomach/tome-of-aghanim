load("./Output/3/GDSC.A.info.RData")

anno <- read.delim("./DATA/GDSC/v17.3_fitted_dose_response.txt", as.is = T)
drugs <- sort(unique(anno$DRUG_NAME))

dr_gdsc_models <- list()
for (kdrug in 1:length(drugs)) {
    drug <- drugs[kdrug]
    tmp <- cbind(idx = c(1:100), all.F1_R2.mat[, drug],
        all.in_sample_R2.mat[, drug], all.avg_CV_R2.mat[, drug])
    tmp <- tmp[order(tmp[, 4], decreasing = T), ]
    best_index <- tmp[1, 1]

    load(paste("./Output/2/", best_index, ".GDSC.model.list.RData", sep = ""))
    model.list[[drug]] -> res_list
    res_list[["best_index"]] <- best_index
    dr_gdsc_models[[drug]] <- res_list
}

save(dr_gdsc_models, file = "./Output/4/dr.GDSC.A.models.RData")
