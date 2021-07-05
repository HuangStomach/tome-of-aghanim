load("../../Output/2/1.CCLE.model.list.RData")

ccle_model_list <- model_list
ccle_obsd <- c()
for (k in seq_len(length(ccle_model_list))) {
    res_list <- ccle_model_list[[k]]
    drug <- names(ccle_model_list)[k]
    ys <- res_list$ys

    ccle_obsd <- cbind(ccle_obsd, ys[, 1])
}
colnames(ccle_obsd) <- names(ccle_model_list)

ccle_pred <- read.table(
    file = "../../Output/3/VAEN_CCLE.A.pred_CCLE.txt",
    header = T, as.is = T, sep = "\t"
)
ccle_pred_full <- read.table(
    file = "../../Output/3/VAEN_CCLE.A.pred_CCLE.full.txt",
    header = T, as.is = T, sep = "\t"
)

colnames(ccle_pred) -> ss
ss[which(ss == "X17.AAG")] <- "17-AAG"
gsub("\\.", "-", ss) -> ss
colnames(ccle_pred) <- ss

colnames(ccle_pred_full) -> ss
ss[which(ss == "X17.AAG")] <- "17-AAG"
gsub("\\.", "-", ss) -> ss
colnames(ccle_pred_full) <- ss

# Prediction
original_ss_pp <- rownames(ccle_obsd)
sapply(original_ss_pp, function(x) {
    strsplit(x, split = "\\.")[[1]][1] -> u
    strsplit(u, split = "_")[[1]] -> v
    v <- v[-1]
    paste(v, collapse = "_")
}) -> tt
names(tt) <- NULL

anno <- read.csv(
    "../../Data/CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv",
    as.is = T
)
drugs <- sort(unique(anno$Compound))

matrix(1, nrow = length(drugs), ncol = length(unique(tt))) ->
    obsd_drug_by_tissue_mat ->
    pred_drug_by_tissue_mat ->
    full_drug_by_tissue_mat

drugs ->
    rownames(full_drug_by_tissue_mat) ->
    rownames(pred_drug_by_tissue_mat) ->
    rownames(obsd_drug_by_tissue_mat)

sort(unique(tt)) ->
    colnames(full_drug_by_tissue_mat) ->
    colnames(pred_drug_by_tissue_mat) ->
    colnames(obsd_drug_by_tissue_mat)

for (kdrug in seq_len(length(drugs))) {
    drug <- drugs[kdrug]
    ys <- cbind(ccle_obsd[, drug], ccle_pred[, drug], ccle_pred_full[, drug])

    sapply(rownames(ys), function(x) {
        strsplit(x, split = "\\.")[[1]][1] -> u
        strsplit(u, split = "_")[[1]] -> v
        v <- v[-1]
        paste(v, collapse = "_")
    }) -> tt
    names(tt) <- NULL

    for (kt in seq_len(ncol(obsd_drug_by_tissue_mat))) {
        x_tissue <- colnames(obsd_drug_by_tissue_mat)[kt]

        x <- ifelse(tt %in% x_tissue, 1, 0)

        which(ys[, 1] != -9) -> ii
        y1 <- ys[ii, 1]
        x1 <- x[ii]
        sum(x1 == 1) -> check
        label <- ifelse(
            mean(y1[which(x1 == 1)]) >
            mean(y1[which(x1 == 0)]), 1, -1
        )
        if (check > 5) {
            label * t.test(y1[which(x1 == 1)], y1[which(x1 == 0)])$p.value ->
            obsd_drug_by_tissue_mat[kdrug, kt]
        }

        which(ys[, 1] != -9) -> ii
        y2 <- ys[ii, 2]
        x2 <- x[ii]
        sum(x2 == 1) -> check
        label <- ifelse(
            mean(y2[which(x2 == 1)]) >
            mean(y2[which(x2 == 0)]), 1, -1
        )
        if (check > 5) {
            pred_drug_by_tissue_mat[kdrug, kt] <- label * t.test(
                y2[which(x2 == 1)], y2[which(x2 == 0)]
            )$p.value
        }

        y3 <- ys[, 3]
        label <- ifelse(
            mean(y3[which(x == 1)]) >
            mean(y3[which(x == 0)]), 1, -1
        )
        full_drug_by_tissue_mat[kdrug, kt] <- label * t.test(
            y3[which(x == 1)], y3[which(x == 0)]
        )$p.value
    }
}

apply(obsd_drug_by_tissue_mat, 2, var) -> obsd_check
apply(pred_drug_by_tissue_mat, 2, var) -> pred_check

print(which(obsd_check == 0))
print(which(pred_check == 0))

obsd_drug_by_tissue_mat <- obsd_drug_by_tissue_mat[, which(obsd_check != 0)]
pred_drug_by_tissue_mat <- pred_drug_by_tissue_mat[, which(obsd_check != 0)]
full_drug_by_tissue_mat <- full_drug_by_tissue_mat[, which(obsd_check != 0)]

drug_ccle <- read.table(
    file = "../../Output/3/VAEN_CCLE.A.pred_TCGA.txt",
    header = T, as.is = T, sep = "\t"
)
colnames(drug_ccle)[3:ncol(drug_ccle)] -> drugs
cancer_types <- unique(drug_ccle[, 2])
sample_type <- substr(drug_ccle[, 1], 14, 15)
ss <- gsub("\\.", "-", drug_ccle[, 1])
drug_ccle[, 1] <- ss

cancer_drug_ccle <- c()
for (ct in 1:length(cancer_types)) {
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
drug_ccle <- cancer_drug_ccle

tcga_drug_by_cancer_mat <- matrix(
    1, nrow = length(drugs), ncol = length(cancer_types)
)
rownames(tcga_drug_by_cancer_mat) <- drugs
colnames(tcga_drug_by_cancer_mat) <- sort(cancer_types)

for (kdrug in seq_len(length(drugs))) {
    y <- drug_ccle[, drugs[kdrug]]

    for (kt in seq_len(ncol(tcga_drug_by_cancer_mat))) {
        x_cancer <- colnames(tcga_drug_by_cancer_mat)[kt]
        x <- ifelse(drug_ccle[, 2] %in% x_cancer, 1, 0)
        if (sum(x == 1) < 5) next
        label <- ifelse(mean(y[which(x == 1)]) > mean(y[which(x == 0)]), 1, -1)
        p <- t.test(y[which(x == 1)], y[which(x == 0)])$p.value
        if (p < 1e-100) p <- 1e-101
        tcga_drug_by_cancer_mat[kdrug, kt] <- label * p
    }
}

drug_ccle <- read.table(
    file = "../../Output/3/VAEN_CCLE.A.pred_TCGA.txt",
    header = T, as.is = T, sep = "\t"
)
immune_cancer <- c("LAML", "DLBC", "THYM")
which(drug_ccle[, 2] %in% immune_cancer) -> ii
drug_ccle <- drug_ccle[-ii, ]
colnames(drug_ccle)[3:ncol(drug_ccle)] -> drugs
cancer_types <- unique(drug_ccle[, 2])
sample_type <- substr(drug_ccle[, 1], 14, 15)
ss <- gsub("\\.", "-", drug_ccle[, 1])
drug_ccle[, 1] <- ss

cancer_drug_ccle <- c()
for (ct in 1:length(cancer_types)) {
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
drug_ccle <- cancer_drug_ccle

tcga_sensitive_mat <- matrix(
    1, nrow = length(drugs), ncol = length(cancer_types)
)
rownames(tcga_sensitive_mat) <- drugs
colnames(tcga_sensitive_mat) <- sort(cancer_types)
sensitive_prop <- c()

for (kdrug in 1:length(drugs)) {
    y <- drug_ccle[, drugs[kdrug]]
    y > quantile(y, probs = .95) -> sensitive_ii
    for (kt in 1:ncol(tcga_sensitive_mat)) {
        x_cancer <- colnames(tcga_sensitive_mat)[kt]
        x <- drug_ccle[, 2] %in% x_cancer
        if (kdrug == 18) {
            sensitive_prop <- c(sensitive_prop, sum(x & sensitive_ii) / sum(x))
        }

        if (sum(x & sensitive_ii) < 5) next

        table(sensitive_ii, x) -> mm
        mm[2:1, 2:1] -> new_mm
        tcga_sensitive_mat[kdrug, kt] <- fisher.test(new_mm)$p.value
    }
}
names(sensitive_prop) <- colnames(tcga_sensitive_mat)

tcga_resistant_mat <- matrix(
    1, nrow = length(drugs), ncol = length(cancer_types)
)
rownames(tcga_resistant_mat) <- drugs
colnames(tcga_resistant_mat) <- sort(cancer_types)
resistant_prop <- c()

for (kdrug in seq_len(length(drugs))) {
    y <- drug_ccle[, drugs[kdrug]]
    y < quantile(y, probs = .05) -> resistant_ii
    for (kt in 1:ncol(tcga_sensitive_mat)) {
        x_cancer <- colnames(tcga_sensitive_mat)[kt]
        x <- drug_ccle[, 2] %in% x_cancer
        if (kdrug == 18) {
            resistant_prop <- c(resistant_prop, sum(x & resistant_ii) / sum(x))
        }
        if (sum(x & resistant_ii) < 5) next

        table(resistant_ii, x) -> mm
        mm[2:1, 2:1] -> new_mm
        tcga_resistant_mat[kdrug, kt] <- fisher.test(new_mm)$p.value
    }
}
names(resistant_prop) <- colnames(tcga_sensitive_mat)

plx4720_pred_dr <- drug_ccle[, 18]

rownames(tcga_sensitive_mat) <- rownames(pred_drug_by_tissue_mat)
rownames(tcga_resistant_mat) <- rownames(pred_drug_by_tissue_mat)

save(
    resistant_prop, sensitive_prop,
    plx4720_pred_dr, obsd_drug_by_tissue_mat,
    pred_drug_by_tissue_mat, full_drug_by_tissue_mat,
    tcga_drug_by_cancer_mat, tcga_sensitive_mat,
    tcga_resistant_mat,
    file = "./BCD.data.RData"
)

write.table(tcga_sensitive_mat, file = "./C.sensitive.txt", quote = F)
write.table(tcga_resistant_mat, file = "./C.resistant.txt", quote = F)