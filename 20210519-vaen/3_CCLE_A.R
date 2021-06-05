library("MASS")
library("magrittr")
library("glmnet")
library("modEvA")
library("vegan")

##
load("./Output/1/tcga_ss_mat.RData")
########
anno <- read.csv("./DATA/CCLE/CCLE_NP24.2009_Drug_data_2015.02.24.csv", as.is = T)
drugs <- sort(unique(anno$Compound))
########

all.sample.size <- all.in_sample_R2.mat <- all.avg_CV_R2.mat <- all.F1_R2.mat <- matrix(0, nrow = 100, ncol = length(drugs))
colnames(all.sample.size) <- colnames(all.in_sample_R2.mat) <- colnames(all.avg_CV_R2.mat) <- colnames(all.F1_R2.mat) <- drugs

all.mat <- c()
for (ksigmoid in 1:100) {
    load(paste("./Output/2/", ksigmoid, ".CCLE.model.list.RData", sep = ""))
    for (kdrug in 1:length(drugs)) {
        drug <- drugs[kdrug]
        model.list[[drugs[kdrug]]] -> res_list
        if (length(res_list) == 0) next
        fit <- res_list$model
        Ys <- res_list$Ys
        which(Ys[, 1] != -9) -> ii
        Ys <- Ys[ii, ]
        if (sd(Ys[, 2]) == 0) next
        all.mat <- rbind(all.mat, c(ksigmoid, res_list$model_summary, cor(Ys[, 1], Ys[, 2])))

        #### way 4, PCC
        recall <- cor(Ys[, 1], Ys[, 2])
        precision <- as.numeric(res_list$model_summary[5])

        all.sample.size[ksigmoid, kdrug] <- nrow(Ys)
        all.in_sample_R2.mat[ksigmoid, kdrug] <- recall
        all.avg_CV_R2.mat[ksigmoid, kdrug] <- precision
        all.F1_R2.mat[ksigmoid, kdrug] <- as.numeric(res_list$model_summary[7])
    }
    cat(ksigmoid, ".", sep = "")
}

save(all.mat, all.sample.size, all.in_sample_R2.mat, all.avg_CV_R2.mat, all.F1_R2.mat, file = "./Output/3/CCLE.A.info.RData")

########

pdf("./Output/3/CCLE.A.ROC.pdf", width = 5, height = 5)
for (k in 1:length(drugs)) {
    drug <- drugs[k]
    plot(x = all.in_sample_R2.mat[, k], y = all.avg_CV_R2.mat[, k], main = drugs[k], xlab = "Self in_sample PCC", ylab = "avg PCC (in_sample)", col = rep("blue", 200), pch = 20, cex = .6)
    tmp <- cbind(idx = c(1:100), all.F1_R2.mat[, drug], all.in_sample_R2.mat[, drug], all.avg_CV_R2.mat[, drug])
    tmp <- tmp[order(tmp[, 4], decreasing = T), ]
    idx <- tmp[1:10, 1]
    points(all.in_sample_R2.mat[idx, k], all.avg_CV_R2.mat[idx, k], pch = 4, col = "red")
}
dev.off()

##########

TCGA.pred_mat <- c()
all.model_summary <- c()
holdout.R2 <- c()
for (k in 1:length(drugs)) {
    drug <- drugs[k]

    tmp <- cbind(idx = c(1:100), all.F1_R2.mat[, drug], all.in_sample_R2.mat[, drug], all.avg_CV_R2.mat[, drug])
    tmp <- tmp[order(tmp[, 4], decreasing = T), ] ### avg_CV_R2
    holdout.R2 <- rbind(holdout.R2, c(drug, tmp[1, 4]))
    best.index <- tmp[1, 1]

    load(paste("./Output/2/", best.index, ".CCLE.model.list.RData", sep = ""))
    model.list[[drug]] -> res_list
    fit <- res_list$model

    TCGA.pred <- read.table(paste("./Output/1/", best.index, ".TCGA_latent.tsv", sep = ""), header = T, sep = "\t", as.is = T)
    TCGA.test.data <- TCGA.pred[, -1]
    TCGA.probabilities <- predict(fit, as.matrix(TCGA.test.data), s = "lambda.min")

    TCGA.pred_mat <- cbind(TCGA.pred_mat, TCGA.probabilities)
    cat("...", drug, ".", sep = "")
}

TCGA.pred_mat <- cbind(TCGA.pred[, 1], "A", TCGA.pred_mat)
gsub("\\.", "-", TCGA.pred_mat[, 1]) -> ss
TCGA.pred_mat[, 1] <- ss
match(TCGA.pred_mat[, 1], tcga_ss_mat[, 1]) -> ii
TCGA.pred_mat[, 2] <- tcga_ss_mat[ii, 2]
colnames(TCGA.pred_mat) <- c("Sample", "Cancer", drugs)
write.table(TCGA.pred_mat, file="./Output/3/VAEN_CCLE.A.pred_TCGA.txt", quote = F, sep = "\t", row.names = FALSE)

############################################################################################

CCLE.PCC <- c()
for (kdrug in 1:length(drugs)) {
    drug <- drugs[kdrug]
    if (drug == "X17.AAG") drug <- "17-AAG"
    gsub("\\.", "-", drug) -> drug

    tmp <- cbind(idx = c(1:100), all.F1_R2.mat[, drug], all.in_sample_R2.mat[, drug], all.avg_CV_R2.mat[, drug])
    tmp <- tmp[order(tmp[, 4], decreasing = T), ]

    pred_mat <- c()
    best.index <- tmp[1, 1]

    load(paste("./Output/2/", best.index, ".CCLE.model.list.RData", sep = ""))
    model.list[[drug]] -> res_list
    Ys <- res_list$Ys

    if (kdrug == 1) {
        self_prediction_mat <- Ys
    } else {
        self_prediction_mat <- cbind(self_prediction_mat, Ys[, 2])
    }
}

self_prediction_mat[, 1] <- rownames(Ys)
colnames(self_prediction_mat) <- c("CELLLINE", drugs)
write.table(self_prediction_mat, file="./Output/3/VAEN_CCLE.A.pred_CCLE.txt", quote = F, sep = "\t", row.names = FALSE)

########

CCLE.pred.full.mat <- c()
for (k in 1:length(drugs)) {
    drug <- drugs[k]
    tmp <- cbind(idx = c(1:100), all.F1_R2.mat[, drug], all.in_sample_R2.mat[, drug], all.avg_CV_R2.mat[, drug])
    tmp <- tmp[order(tmp[, 4], decreasing = T), ]

    pred_mat <- c()
    best.index <- tmp[1, 1]

    load(paste("./Output/2/", best.index, ".CCLE.model.list.RData", sep = ""))
    model.list[[drug]] -> res_list

    CCLE.latent <- read.table(paste("./Output/1/", best.index, ".CCLE_latent.tsv", sep = ""), header = T, sep = "\t", as.is = T)
    CCLE.latent.data <- CCLE.latent[, -1]
    fit <- res_list$model
    CCLE.probabilities <- predict(fit, as.matrix(CCLE.latent.data), s = "lambda.min")

    CCLE.pred.full.mat <- cbind(CCLE.pred.full.mat, CCLE.probabilities)
    cat(drug, ".", sep = "")
}

self.pred_mat <- cbind(CELLINE = CCLE.latent[, 1], CCLE.pred.full.mat)
colnames(self.pred_mat) <- c("CELLINE", drugs)
write.table(self.pred_mat, file = paste("./Output/3/VAEN_CCLE.A.pred_CCLE.full.txt", sep = ""), quote = F, sep = "\t", row.names = FALSE)

########