library(ggplot2)
load("../../Output/4/dr.CCLE.A.models.RData")
drugs <- names(dr_ccle_models)

tcga_pred <- read.table(
    "../../Output/3/VAEN_CCLE.A.pred_TCGA.txt",
    header = T, as.is = T, sep = "\t"
)

dr_ccle_mat <- c()
for (kdrug in seq_len(length(drugs))) {
    drug <- drugs[kdrug]

    dr_ccle_models[[drug]] -> res_list
    ys <- res_list$ys

    which(ys[, 1] != -9) -> ii
    ys <- ys[ii, ]
    dr_ccle_mat <- rbind(dr_ccle_mat,
        cbind(
            drug = drug,
            sample = rownames(ys),
            ActArea = ys[, 1],
            grp = "CCLE Observed DR"
        )
    )

    dr_ccle_mat <- rbind(dr_ccle_mat,
        cbind(drug = drug,
            sample = rownames(ys),
            ActArea = ys[, 2],
            grp = "CCLE Predicted DR"
        )
    )

    tcga_drug_name <- gsub("-", ".", drug)
    if (drug == "17-AAG") tcga_drug_name <- "X17.AAG"
    cur_drug_tcga_pred <- tcga_pred[, tcga_drug_name]
    dr_ccle_mat <- rbind(dr_ccle_mat,
        cbind(
            drug = drug,
            sample = gsub("-", ".", tcga_pred[, 1]
        ),
        ActArea = cur_drug_tcga_pred,
        grp = "TCGA Predicted DR")
    )
}

new_mat <- as.data.frame(dr_ccle_mat)
rownames(new_mat) <- NULL
new_mat[, 3] <- as.numeric(as.character(new_mat[, 3]))
write.table(new_mat, file = "./A.txt", row.names = F, quote = F, sep = "\t")

g <- ggplot(aes(y = ActArea, x = drug, fill = grp), data = new_mat) +
    geom_boxplot() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = c(0.9, 0.9)) +
            xlab("") +
            ylab("Drug response (ActArea)") +
            guides(fill = guide_legend(title = "")
    )
ggsave(g, file = "./A.png", width = 10, height = 5)
