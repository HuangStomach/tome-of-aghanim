drug_ccle <- read.table(
    "../../Output/3/VAEN_CCLE.A.pred_TCGA.txt",
    header = T, as.is = T
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

    blca_ccle <- drug_ccle[which(
        drug_ccle[, 2] == cancer & sample_type == type_code
    ), ]
    cancer_drug_ccle <- rbind(cancer_drug_ccle, blca_ccle)
}
drug_ccle <- cancer_drug_ccle

library("survival")
library("survminer")

response <- read.delim("../../Data/response/drug_response.txt", as.is = T)
response <- response[which(
    response$drug.name == "Paclitaxel" & response$cancers == "BRCA"
), ]

match(response[, 2], substr(drug_ccle[, 1], 1, 12)) -> ii
cbind(response, drug_ccle[ii, ]) -> new2
new2 <- new2[!is.na(ii), ]
dim(new2)

brca_clin_data <- read.delim("BRCA_clinicalMatrix", as.is = T)
match(new2[, 2], substr(brca_clin_data[, 1], 1, 12)) -> ii
brca_clin_data <- brca_clin_data[ii, ]

drug <- "Paclitaxel"

samples <- brca_clin_data[, 1]
match(samples, drug_ccle[, 1]) -> ii
drug_response <- drug_ccle[ii, drug]

surv_data <- brca_clin_data[match(samples, brca_clin_data[, 1]), ]

new3 <- cbind(drug_response, surv_data)
y1 <- Surv(new3[, "X_OS"], new3[, "X_OS_IND"])

xvector <- ifelse(new3[, 1] > median(new3[, 1]), "HR", "LR")
table(xvector)

dat <- data.frame(cbind(new3, xvector))
print(y1)
quit()
fit <- survfit(Surv(X_OS, X_OS_IND) ~ xvector, data = dat)
coxph(y1 ~ xvector)


fit <- survfit(Surv(
    as.numeric(X_OS), as.numeric(X_OS_IND)
) ~ xvector, data = dat)
g1 <- ggsurvplot(
    fit, data = dat,
    risk.table = TRUE, pval = TRUE,
    ggtheme = theme_minimal()
)
ggsave("./C.CCLE.Paclitaxel.BRCA.png", plot = print(g1), width = 5, height = 5)

drug_ccle <- read.table(
    file = "../../Output/3/VAEN_GDSC.A.pred_TCGA.txt",
    header = T, as.is = T, sep = "\t"
)

cancer_types <- unique(drug_ccle[, 2])
sample_type <- substr(drug_ccle[, 1], 14, 15)
cancer_drug_ccle <- c()
for (ct in seq_len(length(cancer_types))) {
    cancer <- cancer_types[ct]

    type_code <- "01"
    if (cancer == "LAML") {
        type_code <- "03"
    }
    if (cancer == "SKCM") {
        type_code <- "06"
    }

    blca_ccle <- drug_ccle[which(
        drug_ccle[, 2] == cancer & sample_type == type_code
    ), ]
    cancer_drug_ccle <- rbind(cancer_drug_ccle, blca_ccle)
}
drug_ccle <- cancer_drug_ccle


response <- read.delim("../../Data/response/drug_response.txt", as.is = T)
response <- response[which(response$drug.name == "Fluorouracil"), ]

match(response[, 2], substr(drug_ccle[, 1], 1, 12)) -> ii
cbind(response, drug_ccle[ii, ]) -> new2
new2 <- new2[!is.na(ii), ]
dim(new2)

new2 <- new2[which(new2[, 1] == "STAD"), ]

brca_clin_data <- read.delim("STAD_clinicalMatrix", as.is = T)
match(new2[, 2], substr(brca_clin_data[, 1], 1, 12)) -> ii
stad_clin_data <- brca_clin_data[ii, ]

library(survival)

drug <- "X5.Fluorouracil"

samples <- stad_clin_data[, 1]
match(samples, drug_ccle[, 1]) -> ii

drug_response <- drug_ccle[ii, drug]
surv_data <- stad_clin_data[match(samples, stad_clin_data[, 1]), ]

new3 <- cbind(drug_response, surv_data)
y1 <- Surv(new3[, "X_OS"], new3[, "X_OS_IND"])

xvector <- ifelse(new3[, 1] > median(new3[, 1]), "HR", "LR")
table(xvector)

dat <- data.frame(cbind(new3, xvector))
fit <- survfit(Surv(X_OS, X_OS_IND) ~ xvector, data = dat)
coxph(y1 ~ xvector)

g1 <- ggsurvplot(
    fit, data = dat,
    risk.table = TRUE, pval = TRUE,
    break.time.by = 500, ggtheme = theme_minimal()
)
ggsave("./B.png", plot = print(g1), width = 5, height = 5)
