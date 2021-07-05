library(ggplot2)

give_n <- function(x) {
    return(c(y = max(x), label = length(x)))
}

drug_ccle <- read.table(
    "../../Output/3/VAEN_CCLE.A.pred_TCGA.txt", header = T, as.is = T
)# 读取根据CCLE对TCGA进行的预测

colnames(drug_ccle)[3:ncol(drug_ccle)] -> drugs # 获取其中药物集合
cancer_types <- unique(drug_ccle[, 2]) # 获取癌症集合
ss <- gsub("\\.", "-", drug_ccle[, 1]) # 样本编号格式化
drug_ccle[, 1] <- ss
sample_type <- substr(drug_ccle[, 1], 14, 15)
# 取14，15位 14-15位为2位数字，01-09表示肿瘤样本，10-16表示正常对照样本

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

    tmp_ccle <- drug_ccle[
        which(drug_ccle[, 2] == cancer & sample_type == type_code),
    ] # 找到符合该癌症且为01类肿瘤样本的预测数据
    cancer_drug_ccle <- rbind(cancer_drug_ccle, tmp_ccle)
}
drug_ccle <- cancer_drug_ccle

clin_data <- read.delim(
    "BRCA_clinicalMatrix", header = T, sep = "\t"
) # 临床信息 也不什么文件我服了
type <- substr(clin_data[, 1], 14, 15) # 第一列也是TCGA样本号

# 也是筛选type_code合规的 但……
clin_data <- clin_data[which(type == type_code), ]
c(1, grep("RFS", colnames(clin_data)), grep("OS", colnames(clin_data))) -> ii
clin_data2 <- clin_data[, ii]
clin_ss <- gsub("\\.", "-", clin_data2[, 1])
rownames(clin_data2) <- clin_ss # 数据筛查赋予行名

which(
    clin_data[, "lab_proc_her2_neu_immunohistochemistry_receptor_status"]
    %in% c("Equivocal", "Negative", "Positive")
) -> ii
clin_data[
    ii, c("sampleID", "lab_proc_her2_neu_immunohistochemistry_receptor_status")
] -> stat
stat <- stat[which(stat[, 1] %in% drug_ccle[, 1]), ] # 筛选一些关键数据

brca_ccle <- drug_ccle[match(stat[, 1], drug_ccle[, 1]), ]
x <- stat[, 2]

dat4plot_ccle <- as.data.frame(
    cbind(Lapatinib = as.numeric(brca_ccle[, "Lapatinib"]),
        x = x)
)
dat4plot_ccle[, 1] <- as.numeric(as.character(dat4plot_ccle[, 1]))

dat4plot_ccle[, 2] <- factor(
    dat4plot_ccle[, 2],
    levels = c("Negative", "Equivocal", "Positive")
)

summary(glm(dat4plot_ccle[, 1] ~ as.numeric(dat4plot_ccle[, 2]))) -> sfit
p <- coef(sfit)[2, 4]

p1 <- ggplot(dat4plot_ccle, aes(x = X, y = Lapatinib, fill = X)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(shape = 21, position = position_jitter(width = 0.3), size = 0.5) +
    guides(colour = FALSE) +
    theme(axis.title.x = element_blank(), legend.title = element_blank()) +
    labs(title = paste("CCLE-based Model\np = ", format(p, digits = 3), sep = ""), x = "", y = "Predicted Response to Lapatinib") +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    stat_summary(fun.data = give_n, geom = "text")

drug.gdsc <- read.delim("../../Output/3/VAEN_GDSC.A.pred_TCGA.txt", header = T, as.is = T)
colnames(drug.gdsc)[3:ncol(drug.gdsc)] -> drugs
cancer_types <- unique(drug.gdsc[, 2])
ss <- gsub("\\.", "-", drug.gdsc[, 1])
drug.gdsc[, 1] <- ss
sample_type <- substr(drug.gdsc[, 1], 14, 15)

cancer.drug.gdsc <- c()
for (ct in 1:length(cancer_types)) {
    cancer <- cancer_types[ct]

    type_code <- "01"
    if (cancer == "LAML") {
        type_code <- "03"
    }
    if (cancer == "SKCM") {
        type_code <- "06"
    }

    tmp.gdsc <- drug.gdsc[which(drug.gdsc[, 2] == cancer & sample_type == type_code), ]
    cancer.drug.gdsc <- rbind(cancer.drug.gdsc, tmp.gdsc)
}
drug.gdsc <- cancer.drug.gdsc

########## CLIN
clin_data <- read.delim("BRCA_clinicalMatrix", header = T, sep = "\t")
type <- substr(clin_data[, 1], 14, 15)
clin_data <- clin_data[which(type == type_code), ]
c(1, grep("RFS", colnames(clin_data)), grep("OS", colnames(clin_data))) -> ii
clin_data2 <- clin_data[, ii]
clin_ss <- gsub("\\.", "-", clin_data2[, 1])
rownames(clin_data2) <- clin_ss

##############

which(clin_data[, "lab_proc_her2_neu_immunohistochemistry_receptor_status"] %in% c("Equivocal", "Negative", "Positive")) -> ii
stat <- clin_data[ii, c("sampleID", "lab_proc_her2_neu_immunohistochemistry_receptor_status")]
stat <- stat[which(stat[, 1] %in% drug.gdsc[, 1]), ]

brca.gdsc <- drug.gdsc[match(stat[, 1], drug.gdsc[, 1]), ]
X <- stat[, 2]

dat4plot.gdsc <- as.data.frame(cbind(Lapatinib = as.numeric(brca.gdsc[, "Lapatinib"]), X = X))
dat4plot.gdsc[, 1] <- as.numeric(as.character(dat4plot.gdsc[, 1]))

dat4plot.gdsc[, 2] <- factor(dat4plot.gdsc[, 2], levels = c("Negative", "Equivocal", "Positive"))

summary(glm(dat4plot.gdsc[, 1] ~ as.numeric(dat4plot.gdsc[, 2]))) -> sfit
p <- coef(sfit)[2, 4]


print(sfit)

p2 <- ggplot(dat4plot.gdsc, aes(x = X, y = Lapatinib, fill = X)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(shape = 21, position = position_jitter(width = 0.3), size = 0.5) +
    guides(colour = FALSE) +
    theme(axis.title.x = element_blank(), legend.title = element_blank()) +
    labs(title = paste("GDSC-based Model\np = ", format(p, digits = 3), sep = ""), x = "", y = "Predicted Response to Lapatinib") +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    stat_summary(fun.data = give_n, geom = "text")


source("../../Lib/multiplot.R")
png("A.ERBB2.Lapatinib.png", height = 800, width = 1000)
multiplot(plotlist = list(p1, p2), layout = matrix(c(1, 2), nrow = 1))
dev.off()
