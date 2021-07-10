library(ggplot2)
source("../../Lib/multiplot.R")

give_n <- function(x) {
    return(c(y = max(x), label = length(x)))
}

drug_ccle <- read.table(
    "../../Output/3/VAEN_CCLE.A.pred_TCGA.txt",
    header = T, as.is = T
)
cancer_types <- unique(drug_ccle[, 2])
sample_type <- substr(drug_ccle[, 1], 14, 15)
ss <- gsub("\\.", "-", drug_ccle[, 1])
drug_ccle[, 1] <- ss

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

    blca_ccle <- drug_ccle[
        which(drug_ccle[, 2] == cancer & sample_type == type_code),
    ]
    cancer_drug_ccle <- rbind(cancer_drug_ccle, blca_ccle)
}
drug_ccle <- cancer_drug_ccle

cancer_ccle <- drug_ccle[which(drug_ccle[, 2] == "LUAD"), ]

load("LUAD.MET_amp.ss.RData")
met_exon14_ii <- match(gsub("\\.", "-", ss), cancer_ccle[, 1])
match(LUAD.MET_amp.ss, cancer_ccle[, 1]) -> met_amp_ii
met_wt_ii <- setdiff(1:nrow(cancer_ccle), met_amp_ii)

###
dat4plot <- data.frame(cbind(
    y = cancer_ccle[, "PF2341066"],
    grp = ifelse(cancer_ccle[, 1] %in% LUAD.MET_amp.ss, "MET Gain", "Other")
))
dat4plot[, 1] <- as.numeric(as.character(dat4plot[, 1]))

p <- t.test(dat4plot[, 1] ~ dat4plot[, 2])$p.value

p1 <- ggplot(dat4plot, aes(y = y, x = grp)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
        shape = 21,
        position = position_jitter(width = 0.3),
        size = 0.5
    ) +
    guides(colour = FALSE) +
    theme(axis.title.x = element_blank(), legend.title = element_blank()) +
    labs(
        title = paste(
            "PF2341066, MET CNV\np = ",
            format(p, digits = 3), sep = ""
        ),
        x = "", y = "Predicted Response to PF2341066"
    ) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    stat_summary(fun.data = give_n, geom = "text")

###
dat4plot <- data.frame(cbind(
    y = cancer_ccle[, "PHA.665752"],
    grp = ifelse(cancer_ccle[, 1] %in% LUAD.MET_amp.ss,
    "MET Gain", "Other")
))
dat4plot[, 1] <- as.numeric(as.character(dat4plot[, 1]))
p <- t.test(dat4plot[, 1] ~ dat4plot[, 2])$p.value

p2 <- ggplot(dat4plot, aes(y = y, x = grp)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
        shape = 21,
        position = position_jitter(width = 0.3),
        size = 0.5
    ) +
    guides(colour = FALSE) +
    theme(axis.title.x = element_blank(), legend.title = element_blank()) +
    labs(
        title = paste(
            "PHA665752, MET CNV\np = ",
            format(p, digits = 3),
            sep = ""
        ), x = "", y = "Predicted Response to PHA665752"
    ) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    stat_summary(fun.data = give_n, geom = "text")

load("LUAD.MET.expr.RData")

x1 <- rep("1", length(gene.expr))
x1[which(gene.expr < quantile(gene.expr, probs = .25))] <- "0"
x1[which(gene.expr > quantile(gene.expr, probs = .75))] <- "2"

dat4plot <- data.frame(cbind(y = cancer_ccle[, "PF2341066"], grp = x1))
dat4plot[, 1] <- as.numeric(as.character(dat4plot[, 1]))

summary(glm(dat4plot[, 1] ~ as.numeric(dat4plot[, 2]))) -> sfit
p <- coef(sfit)[2, 4]

p3 <- ggplot(dat4plot, aes(y = y, x = grp)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
        shape = 21,
        position = position_jitter(width = 0.3),
        size = 0.5
    ) +
    guides(colour = FALSE) +
    theme(axis.title.x = element_blank(), legend.title = element_blank()) +
    labs(
        title = paste(
            "PF2341066, MET expression\np = ",
            format(p, digits = 3),
            sep = ""
        ),
        x = "", y = "Predicted Response to PF2341066"
    ) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = c("0" = "Q25", "1" = "Q25_75", "2" = "Q75")) +
    stat_summary(fun.data = give_n, geom = "text")

x2 <- rep("1", length(gene.expr))
x2[which(gene.expr < quantile(gene.expr, probs = .25))] <- "0"
x2[which(gene.expr > quantile(gene.expr, probs = .75))] <- "2"

dat4plot <- data.frame(cbind(y = cancer_ccle[, "PHA.665752"], grp = x2))
dat4plot[, 1] <- as.numeric(as.character(dat4plot[, 1]))

summary(glm(dat4plot[, 1] ~ as.numeric(dat4plot[, 2]))) -> sfit
p <- coef(sfit)[2, 4]

p4 <- ggplot(dat4plot, aes(y = y, x = grp)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
        shape = 21,
        position = position_jitter(width = 0.3),
        size = 0.5
    ) +
    guides(colour = FALSE) +
    theme(axis.title.x = element_blank(), legend.title = element_blank()) +
    labs(
        title = paste(
            "PHA665752, MET expression\np = ",
            format(p, digits = 3),
            sep = ""
        ),
        x = "", y = "Predicted Response to PHA665752"
    ) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = c("0" = "Q25", "1" = "Q25_75", "2" = "Q75")) +
    stat_summary(fun.data = give_n, geom = "text")

png("./E.met_CCLE.png", height = 1200, width = 1200)
multiplot(plotlist = list(p1, p3, p2, p4), layout = matrix(1:4, nrow = 2))
dev.off()

drug_ccle <- read.table(
    "../../Output/3/VAEN_GDSC.A.pred_TCGA.txt",
    header = T, as.is = T, sep = "\t"
)
cancer_types <- unique(drug_ccle[, 2])
sample_type <- substr(drug_ccle[, 1], 14, 15)
ss <- gsub("\\.", "-", drug_ccle[, 1])
drug_ccle[, 1] <- ss

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

    blca_ccle <- drug_ccle[
        which(drug_ccle[, 2] == cancer & sample_type == type_code),
    ]
    cancer_drug_ccle <- rbind(cancer_drug_ccle, blca_ccle)
}
drug_ccle <- cancer_drug_ccle

cancer_ccle <- drug_ccle[which(drug_ccle[, 2] == "LUAD"), ]

load("LUAD.MET_amp.ss.RData")
met_exon14_ii <- match(gsub("\\.", "-", ss), cancer_ccle[, 1])
match(LUAD.MET_amp.ss, cancer_ccle[, 1]) -> met_amp_ii
met_wt_ii <- setdiff(1:nrow(cancer_ccle), met_amp_ii)

dat4plot <- data.frame(cbind(
    y = cancer_ccle[, "Crizotinib"],
    grp = ifelse(cancer_ccle[, 1] %in% LUAD.MET_amp.ss,
    "MET Gain", "Other")
))
dat4plot[, 1] <- as.numeric(as.character(dat4plot[, 1]))

p <- t.test(dat4plot[, 1] ~ dat4plot[, 2])$p.value

p1 <- ggplot(dat4plot, aes(y = y, x = grp)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
        shape = 21,
        position = position_jitter(width = 0.3),
        size = 0.5
    ) +
    guides(colour = FALSE) +
    theme(axis.title.x = element_blank(), legend.title = element_blank()) +
    labs(
        title = paste(
            "Crizotinib, MET CNV\np = ",
            format(p, digits = 3),
            sep = ""
        ),
        x = "", y = "Predicted Response to Crizotinib"
    ) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    stat_summary(fun.data = give_n, geom = "text")

dat4plot <- data.frame(cbind(
    y = cancer_ccle[, "PHA.665752"],
    grp = ifelse(cancer_ccle[, 1] %in% LUAD.MET_amp.ss,
    "MET Gain", "Other")
))
dat4plot[, 1] <- as.numeric(as.character(dat4plot[, 1]))
p <- t.test(dat4plot[, 1] ~ dat4plot[, 2])$p.value

p2 <- ggplot(dat4plot, aes(y = y, x = grp)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
        shape = 21,
        position = position_jitter(width = 0.3),
        size = 0.5
    ) +
    guides(colour = FALSE) +
    theme(axis.title.x = element_blank(), legend.title = element_blank()) +
    labs(
        title = paste(
            "PHA665752, MET CNV\np = ",
            format(p, digits = 3),
            sep = ""
        ),
        x = "", y = "Predicted Response to PHA665752"
    ) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    stat_summary(fun.data = give_n, geom = "text")

dat4plot <- data.frame(cbind(
    y = cancer_ccle[, "Foretinib"],
    grp = ifelse(cancer_ccle[, 1] %in% LUAD.MET_amp.ss,
    "MET Gain", "Other")
))
dat4plot[, 1] <- as.numeric(as.character(dat4plot[, 1]))
p <- t.test(dat4plot[, 1] ~ dat4plot[, 2])$p.value

p3 <- ggplot(dat4plot, aes(y = y, x = grp)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
        shape = 21,
        position = position_jitter(width = 0.3),
        size = 0.5
    ) +
    guides(colour = FALSE) +
    theme(axis.title.x = element_blank(), legend.title = element_blank()) +
    labs(
        title = paste(
            "Foretinib, MET CNV\np = ",
            format(p, digits = 3),
            sep = ""
        ),
        x = "", y = "Predicted Response to Foretinib"
    ) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    stat_summary(fun.data = give_n, geom = "text")

load("LUAD.MET.expr.RData")

x2 <- rep("1", length(gene.expr))
x2[which(gene.expr < quantile(gene.expr, probs = .25))] <- "0"
x2[which(gene.expr > quantile(gene.expr, probs = .75))] <- "2"

dat4plot <- data.frame(cbind(y = cancer_ccle[, "Crizotinib"], grp = x2))
dat4plot[, 1] <- as.numeric(as.character(dat4plot[, 1]))

summary(glm(dat4plot[, 1] ~ as.numeric(dat4plot[, 2]))) -> sfit
p <- coef(sfit)[2, 4]

p4 <- ggplot(dat4plot, aes(y = y, x = grp)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
        shape = 21,
        position = position_jitter(width = 0.3),
        size = 0.5
    ) +
    guides(colour = FALSE) +
    theme(axis.title.x = element_blank(), legend.title = element_blank()) +
    labs(
        title = paste(
            "Crizotinib, MET expression\np = ",
            format(p, digits = 3),
            sep = ""
        ),
        x = "", y = "Predicted Response to Crizotinib"
    ) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = c("0" = "Q25", "1" = "Q25_75", "2" = "Q75")) +
    stat_summary(fun.data = give_n, geom = "text")


###

x2 <- rep("1", length(gene.expr))
x2[which(gene.expr < quantile(gene.expr, probs = .25))] <- "0"
x2[which(gene.expr > quantile(gene.expr, probs = .75))] <- "2"

dat4plot <- data.frame(cbind(y = cancer_ccle[, "PHA.665752"], grp = x2))
dat4plot[, 1] <- as.numeric(as.character(dat4plot[, 1]))

summary(glm(dat4plot[, 1] ~ as.numeric(dat4plot[, 2]))) -> sfit
p <- coef(sfit)[2, 4]

p5 <- ggplot(dat4plot, aes(y = y, x = grp)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
        shape = 21,
        position = position_jitter(width = 0.3),
        size = 0.5
    ) +
    guides(colour = FALSE) +
    theme(axis.title.x = element_blank(), legend.title = element_blank()) +
    labs(
        title = paste(
            "PHA665752, MET expression\np = ",
            format(p, digits = 3),
            sep = ""
        ),
        x = "", y = "Predicted Response to PHA665752"
    ) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = c("0" = "Q25", "1" = "Q25_75", "2" = "Q75")) +
    stat_summary(fun.data = give_n, geom = "text")

x2 <- rep("1", length(gene.expr))
x2[which(gene.expr < quantile(gene.expr, probs = .25))] <- "0"
x2[which(gene.expr > quantile(gene.expr, probs = .75))] <- "2"

dat4plot <- data.frame(cbind(y = cancer_ccle[, "Foretinib"], grp = x2))
dat4plot[, 1] <- as.numeric(as.character(dat4plot[, 1]))

summary(glm(dat4plot[, 1] ~ as.numeric(dat4plot[, 2]))) -> sfit
p <- coef(sfit)[2, 4]

p6 <- ggplot(dat4plot, aes(y = y, x = grp)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(
        shape = 21,
        position = position_jitter(width = 0.3),
        size = 0.5
    ) +
    guides(colour = FALSE) +
    theme(axis.title.x = element_blank(), legend.title = element_blank()) +
    labs(
        title = paste(
            "Foretinib, MET expression\np = ",
            format(p, digits = 3),
            sep = ""
        ),
        x = "", y = "Predicted Response to Foretinib"
    ) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = c("0" = "Q25", "1" = "Q25_75", "2" = "Q75")) +
    stat_summary(fun.data = give_n, geom = "text")


png("F.met_GDSC.png", height = 1200, width = 1800)
multiplot(
    plotlist = list(p1, p4, p2, p5, p3, p6),
    layout = matrix(1:6, nrow = 2)
)
dev.off()
