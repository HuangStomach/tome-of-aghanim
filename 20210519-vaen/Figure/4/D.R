library(ggplot2)

give_n <- function(x) {
    return(c(y = max(x), label = length(x)))
}

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

response <- read.delim("../../Data/Response/drug_response.txt", as.is = T)
response <- response[which(response$drug.name == "Paclitaxel"), ]

match(response[, 2], substr(drug_ccle[, 1], 1, 12)) -> ii
cbind(response, drug_ccle[ii, ]) -> new2
new2 <- new2[!is.na(ii), ]
dim(new2)

# tapply(new2$Paclitaxel, new2$response, mean)
grep("Disease", new2$response) -> ii
drug <- "Paclitaxel"
p <- t.test(new2[ii, drug], new2[-ii, drug])$p.value


dat4plot <- data.frame(new2[, c("Paclitaxel", "response")])
dat4plot[, 1] <- as.numeric(as.character(dat4plot[, 1]))
dat4plot[, 2] <- factor(dat4plot[, 2], levels = c(
    "Clinical Progressive Disease", "Stable Disease",
    "Partial Response", "Complete Response"
))

p3 <- ggplot(dat4plot, aes(x = response, y = Paclitaxel, fill = response)) +
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
            "CCLE-based Model\np = ",
            format(p, digits = 3),
            sep = ""
        ),
        x = "", y = "Predicted Response to Paclitaxel, TCGA-BRCA"
    ) +
    theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 0.7)
    ) +
    stat_summary(fun.data = give_n, geom = "text") +
    scale_x_discrete(
        labels = c(
            "Clinical Progressive Disease" = "Clinical\nProgressive\nDisease",
            "Stable Disease" = "Stable\nDisease",
            "Partial Response" = "Partial\nResponse",
            "Complete Response" = "Complete\nResponse"
        )
    )
ggsave(p3, file = "./D.Paclitaxel.response.png", width = 5, height = 5)
