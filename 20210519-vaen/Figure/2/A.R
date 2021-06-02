source("../../Lib/unfactor.R")

drugs_match <- read.table("../../Data/Match/durgs_match_2.txt",
    as.is = T, sep = "\t")
load("../../Output/4/dr.CCLE.A.models.RData")
drugs <- names(dr_ccle_models)
x <- c()

for (k in 1:length(drugs)) {
    drug <- drugs[k]
    res_list <- dr_ccle_models[[drug]]
    Ys <- res_list$Ys
    which(Ys[, 1] != -9) -> ii
    Ys <- Ys[ii, ]
    recall <- cor(Ys[, 1], Ys[, 2])

    x <- rbind(x, c(drug, recall))
}
colnames(x) <- c("Drug", "PCC")


load("../../Output/4/dr.GDSC.A.models.RData")
drugs <- names(dr_gdsc_models)
y <- c()
for (k in 1:length(drugs)) {
    drug <- drugs[k]
    res_list <- dr_gdsc_models[[drug]]
    Ys <- res_list$Ys
    which(Ys[, 1] != -9) -> ii
    Ys <- Ys[ii, ]
    recall <- cor(Ys[, 1], Ys[, 2])

    y <- rbind(y, c(drug, recall))
}
colnames(y) <- c("Drug", "PCC")

write.table(x, file = "./Figure2A.CCLE.txt")
write.table(y, file = "./Figure2A.GDSC.txt")

pdf("./A.PCC.pdf", width = 8, height = 4)
match(drugs_match[, 1], x[, 1]) -> sx_ii
match(drugs_match[, 3], y[, 1]) -> sy_ii

a1 <- length(x[-sx_ii, "PCC"])
b1 <- length(y[-sy_ii, "PCC"])
a2 <- length(sx_ii)
b2 <- length(sy_ii)

dat_plot <- rbind(
    cbind(x[-sx_ii, c("Drug", "PCC")], grp = "1"),
    cbind(y[-sy_ii, c("Drug", "PCC")], grp = "2"),
    cbind(x[sx_ii, c("Drug", "PCC")], grp = "3"),
    cbind(y[sy_ii, c("Drug", "PCC")], grp = "4")
)
dat_plot <- unfactor(as.data.frame(dat_plot))
dat_plot <- dat_plot[order(dat_plot[, 3], dat_plot[, 2]), ]

plot(0, 0, xlim = c(-1, a1 + b1 + a2 + b2 + 10), ylim = c(0, 1.05), 
    col = "white", ylab = "", xlab = "")
mtext("Compound index", 1, line = 2)
mtext("In-sample PCC", 2, line = 2)

rect(-1, 0, a1 + 0.5, 1, col = "lightyellow", border = F)
rect(a1 + 0.5, 0, a1 + b1 + 0.5, 1, col = "lightcyan", border = F)
rect(a1 + b1 + 0.5, 0, a1 + b1 + a2 + 0.5, 1, col = "lightgreen", border = F)
rect(a1 + b1 + a2 + 0.5, 0, a1 + b1 + a2 + b2 + 1, 1, 
    col = "lightblue", border = F)

points(dat_plot[, 2], pch = 19, col = "grey", cex = 1)
points(dat_plot[, 2], col = "black", cex = 1)
segments(-1, 0.5, 276, 0.5)

which(dat_plot[, 2] < 0.5 & dat_plot[, 3] == "1") -> ii
for (k in ii) {
    text(k, as.numeric(dat_plot[k, 2]), dat_plot[k, 1], pos = 4, cex = .8)
}

which(dat_plot[, 2] < 0.5 & dat_plot[, 3] == "3") -> ii
text(ii, as.numeric(dat_plot[ii, 2]), dat_plot[ii, 1], pos = 4, cex = .8)

which(dat_plot[, 2] < 0.5 & dat_plot[, 3] == "4") -> ii
text(ii, as.numeric(dat_plot[ii, 2]), dat_plot[ii, 1], pos = 4, cex = .8)

dev.off()
