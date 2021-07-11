latent <- read.table("./Output/1/1.CCLE_latent.tsv", as.is = T, header = T)
# tsne(latent) -> tpc
sapply(rownames(latent), function(x) {
    strsplit(x, split = "\\.")[[1]][1] -> u
    strsplit(u, split = "_")[[1]] -> v
    v <- v[-1]
    paste(v, collapse = "_")
}) -> tt
tissues <- sort(unique(tt))
print(tissues)
