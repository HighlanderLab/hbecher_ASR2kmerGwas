# Trait 1, no epistasis  ---------------------------------------

# str(SP$traits[[1]])
trait = 1
QTL1 = data.frame(
  chr = rep(1:nChr, each = nQtl),
  loc = SP$traits[[trait]]@lociLoc,
  addEff = SP$traits[[trait]]@addEff
)
hist(QTL1$addEff)
dev.copy(png, file = "QTL1GaussAddEff.png")
dev.off()
QTL1$id = paste0(QTL1$chr, "_", QTL1$loc)
#tmp = as.data.frame(SP$traits[[trait]]@epiEff)
#colnames(tmp) = c("qtl1", "qtl2", "epiEff")
#hist(tmp$epiEff)
#dev.copy(png, file = "QTL1GaussEpiEff.png")
#dev.off()
#tmp$chr1 = QTL1$chr[tmp$qtl1]
#tmp$loc1 = QTL1$loc[tmp$qtl1]
#tmp$id1 = paste0(tmp$chr1, "_", tmp$loc1)
#tmp$chr2 = QTL1$chr[tmp$qtl2]
#tmp$loc2 = QTL1$loc[tmp$qtl2]
#tmp$id2 = paste0(tmp$chr2, "_", tmp$loc2)
#head(tmp)
#QTL1 = merge(
#  x = QTL1,
#  y = tmp[, c("id1", "id2", "chr2", "loc2", "epiEff")],
#  by.x = "id",
#  by.y = "id1",
#  all.x = TRUE
#)
# str(QTL1); nrow(QTL1)
#QTL1 = merge(
#  x = QTL1,
#  y = tmp[, c("id2", "id1", "chr1", "loc1", "epiEff")],
#  by.x = "id",
#  by.y = "id2",
#  all.x = TRUE
#)
# str(QTL1); nrow(QTL1)
# head(QTL1)
#sel = is.na(QTL1$id2)
#QTL1$id2[sel] = QTL1$id1[sel]
#QTL1$chr2[sel] = QTL1$chr1[sel]
#QTL1$loc2[sel] = QTL1$loc1[sel]
#QTL1$epiEff.x[sel] = QTL1$epiEff.y[sel]
# str(QTL1); nrow(QTL1)
# View(QTL1)
#QTL1$id1 = NULL
#QTL1$chr1 = NULL
#QTL1$loc1 = NULL
#QTL1$epiEff.y = NULL
# head(QTL1)
#colnames(QTL1) = c(
#  "id",
#  "chr",
#  "loc",
#  "addEff",
#  "id2",
#  "chr2",
#  "loc2",
#  "epiEff"
#)
# head(QTL1)
QTL1 = QTL1[
  order(QTL1$chr, QTL1$loc),
]
# str(QTL1)
# head(QTL1)
write.table(
  x = QTL1,
  file = "QTL1Gauss.txt",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

QTL1$index <- 1:nrow(QTL1)
plot(QTL1$addEff ~ QTL1$index, cex = 0.5)
abline(v = nQtl * c(1:nChr), col = "red", lty = 2)
dev.copy(png, file = "QTL1GaussAddEffMap.png")
dev.off()

# plot(QTL1$epiEff ~ QTL1$index, cex = 0.5)
# abline(v = nQtl * c(1:nChr), col = "red", lty = 2)
# dev.copy(png, file = "QTL1GaussEpiEffMap.png")
# dev.off()

# QTL1 <- QTL1[order(QTL1$chr2, QTL1$loc2), ]
# QTL1$index2 <- 1:nrow(QTL1)
# head(QTL1)
# plot(QTL1$index2 ~ QTL1$index, cex = 10 * QTL1$epiEff)
# abline(v = nQtl * c(1:nChr), col = "red", lty = 2)
# abline(h = nQtl * c(1:nChr), col = "red", lty = 2)
# dev.copy(png, file = "QTL1GaussEpiEffMap2.png")
# dev.off()


# Trait 2, with epistasis ---------------------------------------
# str(SP$traits[[2]])
trait = 2
QTL2 = data.frame(
  chr = rep(1:nChr, each = nQtl),
  loc = SP$traits[[trait]]@lociLoc,
  addEff = SP$traits[[trait]]@addEff
)
hist(QTL2$addEff)
dev.copy(png, file = "QTL2GaussAddEff.png")
dev.off()
QTL2$id = paste0(QTL2$chr, "_", QTL2$loc)
tmp = as.data.frame(SP$traits[[trait]]@epiEff)
colnames(tmp) = c("qtl1", "qtl2", "epiEff")
hist(tmp$epiEff)
dev.copy(png, file = "QTL2GammaEpiEff.png")
dev.off()
tmp$chr1 = QTL2$chr[tmp$qtl1]
tmp$loc1 = QTL2$loc[tmp$qtl1]
tmp$id1 = paste0(tmp$chr1, "_", tmp$loc1)
tmp$chr2 = QTL2$chr[tmp$qtl2]
tmp$loc2 = QTL2$loc[tmp$qtl2]
tmp$id2 = paste0(tmp$chr2, "_", tmp$loc2)
head(tmp)
QTL2 = merge(
  x = QTL2,
  y = tmp[, c("id1", "id2", "chr2", "loc2", "epiEff")],
  by.x = "id",
  by.y = "id1",
  all.x = TRUE
)
# str(QTL2); nrow(QTL2)
QTL2 = merge(
  x = QTL2,
  y = tmp[, c("id2", "id1", "chr1", "loc1", "epiEff")],
  by.x = "id",
  by.y = "id2",
  all.x = TRUE
)
# str(QTL2); nrow(QTL2)
# head(QTL2)
sel = is.na(QTL2$id2)
QTL2$id2[sel] = QTL2$id1[sel]
QTL2$chr2[sel] = QTL2$chr1[sel]
QTL2$loc2[sel] = QTL2$loc1[sel]
QTL2$epiEff.x[sel] = QTL2$epiEff.y[sel]
# str(QTL2); nrow(QTL2)
# View(QTL2)
QTL2$id1 = NULL
QTL2$chr1 = NULL
QTL2$loc1 = NULL
QTL2$epiEff.y = NULL
# head(QTL2)
colnames(QTL2) = c(
  "id",
  "chr",
  "loc",
  "addEff",
  "id2",
  "chr2",
  "loc2",
  "epiEff"
)
# head(QTL2)
QTL2 = QTL2[
  order(QTL2$chr, QTL2$loc),
]
# str(QTL2)
# head(QTL2)
write.table(
  x = QTL2,
  file = "QTL2Gauss.txt",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

QTL2$index <- 1:nrow(QTL2)
plot(QTL2$addEff ~ QTL2$index, cex = 0.5)
abline(v = nQtl * c(1:nChr), col = "red", lty = 2)
dev.copy(png, file = "QTL2GammaAddEffMap.png")
dev.off()

plot(QTL2$epiEff ~ QTL2$index, cex = 0.5)
abline(v = nQtl * c(1:nChr), col = "red", lty = 2)
dev.copy(png, file = "QTL2GammaEpiEffMap.png")
dev.off()

QTL2 <- QTL2[order(QTL2$chr2, QTL2$loc2), ]
QTL2$index2 <- 1:nrow(QTL2)
head(QTL2)
plot(QTL2$index2 ~ QTL2$index, cex = 10 * QTL2$epiEff)
abline(v = nQtl * c(1:nChr), col = "red", lty = 2)
abline(h = nQtl * c(1:nChr), col = "red", lty = 2)
dev.copy(png, file = "QTL2GammaEpiEffMap2.png")
dev.off()


# the same for trait 3 and 4 can be done as needed
# Trait 3 ---------------------------------------
# str(SP$traits[[3]])
trait = 3
QTL3 = data.frame(
  chr = rep(1:nChr, each = nQtl),
  loc = SP$traits[[trait]]@lociLoc,
  addEff = SP$traits[[trait]]@addEff
)
hist(QTL3$addEff)
dev.copy(png, file = "QTL3GammaAddEff.png")
dev.off()
QTL3$id = paste0(QTL3$chr, "_", QTL3$loc)
#tmp = as.data.frame(SP$traits[[trait]]@epiEff)
#colnames(tmp) = c("qtl1", "qtl2", "epiEff")
#hist(tmp$epiEff)
#dev.copy(png, file = "QTL3GaussEpiEff.png")
#dev.off()
#tmp$chr1 = QTL3$chr[tmp$qtl1]
#tmp$loc1 = QTL3$loc[tmp$qtl1]
# tmp$id1 = paste0(tmp$chr1, "_", tmp$loc1)
# tmp$chr2 = QTL3$chr[tmp$qtl2]
# tmp$loc2 = QTL3$loc[tmp$qtl2]
# tmp$id2 = paste0(tmp$chr2, "_", tmp$loc2)
# head(tmp)
# QTL3 = merge(
#   x = QTL3,
#   y = tmp[, c("id1", "id2", "chr2", "loc2", "epiEff")],
#   by.x = "id",
#   by.y = "id1",
#   all.x = TRUE
# )
# # str(QTL3); nrow(QTL3)
# QTL3 = merge(
#   x = QTL3,
#   y = tmp[, c("id2", "id1", "chr1", "loc1", "epiEff")],
#   by.x = "id",
#   by.y = "id2",
#   all.x = TRUE
# )
# # str(QTL3); nrow(QTL3)
# # head(QTL3)
# sel = is.na(QTL3$id2)
# 
# QTL3$id2[sel] = QTL3$id1[sel]
# QTL3$chr2[sel] = QTL3$chr1[sel]
# QTL3$loc2[sel] = QTL3$loc1[sel]
# QTL3$epiEff.x[sel] = QTL3$epiEff.y[sel]
# # str(QTL3); nrow(QTL3)
# # View(QTL3)
# QTL3$id1 = NULL
# QTL3$chr1 = NULL
# QTL3$loc1 = NULL
# QTL3$epiEff.y = NULL
# # head(QTL3)
# colnames(QTL3) = c(
#   "id",
#   "chr",
#   "loc",
#   "addEff",
#   "id2",
#   "chr2",
#   "loc2",
#   "epiEff"
# )
# # head(QTL3)
# QTL3 = QTL3[
#   order(QTL3$chr, QTL3$loc),
# ]
# # str(QTL3)
# head(QTL3)
write.table(
  x = QTL3,
  file = "QTL3Gamma.txt",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
QTL3$index <- 1:nrow(QTL3)
plot(QTL3$addEff ~ QTL3$index, cex = 0.5)
abline(v = nQtl * c(1:nChr), col = "red", lty = 2)
dev.copy(png, file = "QTL3GammaAddEffMap.png")
dev.off()
# plot(QTL3$epiEff ~ QTL3$index, cex = 0.5)
#   
# abline(v = nQtl * c(1:nChr), col = "red", lty = 2)
# dev.copy(png, file = "QTL3GaussEpiEffMap.png")
# dev.off()
# QTL3 <- QTL3[order(QTL3$chr2, QTL3$loc2), ]
# QTL3$index2 <- 1:nrow(QTL3)
# head(QTL3)
# plot(QTL3$index2 ~ QTL3$index, cex = 10 * QTL3$epiEff)
# abline(v = nQtl * c(1:nChr), col = "red", lty = 2)
# abline(h = nQtl * c(1:nChr), col = "red", lty = 2)
# dev.copy(png, file = "QTL3GaussEpiEffMap2.png")
# dev.off()
# Trait 4 ---------------------------------------
# str(SP$traits[[4]])
trait = 4
QTL4 = data.frame(
  chr = rep(1:nChr, each = nQtl),
  loc = SP$traits[[trait]]@lociLoc,
  addEff = SP$traits[[trait]]@addEff
)
hist(QTL4$addEff)
dev.copy(png, file = "QTL4GammaAddEff.png")
dev.off()
QTL4$id = paste0(QTL4$chr, "_", QTL4$loc)
tmp = as.data.frame(SP$traits[[trait]]@epiEff)
colnames(tmp) = c("qtl1", "qtl2", "epiEff")
hist(tmp$epiEff)
dev.copy(png, file = "QTL4GammaEpiEff.png")
dev.off()
tmp$chr1 = QTL4$chr[tmp$qtl1]
tmp$loc1 = QTL4$loc[tmp$qtl1]
tmp$id1 = paste0(tmp$chr1, "_", tmp$loc1)
tmp$chr2 = QTL4$chr[tmp$qtl2]
tmp$loc2 = QTL4$loc[tmp$qtl2]
tmp$id2 = paste0(tmp$chr2, "_", tmp$loc2)
head(tmp)
QTL4 = merge(
  x = QTL4,
  y = tmp[, c("id1", "id2", "chr2", "loc2", "epiEff")],
  by.x = "id",
  by.y = "id1",
  all.x = TRUE
)
# str(QTL4); nrow(QTL4)
QTL4 = merge(
  x = QTL4,
  y = tmp[, c("id2", "id1", "chr1", "loc1", "epiEff")],
  by.x = "id",
  by.y = "id2",
  all.x = TRUE
)
# str(QTL4); nrow(QTL4)
# head(QTL4)
sel = is.na(QTL4$id2)
QTL4$id2[sel] = QTL4$id1[sel]
QTL4$chr2[sel] = QTL4$chr1[sel]
QTL4$loc2[sel] = QTL4$loc1[sel]
QTL4$epiEff.x[sel] = QTL4$epiEff.y[sel]
# str(QTL4); nrow(QTL4)
# View(QTL4)
QTL4$id1 = NULL
QTL4$chr1 = NULL
QTL4$loc1 = NULL
QTL4$epiEff.y = NULL
# head(QTL4)
colnames(QTL4) = c(
  "id",
  "chr",
  "loc",
  "addEff",
  "id2",
  "chr2",
  "loc2",
  "epiEff"
)
# head(QTL4)
QTL4 = QTL4[
  order(QTL4$chr, QTL4$loc),
]
# str(QTL4)
# head(QTL4)
write.table(
  x = QTL4,
  file = "QTL4Gamma.txt",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)
QTL4$index <- 1:nrow(QTL4)
plot(QTL4$addEff ~ QTL4$index, cex = 0.5)
abline(v = nQtl * c(1:nChr), col = "red", lty = 2)
dev.copy(png, file = "QTL4GammaAddEffMap.png")
dev.off()
plot(QTL4$epiEff ~ QTL4$index, cex = 0.5)
abline(v = nQtl * c(1:nChr), col = "red", lty = 2)
dev.copy(png, file = "QTL4GammaEpiEffMap.png")
dev.off()
QTL4 <- QTL4[order(QTL4$chr2, QTL4$loc2), ]
QTL4$index2 <- 1:nrow(QTL4)
head(QTL4)
plot(QTL4$index2 ~ QTL4$index, cex = 10 * QTL4$epiEff)
abline(v = nQtl * c(1:nChr), col = "red", lty = 2)
abline(h = nQtl * c(1:nChr), col = "red", lty = 2)
dev.copy(png, file = "QTL4GammaEpiEffMap2.png")
dev.off()
