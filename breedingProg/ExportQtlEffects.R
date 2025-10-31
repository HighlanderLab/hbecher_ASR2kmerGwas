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
tmp = as.data.frame(SP$traits[[trait]]@epiEff)
colnames(tmp) = c("qtl1", "qtl2", "epiEff")
hist(tmp$epiEff)
dev.copy(png, file = "QTL1GaussEpiEff.png")
dev.off()
tmp$chr1 = QTL1$chr[tmp$qtl1]
tmp$loc1 = QTL1$loc[tmp$qtl1]
tmp$id1 = paste0(tmp$chr1, "_", tmp$loc1)
tmp$chr2 = QTL1$chr[tmp$qtl2]
tmp$loc2 = QTL1$loc[tmp$qtl2]
tmp$id2 = paste0(tmp$chr2, "_", tmp$loc2)
head(tmp)
QTL1 = merge(
  x = QTL1,
  y = tmp[, c("id1", "id2", "chr2", "loc2", "epiEff")],
  by.x = "id",
  by.y = "id1",
  all.x = TRUE
)
# str(QTL1); nrow(QTL1)
QTL1 = merge(
  x = QTL1,
  y = tmp[, c("id2", "id1", "chr1", "loc1", "epiEff")],
  by.x = "id",
  by.y = "id2",
  all.x = TRUE
)
# str(QTL1); nrow(QTL1)
# head(QTL1)
sel = is.na(QTL1$id2)
QTL1$id2[sel] = QTL1$id1[sel]
QTL1$chr2[sel] = QTL1$chr1[sel]
QTL1$loc2[sel] = QTL1$loc1[sel]
QTL1$epiEff.x[sel] = QTL1$epiEff.y[sel]
# str(QTL1); nrow(QTL1)
# View(QTL1)
QTL1$id1 = NULL
QTL1$chr1 = NULL
QTL1$loc1 = NULL
QTL1$epiEff.y = NULL
# head(QTL1)
colnames(QTL1) = c(
  "id",
  "chr",
  "loc",
  "addEff",
  "id2",
  "chr2",
  "loc2",
  "epiEff"
)
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

plot(QTL1$epiEff ~ QTL1$index, cex = 0.5)
abline(v = nQtl * c(1:nChr), col = "red", lty = 2)
dev.copy(png, file = "QTL1GaussEpiEffMap.png")
dev.off()

QTL1 <- QTL1[order(QTL1$chr2, QTL1$loc2), ]
QTL1$index2 <- 1:nrow(QTL1)
head(QTL1)
plot(QTL1$index2 ~ QTL1$index, cex = 10 * QTL1$epiEff)
abline(v = nQtl * c(1:nChr), col = "red", lty = 2)
abline(h = nQtl * c(1:nChr), col = "red", lty = 2)
dev.copy(png, file = "QTL1GaussEpiEffMap2.png")
dev.off()

# str(SP$traits[[2]])
trait = 2
QTL2 = data.frame(
  chr = rep(1:nChr, each = nQtl),
  loc = SP$traits[[trait]]@lociLoc,
  addEff = SP$traits[[trait]]@addEff
)
hist(QTL2$addEff)
dev.copy(png, file = "QTL2GammaAddEff.png")
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
  file = "QTL2Gamma.txt",
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
