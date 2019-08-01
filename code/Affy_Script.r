# only need to run block 1-2
# block 0
#setwd("~/OneDrive - University of Nebraska-Lincoln/BIOC439_TermProject/New")
entrez570 = read.delim(file = "gpl570entrez.csv", sep = ",")
entrez96 = read.delim(file = "gpl96entrez.csv", sep = ",")
entrez97 = read.delim(file = "gpl97entrez.csv", sep = ",")
entrez4685 = read.delim(file = "gpl4685entrez.csv", sep = ",")
entrez8300 = read.delim(file = "gpl8300entrez.csv", sep = ",")

# block 1
#setwd("~/OneDrive - University of Nebraska-Lincoln/BIOC439_TermProject/New/Control/gpl570")
library(affy)
mydata = ReadAffy()
eset = mas5(mydata)
eset_PMA <- mas5calls(mydata)
y <- data.frame(exprs(eset), exprs(eset_PMA), assayDataElement(eset_PMA, "se.exprs"))
y <- y[,sort(names(y))]
y[,25] = rownames(y)
string = colnames(y)
string[25] = "ID"
colnames(y) = string
y_entrez = merge(y, entrez570, by = "ID", all = TRUE)
write.table(y, file="mydata_PMA.xls", quote=F, col.names = NA, sep="\t")
write.table(y_entrez, file="mydata_PMA_entrez.xls", quote=F, sep="\t", col.names = NA)

x = data.frame(exprs(eset))

x[,35] = rownames(x)
string = colnames(x)
string[35] = "ID"
colnames(x) = string
MERGE = merge(x, entrez96, by = "ID", all = TRUE)
write.table(MERGE, file = "MERGE.csv", sep = ",")

# block 2
# organize columns and remove expression data with no entrez ID
#NEW = read.delim(file = "MERGE.csv", sep = ",")
NEW = read.delim(file = "MERGE.csv", sep = ",")
entrez = NEW[,1]
NEW = NEW[,-1]

# block 3
# calculate z-score
NEW_z = scale(log2(NEW))
rownames(NEW_z) = entrez
write.table(NEW_z, file = "data_z.csv", sep = ",")

# boxplot
boxplot_labels = colnames(NEW_z)
boxplot = boxplot.matrix(NEW_z, use.cols = TRUE, outline = TRUE, names = boxplot_labels, main = "GSE2770 (96) Data", xlab = "Samples", ylab = "Normalized Expression Value", col=(c("purple")))

# remove duplicates after sorting highest to lowest average expression
NEW_z = read.delim(file = "data_z.csv", sep = ",")
NEW_z = NEW_z[!duplicated(NEW_z[,1]),]
write.table(NEW_z, file = "data_z_noduplicates.csv", sep = ",", row.names = FALSE)
