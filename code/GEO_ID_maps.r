#!/usr/bin/env Rscript
prjwd <- '~/Dropbox/Helikar/pipelines/data/'
setwd(prjwd)
library(GEOquery)
gpl96 <- getGEO('GPL96', destdir=".")
gpl97 <- getGEO('GPL97', destdir=".")
gpl570 <- getGEO('GPL570', destdir=".")
gpl4685 <- getGEO('GPL4685', destdir=".")
gpl8300 <- getGEO('GPL8300', destdir=".")

entrez_id_96 = gpl96@dataTable@table$ENTREZ_GENE_ID
entrez_id_97 = gpl97@dataTable@table$ENTREZ_GENE_ID
entrez_id_570 = gpl570@dataTable@table$ENTREZ_GENE_ID
#entrez_id_4685 = gpl4685@dataTable@table$ENTREZ_GENE_ID
entrez_id_8300 = gpl8300@dataTable@table$ENTREZ_GENE_ID
id_96 = gpl96@dataTable@table$ID
id_97 = gpl97@dataTable@table$ID
id_570 = gpl570@dataTable@table$ID
#id_4685 = gpl4685@dataTable@table$ID
id_8300 = gpl8300@dataTable@table$ID

entre96 = data.frame(id_96,entrez_id_96)
names(entre96) <- c('ID', 'ENTREZ_GENE_ID')
entre97 = data.frame(id_97,entrez_id_97)
names(entre97) <- c('ID', 'ENTREZ_GENE_ID')
entre570 = data.frame(id_570,entrez_id_570)
names(entre570) <- c('ID', 'ENTREZ_GENE_ID')
#entre4685 = data.frame(id_4685,entrez_id_4685)
#names(entre4685) <- c('ID', 'ENTREZ_GENE_ID')
entre8300 = data.frame(id_8300,entrez_id_8300)
names(entre8300) <- c('ID', 'ENTREZ_GENE_ID')

write.table(entre96, file = "gpl96entrez.csv", row.names = FALSE, quote = FALSE, sep = ",")
write.table(entre97, file = "gpl97entrez.csv", row.names = FALSE, quote = FALSE, sep = ",")
write.table(entre570, file = "gpl570entrez.csv", row.names = FALSE, quote = FALSE, sep = ",")
#write.table(entre4685, file = "gpl4685entrez.csv", row.names = FALSE, quote = FALSE, sep = ",")
write.table(entre8300, file = "gpl8300entrez.csv", row.names = FALSE, quote = FALSE, sep = ",")

#entrez570 = write.delim(file = "gpl570entrez.csv", sep = ",")
#entrez96 = write.delim(file = "gpl96entrez.csv", sep = ",")
#entrez97 = read.delim(file = "gpl97entrez.csv", sep = ",")
#entrez4685 = read.delim(file = "gpl4685entrez.csv", sep = ",")
#entrez8300 = read.delim(file = "gpl8300entrez.csv", sep = ",")
