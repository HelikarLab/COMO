############################ This script is to process similarity scores for multiple genes downloaded from CMap database ###################################################
# read connectivityMap touchstone profiles and convert them into single data matrix 
fnames=dir()
data_comb={};
for (i in 1:length(fnames))
{
	pert_gene_name=strsplit(fnames[i],'_')
	pert_gene_name=pert_gene_name[[1]][1]

	data=read.delim(fnames[i],header=T)

	pert_gene_col=rep(pert_gene_name,length(data$Name))
	data_mat=cbind(pert_gene_col,as.character(data$Name), data$Score, as.character(data$Type))

	data_comb=rbind(data_comb, data_mat)
}

# Change gene names to Entrez_IDs
library('org.Hs.eg.db')

Entrez_col1=mapIds(org.Hs.eg.db, data_comb[,1], 'ENTREZID', 'SYMBOL')
Entrez_col2=mapIds(org.Hs.eg.db, data_comb[,2], 'ENTREZID', 'SYMBOL')



data_out=cbind(data_comb,Entrez_col1,Entrez_col2)
#save(data_comb, data_out, file="CMap_similarity_data_Entrez.RData")



