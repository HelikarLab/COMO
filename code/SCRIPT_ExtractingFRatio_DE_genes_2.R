######################################################################################################################################################
###### This script will work on the output of previous (that should be run one time only) #######################################################################################################################################
# Compare with perturbed profile
# Read matlab generated file into R (output of Knock_out_simulation.m) 
# To run below chunks directly use the CMap_similarity_data_Entrez.RData

MN_gene_pairs=read.table("Gene_Pairs_Inhi_Fratio_DOWN.txt", sep=",",header=T)

# Extract connectivtyMap_score for MN_gene_Pairs
data_out2={};
for( i in 1:nrow(MN_gene_pairs)){
c1_c4=which(MN_gene_pairs[i,1]==data_out[,5])
c2_c5=which(MN_gene_pairs[i,2]==data_out[,6])
matched_ind=intersect(c1_c4,c2_c5)

if(length(matched_ind>0)) {
matched_d = data_out[matched_ind,]
 if(class(matched_d)=="matrix"){
	matched_kd=which(matched_d[,4]=="kd")
	matched_d=matched_d[matched_kd,]
	comb_row=cbind(as.matrix(MN_gene_pairs[i,]),t(as.data.frame(matched_d)))
				}else {
			comb_row=cbind(as.matrix(MN_gene_pairs[i,]),t(as.data.frame(matched_d)))
			
					}


} else {
comb_row=cbind(as.matrix(MN_gene_pairs[i,]),t(as.data.frame(rep("NA",6))))
}

data_out2=rbind(data_out2,comb_row)
}


write.table(data_out2, "MN_gene_pairs_CMap_ScoresDOWN.txt")

save.image("MN_gene_pairs_CMap_ScoresDOWN.RData")
