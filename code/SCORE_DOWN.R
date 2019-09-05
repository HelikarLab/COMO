d_score1=NULL;

p_model_genes=unique(data_out2[,1])

for (i in 1:length(p_model_genes)) 
{
	ind1=which(data_out2[,1]==p_model_genes[i])
	dat1=data_out2[ind1,]
	total_aff=length(unique(dat1[,2]))
	n_aff_down=which(dat1[,3]<0.99)
	n_aff_down=length(unique(dat1[n_aff_down,2]))
	n_aff_up=which(dat1[,3]=="Inf" | as.numeric(dat1[,3])>1)
	n_aff_up=length(unique(dat1[n_aff_up,2]))
	#d_s=((n_aff_up - n_aff_down)/total_aff)
	#d_s=(n_aff_up - n_aff_down)
	d_s=(n_aff_up/total_aff)
	d_s=paste(p_model_genes[i], d_s)
	d_score1=append(d_score1, d_s)
}

write.table(d_score1, "d_scoreDOWN.txt")



############## If above wont run ##################

d_score1=NULL;

p_model_genes=unique(data_out2[,1])

for (i in 1:length(p_model_genes)) 
{
	ind1=which(data_out2[,1]==p_model_genes[i])
	dat1=data_out2[ind1,]

	if (class(dat1) == "character")
	{
		dat1=as.matrix(dat1)
		dat1=t(dat1)
	} 
	else {}
	total_aff=length(unique(dat1[,2]))
	n_aff_down=which(dat1[,3]<0.99)
	n_aff_down=length(unique(dat1[n_aff_down,2]))
	n_aff_up=which(dat1[,3]=="Inf" | as.numeric(dat1[,3])>1)
	n_aff_up=length(unique(dat1[n_aff_up,2]))
	#d_s=((n_aff_up - n_aff_down)/total_aff)
	#d_s=(n_aff_up - n_aff_down)
	d_s=(n_aff_up/total_aff)
	d_s=paste(p_model_genes[i], d_s)
	d_score1=append(d_score1, d_s)
}

write.table(d_score1, "d_scoreDOWN.txt")
