v_stat=NULL;

p_model_genes=unique(data_out2[,1])

for (i in 1:length(p_model_genes)) 
{
	ind1=which(data_out2[,1]==p_model_genes[i])
	dat1=data_out2[ind1,]
	total_aff=length(unique(dat1[,2]))

	n_aff_down=which(dat1[,3]<0.99)
	n_aff_down_val=which(as.numeric(dat1[,6])>50)
	n_aff_down_val=intersect(n_aff_down,n_aff_down_val)

	n_aff_down=length(unique(dat1[n_aff_down,2]))
	n_aff_down_val=length(unique(dat1[n_aff_down_val,2]))

	n_aff_up=which(dat1[,3]=="Inf" | as.numeric(dat1[,3])>1)
	n_aff_up_val=which(as.numeric(dat1[,6])<  -50)
	n_aff_up_val=intersect(n_aff_down,n_aff_down_val)

	n_aff_up=length(unique(dat1[n_aff_up,2]))
	n_aff_up_val=length(unique(dat1[n_aff_up_val,2]))

	#d_s=((n_aff_down  - n_aff_up)/total_aff)
	v_s=paste(p_model_genes[i], total_aff, n_aff_down, n_aff_down_val, n_aff_up, n_aff_up_val, sep=" ")
	v_stat=append(v_stat, v_s)
}

write.table(v_stat, "V_STAT1_Sim50.txt")
