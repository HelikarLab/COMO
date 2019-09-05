

############################################################# 
# Remove duplicates before merging and after merging ###########
######################################################################


test11=read.csv(dir()[6],header=T)

sum1=rowSums(test11[,-1])

test13=cbind(test11,sum1)
test14=test13[rev(order(test13$sum1)),]


#test15=unique(test14, by = "ENTREZ_GENE_ID")

test16=test14[!duplicated(test14$ENTREZ_GENE_ID),]

#####################################################################
# Use other script for replicating Entrez IDs 

data= test16

nsplits=strsplit(as.character(data[,1])," /// ")

### Replicate the multiEntrez rows and Add new names in the first column##

data3=NULL;
for(i in 1:nrow(data))
{
	if(length(nsplits[[i]])>1)
	{
		data2=data[rep(rownames(data)[i], length(nsplits[[i]])), ] 
		#print(data2)
		data2[,1]=as.character(data2[,1])
		for(j in 1:nrow(data2))
		{
			data2[j,1]=nsplits[[i]][j]
		}
		data3=rbind(data3,data2)
	}
}


## Remove rows with multiple Entrez IDs from the original data frame (create reduces dataframe)##

tt1=NULL;
for(i in 1:nrow(data))
{
	tt=length(nsplits[[i]])
	tt1=append(tt1,tt)
}

tt2=which(tt1>1)


data4=data[-tt2,]

## Create final dataframe (Combine replicated dataframe with reduced dataframe)## 

data5= rbind(data4,data3)


##################Repeat again after duplicating Entrez IDs########

test11=data5

sum1=rowSums(test11[,-1])

test13=cbind(test11,sum1)
test14=test13[rev(order(test13$sum1)),]


#test15=unique(test14, by = "ENTREZ_GENE_ID")

test16=test14[!duplicated(test14$ENTREZ_GENE_ID),]


######################################################################
# Use merge function #################################################

merged1=merge(Ex_gse22886_96, Ex_22886_97, by="ENTREZ_GENE_ID", all=TRUE)
merged2=merge(merged1, Ex_2770_8300, by="ENTREZ_GENE_ID", all=TRUE)
merged3=merge(merged2, Ex_2770_96, by="ENTREZ_GENE_ID", all=TRUE)
merged4=merge(merged3, Ex_2770_97, by="ENTREZ_GENE_ID", all=TRUE)
merged5=merge(merged4, Ex_436769_570, by="ENTREZ_GENE_ID", all=TRUE)



write.table(merged5, "Th2_Merged_Affy.txt", sep="\t")


