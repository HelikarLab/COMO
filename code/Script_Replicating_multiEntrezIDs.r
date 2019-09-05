##########################################################################

###Read data from the text file (Change file path here). Use "read.csv" if data is in CSV format##

data=read.table("C:/Users/BHANWAR/Desktop/TEST12.txt",header=TRUE)

###Split Entrz IDs based on the the underscore (as given in files uploaded by Bailee) ##

nsplits=strsplit(as.character(data[,1]),"_")

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


write.table(data5, "data_gse_replicated_entrez.txt", sep="\t")




