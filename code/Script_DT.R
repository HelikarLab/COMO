# Read connectivity map file that has withdrawn drugs removed
con_map_dt=read.csv("Conmap_DT.csv",header=T)
# Read Symbol vs Entrez IDs file retreived form orthoretreiver
Sym2Entrez=read.table("Symbol_Entrez_Conmap_data.txt",header=T)
# Split Drug target symbols

Target_char=as.character(con_map_dt$Target)
spl_dt=list();
for (i in 1:length(Target_char)) #
{
x=strsplit(Target_char[i],split=", ")
spl_dt=append(spl_dt,x)
}

#### Drugs for each target SYMBOL
SymChar=as.character(Sym2Entrez$Symbol)

DT_Symb=NULL;
for(i in 1:length(spl_dt)) #length(spl_dt)
{
	if(length(spl_dt[[i]])>0)
	{
		for(j in 1:length(spl_dt[[i]]))
			{
			x=which((spl_dt[[i]][j]==SymChar))
			DT_Symb1=cbind(SymChar[x], con_map_dt[i,1:4])
			DT_Symb=rbind(DT_Symb,DT_Symb1)
			}	
	}
}

#### Addind entrez for each row of symbol
DT_Entrez=NULL;
#DT_symbChar=as.character(DT_Symb[,1])
for(i in 1:nrow(DT_Symb)) #
{
	x=which(as.character(DT_Symb[i,1])== as.character(Sym2Entrez[,1]))
	DT_Entrez1= Sym2Entrez[x,2]
	DT_Entrez = append(DT_Entrez,DT_Entrez1)
}

#Combine DT_Symb (Drug info per symbol) with DT_Entrez (Entrez for each Symbol DT)

DT_Symb_Entrez=cbind(DT_Entrez,DT_Symb)
write.table(DT_Symb_Entrez, "DRUGS_DATA_ENTREZWISE.txt", sep="\t")

##end of the script







