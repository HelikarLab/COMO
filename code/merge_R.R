## Merge two data files together.
# Create text file for datasets with HGNC as first column followed by 1's and 0's.
> path1 = "~/Documents/datafile1.txt"
> data1=read.table(path1,header=T)
> path2 = "~/Documents/datafile2.txt"
> data2=read.table(path2,header=T)
> mydata1=data1
> mydata2=data2
> mydataset=merge(mydata1, mydata2, by="HGNC",all=TRUE)
> write.table(mydataset,"~/Desktop/merged_file.txt",F)

# Repeat. Use the "merged_file.txt" as path1 for subsequent merges, dataset to be incorporated as path2.