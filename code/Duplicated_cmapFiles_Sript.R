tx1=NULL;
for (i in 1:190) {

tx=tmp1[[i]][1]
tx1= append(tx1,tx)
}


which(duplicated(tx1)==TRUE)
