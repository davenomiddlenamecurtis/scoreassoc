# setwd("/Users/dave_000/OneDrive/sharedseq/SSS2/scores/SSS2.vrare.1/combinedscores")
setwd("d:/dave/OneDrive/sharedseq/SSS2/scores/SSS2.vrare.1/combinedscores")
rm(list=ls())   
paramsOut=dataframe(stringsAsFactors=FALSE)
paramsOut=read.table("paramsOut.txt",header=TRUE)
params.cropped=read.table("params.cropped.txt",header=TRUE)
Table1=cbind(paramsOut[2],paramsOut[1],paramsOut[5:6],params.cropped[1],params.cropped[5:6])
write.table(Table1,"Table1.txt",quote=FALSE,sep="\t",row.names=FALSE)

# Table2=cbind(sub(".txt","",paramsOut[2]))
Table2=cbind(paramsOut[2])
Table2[,1]=as.character(Table2[,1])
for (i in 1:nrow(Table2)) {
	Table2[i,1]=sub(".txt","",Table2[i,1])
}

for (i in 0:4) {
filename=sprintf("paramsOut.trained.fifth.%d.txt",i)
tab=read.table(filename,header=TRUE)
colname=sprintf("Weight.%d",i)
Table2[colname]=tab[1]
colname=sprintf("Min.%d",i)
Table2[colname]=tab[5]
colname=sprintf("Max.%d",i)
Table2[colname]=tab[6]
}
write.table(Table2,"Table2.txt",quote=FALSE,sep="\t",row.names=FALSE)

remove(Table3)
Table3=cbind(paramsOut[2])
for (i in 0:4) {
filename=sprintf("paramsOut.trained.minimal.fifth.%d.txt",i)
tab=read.table(filename,header=TRUE)
colname=sprintf("Weight.%d",i)
Table3[colname]=tab[1]
colname=sprintf("Min.%d",i)
Table3[colname]=tab[5]
colname=sprintf("Max.%d",i)
Table3[colname]=tab[6]
}
write.table(Table3,"Table3.txt",quote=FALSE,sep="\t",row.names=FALSE)
Table3

formTable=cbind(Table2[1]) 
for (i in 1:nrow(formTable)) {
	if (paramsOut[i,1]!=0)
		formTable[i,2]=sprintf("%.1f (%.1f,%.1f)",paramsOut[i,1],paramsOut[i,5],paramsOut[i,6])
	else
		formTable[i,2]="0"
	if (params.cropped[i,1]!=0)
		formTable[i,3]=sprintf("%.1f (%.1f,%.1f)",params.cropped[i,1],params.cropped[i,5],params.cropped[i,6])
	else
		formTable[i,3]="0"
}
colnames(formTable)=c("Name","Fitted weights (range)","Minimal fitted weights (range)")
i=nrow(formTable)+1
formTable[i,1]="tstat"
filename="ttest.out"
tfile=file(filename)
line=readLines(tfile,4)
words=strsplit(line[4],"\\s+")[[1]]
formTable[i,2]=sprintf("%.1f",as.numeric(words[5]))
filename="cropped.ttest.out"
tfile=file(filename)
line=readLines(tfile,4)
words=strsplit(line[4],"\\s+")[[1]]
formTable[i,3]=sprintf("%.1f",as.numeric(words[5]))
formTable
write.table(formTable,"wordTable3.txt",quote=FALSE,sep="\t",row.names=FALSE)

formfifthTable=cbind(Table2[1]) 
for (j in 0:4) {
filename=sprintf("paramsOut.trained.fifth.%d.txt",j)
tab=read.table(filename,header=TRUE)
colname=sprintf("Fitted weights %d (range)",j+1)
for (i in 1:nrow(formTable)) {
	if (tab[i,1]!=0)
		formfifthTable[i,colname]=sprintf("%.1f (%.1f,%.1f)",tab[i,1],tab[i,5],tab[i,6])
	else
		formfifthTable[i,colname]="0"
}
}
i=nrow(formfifthTable)+1
formfifthTable[i,1]="tstat"
for (j in 0:4) {
filename=sprintf("ttest.all.from.%d.out",j)
tfile=file(filename)
line=readLines(tfile,4)
words=strsplit(line[4],"\\s+")[[1]]
words[4]
formfifthTable[i,j+2]=sprintf("%.1f",as.numeric(words[5]))
}
formfifthTable
write.table(formfifthTable,"wordTable5.txt",quote=FALSE,sep="\t",row.names=FALSE)

formfifthTable
write.table(formfifthTable,"wordTable4.txt",quote=FALSE,sep="\t",row.names=FALSE)
remove(formfifthminimalTable)
formfifthminimalTable=cbind(Table2[1]) 
for (j in 0:4) {
filename=sprintf("paramsOut.trained.minimal.fifth.%d.txt",j)
tab=read.table(filename,header=TRUE)
colname=sprintf("Minimal weights %d (range)",j+1)
for (i in 1:nrow(formfifthminimalTable)) {
	if (tab[i,1]!=0)
		formfifthminimalTable[i,colname]=sprintf("%.1f (%.1f,%.1f)",tab[i,1],tab[i,5],tab[i,6])
	else
		formfifthminimalTable[i,colname]="0"
}
}
formfifthminimalTable
i=nrow(formfifthminimalTable)+1
formfifthminimalTable[i,1]="tstat"
for (j in 0:4) {
filename=sprintf("ttest.all.from.minimal.%d.out",j)
tfile=file(filename)
line=readLines(tfile,4)
words=strsplit(line[4],"\\s+")[[1]]
words[4]
formfifthminimalTable[i,j+2]=sprintf("%.1f",as.numeric(words[5]))
}
formfifthminimalTable
write.table(formfifthminimalTable,"wordTable5.txt",quote=FALSE,sep="\t",row.names=FALSE)



