freqCompute=function(mat=mat,refseq=refseq){
	freqR=list()
	refaa=translate_cm(s2c(refseq))
	allmutations=list()
	newMat=c()
	patternLabel=c()
	for(i in 1:dim(mat)[1]){
		specific_row=mat[i,]
		amino_set=setdiff(unique(specific_row[specific_row!=refaa[i]]),c("X","*"))
		for(am in amino_set){
			NY=sum(specific_row==am)
			newMat=rbind(newMat,as.numeric(specific_row==am))
			targetbp=s2c(refseq)[((i-1)*3+1):(i*3)]
			names_label=paste(translate_cm(targetbp),am,sep=as.character(i))
			patternLabel=append(patternLabel,names_label)
			freqd=round(NY/dim(mat)[2],digits=3)
			#names(freqd)=names_label
			allmutations[[names_label]]=freqd
		}
	}
	rownames(newMat)=patternLabel
	freqR$newMat=newMat
	freqR$allmutations=allmutations
	return(freqR)
}
