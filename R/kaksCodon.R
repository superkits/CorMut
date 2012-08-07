kaksCodon <-
function(seq_formated){
	seq_char=tolower(as.character(seq_formated))
	refseq=seq_char[1]
	seq_list=seq_char[2:length(seq_char)]
	mat=sapply(as.list(seq_list),function(x)translate_cm(s2c(x)))
	refseq=tolower(refseq)
	refaa=translate_cm(s2c(refseq))
	fvft=fvft(seq_list,refseq=refseq)
	#final=list()
	kaks_result=c()
	LOD_result=c()
	q_result=c()
	KA_result=c()
	KS_result=c()
	for(i in 1:dim(mat)[1]){
		NS=0
		specific_row=mat[i,]
		NY=sum(specific_row!=refaa[i])
		KA_result=append(KA_result,NY)
		for(k in 1:dim(mat)[2]){
			codon_nu=substr(seq_list[[k]],(i-1)*3+1,i*3)
			#cat(codon_nu,substr(refseq,(i-1)*3+1,i*3),"\n")
			#cat("\n")
			#cat(codon_nu!=substr(refseq,(i-1)*3+1,i*3))
			if(codon_nu!=substr(refseq,(i-1)*3+1,i*3)&specific_row[k]==refaa[i]){
				#cat(specific_row[k],refaa[i],specific_row[k]==refaa[i],"\n")
				NS=NS+1
			}
		}
		KS_result=append(KS_result,NS)
		#compute NtNv
		targetbp=s2c(refseq)[((i-1)*3+1):(i*3)]
		#print(targetbp)
		NtNv=Ntnv(targetbp)
		random_model=(NtNv$Nat*fvft$ft+NtNv$Nav*fvft$fv)/(NtNv$Nst*fvft$ft+NtNv$Nsv*fvft$fv)
		if((random_model==Inf)|(random_model=="NaN")|(random_model==0))random_model=1
		if(NY==0&NS==0){
			KaKs=1
		}else if(NY!=0&NS==0){
			KaKs=(NY)/random_model
		}else{
			KaKs=(NY/NS)/random_model
		}
		kaks_result=append(kaks_result,KaKs)
		q=(NtNv$Nat*fvft$ft+NtNv$Nav*fvft$fv)/(3*fvft$ft+6*fvft$fv)
		q_result=append(q_result,q)
		LOD=0
		#print(NY)
		Nnum=NS+NY
		for(ii in NY:Nnum){
			LOD=choose(Nnum,ii)*(q^ii)*((1-q)^(Nnum-ii))+LOD
			#print(choose(dim(mat)[2],i)*(q^i)*((1-q)^(dim(mat)[2]-i)))
		}
		#print(LOD)
		LOD_result=append(LOD_result,-log10(LOD))
	}
	final=new("kaksCodon",seq_num=length(seq_list),kaks=kaks_result,lod=LOD_result,q=q_result,ka=KA_result,ks=KS_result)
	return(final)
}
