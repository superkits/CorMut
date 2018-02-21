kaksAA <-
function(seq_formated){
	seq_char=tolower(as.character(seq_formated))
	refseq=seq_char[1]
	seq_list=seq_char[2:length(seq_char)]
	mat=sapply(as.list(seq_list),function(x)translate_cm(s2c(x)))
	refseq=tolower(refseq)
	refaa=translate_cm(s2c(refseq))
	fvft=fvft(seq_list,refseq=refseq)
	#final=list()
	final_result=list()
	LOD_result=list()
	q_result=list()
	KA_result=list()
	KS_result=list()
	for(i in 1:dim(mat)[1]){
		NS=0
		specific_row=mat[i,]
		for(k in 1:dim(mat)[2]){
			codon_nu=substr(seq_list[[k]],(i-1)*3+1,i*3)
			#cat(codon_nu!=substr(refseq,(i-1)*3+1,i*3))
			if(codon_nu!=substr(refseq,(i-1)*3+1,i*3)&specific_row[k]==refaa[i]){
				#cat(specific_row[k],refaa[i],specific_row[k]==refaa[i],"\n")
				NS=NS+1
			}
		}
		KS_result=append(KS_result,list(NS))
		#NY=sum(specific_row!=refaa[i])
		amino_set=setdiff(unique(specific_row[specific_row!=refaa[i]]),c("X","*"))
		if(NS==0)NS=1
		#KaKs list of the specific position
		KaKs_slist=c()
		q_slist=c()
		LOD_slist=c()
		KA_slist=c()
		if(length(amino_set)==0){
			KaKs_slist=0;
			LOD_slist=0;
			q_slist=0;
			KA_slist=0;
		}
		for(am in amino_set){
			#NY=sum(specific_row!=refaa[i])
			#NY should not include "X" amino, i.e. gap
			NY=sum(((specific_row!=refaa[i])+(specific_row!="X"))==2)
			NYs=sum(specific_row==am)
			#compute NtNv
			targetbp=s2c(refseq)[((i-1)*3+1):(i*3)]
			names_label=paste(translate_cm(targetbp),am,sep=as.character(i))
			#print(targetbp)
			NtNv=Ntnv_aa(targetbp,aa_m=am)
			#ensure ka/ka has meaning NO zero
			model_fenzi=NtNv$Nat*fvft$ft+NtNv$Nav*fvft$fv
			model_fenmu=NtNv$Nst*fvft$ft+NtNv$Nsv*fvft$fv
			if(model_fenmu==0|model_fenzi==0){
				random_model=((NtNv$Nat+1)*fvft$ft+(NtNv$Nav+1)*fvft$fv)/((NtNv$Nst+1)*fvft$ft+(NtNv$Nsv+1)*fvft$fv)
			}else{
				random_model=model_fenzi/model_fenmu
			}
			#print(c(NtNv$Nst,fvft$ft,NtNv$Nsv,fvft$fv,random_model,NYs/NS))
			#if((random_model==Inf)|(random_model=="NaN")|(random_model==0))random_model=1
			KaKs=(NYs/NS)/random_model
			names(KaKs)=names_label
			KaKs_slist=append(KaKs_slist,KaKs)
			#KA alone
			names(NYs)=names_label
			KA_slist=append(KA_slist,NYs)
			q_fenzi=(NtNv$Nat*fvft$ft+NtNv$Nav*fvft$fv)
			##########################################################Warning
			#if(q_fenzi==0)q_fenzi=0.0001
			q=q_fenzi/(3*fvft$ft+6*fvft$fv)
			names(q)=names_label
			q_slist=append(q_slist,q)
			LOD=0
			Nnum=NS+NY
			for(ii in NYs:Nnum){
				#per_item=choose(Nnum,ii)*(q^ii)*((1-q)^(Nnum-ii))
				per_item=dbinom(ii,Nnum,q)
				#if(per_item=="NaN") print(per_item)
				#if(per_item!="NaN"){
				LOD=per_item+LOD
				#}
				#print(choose(dim(mat)[2],i)*(q^i)*((1-q)^(dim(mat)[2]-i)))
			}
			LOD=-log10(LOD)
			names(LOD)=names_label
			LOD_slist=append(LOD_slist,LOD)
		}
			LOD_result=append(LOD_result,list(LOD_slist))
			final_result=append(final_result,list(KaKs_slist))
			q_result=append(q_result,list(q_slist))
			KA_result=append(KA_result,list(KA_slist))
	}
	final=new("kaksAA",seq_num=length(seq_list),kaks=final_result,lod=LOD_result,q=q_result,ka=KA_result,ks=KS_result)
	return(final)
}
