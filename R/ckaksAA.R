ckaksAA <-
function(seq_formated,kaks=T,lod_cut=2,setPosition=c()){
	#seq_formated=DataFormatCorMut(seq_list)
	if(kaks==TRUE&!is.null(setPosition)) stop("kaks must be FALSE while specify the setPosition")
	seq_char=tolower(as.character(seq_formated))
	refseq=seq_char[1]
	refaa=translate_cm(s2c(refseq))
	seq_list=seq_char[2:length(seq_char)]
	#final_result_list=list()
	mat=sapply(as.list(seq_list),function(x)translate_cm(s2c(x)))
	#producing a list containing
	amino_set_list=list()
	codoni_list=list()
	codonio_list=list()
	Si_list=list()
	Ni_list=list()
	rownames_result=c()
	#retrieve positions##################################
	if(kaks==T){
		kaksfilter=kaksAA(seq_formated)
		positions=1:length(kaksfilter@kaks)
		positionf=c()
		positionaa=list()
		for(i in 1:length(kaksfilter@kaks)){
			positionnames=names(kaksfilter@kaks[[i]])
			#positionaa[[i]]=c()
			flag=0
			if(!is.null(positionnames)){
				positiontag=c()
				pkaks=kaksfilter@kaks[[i]]
				for(j in 1:length(pkaks)){
					if(pkaks[j]>1 & kaksfilter@lod[[i]][names(pkaks[j])]>=lod_cut){
						#print(j)
						monoaa=s2c(names(pkaks[j]))[nchar(names(pkaks[j]))]
						positiontag=append(positiontag,monoaa)
					flag=1
					}
				}
				positionaa[[as.character(i)]]=positiontag
			}
			if(flag==1) positionf=append(positionf,i)
		}
		prange=positionf
	}else if(kaks==FALSE & !is.null(setPosition)){
		prange=as.integer(setPosition)
	}else{
		prange=1:dim(mat)[1]
	}


########################
	for(i in prange){
		specific_row=mat[i,]
		if(kaks==T){
			amino_set=positionaa[[as.character(i)]]
			
		}else{
			amino_set=setdiff(unique(specific_row[specific_row!=refaa[i]]),c("X","*"))
		}
		codoni=sapply(as.list(seq_list),function(x)substr(x,(i-1)*3+1,i*3))
		codonio=substr(refseq,(i-1)*3+1,i*3)
		Si=(((mat[i,]==refaa[i])+(codoni!=codonio))==2)
		Ni=codoni==codonio
		amino_set_list[[as.character(i)]]=amino_set
		codoni_list[[as.character(i)]]=codoni
		codonio_list[[as.character(i)]]=codonio
		Si_list[[as.character(i)]]=Si
		Ni_list[[as.character(i)]]=Ni
		if(length(amino_set)>=1){
			inames=sapply(amino_set,function(xx)paste(c(refaa[i],i,xx),collapse=""),USE.NAMES=FALSE)
			rownames_result=append(rownames_result,inames)
		}
	}
	result_mat=matrix(NA,length(rownames_result),length(rownames_result))
	lod_mat=matrix(NA,length(rownames_result),length(rownames_result))
	rownames(result_mat)=rownames_result
	colnames(result_mat)=rownames_result
	rownames(lod_mat)=rownames_result
	colnames(lod_mat)=rownames_result
	#q values computation
	q_result=list()
	fvft=fvft(seq_list,refseq=refseq)
	for(i in prange){
		targetbp=s2c(refseq)[((i-1)*3+1):(i*3)]
		q_slist=c()
		amino_set=amino_set_list[[as.character(i)]]
		for(am in amino_set){
			names_label=paste(translate_cm(targetbp),am,sep=as.character(i))
			NtNv=Ntnv_aa(targetbp,aa_m=am)
			q_fenzi=(NtNv$Nat*fvft$ft+NtNv$Nav*fvft$fv)
			q=q_fenzi/(3*fvft$ft+6*fvft$fv)
			names(q)=names_label
			q_slist=append(q_slist,q)
		}
		q_result[[as.character(i)]]=q_slist
	}
	q_value=q_result
	#####################
	for(i in prange[2:length(prange)]){
		#fenzi
		amino_set=amino_set_list[[as.character(i)]]
		codoni=codoni_list[[as.character(i)]]
		codonio=codonio_list[[as.character(i)]]
		Si=Si_list[[as.character(i)]]
		Ni=Ni_list[[as.character(i)]]
		indexj02=match(i,prange)-1
		for(j in prange[1:indexj02]){
			amino_set02=amino_set_list[[as.character(j)]]
			codonjo=codonio_list[[as.character(j)]]
			codonj=codoni_list[[as.character(j)]]
			Sj=Si_list[[as.character(j)]]
			Nj=Ni_list[[as.character(j)]]
			for(am in amino_set){
				Mix=mat[i,]==am
				names_label01=paste(translate_cm(s2c(codonio)),am,sep=as.character(i))
				for(am02 in amino_set02){
					names_label02=paste(translate_cm(s2c(codonjo)),am02,sep=as.character(j))
					#fenzi
					Mjx=mat[j,]==am02
					Mij=sum((Mix+Mjx)==2)
					SjMi=sum((Mix+Sj)==2)
					SiMj=sum((Mjx+Si)==2)
					#add
					if(Mij==0){
						fenzi_ij=0
					}else if(Mij>0&SiMj==0){
						fenzi_ij=Mij
					}else{
						fenzi_ij=Mij/SiMj
					}
					if(Mij==0){
						fenzi_ji=0
					}else if(Mij>0&SjMi==0){
						fenzi_ji=Mij
					}else{
						fenzi_ji=Mij/SjMi
					}
					#fenmu
					MjNi=sum((Mjx+Ni)==2)
					MiNj=sum((Mix+Nj)==2)
					SjNi=sum((Ni+Sj)==2)
					SiNj=sum((Nj+Si)==2)
					#added
					#if(MjNi==0)MjNi=1
					#if(MiNj==0)MiNj=1
					#if(SjNi==0)SjNi=1
					#if(SiNj==0)SiNj=1
					#
					if(MiNj==0&SiNj==0){
						fenmu_ij=1
					}else if(MiNj!=0&SiNj==0){
						fenmu_ij=MiNj
					}else if(MiNj==0&SiNj!=0){
						fenmu_ij=1/SiNj
					}else{
						fenmu_ij=MiNj/SiNj
					}
					if(MjNi==0&SjNi==0){
						fenmu_ji=1
					}else if(MjNi!=0&SjNi==0){
						fenmu_ji=MjNi
					}else if(MjNi==0&SjNi!=0){
						fenmu_ji=1/SjNi
					}else{
						fenmu_ji=MjNi/SjNi
					}
					#fenmu_ij=MiNj/SiNj
					#fenmu_ji=MjNi/SjNi
					#if(!is.finite(fenmu_ij))fenmu_ij=1
					#if(!is.finite(fenmu_ji))fenmu_ji=1
					#influ01=paste(names_label01,names_label02,as.character(fenzi_ji/fenmu_ji),sep="\t")
					#influ02=paste(names_label02,names_label01,as.character(fenzi_ij/fenmu_ij),sep="\t")
					##matrix_c[i,j]=fenzi_ji/fenmu_ji
					##matrix_c[j,i]=fenzi_ij/fenmu_ij
					#LOD compute
					Nnumi=SiMj+Mij
					Nnumj=SjMi+Mij
					#print(c(Nnumi,Nnumj,Mij,q_value[i],q_value[j]))
					LODi=0
					LODj=0
					#q value is changed
					q=q_value[[as.character(i)]][[names_label01]]
					for(ii in Mij:Nnumi){
						LODi=choose(Nnumi,ii)*(q^ii)*((1-q)^(Nnumi-ii))+LODi
					}
					q=q_value[[as.character(j)]][[names_label02]]
					for(jj in Mij:Nnumj){
						LODj=choose(Nnumj,jj)*(q^jj)*((1-q)^(Nnumj-jj))+LODj
					}
					coltag01=paste(c(refaa[i],i,am),collapse="")
					coltag02=paste(c(refaa[j],j,am02),collapse="")
					lod_mat[coltag02,coltag01]=round(-log10(LODi),digits=2)
					lod_mat[coltag01,coltag02]=round(-log10(LODj),digits=2)
					result_mat[coltag02,coltag01]=round(fenzi_ij/fenmu_ij,digits=2)
					result_mat[coltag01,coltag02]=round(fenzi_ji/fenmu_ji,digits=2)
				}
			}
		}
	}
	final_result_list=new("ckaks",ckaks=result_mat,lod=lod_mat)
	return(final_result_list)
}
