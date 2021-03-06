ckaksCodon <-
function(seq_formated,kaks=TRUE,lod_cut=2,setPosition=c()){
	#seq_formated=DataFormatCorMut(seq_list)
	if(kaks==TRUE&!is.null(setPosition)) stop("kaks must be FALSE while specify the setPosition")
	seq_char=tolower(as.character(seq_formated))
	refseq=seq_char[1]
	seq_list=seq_char[2:length(seq_char)]
	refaa=translate_cm(s2c(refseq))
	mat=sapply(seq_list,function(x)translate_cm(s2c(x)))
	#matrix_result=list()
	#position filter#################
	kaksfilter=kaksCodon(seq_formated)
	#obtain site positions
	positions=1:length(kaksfilter@kaks)
	kakscheck=((kaksfilter@kaks>1)+(kaksfilter@lod>=lod_cut))==2
	positionf=positions[kakscheck]
	if(kaks==T){
		prange=positionf
	}else if(kaks==FALSE & !is.null(setPosition)){
		prange=as.integer(setPosition)
	}else{
		prange=1:dim(mat)[1]
	}
	#################################
	matrix_c=matrix(NA,length(prange),length(prange))
	matrix_LOD=matrix(NA,length(prange),length(prange))
	rownames(matrix_c)=prange
	colnames(matrix_c)=prange
	rownames(matrix_LOD)=prange
	colnames(matrix_LOD)=prange
	#q values computation
	q_result=list()
	codoni_list=list()
	codonio_list=list()
	Mi_list=list()
	Si_list=list()
	Ni_list=list()
	fvft=fvft(seq_list,refseq=refseq)
	for(i in prange){
		targetbp=s2c(refseq)[((i-1)*3+1):(i*3)]
		NtNv=Ntnv(targetbp)
		q=(NtNv$Nat*fvft$ft+NtNv$Nav*fvft$fv)/(3*fvft$ft+6*fvft$fv)
		q_result[[as.character(i)]]=q
		codoni=sapply(seq_list,function(x)substr(x,(i-1)*3+1,i*3))
		codonio=substr(refseq,(i-1)*3+1,i*3)
		codoni_list[[as.character(i)]]=codoni
		codonio_list[[as.character(i)]]=codonio
		Mi_list[[as.character(i)]]=mat[i,]!=refaa[i]
		Si_list[[as.character(i)]]=(((mat[i,]==refaa[i])+(codoni!=codonio))==2)
		Ni_list[[as.character(i)]]=codoni==codonio
	}
	######################
	for(i in prange[2:length(prange)]){
		#iindex=match(i,prange)
		#fenzi
		codoni=codoni_list[[as.character(i)]]
		codonio=codonio_list[[as.character(i)]]
		Mi=Mi_list[[as.character(i)]]
		Si=Si_list[[as.character(i)]]
		#fenmu
		Ni=Ni_list[[as.character(i)]]
		indexj02=match(i,prange)-1
		for(j in prange[1:indexj02]){
			jindex=match(j,prange)
			#fenzi
			codonjo=codonio_list[[as.character(j)]]
			codonj=codoni_list[[as.character(j)]]
			Mj=Mi_list[[as.character(j)]]
			Sj=Si_list[[as.character(j)]]
			Mij=sum((Mi+Mj)==2)
			SjMi=sum((Mi+Sj)==2)
			SiMj=sum((Mj+Si)==2)
			#smoothing
			if(SiMj==0){
				fenzi_ij=(Mij+1)/(SiMj+1)
			}else{
				fenzi_ij=Mij/SiMj
			}
			if(SjMi==0){
				fenzi_ji=(Mij+1)/(SjMi+1)
			}else{
				fenzi_ji=Mij/SjMi
				}
			#fenmu
			Nj=codonj==codonjo
			MjNi=sum((Mj+Ni)==2)
			MiNj=sum((Mi+Nj)==2)
			SjNi=sum((Ni+Sj)==2)
			SiNj=sum((Nj+Si)==2)
			#added
			if(MjNi==0)MjNi=1
			if(MiNj==0)MiNj=1
			if(SjNi==0)SjNi=1
			if(SiNj==0)SiNj=1
			#
			if(MiNj==0|SiNj==0){
				fenmu_ij=(MiNj+1)/(SiNj+1)
			}else{
				fenmu_ij=MiNj/SiNj
			}
			if(MjNi==0|SjNi==0){
				fenmu_ji=(MjNi+1)/(SjNi+1)
			}else{
				fenmu_ji=MjNi/SjNi
			}
			matrix_c[as.character(i),as.character(j)]=fenzi_ji/fenmu_ji
			matrix_c[as.character(j),as.character(i)]=fenzi_ij/fenmu_ij
			#LOD compute
			Nnumi=SiMj+Mij
			Nnumj=SjMi+Mij
			#print(c(Nnumi,Nnumj,Mij,q_value[i],q_value[j]))
			LODi=0
			LODj=0
			q=q_result[[as.character(i)]]
			for(ii in Mij:Nnumi){
				#LODi=choose(Nnumi,ii)*(q^ii)*((1-q)^(Nnumi-ii))+LODi
				LODi=dbinom(ii,Nnumi,q)+LODi
			}
			q=q_result[[as.character(j)]]
			for(jj in Mij:Nnumj){
				#LODj=choose(Nnumj,jj)*(q^jj)*((1-q)^(Nnumj-jj))+LODj
				LODj=dbinom(jj,Nnumj,q)+LODj
			}
			#print(c(LODi,LODj))
			matrix_LOD[as.character(i),as.character(j)]=-log10(LODj)
			matrix_LOD[as.character(j),as.character(i)]=-log10(LODi)
		}
	}
	matrix_result=new("ckaks",ckaks=matrix_c,lod=matrix_LOD)
	return(matrix_result)
}
