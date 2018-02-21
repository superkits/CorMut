jiAA <-function(seq_formated,kaks=TRUE,lod_cut=2,setPosition=c(),fdr=FALSE){
	#MI_results=list()
	#kaks=F
	#lod_cut=2
	#setPosition=allpcodon
	if(kaks==TRUE&!is.null(setPosition)) stop("kaks must be FALSE while specify the setPosition")
	seq_char=tolower(as.character(seq_formated))
	seq_list=seq_char[2:length(seq_char)]
	mat=sapply(seq_list,function(x)translate_cm(s2c(x)),USE.NAMES=FALSE)
	#reference sequence
	refseq=seq_char[1]
	refaa=translate_cm(s2c(refseq))
	freqComptRsut=freqCompute(mat=mat,refseq=refseq)
	matmat=freqComptRsut$newMat
	#################################
	#retrieve positions##################################
	if(kaks==T){
		kaksfilter=kaksAA(seq_formated)
		positions=filterSites(kaksfilter)$mutation
		matmat=matmat[positions,]
	}else if(kaks==FALSE & !is.null(setPosition)){
		allmutations_org=rownames(matmat)
		allmutations_pos=sapply(allmutations_org,function(xxx)substr(xxx,2,(nchar(xxx)-1)),USE.NAMES=F)
		matmat=matmat[allmutations_pos%in%setPosition,]
	}else{
		matmat=matmat
	}
	allmutations=rownames(matmat)
	result_mat=matrix(NA,length(allmutations),length(allmutations),dimnames=list(allmutations,allmutations))
	pvalue_mat=matrix(NA,length(allmutations),length(allmutations),dimnames=list(allmutations,allmutations))
	OR_mat=matrix(NA,length(allmutations),length(allmutations),dimnames=list(allmutations,allmutations))
	for(kk in 2:length(allmutations)){
		for(mm in 1:(kk-1)){
			AAa=allmutations[kk]
			AAb=allmutations[mm]
			onev=matmat[AAa,]
			twov=matmat[AAb,]
			overlapv=sum((onev+twov)==2)
			v22=sum(twov)-overlapv
			v33=sum(onev)-overlapv
			N_overlapv=sum((onev+twov)==0)
			fmat=matrix(c(overlapv,v22,v33,N_overlapv),byrow=T,nrow=2)
			leastone=sum((onev+twov)>=1)
			jaindex=overlapv/leastone
			result_mat[AAa,AAb]=jaindex
			result_mat[AAb,AAa]=jaindex
			#if(overlapv==0) next
			per_fisher=fisher.test(fmat)
			pvalue_mat[AAa,AAb]=per_fisher$p.value
			pvalue_mat[AAb,AAa]=pvalue_mat[AAa,AAb]
			OR_mat[AAa,AAb]=per_fisher$estimate
			OR_mat[AAb,AAa]=OR_mat[AAa,AAb]
		}
	}
	fdrvalue=fdr
	if(fdrvalue){
		pvalue_mat_FDR=pvalue_mat
	}else{
		pvalue_mat_FDR=matrix(p.adjust(as.vector(pvalue_mat),method="BH",n=sum(!is.na(pvalue_mat))),nrow=dim(pvalue_mat)[1],byrow=F,dimnames=list(allmutations,allmutations))
	}
	JI_results=new("JI",JI=result_mat,p.value=pvalue_mat_FDR,OR=OR_mat)
	return(JI_results)
}
