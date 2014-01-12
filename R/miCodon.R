miCodon <-
function(seq_formated,kaks=TRUE,lod_cut=2,setPosition=c(),fdr=FALSE){
	#MI_results=list()
	if(kaks==TRUE&!is.null(setPosition)) stop("kaks must be FALSE while specify the setPosition")
	seq_char=tolower(as.character(seq_formated))
	seq_list=seq_char[2:length(seq_char)]
	refseq=seq_char[1]
	mat=sapply(seq_list,function(x)translate_cm(s2c(x)),USE.NAMES=FALSE)
	#filter site positions for positive selection sites
	if(kaks==T){
		kaksfilter=kaksCodon(seq_formated)
		positions=1:length(kaksfilter@kaks)
		kakscheck=((kaksfilter@kaks>1)+(kaksfilter@lod>=lod_cut))==2
		positionf=positions[kakscheck]
		prange=positionf
	}else if(kaks==FALSE & !is.null(setPosition)){
		prange=as.integer(setPosition)
	}else{
		prange=1:dim(mat)[1]
	}
	#creating result matrix
	result_mat=matrix(NA,length(prange),length(prange))
	pvalue_mat=matrix(NA,length(prange),length(prange))
	rownames(result_mat)=prange
	colnames(result_mat)=prange
	rownames(pvalue_mat)=prange
	colnames(pvalue_mat)=prange
	#############################################################
	for(i in prange[2:length(prange)]){
		indexj02=match(i,prange)-1
		for(j in prange[1:indexj02]){
			col1=mat[i,]
			col2=mat[j,]
			unique.col1=unique(col1)
			unique.col2=unique(col2)
			if(length(unique.col1)>=2 & length(unique.col2)>=2){
				unique_pairs=c()
				for(k in 1:length(col1)){
					unique_pairs=append(unique_pairs,paste(col1[k],col2[k],sep="#"))
				}
				pairs_count=table(unique_pairs)
				pairs_props=prop.table(pairs_count)
				fmat=matrix(0,nrow=length(unique.col1),ncol=length(unique.col2))
				fmatP=matrix(0,nrow=length(unique.col1),ncol=length(unique.col2))
				rownames(fmat)=unique.col1
				colnames(fmat)=unique.col2
				for(m in 1:length(pairs_props)){
					pair_split=unlist(strsplit(names(pairs_props[m]),split="#"))
					postag1=unique.col1==pair_split[1]
					postag2=unique.col2==pair_split[2]
					fmat[postag1,postag2]=pairs_props[m]
					fmatP[postag1,postag2]=pairs_count[m]
				}
				#compute mutual information
				freqs2d = as.matrix(fmat/sum(fmat))
				freqs.x = apply(fmat, 1, sum)
				freqs.y = apply(fmat, 2, sum)
				H.xy =  -sum(ifelse(freqs2d > 0, freqs2d * log(freqs2d), 0))
				H.x =  -sum(ifelse(freqs.x > 0, freqs.x * log(freqs.x), 0))
				H.y =  -sum(ifelse(freqs.y > 0, freqs.y * log(freqs.y), 0))
				MI = H.x + H.y - H.xy
				result_mat[as.character(i),as.character(j)]=MI
				result_mat[as.character(j),as.character(i)]=MI
				#compute p value
				#print(round(fmat*dim(mat)[1]))
				#pvalue_mat[as.character(i),as.character(j)]=fisher.test(round(fmat*dim(mat)[1]),simulate.p.value=T)$p.value
				pvalue=fisher.test(fmatP,simulate.p.value=T)$p.value
				pvalue_mat[as.character(i),as.character(j)]=pvalue
				pvalue_mat[as.character(j),as.character(i)]=pvalue
			}else{
				result_mat[as.character(i),as.character(j)]=0
				result_mat[as.character(j),as.character(i)]=0
				pvalue_mat[as.character(i),as.character(j)]=NA
				pvalue_mat[as.character(j),as.character(i)]=NA
			}
		}
	}
	fdrvalue=fdr
	if(fdrvalue){
		pvalue_mat_FDR=pvalue_mat
	}else{
		#p value adjustment by FDR
		pvalue_mat_FDR=matrix(p.adjust(as.vector(pvalue_mat),method="BH",n=sum(!is.na(pvalue_mat))),nrow=dim(pvalue_mat)[1],byrow=FALSE)
	}
	rownames(pvalue_mat_FDR)=prange
	colnames(pvalue_mat_FDR)=prange
	MI_results=new("MI",mi=result_mat,p.value=pvalue_mat_FDR)
	return(MI_results)
}
