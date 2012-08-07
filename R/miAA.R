miAA <-
function(seq_formated,kaks=TRUE,lod_cut=2,setPosition=c()){
	#MI_results=list()
	if(kaks==TRUE&!is.null(setPosition)) stop("kaks must be FALSE while specify the setPosition")
	seq_char=tolower(as.character(seq_formated))
	seq_list=seq_char[2:length(seq_char)]
	mat=sapply(seq_list,function(x)translate_cm(s2c(x)),USE.NAMES=FALSE)
	#reference sequence
	refseq=seq_char[1]
	refaa=translate_cm(s2c(refseq))
	#################################
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
	#producing a list containing
	amino_set_list=list()
	rownames_result=c()
	for(i in prange){
		specific_row=mat[i,]
		if(kaks==T){
			amino_set=positionaa[[as.character(i)]]
		}else{
			amino_set=setdiff(unique(specific_row[specific_row!=refaa[i]]),c("X","*"))
		}
		amino_set_list[[as.character(i)]]=amino_set
		if(length(amino_set)>=1){
			inames=sapply(amino_set,function(xx)paste(c(refaa[i],i,xx),collapse=""),USE.NAMES=F)
			rownames_result=append(rownames_result,inames)
		}
	}
	result_mat=matrix(NA,length(rownames_result),length(rownames_result))
	pvalue_mat=matrix(NA,length(rownames_result),length(rownames_result))
	rownames(result_mat)=rownames_result
	colnames(result_mat)=rownames_result
	rownames(pvalue_mat)=rownames_result
	colnames(pvalue_mat)=rownames_result
	for(i in prange[2:length(prange)]){
		amino_set=amino_set_list[[as.character(i)]]
		indexj02=match(i,prange)-1
		for(j in prange[1:indexj02]){
			amino_set02=amino_set_list[[as.character(j)]]
			for(am in amino_set){
				AAi=mat[i,]==am
				for(am02 in amino_set02){
					AAj=mat[j,]==am02
					overlap=((AAi+AAj)==2)
					side1=(as.numeric(AAi)-as.numeric(overlap))==1
					side2=(as.numeric(AAj)-as.numeric(overlap))==1
					NN=dim(mat)[2]-sum(overlap)-sum(side1)-sum(side2)
					fmat=matrix(c(sum(overlap),sum(side1),sum(side2),NN),nrow=2)
					rownames(fmat)=c(am,paste("Non-",am,sep=""))
					colnames(fmat)=c(am02,paste("Non-",am02,sep=""))
					#compute mutual information
					freqs2d = as.matrix(fmat/sum(fmat))
					freqs.x = apply(freqs2d, 1, sum)
					freqs.y = apply(freqs2d, 2, sum)
					H.xy =  -sum(ifelse(freqs2d > 0, freqs2d * log(freqs2d), 0))
					H.x =  -sum(ifelse(freqs.x > 0, freqs.x * log(freqs.x), 0))
					H.y =  -sum(ifelse(freqs.y > 0, freqs.y * log(freqs.y), 0))
					MI = H.x + H.y - H.xy
					coltag01=paste(c(refaa[i],i,am),collapse="")
					coltag02=paste(c(refaa[j],j,am02),collapse="")
					result_mat[coltag01,coltag02]=MI
					result_mat[coltag02,coltag01]=MI
					#compute p value
					#print(round(fmat*dim(mat)[1]))
					pvalue=fisher.test(round(fmat*dim(mat)[1]),simulate.p.value=T,alternative="g")$p.value
					pvalue_mat[coltag01,coltag02]=pvalue
					pvalue_mat[coltag02,coltag01]=pvalue
				}
			}
		}
	}
	#p value adjustment by FDR
	pvalue_mat_FDR=matrix(p.adjust(as.vector(pvalue_mat),method="BH",n=sum(!is.na(pvalue_mat))),nrow=dim(pvalue_mat)[1],byrow=F)
	rownames(pvalue_mat_FDR)=rownames_result
	colnames(pvalue_mat_FDR)=rownames_result
	MI_results=new("MI",mi=result_mat,p.value=pvalue_mat_FDR)
	return(MI_results)
}
