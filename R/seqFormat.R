setGeneric("seqFormat", 
	function(x,format=c("clustal","fasta","DNAStringSet")){
		formated = match.arg(format)
		if(formated != "DNAStringSet"){
			seqAligned <- readDNAMultipleAlignment(filepath=x,format = formated)
		}else{
			seqAligned <- x
		}
		#remove the gaps against the reference sequence.
		allSeq=sapply(as.character(seqAligned,use.names=T),tolower)
		ref01=s2c(allSeq[1])
		seq_f01=as.list(sapply(allSeq,function(x)c2s(s2c(x)[ref01!="-"])))
		seqff01=rawreplace(seq_f01,refseq=seq_f01[[1]])
		seq_f03=lapply(seqff01,c2s)
		ReturnSeq=unlist(seq_f03,use.names = T)
		ReturnSeqf=DNAStringSet(ReturnSeq)
		return(ReturnSeqf)
	}
)