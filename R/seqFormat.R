setGeneric("seqFormat", 
	function(x,format=c("clustal","fasta","mase","phylip","msf")){
		formated = match.arg(format)
		seqAligned <- read.alignment(x,format=formated)
		#remove the gaps against the reference sequence.
		allSeq=as.character(seqAligned$seq,use.names=F)
		ref01=s2c(allSeq[1])
		seq_f01=as.list(sapply(allSeq,function(x)c2s(s2c(x)[ref01!="-"])))
		seqff01=rawreplace(seq_f01,refseq=seq_f01[[1]])
		seq_f03=lapply(seqff01,c2s)
		ReturnSeq=unlist(seq_f03,use.names = F)
		names(ReturnSeq)=seqAligned$nam
		return(ReturnSeq)
	}
)