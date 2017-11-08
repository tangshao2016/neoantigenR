
	##-----------------------------------
	## main function module 2 : alignment
	##-----------------------------------

    ## compute neoantigens using the novel algorithm after the preparation steps
	  #' compute neoantigens using the novel algorithm after the preparation steps
	  #' @return no values will be returned
		neoantigenR.get.peptides<-function(){
	  alignment.indel.detailed.info.table	=	data.frame()
	  for(chosen.gene in target.gene.list){
  		input.file.name=file.path(dataDir, paste(chosen.gene, ".pacbio.dnaseq.protein.sequence.txt", sep=""))
  		if(!is.na(file.size(input.file.name)) & file.size(input.file.name)>0){
  			database.name=chosen.gene
  			matched.database.index=which(uniprot.protein.names.new==chosen.gene & uniprot.protein.swissprot.index==1) ## only swiss-prot ?
  			if(length(matched.database.index)>0){
  				mySequences.all <- readAAStringSet(input.file.name) #Biostrings::readAAStringSet
  				num.unique.isoforms=round(length(mySequences.all)/6)
  				#forward.strand.isoform.idx=c(seq(from=1,to=num.unique.isoforms*6,by=6), seq(from=2,to=num.unique.isoforms*6,by=6),seq(from=3,to=num.unique.isoforms*6,by=6))	#1:(num.unique.isoforms*6)
  				forward.strand.isoform.idx=1:(num.unique.isoforms*6) # 11/7/2017
  				for(n.isof in forward.strand.isoform.idx){
  					ref.sequence.name=names(attributes(mySequences.all[n.isof])$ranges)
  					if(length(grep("_R_", ref.sequence.name))>0 | length(grep("_F_", ref.sequence.name))>0){
  						ref.sequence.name=substr(ref.sequence.name, 1, nchar(ref.sequence.name)-4)
  					}
  					database.seq.isoforms=unlist(lapply(matched.database.index, function(idx){uniprot.human.database[[idx]][1]}))
  					for(ds.idx in 1:length(database.seq.isoforms)){
  						database.isoform.name=attributes(uniprot.human.database[[matched.database.index[ds.idx]]])$Annot
  						database.isoform.uniprotID=strsplit(database.isoform.name, "|",fixed=TRUE)[[1]][2]
  						reference.sequence=toString(mySequences.all[n.isof])
  						database.sequence=database.seq.isoforms[ds.idx]
  						mySequences.new=union(reference.sequence, database.sequence) # database.seq.isoforms)
  						myFirstAlignment=msa::msa(mySequences.new, order="input", type="protein")
  						reference.sequence=gsub("-","z", reference.sequence) # z for stop
  						pacbio.protein.sequences=reference.sequence
  						alignment.full.details.strings=msa::msaConvert(myFirstAlignment, type="bios2mds::align")
  						alignment.full.details.matrix=matrix("", length(alignment.full.details.strings), length(alignment.full.details.strings[[1]]))
  						for(i in 1:length(alignment.full.details.strings)){
  							alignment.full.details.matrix[i,]=alignment.full.details.strings[[i]]
  						}
  						alignment.matrix=alignment.full.details.matrix
  						identity.not.equal.vector=rep(0, ncol(alignment.matrix))
  						identity.equal.matrix=matrix(0, nrow(alignment.matrix)-1, ncol(alignment.matrix))
  						for(i in 1:ncol(alignment.matrix)){
  							if(length(unique(alignment.matrix[1:nrow(alignment.matrix),i]))==1){
  								identity.not.equal.vector[i]=0
  								identity.equal.matrix[1:(nrow(alignment.matrix)-1),i]=1
  							}else if(alignment.matrix[1,i]!="-" & length(unique(alignment.matrix[2:nrow(alignment.matrix),i]))==1 &
  								unique(alignment.matrix[2:nrow(alignment.matrix),i])[1]=="-"){
  								identity.not.equal.vector[i]=1
  							}else if(alignment.matrix[1,i]=="-" & length(unique(alignment.matrix[2:nrow(alignment.matrix),i]))==1 &
  								unique(alignment.matrix[2:nrow(alignment.matrix),i])[1]!="-"){
  								identity.not.equal.vector[i]=-1
  							}else{
  								identity.not.equal.vector[i]=2
  								identity.equal.matrix[which(alignment.matrix[2:nrow(alignment.matrix),i]==alignment.matrix[1,i]),i]=1
  							}
  						}
  						if(length(which(identity.not.equal.vector==0))>minimum.sequence.similarity*(min(c(nchar(reference.sequence), nchar(database.sequence))))){
  							k=1
  							mismatch.indel.regions=list()
  							prev.end=0;curr.start=0
  							cnt=0
  							new.common.seq=TRUE
  							for(j in 1:length(identity.not.equal.vector)){
  								if(identity.not.equal.vector[j]==0){
  									if(new.common.seq==TRUE){ # leading 0
  										cnt=1
  										curr.start=j-1  # a new start
  										new.common.seq=FALSE
  									}else{
  										cnt=cnt+1
  									}
  									if(prev.end>0 & cnt>=minimum.flanking.region.size &
  											(curr.start-prev.end)>=min.mismatch.size & (curr.start-prev.end)<=max.mismatch.size){  # delete the following condition 5/24/2017
  												## & (length(which(identity.not.equal.vector==0))/length(identity.not.equal.vector))>=minimum.sequence.similarity){ # minimum.alignment.match.ratio){
  										this.mismatch.indel.region=c(prev.end, curr.start)
  										left.flanking.seq=paste(alignment.matrix[1, (prev.end-minimum.flanking.region.size):(prev.end-1)], collapse="")
  										right.flanking.seq=paste(alignment.matrix[1, (curr.start+1):(curr.start+minimum.flanking.region.size)], collapse="")
  										mismatch.indel.region.seq=paste(alignment.matrix[1, prev.end:curr.start], collapse="")
  										mismatch.indel.region.seq=gsub("\\-", "", mismatch.indel.region.seq)  # could be all "---" if deletions
  										if(regexpr(mismatch.indel.region.seq, reference.sequence)[1]!=(-1) &
  											regexpr("-", left.flanking.seq)[1]<0 & regexpr("-", right.flanking.seq)[1]<0){
  											mismatch.indel.regions[[k]]=this.mismatch.indel.region;k=k+1
  										}
  										prev.end=0 # reset the prev.end, but keep the counter
  									}
  								}
  								if(identity.not.equal.vector[j]!=0 & cnt>0){
  									new.common.seq=TRUE
  									if(cnt>=minimum.flanking.region.size){  # a valid left flanking sequence, so a begining of mismatch/indel
  										prev.end=j
  										cnt=0     # reset counter
  									}
  								}
  							}
  							if(length(mismatch.indel.regions)>0){
  								for(pr in 1:length(mismatch.indel.regions)){
  									prev.end  =mismatch.indel.regions[[pr]][1]
  									curr.start=mismatch.indel.regions[[pr]][2]
  									indel.region.vals=identity.not.equal.vector[prev.end:curr.start]
  									type.indel.region=
  									if(length(unique(indel.region.vals))==1){
  										if(unique(indel.region.vals)==(-1)){
  											"deletion"
  										}else if(unique(indel.region.vals)==1){
  											"insertion"
  										}else{
  											"mismatch"
  										}
  									}else{
  										"mismatch"
  									}
  									indel.seq=
    									if(type.indel.region=="deletion"){
    									  paste(alignment.matrix[1, c((prev.end-minimum.flanking.region.size):(prev.end-1),(curr.start+1):(curr.start+minimum.flanking.region.size))], collapse="")
    									}else{
    										paste(alignment.matrix[1, prev.end:curr.start], collapse="")
    									}
  									indel.seq.larger=indel.seq  # extend both size for 10 AA for neoantigen prediction
  									if(type.indel.region!="deletion"){
  									  left.end=prev.end-10
  									  right.end=curr.start+10
  									  if(left.end<1){left.end=1}
  									  if(right.end>dim(alignment.matrix)[2]){right.end=dim(alignment.matrix)[2]}
  									  indel.seq.larger=paste(alignment.matrix[1, left.end:right.end], collapse="")
  									}
  									indel.region.seq.overlap.ratio=0
  									reference.indel.region.seq=""
  									if(type.indel.region=="mismatch"){
  										reference.indel.region.seq=paste(alignment.matrix[2, prev.end:curr.start], collapse="")
  										reference.indel.region.seq=gsub("\\-", "", reference.indel.region.seq)
  										indel.region.seq.overlap.ratio=(1-round(length(which(indel.region.vals==2))/length(indel.region.vals),2))
  									}
  									if(type.indel.region=="deletion"){
  										reference.indel.region.seq=paste(alignment.matrix[2, prev.end:curr.start], collapse="")
  									}
  									indel.seq.no.hypen=gsub("\\-", "", indel.seq)
  									if(nchar(indel.seq.no.hypen)/nchar(indel.seq)<0.200000){
  										indel.seq=paste(alignment.matrix[2, prev.end:curr.start], collapse="")  # must be a deletion, any database isoform is good
  									}
  									indel.seq=gsub("\\-", "", indel.seq)
  									match.all.database.isoforms.index=
  									if(nchar(indel.seq)<1460){
  										attributes(regexpr(indel.seq, database.seq.isoforms))$match.length
  									}else{
  										1
  									}
  									if(validate.in.database.boolean==FALSE | (validate.in.database.boolean==TRUE & max(match.all.database.isoforms.index)<0)){
  									  chose.pacbio.transcript.id=as.vector(sample.pacbio.gencode.id.mapping[which(sample.pacbio.gencode.id.mapping[,2]==chosen.gene), 1])
  									  pacbio.transcript.index=ceiling(n.isof/6) # if(n.isof%%6==0){n.isof/6}else{round(n.isof/6)+1}
  									  chose.pacbio.transcript.id=strsplit(chose.pacbio.transcript.id[pacbio.transcript.index], ";")[[1]][1] # need global variable from function above Line-1246
  									  output.file.name=file.path(dataDir, paste(chosen.gene, ".", chose.pacbio.transcript.id, ".F", n.isof, ".Ref", ds.idx, ".pacbio.dnaseq.protein.msa.txt", sep=""))
  										sink(output.file.name)
  										print(myFirstAlignment, show="complete")
  										cat("\n")
  										sink()
  										cat(prev.end, " ", curr.start, " ", indel.seq, "\n")
  										if(type.indel.region=="deletion"){
  											curr.start=prev.end-1
  											prev.end=prev.end-minimum.flanking.region.size
  										}
  										if(type.indel.region=="mismatch"){
  											curr.start=prev.end+nchar(indel.seq)-1
  										}
  										indel.pacbio.info=as.vector(unlist(find.indels.location.by.location(chosen.gene, chose.pacbio.transcript.id, indel.seq, pacbio.protein.sequences,
  											gene.df.pacbio, reference.gff.info.chosen, pacbio.isof.ids, compute.query.exon.details, compute.database.exon.details)))
  										average.coverage=0

  										# prev.end is the location in the alignment matrix, not from its original sequence
  										indel.seq.pos=regexpr(indel.seq, reference.sequence)[1]
  										prev.end.new=indel.seq.pos
  										curr.start.new=prev.end.new + (curr.start-prev.end)

  										indel.seq.detailed.info=c(chosen.gene, ref.sequence.name, database.isoform.name, prev.end.new, curr.start.new, prev.end.new*3-2, curr.start.new*3, indel.seq, average.coverage, n.isof, reference.indel.region.seq, indel.region.seq.overlap.ratio, type.indel.region, database.isoform.uniprotID, indel.seq.larger, indel.pacbio.info)
  										if(length(indel.seq.detailed.info)<36){
  											indel.seq.detailed.info=c(indel.seq.detailed.info, rep("", 36-length(indel.seq.detailed.info)))
  										}

  										indel.seq.detailed.info=t(data.frame(indel.seq.detailed.info))
  										colnames(indel.seq.detailed.info)=c("Gene","PacBioSequenceID","UniprotProteinID","ProteinStart","ProteinEnd","DNAStart","DNAEnd","IndelSequence","Coverage","Isof.Index", "reference.indel.region.seq", "indel.region.seq.overlap.ratio", "type", "uniprotID", "IndelSequenceAnchorSeq", "IndepPositioninGenome","Chrom","PacBio","region","PacBioStart","PacbioEnd","dot","strand","dot2","PacbioTranscriptID","ReferenceChrom","Reference_exon_start","Reference_exon_end", "dot3", "geneName", "strand2", "exontype","geneID","exonID","transcriptID","unkonwn")

  										if(nrow(alignment.indel.detailed.info.table)==0){
  										  alignment.indel.detailed.info.table=indel.seq.detailed.info
  										}else{
  										  alignment.indel.detailed.info.table=rbind(alignment.indel.detailed.info.table, indel.seq.detailed.info) ## Error in match.names(clabs, names(xi)) :	names do not match previous names
  										}
  									}
  								}
  							}
  						}
  					}
  				}
  			}
  		}
	  }
	  if(nrow(alignment.indel.detailed.info.table)>0){
	    colnames(alignment.indel.detailed.info.table)=c("Gene","PacBioSequenceID","UniprotProteinID","ProteinStart","ProteinEnd","DNAStart","DNAEnd","IndelSequence","Coverage","Isof.Index", "reference.indel.region.seq", "indel.region.seq.overlap.ratio", "type", "uniprotID", "IndelSequenceAnchorSeq", "IndepPositioninGenome","Chrom","PacBio","region","PacBioStart","PacbioEnd","dot","strand","dot2","PacbioTranscriptID","ReferenceChrom","Reference_exon_start","Reference_exon_end", "dot3", "geneName", "strand2", "exontype","geneID","exonID","transcriptID","unkonwn")
	  }
	  #print(alignment.indel.detailed.info.table)
	  alignment.indel.detailed.info.table <<- alignment.indel.detailed.info.table
	  #print(alignment.indel.detailed.info.table)
	}


