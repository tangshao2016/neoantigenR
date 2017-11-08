	##--------------------------------
	## utility functions
	##--------------------------------


	## use R to map the gene annotation name and prediction name
	#' read the gff gene prediction data frame
	#' @param gff.file the name of the pairs intersecting model and reference
	#' @return a data frame with annotated namme to gene prediction mapping
	retrieve.geneName.pacbioID.mapping.v2<-function(gff.file){
		#geneName.info=unlist(lapply(as.vector(gff.file[,18]), function(n){strsplit(n, ";")[[1]][7]}))
		#gencode.name=unlist(lapply(geneName.info, function(p){substr(p, 11, nchar(p))}))
		#sample.pacbio.gencode.id.mapping=cbind(gff.file$attributes, gencode.name)
		#sample.pacbio.gencode.id.mapping=sample.pacbio.gencode.id.mapping[!duplicated(sample.pacbio.gencode.id.mapping),]
	  sample.pacbio.gencode.id.mapping=unique(cbind(first(gff.file)$transcript_id, second(gff.file)$gene_name))
	  sample.pacbio.gencode.id.mapping=as.data.frame(sample.pacbio.gencode.id.mapping)
	  colnames(sample.pacbio.gencode.id.mapping)=c("model.name", "ref_name")
	  return(sample.pacbio.gencode.id.mapping)
	}

	## write the dna sequence of a predicted gene based on reference genome object
	#' write the dna sequence of a predicted gene based on reference genome object
	#' @param Hsapiens the name of input file
	#' @param gene.df a data frame (gff) of predicted gene's coordinates
	#' @param file.name the output file name to be written
	#' @param gene.name the gene of interest for dna sequence extraction
	#' @return a dna sequence file will be written to output file
	write.dna.sequence.reference<-function(Hsapiens, gene.df, file.name, gene.name){
		sink(file.name, append = TRUE)
		for(transcript in unique(gene.df$pacBio)){
			isof.seq=""
			for(index in which(gene.df$pacBio==transcript)){
				this.seq=getSeq(Hsapiens, gene.df[index, 1], gene.df[index, 2], gene.df[index, 3])
				isof.seq=paste(isof.seq, this.seq, sep="")
			}
			for(strand in c("+", "-")){
				strand.name=if(strand=="+"){"F"}else{"R"}
				for(j in c(0,1,2)){
					if(nchar(isof.seq)>6){
						seqv=seqinr::translate(s2c(isof.seq),sens=strand.name, frame=j)
						stop.count=length(which(seqv=="*"))
						seqv[which(seqv=="*")]="-"  # 10/19/2016, tried to disable it, or change it to "+" doesn't work
						if(stop.count<1000){
							pseq=paste(seqv, collapse="")
							cat(paste(">", gene.name, "_", strand.name, "_", j, sep=""))
							cat("\n")
							cat(pseq)
							cat("\n")
						}
					}
				}
			}
		}
		sink()
	}

	## find the place (exon) of the indel sequence in annotation and predicted genes
	#' find the place (exon) of the indel sequence in annotation and predicted genes
	#' @param chosen.gene the gene of interest
	#' @param chose.pacbio.transcript.id chosen gene's prediction name
	#' @param indel.seq the chosen gene's alternative sequences
	#' @param pacbio.protein.sequences the protein sequence of the predicted gene
	#' @param gene.df.pacbio the data frame of the predicted genes (Gff)
	#' @param reference.gff.info.chosen the data frame of the annnotated genes
	#' @param pacbio.isof.ids the predicted gene's prediction id
	#' @param compute.query.exon.details whether to find the neoantigen's position in predictions
	#' @param compute.database.exon.details whether to find the neoantigen's position in annotations
	#' @return a vector of the predicted gene's genomic locations and locations in annotations
	find.indels.location.by.location<-function(chosen.gene, chose.pacbio.transcript.id, indel.seq,
		  pacbio.protein.sequences, gene.df.pacbio, reference.gff.info.chosen, pacbio.isof.ids,
		  compute.query.exon.details=FALSE,
		  compute.database.exon.details=FALSE){

		dna.target.pos=gregexpr(indel.seq, pacbio.protein.sequences)[[1]][1]*3-2
		this.gene.gff=gene.df.pacbio[which(gene.df.pacbio$feature=="exon" & pacbio.isof.ids==chose.pacbio.transcript.id),]
		reference.gene.gff=reference.gff.info.chosen[which(reference.gff.info.chosen$pacBio == chosen.gene),]
		strand=unique(this.gene.gff[,'strand'])[1]
		target.pos=0
		pacbio.info=c()
		refseq.info=c()

		if(dna.target.pos>0){
			if(compute.query.exon.details==TRUE & nrow(this.gene.gff)>0){
				prev.exon.length.sum=0
				if(strand=="+"){
  				for(i in 1:nrow(this.gene.gff)){
  					sstart=this.gene.gff[i, 'start']
  					send=this.gene.gff[i, 'end']
  					srange=send - sstart + 1 # should add 1 ?
  					if(prev.exon.length.sum<=dna.target.pos & dna.target.pos<=(prev.exon.length.sum+srange)){
  						target.pos=sstart+(dna.target.pos-prev.exon.length.sum)-1
  						pacbio.info=as.character(unlist(this.gene.gff[i,]))
  						break
  					}
  					prev.exon.length.sum=prev.exon.length.sum+srange
  				}
				}else{
				  for(i in nrow(this.gene.gff):1){
				    sstart=this.gene.gff[i, 'start']
				    send=this.gene.gff[i, 'end']
				    srange=send - sstart + 1 # should add 1 ?
				    if(prev.exon.length.sum<=dna.target.pos & dna.target.pos<=(prev.exon.length.sum+srange)){
				      target.pos=send-(dna.target.pos-prev.exon.length.sum)+1
				      pacbio.info=as.character(unlist(this.gene.gff[i,]))
				      break
				    }
				    prev.exon.length.sum=prev.exon.length.sum+srange
				  }
	  		}
			}

			if(compute.database.exon.details==TRUE & nrow(reference.gene.gff)){
				for(i in 1:nrow(reference.gene.gff)){
					sstart=reference.gene.gff[i, 'start']
					send=reference.gene.gff[i, 'end']
					srange=send - sstart +  1 # should I add 1 ?
					if(sstart<=target.pos & target.pos<=send){
						refseq.info=as.character(unlist(reference.gene.gff[i,]))
						break
					}
				}
			}
		}
		return(c(target.pos, pacbio.info, refseq.info))
	}



