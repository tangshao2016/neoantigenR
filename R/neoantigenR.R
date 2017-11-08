
	##--------------------------------
	## load required packages
	##--------------------------------


	library(seqinr)
	library(msa)
	library(GenomicRanges)
	library(Gviz)
	library(Biostrings)
	library(BSgenome)
	library(BSgenome.Hsapiens.UCSC.hg19)
	library(rtracklayer)


	##-----------------------------------
	## data preprossing
	##-----------------------------------

	## initialize the parameters and functions for neoantigen prediction
	#' initialize the parameters and functions for neoantigen prediction
	#' @return global variables will be set after initialization
	neoantigenR.initialize<-function(){

	  mainDir		<<-	output.folder
	  #subDir		<<-	"data"
	  dataDir		<<-	file.path(mainDir) #, subDir)
	  #if (!file.exists(subDir)){dir.create(dataDir)}
	  if (!file.exists(dataDir)){dir.create(dataDir)}


	  #gene.df.pacbio					=	gffRead(pacbio.gff)
	  #pacbio.isof.ids					=	as.vector(gene.df.pacbio[, 'attributes'])
	  #pacbio.isof.ids					=	unlist(lapply(pacbio.isof.ids, function(id){id=strsplit(id,  " ")[[1]][4]; substr(id, 2, nchar(id)-2)}))
	  #gene.df.pacbio[, 'attributes']	=	pacbio.isof.ids
	  #pacbio.isof.ids <<- pacbio.isof.ids
	  #gene.df.pacbio <<- gene.df.pacbio


	  genome <- Seqinfo(genome = "hg19")
	  gr_a <- import(pacbio.gff)
	  gr_b <- import(reference.gff.file)#"d:/coding/bioconductor/neoantigenR/neoantigenR/inst/extdata/gencode.v19.annotation.gff3")
	  gff.pairs <- findOverlapPairs(gr_a, gr_b, ignore.strand = TRUE)
	  gff.model <- first(gff.pairs)
	  gff.ref   <- second(gff.pairs)


	  model.table=as.data.frame(elementMetadata(gff.model))
	  gene.df.pacbio.tmp=cbind(as.vector(seqnames(gff.model)), model.table[, c(1,2)], as.data.frame(gff.model@ranges)[,1:2], model.table[, c(3)], as.vector(strand(gff.model)), rep(".", nrow(model.table)), model.table[, 6])
	  colnames(gene.df.pacbio.tmp)=c("seqname","source","feature","start","end","score","strand","frame","attributes")
	  gene.df.pacbio<<-gene.df.pacbio.tmp


	  pacbio.transcript.ids <<- gff.model$transcript_id
	  gff.names <<- gff.ref$gene_name


	  gff.file=cbind(as.vector(seqnames(gff.model)), as.data.frame(gff.model@ranges)[,1:2], gff.model$transcript_id, as.data.frame(gff.model@ranges)[,3], as.vector(strand(gff.model)), elementMetadata(gff.ref)[, c('type','gene_id', 'transcript_id', 'protein_id', 'gene_name')])
	  colnames(gff.file)=c("chromosome","start","end", "pacBio", "width","strand","feature","gene","exon","transcript","symbol")
	  gff.file<<-as.data.frame(gff.file)

	  # ============> delete 11/9/2017

	  #gff.file = gffRead(pacbio.gencode.overlapping.file)
	  #gff.file = gff.file[which(gff.file[,12]=="exon"),]
	  #gff.attributes = as.vector(gff.file[, 18])
	  #gff.names = unlist(lapply(gff.attributes, function(g){v=strsplit(g, ";")[[1]][7]; substr(v, 11, nchar(v))}))
	  #pacbio.transcript.ids = unlist(lapply(as.vector(gff.file$attributes), function(g){v=strsplit(g, ";")[[1]]; if(length(v)>1){v=v[2]}else{v=v[1]}; v2=strsplit(v, '"')[[1]];if(length(v2)>1){v2=v2[2]}else{v2=v2[1]}; v2}))
	  #gff.file$attributes = pacbio.transcript.ids # change "gene_id \"PB.1\"; transcript_id \"PB.1.1\";" ==> "PB.1.1"
	  #pacbio.transcript.ids <<- pacbio.transcript.ids
	  #gff.file <<- gff.file
	  #gff.names <<- gff.names

	  #11/09sample.pacbio.gencode.id.mapping <<- retrieve.geneName.pacbioID.mapping.v2(gff.file)
	  sample.pacbio.gencode.id.mapping <<- retrieve.geneName.pacbioID.mapping.v2(gff.pairs)
	  target.gene.list<<-unique(gff.names)


	  ref.table=as.data.frame(elementMetadata(gff.ref))
	  reference.gff.info.chosen.tmp=cbind(as.vector(seqnames(gff.ref)), as.data.frame(gff.ref@ranges)[,1:2], ref.table$gene_name, rep(".", nrow(ref.table)), as.vector(strand(gff.ref)), ref.table[, c('type', 'gene_id', 'ID', 'transcript_id', 'protein_id')])
	  colnames(reference.gff.info.chosen.tmp)=c("chromosome","start","end", "pacBio", "width","strand","feature","gene","exon","transcript","symbol")
	  reference.gff.info.chosen<<-reference.gff.info.chosen.tmp

	  #reference.gff					<<-	gffRead(reference.gff.file)
	  #att.info						<<-	data.frame(do.call('rbind', strsplit(as.character(reference.gff[,'attributes']),';',fixed=TRUE)))
	  #reference.gff.info				<<-	data.frame(reference.gff, att.info)
	  #reference.gff.info				<<-	reference.gff.info[which(reference.gff.info$feature=="exon"),]
	  #chosen.column.indicies			<<-	c("seqname", "start", "end", "X7", "score", "strand", "feature", "X3", "X1", "X4", "X14")
	  #reference.gff.info.chosen		<<-	reference.gff.info[, chosen.column.indicies]
	  #reference.gff.info.chosen[,4]	=	unlist(lapply(as.vector(reference.gff.info.chosen[,4]), function(r){strsplit(r, "=")[[1]][2]}))
	  #colnames(reference.gff.info.chosen)	=	c("chromosome","start","end", "pacBio", "width","strand","feature","gene","exon","transcript","symbol")
	  #reference.gff.info.chosen		<<-	data.frame(reference.gff.info.chosen, stringsAsFactors=FALSE)
	  #reference.gff.info.chosen[,8]	=	as.character(reference.gff.info.chosen[,8])
	  #reference.gff.info.chosen[,9]	=	as.character(reference.gff.info.chosen[,9])
	  #reference.gff.info.chosen[,10]	=	as.character(reference.gff.info.chosen[,10])


	  uniprot.human.database			<<-	read.fasta(protein.database.file.name, seqtype = "AA" ,  as.string = TRUE)
	  uniprot.protein.names.new		<<-	unlist(lapply(1:length(uniprot.human.database), function(s){name=attributes(uniprot.human.database[[s]])$Annot
	  name.part.pos=regexpr("GN=", name)[1]; pref.name=strsplit(substr(name, name.part.pos, nchar(name)), " ")[[1]][1]; substr(pref.name, 4, nchar(pref.name))}))
	  uniprot.protein.swissprot.index	<<-	unlist(lapply(1:length(uniprot.human.database), function(s){name=attributes(uniprot.human.database[[s]])$name
	  if(substr(name, 1, 2)=="sp"){1}else{0}}))
	}



	##-----------------------------------
	## main function module 1 : isoforms
	##-----------------------------------

	## retrieve the predicted gene's exon coordinates and prepare them for visualization
	#' retrieve the predicted gene's exon coordinates and prepare them for visualization
  #' @return a list of gene models will be produced from the GFF file
	neoantigenR.get.Model<-function(){
  	for(suthee.gene in target.gene.list){

  		gene.intersection=gff.file[which(gff.names==suthee.gene),]
  		if(nrow(gene.intersection)>2){
  			#colnames(gene.intersection)=paste("V", 1:dim(gene.intersection)[2], sep="")
  			#gene.intersection=gene.intersection[which(gene.intersection[,12]=="exon" & gene.intersection[,3]=="exon"),]
  			#vcf.info=data.frame(do.call('rbind', strsplit(as.character(gene.intersection[,18]),';',fixed=TRUE)))
  			#gene.intersection=data.frame(gene.intersection, vcf.info)
  			if(length(unique(gene.intersection$V1))>1){
  				gene.intersection=gene.intersection[which(gene.intersection$V1==unique(gene.intersection$V1)[1]),]
  			}
  			grtrack.list=list()
  			m=1
  			#chosen.column.indicies=c('V1','V4','V5',  'V9', 'V19', 'V7','V3', 'X3','X1','X4','X15')

  			if(dim(gene.intersection)[1]>2){# & length(intersect(chosen.column.indicies, colnames(gene.intersection)))==11){

  				#gene.df=gene.intersection[,chosen.column.indicies]

  			  gene.df=gene.intersection
  			  colnames(gene.df)=c("chromosome","start","end", "pacBio", "width","strand","feature","gene","exon","transcript","symbol")
  			  gene.df=gene.df[!duplicated(gene.df),]
  				gene.df=gene.df[which(gene.df$feature=="exon"),]
  				#gene.df$width=0

  				gene.df=gene.df[order(gene.df$pacBio, gene.df$start),]
  				chr.start.end.list=unlist(lapply(1:dim(gene.df)[1], function(n){paste(gene.df[n,1:3], collapse="-")}))
  				gene.df=gene.df[which(duplicated(chr.start.end.list)==FALSE),]
  				if(length(unique(gene.df$pacBio))>10){
  					gene.df=gene.df[which(gene.df$pacBio%in%unique(gene.df$pacBio)[1:10]),]
  					duplication.index=duplicated(gene.df[,c("chromosome","start","end", "pacBio")])
  					gene.df=gene.df[which(duplication.index==FALSE),]
  					gene.df=data.frame(gene.df, stringsAsFactors=FALSE)
  					gene.df$transcript=as.vector(gene.df$transcript)
  				}

  				if(use.reference.annotation.gff==TRUE){
  					reference.gff.info.chosen.gene=reference.gff.info.chosen[which(reference.gff.info.chosen$pacBio==suthee.gene),]
  					gene.df.gff=reference.gff.info.chosen.gene[order(reference.gff.info.chosen.gene$transcript, reference.gff.info.chosen.gene$start),]
  					gene.df=rbind(gene.df, gene.df.gff)
  				}
  				granges.df=makeGRangesFromDataFrame(gene.df)
  				granges.df@seqinfo@genome=org
  				my.gen <- genome(granges.df)
  				my.chr <- as.character(unique(seqnames(granges.df)))
  				my.itrack <- IdeogramTrack(genome = my.gen, chromosome = my.chr)
  				my.gtrack <- GenomeAxisTrack()

    			grtrack.list[[m]]=my.itrack; m=m+1
  				grtrack.list[[m]]=my.gtrack; m=m+1

  				isoform.names=paste(gene.df$pacBio, gene.df$transcript, sep=":")
  				for(isoform in unique(isoform.names)){
  					this.gene.df=gene.df[which(isoform.names==isoform),]	# gene.df[which(gene.df$pacBio==isoform),]
  					granges.df=makeGRangesFromDataFrame(this.gene.df)
  					granges.df@seqinfo@genome=org
  					if(write.dna.seq.by.reference.genome==TRUE){
  						file.name=file.path(dataDir, paste(suthee.gene, ".pacbio.dnaseq.protein.sequence.txt", sep=""))
  						write.dna.sequence.reference(Hsapiens, this.gene.df[, 1:4], file.name, suthee.gene)
  					}
  					my.gen <- genome(granges.df)
  					my.chr <- as.character(unique(seqnames(granges.df)))
  					my.grtrack <- GeneRegionTrack(granges.df, genome = my.gen, chromosome = my.chr, name = isoform) #, transcriptAnnotation = "symbol")
  					my.itrack <- IdeogramTrack(genome = my.gen, chromosome = my.chr)
  					my.gtrack <- GenomeAxisTrack()
  					my.atrack <- AnnotationTrack(granges.df, name = "CpG")
  					grtrack.list[[m]]=my.grtrack
  					m=m+1
  				}

  				transcript.type="PacBio"

  			}
  			if(length(grtrack.list)>0){
  				pdf(file.path(dataDir, paste(suthee.gene, "_", transcript.type, ".isoforms.pdf", sep="")))
  				plotTracks(grtrack.list)
  				dev.off()
  			}
  		}
  	}
	}


	##-----------------------------------
	## Postprocessing, outputing results
	##-----------------------------------

	## write the predicted neoantigens for final results
	#' write the predicted neoantigens for final results
	#' @return a predicted neoantigen file will be written
	neoantigenR.write<-function(){
  	alignment.indel.detailed.info.table=data.frame(alignment.indel.detailed.info.table)
  	alignment.indel.detailed.info.table=alignment.indel.detailed.info.table[!duplicated(alignment.indel.detailed.info.table),]
  	write.table(alignment.indel.detailed.info.table, file.path(dataDir, "putative_neoantigen_candidates.txt"), sep=";", quote=FALSE, row.names=FALSE, col.names=TRUE)

  	sink(file.path(dataDir, "putative_neoantigen_candidate_sequences.fasta"))
  	for(ik in 1:dim(alignment.indel.detailed.info.table)[1]){
  	  header=paste(">", as.vector(alignment.indel.detailed.info.table[ik, 1]), sep="")
  	  cat(header, "\n")
  	  indel.w.anchor.seq=as.vector(alignment.indel.detailed.info.table[ik, 'IndelSequenceAnchorSeq'])
  	  indel.w.anchor.seq=gsub("-", "", indel.w.anchor.seq)
  	  cat(indel.w.anchor.seq, "\n")
  	}
  	sink()
	}






