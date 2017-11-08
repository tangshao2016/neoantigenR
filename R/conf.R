	##--------------------------------
	## required input file names
	##--------------------------------

	#' the organism used to search protein database
	#' @examples
	#' org="hg19"
	#' @return the organism used to search protein database
	org="hg19"
	#' the path for reference protein database file (Swiss-uniprot)
	#' @examples
	#' protein.database.file.name="swissuniprots.fasta"
	#' @return the path for reference protein database file (Swiss-uniprot)
	protein.database.file.name		   =	system.file("extdata", "swissuniprots.fasta", package="neoantigenR")
	#' the path for reference annotated gene database file
	#' @examples
	#' reference.gff.file="gencode.v19.annotation.gff3"
	#' @return the path for reference annotated gene database file
	reference.gff.file					     =	system.file("extdata", "gencode.v19.annotation.gff3", package="neoantigenR")
	#' the predicted gff file
	#' @examples
	#' pacbio.gff="cufflinks.gff"
	#' @return the predicted gff file
	pacbio.gff							         =	system.file("extdata", "model.gff", package="neoantigenR")
	#' the predicted gff file intersecting with reference gene annotation file by BEDTools
	#' @examples
	#' pacbio.gencode.overlapping.file	 = "bedtool.intersect.overlaps.txt"
	#' @return the predicted gff file intersecting with reference gene annotation file by BEDTools
	pacbio.gencode.overlapping.file	 =	system.file("extdata", "bedtool.intersect.overlaps.txt", package="neoantigenR")
	#' the folder of the output files
	#' @examples
	#' output.folder = "analysis"
	#' @return the folder of the output files
	output.folder					         	 =	system.file("extdata", package="neoantigenR")


	
	
	##-----------------------------------
	## setup parameters
	##-----------------------------------

	#' indicate whether we will use reference gene annotation or now
	#' @examples
	#' use.reference.annotation.gff=TRUE
	#' @return indicate whether we will use reference gene annotation or now
	use.reference.annotation.gff      = TRUE
	#' indicate whether we will write dna sequence for protein translation
	#' @examples
	#' write.dna.seq.by.reference.genome	=	TRUE
	#' @return indicate whether we will write dna sequence for protein translation
	write.dna.seq.by.reference.genome	=	TRUE 	# if true, name Hsapiens instance
	#' the alternative isoform and reference protein sequence must be 80 percentage or more overlap
	#' @examples
	#' percOverlap	= 0.8
	#' @return the alternative isoform and reference protein sequence must be 80 percentage or more overlap
	percOverlap							=	0.8 	        # used in determing the overlap in protein MSA analysis, disabled !!
	#' mimimum anchor size in the indel region
	#' @examples
	#' minimum.flanking.region.size = 10
	#' @return mimimum anchor size in the indel region
	minimum.flanking.region.size		=	10
	#' at least 3 amino acids must be different from the annotation
	#' @examples
	#' min.mismatch.size=3
	#' @return at least 3 amino acids must be different from the annotation
	min.mismatch.size					=	3  		      # need at least a region >= 3 AAs
	#' at most 1000 amino is different from the annotation
	#' @examples
	#' max.mismatch.size=1000
	#' @return at most 1000 amino is different from the annotation
	max.mismatch.size					=	1000
	#' the alternative isoform and reference protein sequence must be 50 percentage or more overlap
	#' @examples
	#' minimum.sequence.similarity     = 0.5
	#' @return the alternative isoform and reference protein sequence must be 50 percentage or more overlap
	minimum.sequence.similarity     = 0.5   # at least 50 percentage sequence overlap
	#' look for the genomic location of the indel sequence in the predicted gff gene models
	#' @examples
	#' compute.query.exon.details=TRUE
	#' @return look for the genomic location of the indel sequence in the predicted gff gene models
	compute.query.exon.details			=	TRUE
	#' look for the genomic location of the indel sequence in the annotated gff gene models
	#' @examples
	#' compute.database.exon.details=TRUE
	#' @return look for the genomic location of the indel sequence in the annotated gff gene models
	compute.database.exon.details		=	TRUE
	#' whether to validate indels in the annotated protein database
	#' @examples
	#' validate.in.database.boolean=TRUE
	#' @return whether to validate indels in the annotated protein database
	validate.in.database.boolean		=	TRUE

	#' the annotated gene models for chosen genes
	#' @examples
	#' reference.gff.info.chosen=data.frame()
	#' @return the annotated gene models for chosen genes
	reference.gff.info.chosen			=	data.frame()

	#' coverage file to calculate sequence coverage for indel regions with BAM files
	#' @examples
	#' coverage.file= "alignment.bam"
	#' @return coverage file to calculate sequence coverage for indel regions with BAM files
	coverage.file						= ""			# might need to be removed
	#' the predicted gene model's gene ids (prediction names)
	#' @examples
	#' refseq.pacbio.name.ids=c()
	#' @return the predicted gene model's gene ids (prediction names)
	refseq.pacbio.name.ids				= c()  		# might need to be removed


	##-----------------------------------
	## define global variables
	##-----------------------------------

	#' the gff file in the dat frame mode
	#' @examples
	#' gff.file=data.frame()
	#' @return the gff file in the dat frame mode
	gff.file=data.frame()
	#' the gff file's gene names
	#' @examples
	#' gff.names="prediction.gff"
	#' @return the gff file's gene names
	gff.names=""
	#' the directly where code and example data are stored
	#' @examples
	#' mainDir="data"
	#' @return the directly where code and example data are stored
	mainDir=""
	#' the directly where data is stored
	#' @examples
	#' dataDir=""
	#' @return the directly where data is stored
	dataDir=""
	#' the human protein names from uniprot
	#' @examples
	#' uniprot.protein.names.new=""
	#' @return the human protein names from uniprot
	uniprot.protein.names.new=""
	#' the human protein names from uniprot swissport (only the rows or indexes out of the entire uniprot file)
	#' @examples
	#' uniprot.protein.swissprot.index=""
	#' @return the human protein names from uniprot swissport (only the rows or indexes out of the entire uniprot file)
	uniprot.protein.swissprot.index=""
	#' the annotated gene in data frame mode
	#' @examples
	#' reference.gff=data.frame()
	#' @return the annotated gene in data frame mode
	reference.gff=data.frame()
	#' the annotated gene in data frame mode with attributes splitted
	#' @examples
	#' reference.gff.info=data.frame()
	#' @return the annotated gene in data frame mode with attributes splitted
	reference.gff.info=data.frame()
	#' the chosen column of gff file (to remove unused columns from analysis)
	#' @examples
	#' chosen.column.indicies=c()
	#' @return the chosen column of gff file (to remove unused columns from analysis)
	chosen.column.indicies=c()
	#' the uniprot human dataset in data frame
	#' @examples
	#' uniprot.human.database=data.frame()
	#' @return the uniprot human dataset in data frame
	uniprot.human.database=data.frame()
	#' the gene list finally chosen for neoantigen analysis
	#' @examples
	#' target.gene.list=c("MET")
	#' @return the gene list finally chosen for neoantigen analysis
	target.gene.list=c()
	#' the mapping between predicted gene names and reference annotated gene names
	#' @examples
	#' sample.pacbio.gencode.id.mapping=data.frame()
	#' @return the mapping between predicted gene names and reference annotated gene names
	sample.pacbio.gencode.id.mapping=data.frame()
	#' the data frame of predicted gene models
	#' @examples
	#' gene.df.pacbio=data.frame()
	#' @return the data frame of predicted gene models
	gene.df.pacbio=data.frame()
	#' the isoform or gene names from predicted gff model
	#' @examples
	#' pacbio.isof.ids=c()
	#' @return the isoform or gene names from predicted gff model
	pacbio.isof.ids=c()
	#' the attribute information of the gff file
	#' @examples
	#' att.info=c()
	#' @return the attribute information of the gff file
	att.info=c()

