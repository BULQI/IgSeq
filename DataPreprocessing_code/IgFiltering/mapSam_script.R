# R code to process BAM/SAM files
#   Feng 2020 June
#  start this to process the BAM file after bowtie2 filtering 
#	and write out fasta file  for the next step of processing, either umi or isotype demux.
#
#------June 9th, for now we need to do more experiments for making decisions about paramters/conditions
#		to pick good sequens
# 	see details about the conditions in ../pipeline.docx bullet point #3.
# Note: this now starts reading the quality trimmed data.
#this is the R code to do interactive analysis.
# see the R script in mapSam_script.R
#-----update June 12 2020
#this is from mapSam and now it is turned into R CMD BATCH script.
#
library(Rsamtools)

###############for command line input
## defaults should be either NULL or a named list
parseCommandArgs <- function(defaults=NULL, envir=globalenv()) {
  for (arg in commandArgs(TRUE))
    eval(parse(text=arg), envir=envir)

  for (nm in names(defaults))
    assign(nm, mget(nm, ifnotfound=list(defaults[[nm]]), envir=envir)[[1]], pos=envir)
}

defaults <- list(bamfiles=c()) ## default values of any arguments we might pass

## example usage:
parseCommandArgs(defaults)

#load bam files
#setwd("~/Windows/windowsD/feng/LAB/hg/IgSeq_MS/pipeline/IgFilter")
print(getwd())
#check for input files
if(length(bamfiles)==0)
{
	cat("Info:no file specified. please check\n");
	quit(save="no")
} else {
	cat("Info: total ", length(bamfiles), " specified\n");
}

#define a function to get the softclip numbers 
getSoftClipNum<-function( cigar, mode=c("start", "end"))
{
	#cigar<-aln.conc.read1$cigar
	if(missing(cigar))
	{
		stop("no cigar string input, please check")
	}
	mode=match.arg(mode);
	#remove soft clipping from both end
	clip<-gregexpr("^[0-9]+S", cigar)

	if(mode=="end")
	{
			clip<-gregexpr("[0-9]+S$", cigar)
	}
	clip.start<-unlist(clip)
	clip.match.length<-unlist(lapply(clip, FUN=function(x){attr(x, "match.length")}))

	clip.num<-rep(0, length(clip.start))

	#start.clip.start ==-1 mean ,no soft clip

	temp.dt<-data.frame(cigar[clip.start!=-1], clip.start[clip.start!=-1], clip.match.length[clip.start!=-1]) 
	clip.num[clip.start!=-1]<-apply(temp.dt, 1, FUN=function(x){as.numeric(substr(x[1], as.numeric(x[2]), as.numeric(x[2])+as.numeric(x[3])-1-1))})

	return (clip.num)
}
for(samFile in bamfiles)
{
		cat("Processing ", samFile, ".............\n");
		if(! file.exists(samFile))
		{
			cat("the file \"", samFile, "\" not exist , please check!\n")
			next
		}
		scamParam<-ScanBamParam(what=scanBamWhat(),tag=c("AS","XS"))
		samf<-scanBam(samFile, param=scamParam)
		
		#samf<-read.table(file=samFile, sep="\t", header=T)
		aln<-samf[[1]]

		#for each aln, 
		aln<-(do.call("DataFrame", aln))

		#first get only good pairs, concordantly, flag==99 (R1) and flag==147 (R2)
		aln.conc<-aln[aln$flag==99 | aln$flag==147,]; #do not include 83 and 163. another set of proper alignment. R1 and R2 reversed. very small number anyway

		#turn it in to a data frame to have two reads for each sequence on the same row
		aln.conc.read1<-aln.conc[aln.conc$strand=="+",]
		rownames(aln.conc.read1)<-aln.conc.read1$qname
		aln.conc.read2<-aln.conc[aln.conc$strand=="-",]
		rownames(aln.conc.read2)<-aln.conc.read2$qname
		#colnames(aln.conc.read2)<-paste0(colnames(aln.conc.read2), "_R2")
		aln.conc.read2<-aln.conc.read2[aln.conc.read1$qname,]
		#aln.concs<-cbind(aln.conc.read1, aln.conc.read2[aln.conc.read1$qname, c(-1)])

		#check for read 2 for the mapping, first mapping length, at least 100nt
		#start getting the softClipNum
		
		#get the matching length from the CIGAR
		cigar<-aln.conc.read2$cigar
		
		start_clip.num<-getSoftClipNum(cigar, mode="start")
		end_clip.num<-getSoftClipNum(cigar, mode="end")

		#remove soft clipping from both end
		vlen.R2<-aln.conc.read2$qwidth-end_clip.num-start_clip.num
		
		aln.conc.read2$vlen_R2<-vlen.R2

		aln.conc.read2<-aln.conc.read2[vlen.R2>=100, ]
		aln.conc.read1<-aln.conc.read1[vlen.R2>=100, ]

		#now working on read1, 
		cigar<-aln.conc.read1$cigar
		
		start_clip.num<-getSoftClipNum(cigar, mode="start")
		end_clip.num<-getSoftClipNum(cigar, mode="end")

		vlen.R1<-aln.conc.read1$qwidth-end_clip.num-start_clip.num
		aln.conc.read1$vlen_R1<-vlen.R1

		aln.conc.read1<-aln.conc.read1[vlen.R1>50,]
		aln.conc.read2<-aln.conc.read2[vlen.R1>50,]

		#now testing  of gaps
		gap<-aln.conc.read1$mpos-(aln.conc.read1$pos+aln.conc.read1$vlen_R1-1)
		sum(gap< 0)/length(gap)
		sum(gap< -5)/length(gap)
		sum(gap< -10)/length(gap)
		aln.conc.read1$gap<-gap   #no pick on gap for now <------------------------------

		#now let's do alignment score/quality. note not mapq
		# mapq is more about unique mapping
		score.threshold<-1.8
		sum(aln.conc.read2$tag.AS/aln.conc.read2$vlen_R2 < score.threshold)/dim(aln.conc.read2)[1]
		sum(aln.conc.read1$tag.AS/aln.conc.read1$vlen_R1 < score.threshold)/dim(aln.conc.read1)[1]

		sum(aln.conc.read1$tag.AS/aln.conc.read1$vlen_R1 > score.threshold& aln.conc.read2$tag.AS/aln.conc.read2$vlen_R2 > score.threshold)/dim(aln.conc.read1)[1]
		index<-which(aln.conc.read1$tag.AS/aln.conc.read1$vlen_R1 > score.threshold& aln.conc.read2$tag.AS/aln.conc.read2$vlen_R2 > score.threshold)

		aln.conc.read1<-aln.conc.read1[index,]
		aln.conc.read2<-aln.conc.read2[index,]
		
		#start writing the output.
		r1<-DNAStringSet(aln.conc.read1$seq)
		names(r1)<-aln.conc.read1$qname
		r2<-DNAStringSet(aln.conc.read2$seq)
		names(r2)<-aln.conc.read2$qname

		r1.qual<-PhredQuality(aln.conc.read1$qual)
		r2.qual<-PhredQuality(aln.conc.read2$qual)

		r1<-QualityScaledDNAStringSet(r1, r1.qual)
		r2<-QualityScaledDNAStringSet(r2, r2.qual)
		r2<-reverseComplement(r2)

		#writeQualityScaledXStringSet(r1, "Sample1_R1_filtered.fastq", #format="fastq",
		#                 compress=F)
		#writeXStringSet(r2, "Sample1_R2_filtered.fastq", format="fastq",
		 #               qualities=r2.qual, compress=T)
		
		#get file names.
		outFile<-sub("\\..*$", "", basename(samFile))
		#write output files 
		cat("writing output file for ", outFile, "\n")
		#next;
		#read 1
		r1.str<-paste0("@",names(r1),"\n", r1, "\n+\n",quality(r1))
		r1.dt<-data.frame(r1.str)
		r1.dt[,1]<-as.character(r1.dt[,1])
		write.table(file=paste0(outFile, "_R1_filtered.fastq"), x=r1.dt, row.names=F, col.names=F, quote=F) 

		#read 2
		r2.str<-paste0("@",names(r2),"\n", r2, "\n+\n",quality(r2))
		r2.dt<-data.frame(r2.str)
		r2.dt[,1]<-as.character(r2.dt[,1])
		write.table(file=paste0(outFile, "_R2_filtered.fastq"), x=r2.dt, row.names=F, col.names=F, quote=F) 
		cat("....Done\n")
}