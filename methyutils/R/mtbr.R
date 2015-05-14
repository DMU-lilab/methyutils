
####required libraries#######
library(Biostrings)
library(reshape2)
library(data.table)
###############CONSTANTS########################
kCcontext <- list(
	CG=c("CG"),
	CHG=c("CAG","CCG","CTG"),
	CHH=c("CAA","CAC","CAT","CCA","CCC","CCT","CTA","CTC","CTT")
)
# get C context from a DNAString object.
# Arguments:
# Return:

GetCcontextPosition <- function(dna.seq, context="CG", return.type="logicalMarker"){
	
	# check input parameters

	if(class(dna.seq) != "DNAString"){
		stop("dna.seq must to be DNAString")
	}
	
	if(!(context %in% names(kCcontext))){
		 stop("context must be: ", paste("\"", names(kCcontext), "\"", collapse=", ", sep=""))  
	}	
	
	# search posi by matchPattern function
	posi.views <- c()
	con <- kCcontext[[context]]

	for (i in con ){
		posi <- matchPattern(i,dna.seq)
		posi.views <- c(posi.views,start(posi))
	}	
		
 	# return C context according to return.type
	
	if(return.type == "logicalMarker"){
		rt <- logical(length(dna.seq))
		rt[posi.views] <- TRUE
		return(rt)
	} else if(return.type == "integerPosi"){
		return(posi.views)
	} else {
		stop("undefined return type. Only logicalMarker and integerPosi are available now.")
	}
} 

# get the records at CG sites from the mtbr that have informations for all C sites

GetCGmtbr <- function(mtbr, cg.position){

	# check the col count
	mtbr.colnames <- c("chrom", "posi", "strand", "rC", "rT")
	
	if(ncol(mtbr) != 5){
		stop("colcount must be 5, colnames must be: c(", paste("\"", mtbr.colnames, "\"", collapse=", ", sep=""), ")" )
	} 

	# get mtbr at CG site
	
	setnames(mtbr,names(mtbr),mtbr.colnames)
	
	mtbr$strand[mtbr$strand == "+"] <- "p"
	mtbr$strand[mtbr$strand == "-"] <- "n"
	cg.position[which(cg.position) + 1] <- TRUE
	cg.mtbr <- mtbr[cg.position[mtbr$posi],]
	
	m <- melt(cg.mtbr, id=c("chrom","posi","strand"), variable.name="read", value.name="cv") 
	
	m$posi <- ifelse(m$strand== "n",m$posi-1,m$posi)
	
	
	cst <- dcast(m, chrom + posi ~ read + strand, value.var="cv")
	cst$posi <- as.integer(cst$posi)
	cst[is.na(cst)] <- 0
			
	return(cst)
}



# save cg.mtbr to Rdata file
SaveMtbrRdata <- function(cg.mtbr,file.name){
	save(cg.mtbr,file = file.name)	
}

	
WriteSeparatedWig <- function(cgmtbr,wigfilePath,chr){

	head <- paste("variableStep chrom=", chr, " span=1\n", sep="")
	
	# generate the score and cover 
	
	Sp <- data.frame(posi=cgmtbr$posi,score=cgmtbr$rC_p/(cgmtbr$rC_p + cgmtbr$rT_p))
	Sp <- Sp[!is.na(Sp$score),]
	Sp <- Sp[Sp$score>0,]
	Rmp <- data.frame(posi=cgmtbr$posi,score=cgmtbr$rC_p)
	Rmp <- Rmp[!is.na(Rmp$score),]
	Rmp <- Rmp[Rmp$score>0,]
	Rtotalp <- data.frame(posi=cgmtbr$posi,score=(cgmtbr$rC_p + cgmtbr$rT_p))
	Rtotalp <- Rtotalp[!is.na(Rtotalp$score),]
	Rtotalp <- Rtotalp[Rtotalp$score>0,]
	Sn <- data.frame(posi=cgmtbr$posi+1,score=(-cgmtbr$rC_n/(cgmtbr$rC_n + cgmtbr$rT_n)))
	Sn <- Sn[!is.na(Sn$score),]
	Sn <- Sn[Sn$score<0,]
	Rmn <- data.frame(posi=cgmtbr$posi+1,score=-cgmtbr$rC_n)
	Rmn <- Rmn[!is.na(Rmn$score),]
	Rmn <- Rmn[Rmn$score<0,]
	Rtotaln <- data.frame(posi=cgmtbr$posi+1,score=-(cgmtbr$rC_n + cgmtbr$rT_n))
	Rtotaln <- Rtotaln[!is.na(Rtotaln$score),]
	Rtotaln <- Rtotaln[Rtotaln$score<0,]
	
	## save the wig file 
	
	Spf <- file.path(wigfilePath,"methy.score_p.wig")
	Rmpf <- file.path(wigfilePath,"methy.cover_p.wig")
	Rtotalpf <- file.path(wigfilePath,"total.cover_p.wig")
	Snf <- file.path(wigfilePath,"methy.score_n.wig")
	Rmnf <- file.path(wigfilePath,"methy.cover_n.wig")
	Rtotalnf <- file.path(wigfilePath,"total.cover_n.wig")

	cat(head, file=Spf,append=T)
	write.table(Sp, Spf, row.names=F, col.names=F, quote=F, append=T)
	cat(head, file=Rmpf,append=T)
	write.table(Rmp, Rmpf, row.names=F, col.names=F, quote=F, append=T)
	cat(head, file=Rtotalpf,append=T)
	write.table(Rtotalp, Rtotalpf, row.names=F, col.names=F, quote=F, append=T)
	cat(head, file=Snf,append=T)
	write.table(Sn, Snf, row.names=F, col.names=F, quote=F, append=T)
	cat(head, file=Rmnf,append=T)
	write.table(Rmn, Rmnf, row.names=F, col.names=F, quote=F, append=T)
	cat(head, file=Rtotalnf,append=T)
	write.table(Rtotaln, Rtotalnf, row.names=F, col.names=F, quote=F, append=T)	
}

WriteMixedWig <- function(cgmtbr,wigfilePath,chr){

	head <- paste("variableStep chrom=", chr, " span=1\n", sep="")
		
	## generate score and cover

	cgscore <- data.frame(posi=cgmtbr$posi,score=(cgmtbr$rC_p + cgmtbr$rC_n)/(cgmtbr$rC_p + cgmtbr$rT_p + cgmtbr$rC_n + cgmtbr$rT_n))
	cgscore <- cgscore[!is.na(cgscore$score),]
	cgscore <- cgscore[cgscore$score>0,]
	Rm <- data.frame(posi=cgmtbr$posi,score=(cgmtbr$rC_p + cgmtbr$rC_n))
	Rm <- Rm[!is.na(Rm$score),]
	Rm <- Rm[Rm$score>0,]
	Rtotal <- data.frame(posi=cgmtbr$posi,score=(cgmtbr$rC_p + cgmtbr$rT_p + cgmtbr$rC_n + cgmtbr$rT_n))
	Rtotal <- Rtotal[!is.na(Rtotal$score),]
	Rtotal <- Rtotal[Rtotal$score>0,]
		
	## save the wig file
		
	cgscoref <- file.path(wigfilePath,"methy.score.wig")
	Rmf <- file.path(wigfilePath,"methy.cover.wig")
	Rtotalf <- file.path(wigfilePath,"total.cover.wig")
	cat(head,file=cgscoref,append=T)
	write.table(cgscore, cgscoref, row.names=F, col.names=F, quote=F, append=T)
	cat(head, file=Rmf,append=T)
	write.table(Rm, Rmf, row.names=F, col.names=F, quote=F, append=T)
	cat(head, file=Rtotalf,append=T)
	write.table(Rtotal, Rtotalf, row.names=F, col.names=F, quote=F, append=T)
}


GetValueFromMtbrByRegion <- function(cg.mtbr,region,geonome="mm9"){
		################################
		##region is a data.frame with 3 columns of chrom, start, and end
		###############################
		#find mtbr in region
		colnames(region) <- c("chrom","start","end")
		cg.mtbr$rC <- cg.mtbr$rC_p + cg.mtbr$rC_n
		cg.mtbr$rT <- cg.mtbr$rT_p + cg.mtbr$rT_n
		cg.mtbr$coverNum <- cg.mtbr$rC + cg.mtbr$rT
		cg.mtbr$st <- findInterval(cg.mtbr$posi,region$start) 
		cg.mtbr$en <- findInterval(cg.mtbr$posi,region$end+1)+1 
		cgmtbrInregion <- cg.mtbr[cg.mtbr$st==cg.mtbr$en,]
		# count the region rC rT and coverNum
		rC <- c()
		rT <- c()
		coverNum <- c()
			for (j in c(1:length(region$chrom))){
			a <- sum(cgmtbrInregion[cgmtbrInregion$st == j,7])
			b <- sum(cgmtbrInregion[cgmtbrInregion$st == j,8])
			c <- sum(cgmtbrInregion[cgmtbrInregion$st == j,9])
			rC <- c(rC,a)
			rT <- c(rT,b)
			coverNum <- c(coverNum,c)
		}
		region$rC <- rC
		region$rT <- rT
		region$coverNum <- coverNum
		region$score <- region$rC/region$coverNum
		#region <- region[!is.na(region$score),]
		# count region cgNum
		chr <- unique(cg.mtbr$chrom)
		library("methyutils")
		if(geonome == "mm9"){
		
			library(BSgenome.Mmusculus.UCSC.mm9)
			dna.seq <- Mmusculus[[chr]]
			cposition <- GetCcontextPosition(dna.seq)
			region$cgNum <- apply(region[,c("start","end")], 1, function(rg)sum(cposition[rg["start"]:rg["end"]]))
			return(region)
		}
		else if (type == "hg19"){
		
			library(BSgenome.Hsapiens.UCSC.hg19)
			dna.seq <- Hsapiens[[chr]]
			cposition <- GetCcontextPosition(dna.seq)
			region$cgNum <- apply(region[,c("start","end")], 1, function(rg)sum(cposition[rg["start"]:rg["end"]]))
			return(region)
		}
		else {
			 stop("undefined type. Only mm9 and hg19 are available now.")
		}
} 


