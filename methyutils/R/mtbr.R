
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
	##add tow row 
	addtworows <- cg.mtbr[1:2,]
	addtworows$posi <- c(-2:-1)
	addtworows$strand <- c("p","n")
	cg.mtbr <- rbind(cg.mtbr,addtworows)
	##	

	m <- melt(cg.mtbr, id=c("chrom","posi","strand"), variable.name="read", value.name="cv") 
	
	m$posi <- ifelse(m$strand== "n",m$posi-1,m$posi)
	
	
	cst <- dcast(m, chrom + posi ~ read + strand, value.var="cv")
	cst$posi <- as.integer(cst$posi)
	cst[is.na(cst)] <- 0
	
	##remove one row
	cst <- cst[-1,]
	##		
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
	#
	Sn$posi <- as.integer(Sn$posi)
	#
	Rmn <- data.frame(posi=cgmtbr$posi+1,score=-cgmtbr$rC_n)
	Rmn <- Rmn[!is.na(Rmn$score),]
	Rmn <- Rmn[Rmn$score<0,]
	#
	Rmn$posi <- as.integer(Rmn$posi)
	#
	Rtotaln <- data.frame(posi=cgmtbr$posi+1,score=-(cgmtbr$rC_n + cgmtbr$rT_n))
	Rtotaln <- Rtotaln[!is.na(Rtotaln$score),]
	Rtotaln <- Rtotaln[Rtotaln$score<0,]
	#
	Rtotaln$posi <- as.integer(Rtotaln$posi)
	#
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



GetValueFromMtbrByRegion <- function(cg.mtbr,region,genome="mm9"){
		################################
		##region is a data.frame with 3 columns of chrom, start, and end
		###############################
		#find mtbr in region
		colnames(region) <- c("chrom","start","end")
		cg.mtbr$rC <- cg.mtbr$rC_p + cg.mtbr$rC_n
		cg.mtbr$rT <- cg.mtbr$rT_p + cg.mtbr$rT_n
		
		# count region's rC rT coverNum score and cgNum
		chr <- unique(cg.mtbr$chrom)
			
		if(genome == "mm9"){
			library("BSgenome.Mmusculus.UCSC.mm9")
			dna.seq <- Mmusculus[[chr]]
		}else if (genome == "hg19"){
			library("BSgenome.Hsapiens.UCSC.hg19")
			dna.seq <- Hsapiens[[chr]]
		}else {
			 stop("undefined genome. Only mm9 and hg19 are available now.")
		}
		
		cposition <- GetCcontextPosition(dna.seq)
		region$cgNum <- apply(region[,c("start","end")], 1, function(rg)sum(cposition[rg["start"]:rg["end"]]))
		## 
		rC <- integer(length(dna.seq))
		rC[cg.mtbr$posi] <- cg.mtbr$rC
		rT <- integer(length(dna.seq))
		rT[cg.mtbr$posi] <- cg.mtbr$rT
		coverNum <- rC + rT
		##
		cgCovermtbr <- cg.mtbr[cg.mtbr$rC>0,]
		cgCover <- logical(length(dna.seq))
		cgCover[cgCovermtbr$posi] <- TRUE
		##
		region$rC <- apply(region[,c("start","end")], 1, function(rg)sum(rC[rg["start"]:rg["end"]]))
		region$rT <- apply(region[,c("start","end")], 1, function(rg)sum(rT[rg["start"]:rg["end"]]))
		region$cgCover <- apply(region[,c("start","end")], 1, function(rg)sum(cgCover[rg["start"]:rg["end"]]))
		region$coverNum <- region$rC + region$rT
		region$score <- region$rC/region$coverNum
			
		return(region)
}
		
##sliding windows
swsCalc<-function(x, win=list(L=75, R=75)){
    library(zoo)
    xLen <- length(x)
    winLen<-win$L + win$R + 1
    sws <- rollsum(x,winLen)
    sws_head<-tail(cumsum(x[1:(winLen-1)]),win$L)
    sws_tail<-rev(tail(cumsum(rev(x[(xLen - winLen + 2):xLen])),win$R))
    sws <- append(sws_head, sws)
    sws <- append(sws, sws_tail)
    return(sws)
}

##regionAsBed
regionAsBed <- function(marker, cf.length=5, tolerance.length=NULL,chrom){
    #library(zoo)
    r <- rle(marker)
    if(is.null(tolerance.length)){
        ##The tolerance length is NUL, which means any number of false
        ##values, even 1, can break a regionn into two.
        end <- with(r,cumsum(lengths)[values & lengths>cf.length])
        start <- end - with(r,lengths[values & lengths>cf.length]) + 1
        #df <- data.frame(chrom,start,end)
        if(length(start) == 0 || length(end) == 0){
            df <- data.frame(chrom=character(), start=integer(), end=integer())
        } else {
            df <- data.frame(chrom, start, end)
        }
    } else {
        ##Tolerance length is not null, which means if the number of false
        ##is less than or equal to tolerance length, the two regions will
        ##be combined together
        ##Get all the starts and ends
        end<-cumsum(r$lengths)
        start<-end - r$lengths + 1
        revision.mk <- with(r, lengths<=tolerance.length & !values)
        ##The position that need to be revised
        revision.posi <- data.frame(cbind(start,end))[revision.mk,]
        ##revise the orignal marker 
        invisible(apply(revision.posi, 1, function(x) marker[x[1]:x[2]] <<- T))
        r <- rle(marker)
        end <- with(r,cumsum(lengths)[values & lengths>cf.length])
        start <- end - with(r,lengths[values & lengths>cf.length]) + 1
        if(length(start) == 0 || length(end) == 0){
            df <- data.frame(chrom=character(), start=integer(), end=integer())
        } else {
            df <- data.frame(chrom, start, end)
        }

    }
    return(df)
}

##cgDensity

cgDensity <- function(cg.mtbr,genome="mm9",window=300,maxpercent=0.005,overlap=40,gap=-40,cf.length=100, tolerance.length=2){
		
	cg.mtbr$rC <- cg.mtbr$rC_p + cg.mtbr$rC_n
	cg.mtbr$rT <- cg.mtbr$rT_p + cg.mtbr$rT_n
		
	# count region's rC rT coverNum score and cgNum
	chr <- unique(cg.mtbr$chrom)
			
	if(genome == "mm9"){
		library("BSgenome.Mmusculus.UCSC.mm9")
		dna.seq <- Mmusculus[[chr]]
	}else if (genome == "hg19"){
		library("BSgenome.Hsapiens.UCSC.hg19")
		dna.seq <- Hsapiens[[chr]]
	}else if (genome == "hg38"){
                library("BSgenome.Hsapiens.UCSC.hg38")
                dna.seq <- Hsapiens[[chr]]

	}else {
		 stop("undefined genome. Only mm9 and hg19,hg38 are available now.")
	}
		
	## get the CG positon
	cposition <- GetCcontextPosition(dna.seq)
	## get the rC and rT 
	rC <- integer(length(dna.seq))
	rC[cg.mtbr$posi] <- cg.mtbr$rC
	rT <- integer(length(dna.seq))
	rT[cg.mtbr$posi] <- cg.mtbr$rT
	## slidingwindow the cg ,rC and rT
	win <- list(L=window,R=window)
	cgs <- swsCalc(cposition,win)
	rCs <- swsCalc(rC,win)
        rTs <- swsCalc(rT,win)
        ## get methylation score
        score <- rCs/(rCs + rTs)
        score[is.na(score[])] <- 0
	##unity the cg number
        ordcgs <- order(cgs)
        maxNum <- cgs[ordcgs[length(dna.seq)*(1-maxpercent)]]
        cgsu <- cgs*(100/maxNum)
        cgsu[cgsu[]>=100] <- 100
	##unity the methylation
	scoreu <- score * 100
	## make one data.frame for cg and methylation
        cs <- data.frame(posi=1:length(cgsu), cg=cgsu[], methy=scoreu[])
        cs$icmp <- cs$cg + cs$methy - 100
        ##cutoff for incomplemtarity
        cutoff <- list(gap=gap, overlap=overlap, len=cf.length, tolerate=tolerance.length)
        ##marker for overlap and gap incomplementarity
        mlap <- cs$icmp > cutoff$overlap
        mgap <- cs$icmp < cutoff$gap
	##find incomplementary overlap and gap regions
	rglap <- regionAsBed(marker=mlap, cf.length=cutoff$len, tolerance.length=cutoff$tolerate,chrom=chr)
	rggap <- regionAsBed(marker=mgap, cf.length=cutoff$len, tolerance.length=cutoff$tolerate,chrom=chr)
	##covernt index to chrom coordinate
	rglapc <- data.frame(chrom=rglap$chrom, start=cs$posi[rglap$start], end=cs$posi[rglap$end])
	rggapc <- data.frame(chrom=rggap$chrom, start=cs$posi[rggap$start], end=cs$posi[rggap$end])
	## gap and overlap content
	icmpValue <- list()
	icmpValue$overlap <- rglapc
	icmpValue$gap <- rggapc
		
	return(icmpValue)
				
}


