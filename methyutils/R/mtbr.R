
####required libraries#######
require(DNAStrings)
require(reshape2)
###############CONSTANTS########################
####C context pwm matrix####
kCcontextPWM <- list(
	CG=matrix(c(A=0,C=100,G=0,T=0,A=0, C=0, G=100, T=0), nrow=4, dimnames=list(c("A","C","G","T"), 1:2)),
	CHG=matrix(c(A=0,C=100,G=0,T=0,A=33, C=33, G=0, T=33,A=0, C=0, G=100, T=0), nrow=4, dimnames=list(c("A","C","G","T"), 1:3)),
	CHH=matrix(c(A=0,C=100,G=0,T=0,A=33, C=33, G=0, T=33,A=33, C=33, G=0, T=33), nrow=4, dimnames=list(c("A","C","G","T"), 1:3))
)

# get C context from a DNAString object.
# Arguments:
# Return:

GetCcontextPosition <- function(dna.seq, context="CG", return.type="logicalMarker"){
	
	# check input parameters

	if(class(dna.seq) != "DNAString"){
		stop("dna.seq must to be DNAString")
	}
	
	if(!(context %in% names(kCcontextPWM))){
		 stop("context must be: ", paste("\"", names(kCcontextPWM), "\"", collapse=", ", sep=""))  
	}	
	
	# search posi by matchPWM function
	
	posi.views <- matchPWM(kCcontextPWM[[context]], subject=dna.seq, min.score="99%")
	
 	# return C context according to return.type
	
	if(return.type == "logicalMarker"){
		rt <- logical(length(dna.seq))
		rt[start(posi.views)] <- TRUE
		return(rt)
	} else if(return.type == "integerPosi"){
		return(start(posi.views))
	} else {
		stop("undefined return type. Only logicalMarker and integerPosi are available now.")
	}
} 

# get the records at CG sites from the mtbr that have informations for all C sites

GetCGmtbr <- function(mtbr, cg.position){
	# check the col count
	mtbr.colnames <- c("chrom", "posi", "gBase", "rC", "rT")
	if(ncol(mtbr) != 5){
		stop("colcount must be 5, colnames must be: c(", paste("\"", mtbr.colnames, "\"", collapse=", ", sep=""), ")" )
	} 

	# get mtbr at CG site

	colnames(mtbr) <- mtbr.colnames	
 	cg.position[which(cg.position) + 1] <- TRUE
	cg.mtbr <- mtbr[cg.position[mtbr$posi],]
	
	# put rC and rT for a CG site in one row by melting first and then casting
	
	m <- melt(cg.mtbr, id=c("chrom","posi","gBase"), variable.name="read", value.name="cv") 
	cst <- dcast(m, chrom + posi ~ read + gBase, value.var="cv")
	colnames(cst) <- c("chrom", "posi", "rC_p", "rC_n", "rT_p", "rT_n")
	
	return(cst)
}



# save mtbr to Rdata file
SaveMtbrRdata <- function(mtbr,file.name){
	save(mtbr,file = file.name)	
}
