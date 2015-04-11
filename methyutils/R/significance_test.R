ChisqTest <- function(mat) {
	.Call("_ChisqTest", as.vector(t(mat)), nrow(mat), ncol(mat))
}

MultiChisqTest <- function(data, nrow, ncol) {
	
	# Totoal data size

	data.size <- length(data) 
	
	# Matrix size

	table.size <- nrow * ncol

	if(!(data.size %% table.size == 0)) {
		stop("ncount is not a multiple of contingency table size ( = nrow * ncol).")
	}

	.Call("_MultiChisqTest", as.numeric(data), data.size, as.integer(nrow), as.integer(ncol))
}

FisherTest <- function(vec) {
	
	if(!is.vector(vec) || length(vec) != 4) {
		stop("vec should be a vector of a 2x2 contingency table by row")
	}

	# Total data size

	data.size <- length(vec) 

	.Call("_FisherTest", 
		  as.integer(vec[1]), as.integer(vec[2]), 
		  as.integer(vec[3]), as.integer(vec[4]))
}

MultiFisherTest <- function(vec) {
	if(!is.vector(vec)) {
		stop("vec should be a vector")
	}

	data.size <- length(vec)

	if(data.size %% 4 != 0) {
		stop("the size of vec should be a multiple of 4 which represents a 2x2 contingency table")
	}

	.Call("_MultiFisherTest", as.integer(vec), data.size)
}