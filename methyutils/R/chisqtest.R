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

	.Call("_MultiChisqTest", data, data.size, as.integer(nrow), as.integer(ncol))
}