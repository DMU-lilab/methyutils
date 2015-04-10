#include <R.h>
#include <Rinternals.h>

#include "chisq.h"

SEXP _ChisqTest(SEXP vect, SEXP nrow, SEXP ncol)
{
	SEXP pvalue = PROTECT(allocVector(REALSXP, 1));
	
	REAL(pvalue)[0] = tchisqtest(REAL(vect), INTEGER(nrow)[0], INTEGER(ncol)[0]);
	
	UNPROTECT(1);

	return pvalue;
}

SEXP _MultiChisqTest(SEXP vect, SEXP ncount, SEXP nrow, SEXP ncol)
{
	int i, j;
	double *pData = NULL;
	double *pPvalue = NULL;

	int nt = INTEGER(ncount)[0];
	int nr = INTEGER(nrow)[0];
	int nc = INTEGER(ncol)[0];
	int n = (nt / (nr * nc));

	SEXP pvalues = PROTECT(allocVector(REALSXP, n));

	for(i = 0, pData = REAL(vect), pPvalue = REAL(pvalues); 
		i < n; 
		i++, pData += nr * nc, pPvalue++)
	{
		*pPvalue = tchisqtest(pData, nr, nc);
	}

	UNPROTECT(1);

	return pvalues;
}