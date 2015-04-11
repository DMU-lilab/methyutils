#include <R.h>
#include <Rinternals.h>
#include <R_ext/Print.h>

#include "chisq.h"
#include "fisher.h"

SEXP _ChisqTest(SEXP vect, SEXP nrow, SEXP ncol)
{
	SEXP pvalue = PROTECT(allocVector(REALSXP, 1));
	
	REAL(pvalue)[0] = tchisqtest(REAL(vect), INTEGER(nrow)[0], INTEGER(ncol)[0]);
	
	UNPROTECT(1);

	return pvalue;
}

SEXP _MultiChisqTest(SEXP vect, SEXP ncount, SEXP nrow, SEXP ncol)
{
	int i;
	double *pData = NULL;
	double *pPvalue = NULL;

	int nt = INTEGER(ncount)[0];
	int nr = INTEGER(nrow)[0];
	int nc = INTEGER(ncol)[0];
	int n = (int)(nt / (nr * nc));

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

SEXP _FisherTest(SEXP d11, SEXP d12, SEXP d21, SEXP d22)
{
	double left, right, twotail, prob;
	SEXP pvalue = PROTECT(allocVector(REALSXP, 1));
	
	prob = fishertest(INTEGER(d11)[0], INTEGER(d12)[0], 
					  INTEGER(d21)[0], INTEGER(d22)[0],
					  &left, &right, &twotail);
	REAL(pvalue)[0] = twotail;

	UNPROTECT(1);

	return pvalue;
}

SEXP _MultiFisherTest(SEXP vect, SEXP ncount)
{
	int i;
	int *pData = NULL;
	double *pPvalue = NULL;
	double left, right, twotail, prob;
	int n = (int)(INTEGER(ncount)[0] / 4);

	SEXP pvalues = PROTECT(allocVector(REALSXP, n));
	for(i = 0, pData = INTEGER(vect), pPvalue = REAL(pvalues);
		i < n;
		i++, pData += 4, pPvalue++)
	{
		prob = fishertest(*pData, *(pData + 1), *(pData + 2), *(pData + 3),
				          &left, &right, &twotail);
		*pPvalue = twotail;
	}

	UNPROTECT(1);

	return pvalues;
}