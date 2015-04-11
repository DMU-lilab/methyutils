#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mathutils.h"
#include "chisq.h"

double pchisq(double z, int df)
{
	return gammp(df / 2.0, z / 2.0);
}

double qchisq(double z, int df)
{
	return gammq(df / 2.0, z / 2.0);
}

/* 
 * Chi-Square test for table , return P-value
 * Whe P-value < 0, errors encountered.
 * 
 */

double tchisqtest(double *table, int nrow, int ncol)
{
	int i, j;
	int ncount = 0; 
	int df = 0;
	double sum = 0.0;
	double *rowsums = NULL;
	double *colsums = NULL;
	double *extable = NULL;
	double chisq = 0.0;
	double pvalue = -1.0;
	double value = 0;

	if(table == NULL)
	{
		return -1.0;
	}

	ncount = nrow * ncol;
	df = (nrow - 1) * (ncol - 1);

	if((rowsums = (double *)malloc(sizeof(double) * nrow)) == NULL)
	{
		goto error;
	}

  	if((colsums = (double *)malloc(sizeof(double) * ncol)) == NULL)
  	{
  		goto error;
  	}

  	if((extable = (double *)malloc(sizeof(double) * ncount)) == NULL)
  	{
  		goto error;
  	}

  	memset(rowsums, 0, sizeof(double) * nrow);
  	memset(colsums, 0, sizeof(double) * ncol);
  	memset(extable, 0, sizeof(double) * ncount);

  	/* calculate row & col sums */

	for(i = 0; i < nrow; i++)
	{
		for(j = 0; j < ncol; j++)
		{
			value = *(table + i * ncol + j);

			/* sum rows */
  	
			*(rowsums + i) += value;
			
			/* sum cols */

			*(colsums + j) += value;

			sum += value;
		}
	}

	if(sum == 0) 
	{
		goto error;
	}

	/* calculate expected value table */

	for(i = 0; i < nrow; i++)
	{
		for(j = 0; j < ncol; j++)
		{
			*(extable + i * ncol + j) = (*(rowsums + i)) * (*(colsums + j)) / sum;
		}
	}

	/* calculate Chi^2 value */

	for(i = 0; i < ncount; i++)
	{
		if(*(extable + i) == 0)
		{
			goto error;
		}

		chisq += (*(table + i) - *(extable + i)) * (*(table + i) - *(extable + i)) / *(extable + i);
	}

	/* calculate p value */

	pvalue = 1 - pchisq(chisq, df);

error:

	/* clean up memory */
	
	free(rowsums);
	free(colsums);
	free(extable);

	return pvalue;
}