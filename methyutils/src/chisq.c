#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "chisq.h"

/* 
 * Returns  the  value ln[Gamma(xx)]for xx>0. 
 *
 */	

double gammln(double a)
{
	int i;
	double x, y, tmp, ser;
	const double cof[6] = 
	{ 
		76.18009172947146,
		-86.50532032941677,
		24.01409824083091,
		-1.231739572450155,
		0.1208650973866179e-2,
		-0.5395239384953e-5
	};

	y = x = a;
	tmp = x + 5.5;
	tmp -= (x + 0.5) * log(tmp);
	ser = 1.000000000190015;
	for (i = 0; i <= 5; i++)
	{
		y += 1;
		ser += cof[i] / y;	
	} 

	return -tmp + log(2.5066282746310005 * ser / x);
}

/*
 * Returns the incomplete gamma function Q(a, x)
 * evaluated by its  continued fraction representation as
 * gammcf. Also returns ln Gamma(a) as gln.
 *
 */

void gcf(double *gammcf, double a, double x, double *gln)
{
	int i;
	double an, b, c, d, del, h;

	/* Set up for evaluating continued fraction by modified Lentz's method with b0 = 0. */

	*gln = gammln(a);
	b = x + 1.0 - a;
	c = 1.0 / FPMIN;
	d = 1.0 / b;
	h = d;
	
	for (i = 1; i <= ITMAX; i++) 
	{
		/* Iterate to convergence. */
		
		an = -i * (i - a);
		b += 2.0;
		d = an * d + b;
		if (fabs(d) < FPMIN) 
		{
			d = FPMIN;
		}
		c = b + an / c;

		if (fabs(c) < FPMIN) 
		{
			c = FPMIN;
		}
		
		d = 1.0 / d;
		del = d * c;
		h *= del;
		if (fabs(del-1.0) < EPS)
		{
			break;	
		} 
	}
	
	if (i > ITMAX)
	{
		fprintf(stderr, "[Error]: a too large, ITMAX too small in gcf.\n");	
	} 
	
	*gammcf = exp(-x + a * log(x) - (*gln)) * h; /* Put  factors in  front. */
}


/*
 * Returns the incomplete gamma function P(a, x) evaluated by its series 
 * representation as gamser. Also  returns lnGamma(a) as gln.
 *
 */ 

void gser(double *gamser, double a, double x, double *gln)
{
	int n;
	double sum, del, ap;
	*gln = gammln(a);
	if (x <= 0.0) 
	{
		if (x < 0.0) 
		{
			fprintf(stderr, "[Error]: x less than 0 in routine gser.\n");
		}
		
		*gamser=0.0;
		
		return;
	} 
	else 
	{
		ap = a;
		del = sum = 1.0 / a;
		for (n = 1; n <= ITMAX; n++) 
		{
			ap += 1.0;
			del *= x / ap;
			sum += del;
			if (fabs(del) < fabs(sum) * EPS) 
			{
				*gamser = sum * exp(-x + a * log(x) - (*gln));
				return;
			}
		}

		fprintf(stderr, "[Error]: a too large, ITMAX too small in routine gser.\n");
		return;
	}	
}

/* 
 * Returns  the  incomplete  gamma  function P(a,x)
 */

double gammp(double a, double x)
{
	double gamser, gammcf, gln;

	if (x < 0.0 || a <= 0.0)
	{
		fprintf(stderr, "[Error]: Invalid arguments in routine gammp.\n");
		return -1.0;	
	} 

	if (x < (a + 1.0)) 
	{
		gser(&gamser, a, x, &gln); /* Use  the  series representation. */
		return gamser;
	} 
	else 
	{
		gcf(&gammcf, a, x, &gln); /* Use the continued  fraction representation */
		return 1.0 - gammcf; /* and  take  its  complement. */
		
	}
}

/* 
 * Returns  the  incomplete  gamma  function Q(a, x) = 1 − P(a, x).
 */

double gammq(double a, double x)
{
	double gamser, gammcf, gln;
	
	if (x < 0.0 || a <= 0.0)
	{
		fprintf(stderr, "[Error]: Invalid arguments in routine gammq.\n");
		return -1.0;
	}

	if (x < (a+1.0)) 
	{
		gser(&gamser, a, x, &gln); /* Use  the  series representation */
		return 1.0 - gamser; /* and  take  its  complement. */
	} 
	else 
	{
		gcf(&gammcf, a, x, &gln); /* Use the continued fraction representation. */
		return gammcf;
	}
}


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