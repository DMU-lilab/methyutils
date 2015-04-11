#ifndef	FISHER_H
#define FISHER_H

typedef struct 
{
	int n11, n1_, n_1, n;
	double p;
} hgacc_t;

#ifdef __cplusplus
	extern "C" {
#endif

double lbinom(int n, int k);
double hypergeo(int n11, int n1_, int n_1, int n);
double hypergeo_acc(int n11, int n1_, int n_1, int n, hgacc_t *aux);

double fishertest(int n11, int n12, int n21, int n22, double *_left, double *_right, double *two);

#ifdef __cplusplus
	}
#endif

#endif
