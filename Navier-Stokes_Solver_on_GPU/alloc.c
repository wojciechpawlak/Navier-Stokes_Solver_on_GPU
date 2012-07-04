#include <stdio.h>
#include <stdlib.h>

#include "datadef.h"

/*
 * Allocates memory for a [nrl,nrh]x[ncl,nch]-array of REAL-type
 */
REAL **RMATRIX(int nrl,int nrh,int ncl,int nch)
{
	int i;
	REAL **m;
	if((m = (REAL**) malloc((unsigned) (nrh-nrl+1)*sizeof(double*))) == NULL)
	{	
		printf("no memory\n");
		exit(0);
	}
	m -= nrl;
	for(i=nrl;i<=nrh;i++)
	{
		if((m[i] = (REAL*) malloc((unsigned) (nch-ncl+1)*sizeof(double)))==NULL)
		{	
			printf("no memory\n");
			exit(0);
		}
		m[i] -= ncl;
	}
	return m;
} 

/*
 * Frees the memory of an array allocated with RMATRIX
 */
void FREE_RMATRIX(REAL** m,int nrl,int nrh,int ncl,int nch)
{
	int i;
	for (i=nrh;i>=nrl;i--) free((void*) (m[i]+ncl));
	free((char*) (m+nrl));
}

/*
 * Allocates memory for a [nrl,nrh]x[ncl,nch]-array of integer-type
 */
int **IMATRIX(int nrl,int nrh,int ncl,int nch)
{
	int i;
	int **m;
	if((m = (int**) malloc((unsigned) (nrh-nrl+1)*sizeof(int*))) == NULL)
	{	
		printf("no memory\n");
		exit(0);
	}
	m -= nrl;
	for(i=nrl;i<=nrh;i++)
	{
		if((m[i] = (int*) malloc((unsigned) (nch-ncl+1)*sizeof(int)))==NULL)
		{	
			printf("no memory\n");
			exit(0);
		}
		m[i] -= ncl;
	}
	return m;
} 

/*
 * Frees the memory of an array allocated with IMATRIX
 */
void FREE_IMATRIX(int** m,int nrl,int nrh,int ncl,int nch)
{
	int i;
	for (i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}