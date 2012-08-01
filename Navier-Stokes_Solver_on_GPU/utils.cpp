#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "datadef.h"

/*
 * Write the values of U, V, P, TEMP, FLAG into a binary file
 * for subsequent calculations
 */
void WRITE_bin(REAL **U,REAL **V,REAL **P,REAL **TEMP,int **FLAG,
	int imax,int jmax,char* file)
{
	int i;
	FILE *fp;

	fp = fopen(file, "wb"); 

	fwrite(&imax, sizeof(int), 1, fp);
	fwrite(&jmax, sizeof(int), 1, fp);

	for(i=0;i<=imax+1;i+=1)
		fwrite(U[i], sizeof(REAL), jmax+2, fp);
	for(i=0;i<=imax+1;i+=1)
		fwrite(V[i], sizeof(REAL), jmax+2, fp);
	for(i=0;i<=imax+1;i+=1)
		fwrite(P[i], sizeof(REAL), jmax+2, fp);
	for(i=0;i<=imax+1;i+=1)
		fwrite(TEMP[i], sizeof(REAL), jmax+2, fp);
	for(i=0;i<=imax+1;i+=1)
		fwrite(FLAG[i], sizeof(int), jmax+2, fp);
	fclose(fp);
}

/*
 * Write the values of U, V, P, TEMP, FLAG into a text file
 */
void WRITE_txt(REAL **U,REAL **V,REAL **P,REAL **TEMP,int **FLAG,
	int imax,int jmax,char* file)
{
	int i, j;
	FILE *fp;

	fp = fopen(file, "a");	// append or overwrite

	fprintf(fp, "%d %d\n", imax, jmax);

	for (j=1; j<=jmax; j+=1) {
		for(i=0;i<=imax+1;i+=1) {
			fprintf(fp, "%f ", U[i][j]);
		}
	}
	fprintf(fp, "\n");

	for (j=1; j<=jmax; j+=1) {
		for(i=0;i<=imax+1;i+=1) {
			fprintf(fp, "%f ", V[i][j]);
		}
	}
	fprintf(fp, "\n");

	for (j=1; j<=jmax; j+=1) {
		for(i=0;i<=imax+1;i+=1) {
			fprintf(fp, "%f ", P[i][j]);
		}
	}
	fprintf(fp, "\n");

	for (j=1; j<=jmax; j+=1) {
		for(i=0;i<=imax+1;i+=1) {
			fprintf(fp, "%f ", TEMP[i][j]);
		}
	}
	fprintf(fp, "\n");

	for (j=1; j<=jmax; j+=1) {
		for(i=0;i<=imax+1;i+=1) {
			fprintf(fp, "%d ", FLAG[i][j]);
		}
	}
	fprintf(fp, "\n");
	fprintf(fp, "\n");

	fclose(fp);
}

/*
 * Read initial values from a binary file
 */
int READ_bin(REAL **U,REAL **V,REAL **P,REAL **TEMP,int **FLAG,
	int imax,int jmax,char* file)
{
	int i,j;
	FILE *fp;

	if(strcmp(file, "none") == 0) return(-1);

	if( (fp = fopen(file,"rb")) == NULL){
		printf("Error in READ_bin: File %s cannot be opened!\n", file);
		return(1);
	}

	fread(&i, sizeof(int), 1, fp);
	fread(&j, sizeof(int), 1, fp);

	if(i!=imax || j!=jmax){
		printf("ATTENTION: imax and jmax have wrong values in %s\n",file);
		printf("imax = %d  jmax = %d\n", i, j);
		return(2);
	}

	for(i=0;i<=imax+1;i+=1)
		fread(U[i], sizeof(REAL), jmax+2, fp);
	for(i=0;i<=imax+1;i+=1)
		fread(V[i], sizeof(REAL), jmax+2, fp);
	for(i=0;i<=imax+1;i+=1)
		fread(P[i], sizeof(REAL), jmax+2, fp);
	for(i=0;i<=imax+1;i+=1)
		fread(TEMP[i], sizeof(REAL), jmax+2, fp);
	for(i=0;i<=imax+1;i+=1)
		fread(FLAG[i], sizeof(int), jmax+2, fp);
	fclose(fp);

	return(0);
}

/*
 * Read in an OpenCL kernel file and stores it as a char pointer
 */
char* READ_kernelSource(const char *sourceFilename) {

	FILE *fp;
	int err;
	int size;

	char *source;

	fp = fopen(sourceFilename, "rb");
	if(fp == NULL) {
		printf("Could not open kernel file: %s\n", sourceFilename);
		exit(-1);
	}

	err = fseek(fp, 0, SEEK_END);
	if(err != 0) {
		printf("Error seeking to end of file\n");
		exit(-1);
	}

	size = ftell(fp);
	if(size < 0) {
		printf("Error getting file position\n");
		exit(-1);
	}

	err = fseek(fp, 0, SEEK_SET);
	if(err != 0) {
		printf("Error seeking to start of file\n");
		exit(-1);
	}

	source = (char*)malloc(size+1);
	if(source == NULL) {
		printf("Error allocating %d bytes for the program source\n", size+1);
		exit(-1);
	}

	err = fread(source, 1, size, fp);
	if(err != size) {
		printf("only read %d bytes\n", err);
		exit(0);
	}

	source[size] = '\0';

	return source;
}

void copy_array_real(REAL** src, REAL** dst, int m, int n)
{
	for (int r = 0; r < m; r++) {
		for (int s = 0; s < n; s++) {
			dst[r][s] = src[r][s];
		}
	}  
}

void copy_array_real_2d_to_1d(REAL** src, REAL* dst, int m, int n)
{
	for (int r = 0; r < m; r++) {
		for (int s = 0; s < n; s++) {
			dst[r*n+s] = src[r][s];
		}
	}  
}

void copy_array_int(int** src, int** dst, int m, int n)
{
	for (int r = 0; r < m; r++) {
		for (int s = 0; s < n; s++) {
			dst[r][s] = src[r][s];
		}
	}
}

void copy_array_int_2d_to_1d(int** src, int* dst, int m, int n)
{
	for (int r = 0; r < m; r++) {
		for (int s = 0; s < n; s++) {
			dst[r*n+s] = src[r][s];
		}
	}
}

void print_array(REAL** A, int m, int n)
{
	for (int r = 0; r < m; r++) {
		for (int s = 0; s < n; s++) {
			printf("%6.2f\t", A[r][s]);
		}
		printf("\n");
	}  
}

void print_array_to_file(REAL** A, int m, int n, char* filename)
{
	FILE *fp;
	
	fp = fopen(filename, "w");

	for (int r = 0; r < m; r++) {
		for (int s = 0; s < n; s++) {
			fprintf(fp, "%f\t", A[r][s]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}

void print_1darray_to_file(REAL* A, int m, int n, char* filename)
{
	FILE *fp;

	fp = fopen(filename, "w");

	for (int r = 0; r < m; r++) {
		for (int s = 0; s < n; s++) {
			fprintf(fp, "%f\t", A[r*n+s]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}

void print_array_int_to_file(int** A, int m, int n, char* filename)
{
	FILE *fp;

	fp = fopen(filename, "w");

	for (int r = 0; r < m; r++) {
		for (int s = 0; s < n; s++) {
			fprintf(fp, "%d\t", A[r][s]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}

void print_1darray_int_to_file(int* A, int m, int n, char* filename)
{
	FILE *fp;

	fp = fopen(filename, "w");

	for (int r = 0; r < m; r++) {
		for (int s = 0; s < n; s++) {
			fprintf(fp, "%d\t", A[r*n+s]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}

bool compare_array(REAL** A1, REAL* A2, int m, int n)
{
	REAL sum1, sum2;
	sum1 = 0;
	sum2 = 0;
	for (int r = 0; r < m; r++) {
		for (int s = 0; s < n; s++) {
			sum1 += A1[r][s];
			sum2 += A2[r*n+s];
		}
	}  

	printf("CPU: %10.10f\tGPU: %10.10f\t", sum1, sum2);
	if (fabs(sum1 - sum2) < 10e-5)
		return true;
	else
		return false;
}
