#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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