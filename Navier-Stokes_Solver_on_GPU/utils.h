#pragma once
#ifndef UTILS_H_
#define UTILS_H_

void WRITE_bin(REAL **U, REAL **V, REAL **P, REAL **TEMP, int **FLAG,
	int imax, int jmax, char* file);

void WRITE_txt(REAL **U, REAL **V, REAL **P, REAL **TEMP, int **FLAG,
	int imax, int jmax, char* file);

int READ_bin(REAL **U, REAL **V, REAL **P, REAL **TEMP, int **FLAG,
	int imax, int jmax, char* file);

char* READ_kernelSource(const char *sourceFilename);

#endif /* UTILS_H_ */