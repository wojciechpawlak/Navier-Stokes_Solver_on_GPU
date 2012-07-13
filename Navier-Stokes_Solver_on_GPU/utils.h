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

void copy_array_real(REAL** src, REAL** dst, int m, int n);

void copy_array_real_2d_to_1d(REAL** src, REAL* dst, int m, int n);

void copy_array_int(int** src, int** dst, int m, int n);

void copy_array_int_2d_to_1d(int** src, int* dst, int m, int n);

void print_array(REAL** A, int m, int n);

void print_array_to_file(REAL** A, int m, int n, char* filename);

void print_1darray_to_file(REAL* A, int m, int n, char* filename);

void print_array_int_to_file(int** A, int m, int n, char* filename);

void print_1darray_int_to_file(int* A, int m, int n, char* filename);

bool compare_array(REAL** A1, REAL* A2, int m, int n);

#endif /* UTILS_H_ */