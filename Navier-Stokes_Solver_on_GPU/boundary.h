#pragma once
#ifndef BOUNDARY_H_
#define BOUNDARY_H_

void SETBCOND(REAL **U, REAL **V, REAL **P, REAL **TEMP, int **FLAG,
	int imax, int jmax, int wW, int wE, int wN, int wS);

void SETSPECBCOND(char* problem, REAL **U, REAL **V, REAL **P, REAL **TEMP,
	int imax, int jmax, REAL UI, REAL VI);

void SETBCOND_1d(REAL *U, REAL *V, REAL *P, REAL *TEMP, int *FLAG,
	int imax, int jmax, int wW, int wE, int wN, int wS);

void SETSPECBCOND_1d(char* problem, REAL *U, REAL *V, REAL *P,REAL *TEMP,
	int imax, int jmax, REAL UI, REAL VI);

#endif /* BOUNDARY_H_ */