#pragma once
#ifndef INIT_H_
#define INIT_H_

int READ_PARAMETER(char *Inputfile, char *problem,
	REAL *xlength, REAL *ylength, int *imax, int *jmax,
	REAL *delx, REAL *dely,
	REAL *t_end, REAL *delt, REAL *tau,
	REAL *del_trace, REAL *del_inj, REAL *del_streak, REAL *del_vec,
	char *outputfile, char *tracefile, char *streakfile,
	char *infile, char *outfile,
	int *N, REAL *pos1x, REAL *pos1y, REAL *pos2x, REAL *pos2y,
	int *itermax, REAL *eps, REAL *omg, REAL *alpha, int *p_bound,
	REAL *Re, REAL *Pr, REAL *beta, REAL *GX, REAL *GY,
	REAL *UI, REAL *VI ,REAL *TI,
	int *wW, int *wE, int *wN, int *wS);

void INIT_UVP(char *problem,
	REAL **U, REAL **V, REAL **P, REAL **TEMP, int imax, int jmax,
	REAL UI, REAL VI, REAL TI);

void INIT_FLAG(char *problem, int **FLAG, int imax, int jmax,
	REAL delx, REAL dely, int *ibound);

#endif /* INIT_H_ */
