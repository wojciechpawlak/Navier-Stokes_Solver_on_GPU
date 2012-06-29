#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "datadef.h"
#include "init.h"
#include "boundary.h"
#include "uvp.h"
#include "visual.h"
#include "surface.h"

/*----------------------------------------------------------------- */
/*                 M A I N     P R O G R A M                        */
/*------------------------------------------------------------------*/
int main(int argc, char *Inputfile[])
{
	char problem[30];
	char infile[30], outfile[30];
	REAL xlength, ylength;
	int  imax, jmax;
	REAL delx, dely;
	REAL t_end, delt, tau;
	REAL del_trace, del_inj, del_streak, del_vec;
	char vecfile[30], tracefile[30], streakfile[30];
	int N;
	REAL pos1x, pos2x, pos1y, pos2y;
	int itermax;
	REAL eps, omg, gamma;
	int  p_bound;
	REAL Re, Pr, beta, GX, GY, UI, VI, TI;
	int wW, wE, wN, wS;
	int itersor=0, write;

	REAL t;
	REAL res;
	REAL **U, **V, **P, **PSI, **ZETA, **RHS, **F, **G, **TEMP, **HEAT;
	int  **FLAG;
	int  ppc, ifull=0, isurf=0, ibound=0;
	struct particleline *Particlelines;
	int init_case, cycle;

	//REAL xlength, ylength, delx, dely, t_end, delt, tau, t;
	//REAL del_trace, del_inj, del_streak, del_vec;
	//REAL pos1x, pos2x, pos1y, pos2y;
	//int  imax, jmax, wW, wE, wN, wS, itermax, itersor=0, write, N;
	//REAL Re, Pr, GX, GY, UI, VI, TI, beta;
	//REAL eps,omg, gamma, res;
	//int  p_bound;
	//REAL **U, **V, **P, **PSI, **ZETA, **RHS, **F, **G, **TEMP, **HEAT;
	//int  **FLAG;
	//int  ppc, ifull=0, isurf=0, ibound=0;
	//char problem[30];
	//char vecfile[30], tracefile[30], streakfile[30];
	//char infile[30], outfile[30];
	//struct particleline *Particlelines;
	//int init_case, cycle;

	/* READ the parameters of the problem.                */
	/* Stop if problem type or inputfile are not defined  */       
	/*----------------------------------------------------*/
	if (READ_PARAMETER(Inputfile[1],problem,
		&xlength, &ylength, &imax, &jmax, &delx, &dely,
		&t_end, &delt, &tau, 
		&del_trace, &del_inj, &del_streak, &del_vec,
		vecfile,tracefile,streakfile,
		infile, outfile,
		&N, &pos1x, &pos1y, &pos2x, &pos2y,
		&itermax,&eps,&omg,&gamma,&p_bound,
		&Re, &Pr, &beta, &GX, &GY, &UI, &VI, &TI,
		&wW, &wE, &wN, &wS) != 0 ) {
		
		return(1); 
	}

	/* Allocate memory for the arrays */
	/*--------------------------------*/
	U    = RMATRIX(0, imax+1, 0, jmax+1);
	V    = RMATRIX(0, imax+1, 0, jmax+1);
	F    = RMATRIX(0, imax+1, 0, jmax+1);
	G    = RMATRIX(0, imax+1, 0, jmax+1);
	P    = RMATRIX(0, imax+1, 0, jmax+1);
	TEMP = RMATRIX(0, imax+1, 0, jmax+1);
	PSI  = RMATRIX(0, imax,	  0, jmax);
	ZETA = RMATRIX(1, imax-1, 1, jmax-1);
	HEAT = RMATRIX(0, imax,   0, jmax);
	RHS  = RMATRIX(0, imax+1, 0, jmax+1); 
	FLAG = IMATRIX(0, imax+1, 0, jmax+1);
	ppc  = 4;                             

	/* Read initial values from file "infile" */
	/*----------------------------------------*/
	init_case = READ_bin(U, V, P, TEMP, FLAG, imax, jmax, infile); 

	if (init_case > 0) return(1);	/* Error while reading "infile" */
	if (init_case < 0) {            /* Set initial values if        */
									/* "infile" is not specified    */
		INIT_UVP(problem, U, V, P, TEMP, imax, jmax, UI, VI, TI); 
		INIT_FLAG(problem, FLAG, imax, jmax, delx, dely, &ibound);
	}

	/* Initialize particles for streaklines or particle tracing */
	/*----------------------------------------------------------*/
	if (strcmp(streakfile, "none") || strcmp(tracefile, "none")) {
		Particlelines = SET_PARTICLES(N, pos1x, pos1y, pos2x, pos2y);
	}

	/* Initialize particles for free boundary problems */      
	/*-------------------------------------------------*/
	if (!strcmp(problem, "drop") || !strcmp(problem, "dam")) {
		Particlelines = INIT_PARTICLES(&N, imax, jmax, delx, dely,
			ppc, problem, U, V);
	}

	SETBCOND(U, V, P, TEMP, FLAG, imax, jmax, wW, wE, wN, wS);
	SETSPECBCOND(problem, U, V, P, TEMP, imax, jmax, UI, VI);


	/* t i m e    l o o p */
	/*--------------------*/
	for (t=0.0, cycle=0; t < t_end; t+=delt, cycle++) {
		COMP_delt(&delt, t, imax, jmax, delx, dely, U, V, Re, Pr, tau, &write,
			del_trace, del_inj, del_streak, del_vec);    

		/* Determine fluid cells for free boundary problems */
		/* and set boundary values at free surface          */
		/*--------------------------------------------------*/
		if (!strcmp(problem, "drop") || !strcmp(problem, "dam") ||
			!strcmp(problem, "molding") || !strcmp(problem, "wave")) {
			MARK_CELLS(FLAG, imax, jmax, delx, dely, &ifull, &isurf,
				N, Particlelines);
			SET_UVP_SURFACE(U, V, P, FLAG, GX, GY,
				imax, jmax, Re, delx, dely, delt);
		} else {
			ifull = imax*jmax-ibound;
		}

		/* Compute new temperature */
		/*-------------------------*/
		COMP_TEMP(U, V, TEMP, FLAG, imax, jmax,
			delt, delx, dely, gamma, Re, Pr);

		/* Compute tentative velocity field (F, G) */
		/*----------------------------------------*/
		COMP_FG(U, V, TEMP, F, G, FLAG, imax, jmax,
			delt, delx, dely, GX, GY, gamma, Re, beta);

		/* Compute right hand side for pressure equation */
		/*-----------------------------------------------*/
		COMP_RHS(F, G, RHS, FLAG, imax, jmax, delt, delx, dely);

		/* Solve the pressure equation by successive over relaxation */
		/*-----------------------------------------------------------*/
		if (ifull > 0) {
			itersor = POISSON(P, RHS, FLAG, imax, jmax, delx, dely,
			eps, itermax, omg, &res, ifull, p_bound);
		}

		printf("t_end= %1.5g, t= %1.3e, delt= %1.1e, iterations %3d, res: %e, F-cells: %d, S-cells: %d, B-cells: %d\n",
			t_end, t+delt, delt, itersor, res, ifull, isurf, ibound);  

		/* Compute the new velocity field */
		/*--------------------------------*/
		ADAP_UV(U, V, F, G, P, FLAG, imax, jmax, delt, delx, dely);

		/* Set boundary conditions */
		/*-------------------------*/
		SETBCOND(U, V, P, TEMP, FLAG, imax, jmax, wW, wE, wN, wS);
		/* Set special boundary conditions */
		/* Overwrite preset default values */
		/*---------------------------------*/
		SETSPECBCOND(problem, U, V, P, TEMP, imax, jmax, UI, VI);

		if (!strcmp(problem, "drop") || !strcmp(problem, "dam") ||
			!strcmp(problem, "molding") || !strcmp(problem, "wave")) {
			SET_UVP_SURFACE(U, V, P, FLAG, GX, GY, imax, jmax, Re, delx, dely, delt);
		}

		/* Write data for visualization */
		/*------------------------------*/
		if ((write & 8) && strcmp(vecfile, "none")) {     
			COMPPSIZETA(U, V, PSI, ZETA, FLAG, imax, jmax, delx, dely);
			COMP_HEAT(U, V, TEMP, HEAT, FLAG, Re, Pr, imax, jmax, delx, dely);
			OUTPUTVEC_bin(U, V, P, TEMP, PSI, ZETA, HEAT, FLAG, xlength, ylength,
				imax, jmax, vecfile);
		}
		if ((write & 8) && strcmp(outfile, "none")) {
			WRITE_bin(U, V, P, TEMP, FLAG, imax, jmax, outfile);
		}
		if (strcmp(tracefile, "none")) {
			PARTICLE_TRACING(tracefile, t, imax, jmax, delx, dely, delt,
				U, V, FLAG,	N, Particlelines, write);
		}
		if (strcmp(streakfile, "none")) {
			STREAKLINES(streakfile, write, imax, jmax, delx, dely, delt, t,
				U, V, FLAG, N, Particlelines);
		}
	}           
	/* e n d   o f   t i m e   l o o p */

	if (strcmp(vecfile,"none"))
	{     
		COMPPSIZETA(U, V, PSI, ZETA, FLAG, imax, jmax, delx, dely);
		COMP_HEAT(U, V, TEMP, HEAT, FLAG, Re, Pr, imax, jmax, delx, dely);
		OUTPUTVEC_bin(U, V, P, TEMP, PSI, ZETA, HEAT, FLAG, xlength, ylength,
			imax, jmax, vecfile);
	}
	if (strcmp(outfile,"none")) {
		WRITE_bin(U, V, P, TEMP, FLAG, imax, jmax, outfile);
	}

	/* free memory */
	/*-------------*/
	FREE_RMATRIX(U,		0, imax+1,	0, jmax+1);
	FREE_RMATRIX(V,		0, imax+1,	0, jmax+1);
	FREE_RMATRIX(F,		0, imax+1,	0, jmax+1);
	FREE_RMATRIX(G,		0, imax+1,	0, jmax+1);
	FREE_RMATRIX(P,		0, imax+1,	0, jmax+1);
	FREE_RMATRIX(TEMP,	0, imax+1,	0, jmax+1);
	FREE_RMATRIX(PSI,	0, imax,	0, jmax);
	FREE_RMATRIX(ZETA,	1, imax-1,	1, jmax-1);
	FREE_RMATRIX(HEAT,	0, imax,	0, jmax);
	FREE_RMATRIX(RHS,	0, imax+1,	0, jmax+1);
	FREE_IMATRIX(FLAG,	0, imax+1,	0, jmax+1);

	printf("End of program\n");
	return(0);
}
