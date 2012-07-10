#pragma OPENCL EXTENSION cl_khr_fp64 : enable 
#include "datadef.h"

__kernel
void COMP_TEMP_kernel(__global REAL *U,
	__global REAL *V,
	__global REAL *TEMP,
	__global REAL *TEMP_new,
	__global int *FLAG,
	int imax,
	int jmax,
	REAL delt,
	REAL delx,
	REAL dely,
	REAL gamma,
	REAL Re,
	REAL Pr)
{
	REAL LAPLT, DUTDX, DVTDY/*, indelx2, indely2*/;

	imax = imax + 2;
	jmax = jmax + 2;

	int i = get_global_id(0);
	int j = get_global_id(1);

	//indelx2 = 1./delx/delx;
	//indely2 = 1./dely/dely;

	if ((i > 0 && i < imax - 1) && (j > 0 && j < jmax - 1)) {
		if( (FLAG[i*jmax + j] & C_F) && (FLAG[i*jmax + j] < C_E) ) {
			LAPLT = (TEMP[(i+1)*jmax + j]-2.0*TEMP[i*jmax + j]+TEMP[(i-1)*jmax + j])*(1./delx/delx) +
				(TEMP[i*jmax + j+1]-2.0*TEMP[i*jmax + j]+TEMP[i*jmax + j-1])*(1./dely/dely);
			DUTDX = ( (U[i*jmax + j]*0.5*(TEMP[i*jmax + j]+TEMP[(i+1)*jmax + j]) -
				U[(i-1)*jmax + j]*0.5*(TEMP[(i-1)*jmax + j]+TEMP[i*jmax + j])) +
				gamma*(fabs(U[i*jmax + j])*0.5*(TEMP[i*jmax + j]-TEMP[(i+1)*jmax + j]) -
				fabs(U[(i-1)*jmax + j])*0.5*(TEMP[(i-1)*jmax + j]-TEMP[i*jmax + j]))
				)/delx;
			DVTDY = ( (V[i*jmax + j]*0.5*(TEMP[i*jmax + j]+TEMP[i*jmax + j+1]) -
				V[i*jmax + j-1]*0.5*(TEMP[i*jmax + j-1]+TEMP[i*jmax + j])) +
				gamma*(fabs(V[i*jmax + j])*0.5*(TEMP[i*jmax + j]-TEMP[i*jmax + j+1]) -
				fabs(V[i*jmax + j-1])*0.5*(TEMP[i*jmax + j-1]-TEMP[i*jmax + j]))
				)/dely;
				
			TEMP_new[i*jmax + j] = TEMP[i*jmax + j]+delt*(LAPLT/Re/Pr - DUTDX - DVTDY);
		}
	}
}

__kernel
void COMP_FG_kernel(__global REAL *U,
					__global REAL *V,
					__global REAL *TEMP,
					__global REAL *F,
					__global REAL *G,
					__global int *FLAG,
					int imax,
					int jmax,
					REAL delt,
					REAL delx,
					REAL dely,
					REAL GX,
					REAL GY,
					REAL gamma,
					REAL Re,
					REAL beta)
{
	REAL DU2DX, DUVDY, DUVDX, DV2DY, LAPLU, LAPLV;

	imax = imax + 2;
	jmax = jmax + 2;

	int i = get_global_id(0);
	int j = get_global_id(1);


	if ((i > 0 && i < imax - 1) && (j > 0 && j < jmax - 1)) {

		/* only if both adjacent cells are fluid cells */
		if ( ((FLAG[i*jmax + j] & C_F) && (FLAG[i*jmax + j] < C_E)) &&
			((FLAG[(i+1)*jmax + j] & C_F) && (FLAG[(i+1)*jmax + j] < C_E)) ) {
				DU2DX = ((U[i*jmax + j]+U[(i+1)*jmax + j])*(U[i*jmax + j]+U[(i+1)*jmax + j])+
					gamma*fabs(U[i*jmax + j]+U[(i+1)*jmax + j])*(U[i*jmax + j]-U[(i+1)*jmax + j])-
					(U[(i-1)*jmax + j]+U[i*jmax + j])*(U[(i-1)*jmax + j]+U[i*jmax + j])-
					gamma*fabs(U[(i-1)*jmax + j]+U[i*jmax + j])*(U[(i-1)*jmax + j]-U[i*jmax + j]))
					/(4.0*delx);
				DUVDY = ((V[i*jmax + j]+V[(i+1)*jmax + j])*(U[i*jmax + j]+U[i*jmax +(j+1)])+
					gamma*fabs(V[i*jmax + j]+V[(i+1)*jmax + j])*(U[i*jmax + j]-U[i*jmax +(j+1)])-
					(V[i*jmax +(j-1)]+V[(i+1)*jmax +(j-1)])*(U[i*jmax +(j-1)]+U[i*jmax + j])-
					gamma*fabs(V[i*jmax +(j-1)]+V[(i+1)*jmax +(j-1)])*(U[i*jmax +(j-1)]-U[i*jmax + j]))
					/(4.0*dely);
				LAPLU = (U[(i+1)*jmax + j]-2.0*U[i*jmax + j]+U[(i-1)*jmax + j])/delx/delx+
					(U[i*jmax +(j+1)]-2.0*U[i*jmax + j]+U[i*jmax +(j-1)])/dely/dely;

				F[i*jmax + j] = U[i*jmax + j]+delt*(LAPLU/Re-DU2DX-DUVDY+GX)
					-delt*beta*GX*(TEMP[i*jmax + j]+TEMP[(i+1)*jmax + j])/2;
		} else {
			F[i*jmax + j] = U[i*jmax + j];
		}

		if( ((FLAG[i*jmax + j] & C_F) && (FLAG[i*jmax + j] < C_E)) &&
			((FLAG[i*jmax +(j+1)] & C_F) && (FLAG[i*jmax +(j+1)] < C_E)) ) {
				DUVDX = ((U[i*jmax + j]+U[i*jmax +(j+1)])*(V[i*jmax + j]+V[(i+1)*jmax + j])+
					gamma*fabs(U[i*jmax + j]+U[i*jmax +(j+1)])*(V[i*jmax + j]-V[(i+1)*jmax + j])-
					(U[(i-1)*jmax + j]+U[(i-1)*jmax +(j+1)])*(V[(i-1)*jmax + j]+V[i*jmax + j])-
					gamma*fabs(U[(i-1)*jmax + j]+U[(i-1)*jmax +(j+1)])*(V[(i-1)*jmax + j]-V[i*jmax + j]))
					/(4.0*delx);
				DV2DY = ((V[i*jmax + j]+V[i*jmax +(j+1)])*(V[i*jmax + j]+V[i*jmax +(j+1)])+
					gamma*fabs(V[i*jmax + j]+V[i*jmax +(j+1)])*(V[i*jmax + j]-V[i*jmax +(j+1)])-
					(V[i*jmax +(j-1)]+V[i*jmax + j])*(V[i*jmax +(j-1)]+V[i*jmax + j])-
					gamma*fabs(V[i*jmax +(j-1)]+V[i*jmax + j])*(V[i*jmax +(j-1)]-V[i*jmax + j]))
					/(4.0*dely);

				LAPLV = (V[(i+1)*jmax + j]-2.0*V[i*jmax + j]+V[(i-1)*jmax + j])/delx/delx+
					(V[i*jmax +(j+1)]-2.0*V[i*jmax + j]+V[i*jmax +(j-1)])/dely/dely;

				G[i*jmax + j] = V[i*jmax + j]+delt*(LAPLV/Re-DUVDX-DV2DY+GY)
					-delt*beta*GY*(TEMP[i*jmax + j]+TEMP[i*jmax +(j+1)])/2;;		      
		} else {
			G[i*jmax + j] = V[i*jmax + j];
		}

	} else if ((i == 0) && (j > 0 && j < jmax-1))  {
		// F at external boundary
		F[j]    =  U[j];

	//} else if ((i == imax-1) && (j > 0 && j < jmax-1))  {
	//	F[(imax-1)*jmax + j] = U[(imax-1)*jmax + j]; //TODO commented out because wrong CPU code

	} else if ((j == 0) && (i > 0 && i < imax-1)) {
		// G at external boundary
		G[i*jmax]    = V[i*jmax];

	//} else if ((j == jmax-1) && (i > 0 && i < imax-1)) {
	//	G[i*jmax + (jmax-1)] = V[i*jmax + (jmax-1)];
	}
}