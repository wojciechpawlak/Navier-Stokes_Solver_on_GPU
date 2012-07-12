#pragma OPENCL EXTENSION cl_khr_fp64 : enable 
#include "datadef.h"

__kernel
void COMP_TEMP_kernel(	__global REAL *U,
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
	REAL LAPLT, DUTDX, DVTDY;

	imax = imax + 2;
	jmax = jmax + 2;

	int i = get_global_id(0);
	int j = get_global_id(1);

	if ((i > 0 && i < imax - 1) && (j > 0 && j < jmax - 1)) {
		if ( (FLAG[i*jmax + j] & C_F) && (FLAG[i*jmax + j] < C_E) ) {

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

		// only if both adjacent cells are fluid cells
		if (	((FLAG[i*jmax + j] & C_F) && (FLAG[i*jmax + j] < C_E)) &&
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

		if (	((FLAG[i*jmax + j] & C_F) && (FLAG[i*jmax + j] < C_E)) &&
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

__kernel
void COMP_RHS_kernel(	__global REAL *F,
						__global REAL *G,
						__global REAL *RHS,
						__global int *FLAG,
						int imax,
						int jmax,
						REAL delt,
						REAL delx,
						REAL dely)
{
	imax = imax + 2;
	jmax = jmax + 2;
	
	int i = get_global_id(0);
	int j = get_global_id(1);

	if ((i > 0 && i < imax - 1) && (j > 0 && j < jmax - 1)) {
		if ((FLAG[i*jmax + j] & C_F) && (FLAG[i*jmax + j] < 0x0100)) {
			// only for fluid and non-surface cells
			RHS[i*jmax + j] = ((F[i*jmax + j]-F[(i-1)*jmax + j])/delx
				+ (G[i*jmax + j]-G[i*jmax + j-1])/dely)/delt;
		}
	}
}

__kernel
void POISSON_p0_kernel(	__global REAL *P,
						__global int *FLAG,
						int imax,
						int jmax,
						__global REAL *p0)
{
	imax = imax + 2;
	jmax = jmax + 2;
	
	int i = get_global_id(0);
	int j = get_global_id(1);

	if ((i > 0 && i < imax - 1) && (j > 0 && j < jmax - 1)) {
		if (FLAG[i*jmax + j] & C_F) {
			*p0 += P[i*jmax + j]*P[i*jmax + j];
		}
	}
}

//__kernel
//void POISSON_1_relaxation_kernel(	__global REAL *P,
//									__global REAL *RHS,
//									__global int *FLAG,
//									int imax,
//									int jmax,
//									REAL delx,
//									REAL dely,
//									REAL omg)
//{
//	REAL rdx2, rdy2;
//	REAL add, beta_2, beta_mod;
//
//	imax = imax + 2;
//	jmax = jmax + 2;
//
//	int i = get_global_id(0);
//	int j = get_global_id(1);
//	
//	rdx2 = 1./delx/delx;
//	rdy2 = 1./dely/dely;
//	beta_2 = -omg/(2.0*(rdx2+rdy2));
//
//	if ((i > 0 && i < imax - 1) && (j > 0 && j < jmax - 1)) { 
//		/* five point star for interior fluid cells */
//		if (FLAG[i*jmax + j] == 0x001f) {
//			P[i*jmax + j] = (1.-omg)*P[i*jmax + j] - 
//				beta_2*((P[(i+1)*jmax + j]+P[(i-1)*jmax + j])*rdx2 +
//				(P[i*jmax + j+1]+P[i*jmax + j-1])*rdy2 - RHS[i*jmax + j]);
//		}
//		/* modified star near boundary */
//		else if ((FLAG[i*jmax + j] & C_F) && (FLAG[i*jmax + j] < 0x0100)) { 
//			beta_mod = -omg/((eps_E+eps_W)*rdx2+(eps_N+eps_S)*rdy2);
//			
//			P[i*jmax + j] = (1.-omg)*P[i*jmax + j] -
//				beta_mod*( (eps_E*P[(i+1)*jmax + j]+eps_W*P[(i-1)*jmax + j])*rdx2 +
//				(eps_N*P[i*jmax + j+1]+eps_S*P[i*jmax + j-1])*rdy2 - RHS[i*jmax + j]);
//		}
//	}
//}

__kernel
void POISSON_1_comp_res_kernel(	__global REAL *P,
								__global REAL *RHS,
								__global int *FLAG,
								int imax,
								int jmax,
								REAL delx,
								REAL dely,
								__global REAL *res)
{
	imax = imax + 2;
	jmax = jmax + 2;

	int i = get_global_id(0);
	int j = get_global_id(1);

	REAL add;
	REAL rdx2, rdy2;

	rdx2 = 1./delx/delx;
	rdy2 = 1./dely/dely;
		
	if ((i > 0 && i < imax - 1) && (j > 0 && j < jmax - 1)) { 
		/* only fluid cells */
		if ((FLAG[i*jmax + j] & C_F) && (FLAG[i*jmax + j] < 0x0100)) {
			add = (	eps_E_d*(P[(i+1)*jmax + j]-P[i*jmax + j]) - 
					eps_W_d*(P[i*jmax + j]-P[(i-1)*jmax + j])) * rdx2
					+ (	eps_N_d*(P[i*jmax + j+1]-P[i*jmax + j]) -
						eps_S_d*(P[i*jmax + j]-P[i*jmax + j-1])) * rdy2
					- RHS[i*jmax + j];

			*res += add*add;
		}
	}
}

//__kernel
//void POISSON_2_copy_boundary_kernel(__global REAL *P,
//									__global int *FLAG,
//									int imax,
//									int jmax)
//{
//	imax = imax + 2;
//	jmax = jmax + 2;
//
//	int i = get_global_id(0);
//	int j = get_global_id(1);
//
//	/* copy values at external boundary */
//	/*----------------------------------*/
//	if ((j == 0) && (i > 0 && i < imax-1)) {
//		P[i*jmax]      = P[i*jmax + 1];
//		P[i*jmax + (jmax-1)] = P[i*jmax + (jmax-2)];
//	} else if ((i == 0) && (j > 0 && j < jmax-1))  {
//		P[j]      = P[jmax + j];
//		P[(imax-1)*jmax + j] = P[(imax-2)*jmax + j];
//	}
//	/* and at interior boundary cells */
//	/*--------------------------------*/
//	else if ((i > 0 && i < imax - 1) && (j > 0 && j < jmax - 1)) {
//		if (FLAG[i*jmax + j] >= B_N && FLAG[i*jmax + j] <= B_SO) {
//			switch (FLAG[i*jmax + j]) {
//			case B_N:	{ P[i*jmax + j] = P[i*jmax + j+1];							break; }
//			case B_O:	{ P[i*jmax + j] = P[(i+1)*jmax + j];						break; }
//			case B_S:	{ P[i*jmax + j] = P[i*jmax + j-1];							break; } 
//			case B_W:	{ P[i*jmax + j] = P[(i-1)*jmax + j];						break; }
//			case B_NO:	{ P[i*jmax + j] = 0.5*(P[i*jmax + j+1]+P[(i+1)*jmax + j]);	break; }
//			case B_SO:	{ P[i*jmax + j] = 0.5*(P[i*jmax + j-1]+P[(i+1)*jmax + j]);	break; }
//			case B_SW:	{ P[i*jmax + j] = 0.5*(P[i*jmax + j-1]+P[(i-1)*jmax + j]);	break; }
//			case B_NW:	{ P[i*jmax + j] = 0.5*(P[i*jmax + j+1]+P[(i-1)*jmax + j]);	break; }
//			default: break;
//			}
//		}
//	}
//			
//}

//__kernel
//void POISSON_2_relaxation_kernel(	__global REAL *P,
//									__global REAL *RHS,
//									__global int *FLAG,
//									int imax,
//									int jmax,
//									REAL delx,
//									REAL dely,
//									REAL omg)
//{
//	imax = imax + 2;
//	jmax = jmax + 2;
//
//	int i = get_global_id(0);
//	int j = get_global_id(1);
//
//	REAL rdx2, rdy2;
//	REAL beta_2;
//
//	rdx2 = 1./delx/delx;
//	rdy2 = 1./dely/dely;
//	beta_2 = -omg/(2.0*(rdx2+rdy2));
//	
//	if ((i > 0 && i < imax - 1) && (j > 0 && j < jmax - 1)) {	
//		if ((FLAG[i*jmax + j] & C_F) && (FLAG[i*jmax + j] < 0x0100)) {
//			P[i*jmax + j] = (1.-omg)*P[i*jmax + j] - 
//				beta_2*((P[(i+1)*jmax + j]+P[(i-1)*jmax + j])*rdx2 +
//				(P[i*jmax + j+1]+P[i*jmax + j-1])*rdy2 - RHS[i*jmax + j]);
//		}
//	}
//}

__kernel
void POISSON_2_comp_res_kernel(	__global REAL *P,
								__global REAL *RHS,
								__global int *FLAG,
								int imax,
								int jmax,
								REAL delx,
								REAL dely,
								__global REAL *res)
{
	imax = imax + 2;
	jmax = jmax + 2;

	int i = get_global_id(0);
	int j = get_global_id(1);

	REAL add;
	REAL rdx2, rdy2;

	rdx2 = 1./delx/delx;
	rdy2 = 1./dely/dely;

	if ((i > 0 && i < imax - 1) && (j > 0 && j < jmax - 1)) {
		/* only fluid cells */
		if ((FLAG[i*jmax + j] & C_F) && (FLAG[i*jmax + j] < 0x0100)) {
			add =	(P[(i+1)*jmax + j]-2*P[i*jmax + j]+P[(i-1)*jmax + j])*rdx2
					+ (P[i*jmax + j+1]-2*P[i*jmax + j]+P[i*jmax + j-1])*rdy2
					- RHS[i*jmax + j];
						
			*res += add*add;
		}
	}
}

__kernel
void ADAP_UV_kernel(__global REAL *U,
					__global REAL *V,
					__global REAL *F,
					__global REAL *G,
					__global REAL *P,
					__global int *FLAG,
					int imax,
					int jmax,
					REAL delt,
					REAL delx,
					REAL dely)
{
	imax = imax + 2;
	jmax = jmax + 2;

	int i = get_global_id(0);
	int j = get_global_id(1);

	if ((i > 0 && i < imax - 2) && (j > 0 && j < jmax - 1)) { // imax - 2
		// only if both adjacent cells are fluid cells
		if (((FLAG[i*jmax + j] & C_F) && (FLAG[i*jmax + j] < C_E)) &&
			((FLAG[(i+1)*jmax + j] & C_F) && (FLAG[(i+1)*jmax + j] < C_E))) {

			U[i*jmax + j] = F[i*jmax + j]-(P[(i+1)*jmax + j]-P[i*jmax + j])*delt/delx;
		}
	} 
	
	if ((i > 0 && i < imax - 1) && (j > 0 && j < jmax - 2)) { // jmax - 2
		// only if both adjacent cells are fluid cells
		if (((FLAG[i*jmax + j] & C_F) && (FLAG[i*jmax + j] < C_E)) &&
			((FLAG[i*jmax + j+1] & C_F) && (FLAG[i*jmax + j+1] < C_E))) {

			V[i*jmax + j] = G[i*jmax + j]-(P[i*jmax + j+1]-P[i*jmax + j])*delt/dely;
		}
	}

}