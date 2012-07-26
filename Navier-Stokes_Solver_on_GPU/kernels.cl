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

//__kernel
//void POISSON_1_comp_res_kernel(	__global REAL *P,
//								__global REAL *RHS,
//								__global int *FLAG,
//								int imax,
//								int jmax,
//								REAL delx,
//								REAL dely,
//								__global REAL *res)
//{
//	imax = imax + 2;
//	jmax = jmax + 2;
//
//	int i = get_global_id(0);
//	int j = get_global_id(1);
//
//	REAL add;
//	REAL rdx2, rdy2;
//
//	rdx2 = 1./delx/delx;
//	rdy2 = 1./dely/dely;
//		
//	if ((i > 0 && i < imax - 1) && (j > 0 && j < jmax - 1)) { 
//		/* only fluid cells */
//		if ((FLAG[i*jmax + j] & C_F) && (FLAG[i*jmax + j] < 0x0100)) {
//			add = (	eps_E_d*(P[(i+1)*jmax + j]-P[i*jmax + j]) - 
//					eps_W_d*(P[i*jmax + j]-P[(i-1)*jmax + j])) * rdx2
//					+ (	eps_N_d*(P[i*jmax + j+1]-P[i*jmax + j]) -
//						eps_S_d*(P[i*jmax + j]-P[i*jmax + j-1])) * rdy2
//					- RHS[i*jmax + j];
//
//			*res += add*add;
//		}
//	}
//}

__kernel
void POISSON_2_copy_boundary_kernel(__global REAL *P,
									__global int *FLAG,
									int imax,
									int jmax)
{
	imax = imax + 2;
	jmax = jmax + 2;

	int i = get_global_id(0);
	int j = get_global_id(1);

	/* copy values at external boundary */
	/*----------------------------------*/
	if ((j == 0) && (i > 0 && i < imax-1)) {
		P[i*jmax]      = P[i*jmax + 1];
		P[i*jmax + (jmax-1)] = P[i*jmax + (jmax-2)];
	} else if ((i == 0) && (j > 0 && j < jmax-1))  {
		P[j]      = P[jmax + j];
		P[(imax-1)*jmax + j] = P[(imax-2)*jmax + j];
	}
	/* and at interior boundary cells */
	/*--------------------------------*/
	else if ((i > 0 && i < imax - 1) && (j > 0 && j < jmax - 1)) {
		if (FLAG[i*jmax + j] >= B_N && FLAG[i*jmax + j] <= B_SO) {
			switch (FLAG[i*jmax + j]) {
			case B_N:	{ P[i*jmax + j] = P[i*jmax + j+1];							break; }
			case B_O:	{ P[i*jmax + j] = P[(i+1)*jmax + j];						break; }
			case B_S:	{ P[i*jmax + j] = P[i*jmax + j-1];							break; } 
			case B_W:	{ P[i*jmax + j] = P[(i-1)*jmax + j];						break; }
			case B_NO:	{ P[i*jmax + j] = 0.5*(P[i*jmax + j+1]+P[(i+1)*jmax + j]);	break; }
			case B_SO:	{ P[i*jmax + j] = 0.5*(P[i*jmax + j-1]+P[(i+1)*jmax + j]);	break; }
			case B_SW:	{ P[i*jmax + j] = 0.5*(P[i*jmax + j-1]+P[(i-1)*jmax + j]);	break; }
			case B_NW:	{ P[i*jmax + j] = 0.5*(P[i*jmax + j+1]+P[(i-1)*jmax + j]);	break; }
			default: break;
			}
		}
	}
			
}

__kernel
void POISSON_2_relaxation_kernel(	__global REAL *P,
									__global REAL *RHS,
									__global int *FLAG,
									int imax,
									int jmax,
									REAL delx,
									REAL dely,
									REAL omg)
{
	imax = imax + 2;
	jmax = jmax + 2;

	int i = get_global_id(0);
	int j = get_global_id(1);

	REAL rdx2, rdy2;
	REAL beta_2;

	rdx2 = 1./delx/delx;
	rdy2 = 1./dely/dely;
	beta_2 = -omg/(2.0*(rdx2+rdy2));
	
	if ((i > 0 && i < imax - 1) && (j > 0 && j < jmax - 1)) {	
		if ((FLAG[i*jmax + j] & C_F) && (FLAG[i*jmax + j] < 0x0100)) {
			P[i*jmax + j] = (1.-omg)*P[i*jmax + j] - 
				beta_2*((P[(i+1)*jmax + j]+P[(i-1)*jmax + j])*rdx2 +
				(P[i*jmax + j+1]+P[i*jmax + j-1])*rdy2 - RHS[i*jmax + j]);
		}
	}
}

__kernel
void POISSON_2_comp_res_kernel(	__global REAL *P,
								__global REAL *RHS,
								__global int *FLAG,
								int imax,
								int jmax,
								REAL delx,
								REAL dely,
								__global REAL *res_result,
								__local REAL *scratch)
{
	imax = imax + 2; // imax == get_global_size(0)
	jmax = jmax + 2; // jmax == get_global_size(1)

	unsigned int gid = get_global_id(0);
	unsigned int tid = get_local_id(0);
	
	int i = gid / imax;
	int j = gid % jmax;

	int ti = tid / imax;
	int tj = tid % jmax;

	//int j = get_global_id(1);

	//int id = i * jmax + j; 

	//__local REAL scratch[16*16];

	//int li = get_local_id(0);
	//int lj = get_local_id(1);

	//int local_size_i = get_local_size(0);
	//int local_size_j = get_local_size(1);

	//int lid = li * get_local_size(1) + lj; 

	REAL add;
	REAL rdx2, rdy2;

	//REAL result = 0.0;

	rdx2 = 1./delx/delx;
	rdy2 = 1./dely/dely;

	if ((i > 0 && i < imax - 1) && (j > 0 && j < jmax - 1)) {
		/* only fluid cells */
		if ((FLAG[i*jmax + j] & C_F) && (FLAG[i*jmax + j] < 0x0100)) {
			add =	(P[(i+1)*jmax + j]-2*P[i*jmax + j]+P[(i-1)*jmax + j])*rdx2
					+ (P[i*jmax + j+1]-2*P[i*jmax + j]+P[i*jmax + j-1])*rdy2
					- RHS[i*jmax + j];
						
			scratch[tid] = add*add;
		} else {
		scratch[tid] = 0.0;
		}
	} else {
		scratch[tid] = 0.0;
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	// do reduction in shared mem
	for (unsigned int s = 1; s < get_local_size(0); s *= 2) {
		// modulo arithmetic is slow!
		if ((tid % (2*s)) == 0) {
			scratch[tid] += scratch[tid + s];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	if (tid == 0) {
		res_result[get_group_id(0)*get_num_groups(1) + get_group_id(1)] = scratch[0];
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

__kernel
void SETBCOND_outer_kernel(__global REAL *U,
					__global REAL *V,
					__global REAL *P,
					__global REAL *TEMP,
					int imax,
					int jmax,
					int wW,
					int wE,
					int wN,
					int wS)
{
	imax = imax + 2;
	jmax = jmax + 2;

	int i = get_global_id(0);
	int j = get_global_id(1);
	
	if ((i >= 0 && i <= imax - 1) && (j >= 0 && j <= jmax - 1)) { // TODO important

		if (wW == 1) {											// slip
			U[0*jmax + j] = 0.0;								// u = 0
			V[0*jmax + j] = V[1*jmax + j];						// dv/dn = 0
		}
		if (wW == 2) {											// no-slip
			U[0*jmax + j] = 0.0;								// u = 0
			V[0*jmax + j] = (-1.0)*V[1*jmax + j];				// v=0 at the boundary by averaging
		}
		if (wW == 3) {											// outflow
			U[0*jmax + j] = U[1*jmax + j];
			V[0*jmax + j] = V[1*jmax + j];
		}
		if (wW == 4) {											// periodic
			U[0*jmax + j] = U[(imax-2)*jmax + j];
			V[0*jmax + j] = V[(imax-2)*jmax + j];				// left and right cells
			V[1*jmax + j] = V[(imax-2)*jmax + j];				// are overlapping
			P[1*jmax + j] = P[(imax-2)*jmax + j]; 
		}

		TEMP[0*jmax + j] = TEMP[1*jmax + j];					// dT/dn = 0

		if (wE == 1) {											// free-slip
			U[(imax-2)*jmax + j] = 0.0;         
			V[(imax-1)*jmax + j] = V[(imax-2)*jmax + j];  
		}
		if (wE == 2) {											// no-slip
			U[(imax-2)*jmax + j] = 0.0;
			V[(imax-1)*jmax + j] = (-1.0)*V[(imax-2)*jmax + j];
		}
		if (wE == 3) {											// outflow
			U[(imax-2)*jmax + j] = U[(imax-3)*jmax + j];
			V[(imax-1)*jmax + j] = V[(imax-2)*jmax + j];
		}
		if (wE == 4) {											// periodic
			U[(imax-2)*jmax + j] = U[1*jmax + j];
			V[(imax-1)*jmax + j] = V[2*jmax + j];
		}

		TEMP[(imax-1)*jmax + j] = TEMP[(imax-2)*jmax + j];

		if (wN == 1) {
			V[i*jmax + (jmax-2)] = 0.0;
			U[i*jmax + (jmax-1)] = U[i*jmax + (jmax-2)];
		}
		if (wN == 2) { 
			V[i*jmax + (jmax-2)] = 0.0;
			U[i*jmax + (jmax-1)] = (-1.0)*U[i*jmax + (jmax-2)];
		}
		if (wN == 3) { 
			V[i*jmax + (jmax-2)] = V[i*jmax + (jmax-3)];
			U[i*jmax + (jmax-1)] = U[i*jmax + (jmax-2)];
		}
		if (wN == 4) { 
			V[i*jmax + (jmax-2)] = V[i*jmax + 1];
			U[i*jmax + (jmax-1)] = U[i*jmax + 2];
		}

		TEMP[i*jmax + 0] = TEMP[i*jmax + 1];

		if (wS == 1) { 
			V[i*jmax + 0] = 0.0;
			U[i*jmax + 0] = U[i*jmax + 1];
		}
		if (wS == 2) { 
			V[i*jmax + 0] = 0.0;
			U[i*jmax + 0] = (-1.0)*U[i*jmax + 1];
		}
		if (wS == 3) { 
			V[i*jmax + 0] = V[i*jmax + 1];
			U[i*jmax + 0] = U[i*jmax + 1];
		}
		if (wS == 4 ) { 
			V[i*jmax + 0] = V[i*jmax + (jmax-3)];
			U[i*jmax + 0] = U[i*jmax + (jmax-3)];
			U[i*jmax + 1] = U[i*jmax + (jmax-2)];
			P[i*jmax + 1] = P[i*jmax + (jmax-2)];
		}

		TEMP[i*jmax + (jmax-1)] = TEMP[i*jmax + (jmax-2)]; 
	}
}

__kernel
void SETBCOND_inner_kernel(__global REAL *U,
					__global REAL *V,
					__global REAL *TEMP,
					__global int *FLAG,
					int imax,
					int jmax)
{
	imax = imax + 2;
	jmax = jmax + 2;

	int i = get_global_id(0);
	int j = get_global_id(1);

	if ((i > 0 && i < imax - 1) && (j > 0 && j < jmax - 1)) {
		if (FLAG[i*jmax + j] & 0x000f) {
			// The mask 0x000f filters the obstacle
			// cells adjacent to fluid cells
			switch (FLAG[i*jmax + j]) {
			case B_N: { 
				V[i*jmax + j]		= 0.0;
				U[i*jmax + j]		= -U[i*jmax + (j+1)];
				U[(i-1)*jmax + j]	= -U[(i-1)*jmax + (j+1)];
				TEMP[i*jmax + j]	= TEMP[i*jmax + (j+1)];
				break;
			}
			case B_O: { 
				U[i*jmax + j]		= 0.0;
				V[i*jmax + j]		= -V[(i+1)*jmax + j];
				V[i*jmax + (j-1)]	= -V[(i+1)*jmax + (j-1)];
				TEMP[i*jmax + j]	= TEMP[(i+1)*jmax + j];
				break;
			}
			case B_S: { 
				V[i*jmax + (j-1)]	= 0.0;
				U[i*jmax + j]		= -U[i*jmax + (j-1)];
				U[(i-1)*jmax + j]	= -U[(i-1)*jmax + (j-1)];
				TEMP[i*jmax + j]	= TEMP[i*jmax + (j-1)];
				break;
			}
			case B_W: { 
				U[(i-1)*jmax + j]	= 0.0;
				V[i*jmax + j]		= -V[(i-1)*jmax + j];
				V[i*jmax + (j-1)]	= -V[(i-1)*jmax + (j-1)];
				TEMP[i*jmax + j]	= TEMP[(i-1)*jmax + j];
				break;
			}
			case B_NO: { 
				V[i*jmax + j]		= 0.0;
				U[i*jmax + j]		= 0.0;
				V[i*jmax + (j-1)]	= -V[(i+1)*jmax + (j-1)];
				U[(i-1)*jmax + j]	= -U[(i-1)*jmax + (j+1)];
				TEMP[i*jmax + j]	= 0.5*(TEMP[i*jmax + (j+1)]+TEMP[(i+1)*jmax + j]);
				break;
			}
			case B_SO: { 
				V[i*jmax + (j-1)]	= 0.0;
				U[i*jmax + j]		= 0.0;
				V[i*jmax + j]		= -V[(i+1)*jmax + j];
				U[(i-1)*jmax + j]	= -U[(i-1)*jmax + (j-1)];
				TEMP[i*jmax + j]	= 0.5*(TEMP[i*jmax + (j-1)]+TEMP[(i+1)*jmax + j]);
				break;
			}
			case B_SW: { 
				V[i*jmax + (j-1)]	= 0.0;
				U[(i-1)*jmax + j]	= 0.0;
				V[i*jmax + j]		= -V[(i-1)*jmax + j];
				U[i*jmax + j]		= -U[i*jmax + (j-1)];
				TEMP[i*jmax + j]	= 0.5*(TEMP[i*jmax + (j-1)]+TEMP[(i-1)*jmax + j]);
				break;
				}
			case B_NW: { 
				V[i*jmax + j]		= 0.0;
				U[(i-1)*jmax + j]	= 0.0;
				V[i*jmax + (j-1)]	= -V[(i-1)*jmax + (j-1)];
				U[i*jmax + j]		= -U[i*jmax + (j+1)];
				TEMP[i*jmax + j]	= 0.5*(TEMP[i*jmax + (j+1)]+TEMP[(i-1)*jmax + j]);
				break;
			}
			default : break;
			}
		}
	}
}

__kernel
void SETSPECBCOND_kernel(	/*char* problem,*/ //TODO should run with every problem
							__global REAL *U,
							__global REAL *V,
							__global REAL *P,
							__global REAL *TEMP,
							int imax,
							int jmax,
							REAL UI,
							REAL VI)
{
	imax = imax + 2;
	jmax = jmax + 2;

	int i = get_global_id(0);
	int j = get_global_id(1);

	//if ((strcmp(problem, "drop") == 0) || (strcmp(problem, "dam") == 0)) {

	//}

	/*-----------------------------------------------------------*/
	/* Driven Cavity: U = 1.0 at the upper boundary              */
	/*-----------------------------------------------------------*/
	//else if (strcmp(problem, "dcavity") == 0) {
		if ((i >= 0 && i <= imax - 2)) {
			U[i*jmax + (jmax-1)] = 2.0 - U[i*jmax + (jmax-2)];
		}
	//}

	//TODO make it work for other 
	/*-----------------------------------------------------------------*/
	/* Flow past a backward facing step, with or without free boundary */
	/*                  U = 1.0 at the left boundary                   */
	/*-----------------------------------------------------------------*/
	//else if (strcmp(problem, "backstep") == 0 || strcmp(problem, "wave") == 0) {
	//	for (j = (jmax/2)+1; j <= jmax-2; j++) {
	//		U[0*jmax + j]    =   1.0;
	//	}
	//}

	/*--------------------------------------------------------------*/
	/* Flow past an obstacle: U = 1.0 at left boundary              */
	/*--------------------------------------------------------------*/
	//else if (strcmp(problem, "plate") == 0 || strcmp(problem, "circle") == 0) {
	//	V[0*jmax + 0] = 2*VI-V[1*jmax + 0];
	//	for (j = 1; j <= jmax-2; j++)
	//	{
	//		U[0*jmax + j] = UI;
	//		V[0*jmax + j] = 2*VI-V[1*jmax + j];
	//	}
	//}

	//TODO may need fix
	/*---------------------------------------------------------------------*/
	/* Inflow for injection molding: U = 1.0 in the mid of left boundary   */
	/*---------------------------------------------------------------------*/
	//else if (strcmp(problem, "molding") == 0) {
	//	for (j = (int)(0.4*jmax)+1; j <= (int)(0.6*jmax); j++) {
	//		U[0*jmax + j] = 1.0;       
	//	}
	//}

	/*------------------------------------------------------------------*/
	/* natural convection or fluidtrap: left T = 0.5 right T = -0.5     */
	/*                          upper and lower wall adiabatic          */
	/*------------------------------------------------------------------*/
	//else if (strcmp(problem, "convection") == 0 || strcmp(problem, "fluidtrap") == 0) {
	//	for (j = 0; j <= jmax-1; j++)
	//	{
	//		TEMP[0*jmax + j] = 2*(0.5)-TEMP[1*jmax + j];				// left wall heated
	//		TEMP[(imax-1)*jmax + j] = 2*(-0.5)-TEMP[(imax-2)*jmax + j]; // right wall cooled
	//	}
	//	for (i=0; i <= imax+1; i++)
	//	{
	//		TEMP[i*jmax + 0] = TEMP[i*jmax + 1];
	//		TEMP[i*jmax + (jmax-1)] = TEMP[i*jmax + (jmax-2)];			// adiabatic walls
	//	}
	//}

	/*----------------------------------------------------*/
	/* Rayleigh-Benard flow: top T = -0.5 bottom T = 0.5  */
	/*                       left and right adiabatic     */
	/*----------------------------------------------------*/
	//else if (strcmp(problem, "rayleigh") == 0) {
	//	for (j = 0; j <= jmax-1; j++)
	//	{
	//		TEMP[0*jmax + j] = TEMP[1*jmax + j];         
	//		TEMP[(imax-1)*jmax + j] = TEMP[(imax-2)*jmax + j];			// adiabatic walls
	//	}
	//	for (i=0; i <= imax+1; i++)
	//	{
	//		TEMP[i*jmax + 0] = 2*(0.5)-TEMP[i*jmax + 1];				// lower wall heated
	//		TEMP[i*jmax + (jmax-1)] = 2*(-0.5)-TEMP[i*jmax + (jmax-2)];	// upper wall cooled
	//	}
	//	return;
	//}

	/*------*/
	/* ELSE */
	/*------*/
	//printf("Problem %s not defined!\n", problem);

}