#pragma OPENCL EXTENSION cl_khr_fp64 : enable 
#pragma extension cl_nv_compiler_options : enable
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
						REAL Pr,
						__local REAL *TEMP_l)
{
	REAL LAPLT, DUTDX, DVTDY;
	// TODO change delx dely
	imax = imax + 2;
	jmax = jmax + 2;

	int i = get_global_id(1);
	int j = get_global_id(0);

	int ti = get_global_id(1);
	int tj = get_global_id(0);

	int local_size_i = get_local_size(1);
	int local_size_j = get_local_size(0);

	//TEMP_l[ti*local_size_j+tj] = TEMP[i*jmax + j];
	

	//barrier(CLK_LOCAL_MEM_FENCE);

	//TEMP_new[i*jmax + j] = TEMP_l[ti*local_size_j + tj];
	TEMP_new[i*jmax + j] = TEMP[i*jmax + j];

	//if ((i > 0 && i < imax - 1) && (j > 0 && j < jmax - 1)) {
	//	if ( (FLAG[i*jmax + j] & C_F) && (FLAG[i*jmax + j] < C_E) ) {
	//		if (ti < local_size_i-1 && tj < local_size_j-1) {
	//		LAPLT = (TEMP_l[(ti+1)*local_size_j + tj]-2.0*TEMP_l[ti*local_size_j + tj]+TEMP_l[(ti-1)*local_size_j + tj])*(1./delx/delx) +
	//			(TEMP_l[ti*local_size_j + tj+1]-2.0*TEMP_l[ti*local_size_j + tj]+TEMP_l[ti*local_size_j + tj-1])*(1./dely/dely);
	//		DUTDX = ( (U[i*jmax + j]*0.5*(TEMP_l[ti*local_size_j + tj]+TEMP_l[(ti+1)*local_size_j + tj]) -
	//			U[(i-1)*jmax + j]*0.5*(TEMP_l[(ti-1)*local_size_j + tj]+TEMP_l[ti*local_size_j + tj])) +
	//			gamma*(fabs(U[i*jmax + j])*0.5*(TEMP_l[ti*local_size_j + tj]-TEMP_l[(ti+1)*local_size_j + tj]) -
	//			fabs(U[(i-1)*jmax + j])*0.5*(TEMP_l[(ti-1)*local_size_j + tj]-TEMP_l[ti*local_size_j + tj]))
	//			)/delx;
	//		DVTDY = ( (V[i*jmax + j]*0.5*(TEMP_l[ti*local_size_j + tj]+TEMP_l[ti*local_size_j + tj+1]) -
	//			V[i*jmax + j-1]*0.5*(TEMP_l[ti*local_size_j + tj-1]+TEMP_l[ti*local_size_j + tj])) +
	//			gamma*(fabs(V[i*jmax + j])*0.5*(TEMP_l[ti*local_size_j + tj]-TEMP_l[ti*local_size_j + tj+1]) -
	//			fabs(V[i*jmax + j-1])*0.5*(TEMP_l[ti*local_size_j + tj-1]-TEMP_l[ti*local_size_j + tj]))
	//			)/dely;
	//			
	//		TEMP_new[i*jmax + j] = TEMP_l[ti*local_size_j + tj]+delt*(LAPLT/Re/Pr - DUTDX - DVTDY);
	//		}
	//	}
	//}
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

	int i = get_global_id(1);
	int j = get_global_id(0);


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
	
	int i = get_global_id(1);
	int j = get_global_id(0);

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
						__global REAL *p0_result,
						__local REAL *scratch)
{
	imax = imax + 2;
	jmax = jmax + 2;

	unsigned int gid = get_global_id(0);
	unsigned int tid = get_local_id(0);
	
	int i = gid / imax;
	int j = gid % jmax;
	
	if ((i > 0 && i < imax - 1) && (j > 0 && j < jmax - 1)) {
		if (FLAG[i*jmax + j] & C_F) {
			scratch[tid] = P[i*jmax + j]*P[i*jmax + j];
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
		p0_result[get_group_id(0)*get_num_groups(1) + get_group_id(1)] = scratch[0];
	}
}

__kernel
void POISSON_1_relaxation_kernel(	__global REAL *P,
									__global REAL *RHS,
									__global int *FLAG,
									int imax,
									int jmax,
									REAL delx,
									REAL dely,
									REAL omg)
{
	  //TODO implement
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
}

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
	 //TODO implement
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
}

__kernel
void POISSON_2_copy_boundary_kernel(__global REAL *P,
									__global int *FLAG,
									int imax,
									int jmax)
{
	imax = imax + 2;
	jmax = jmax + 2;

	int i = get_global_id(1);
	int j = get_global_id(0);

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

//TODO It is impossible to implement this kernel without changing the algorithm.
// No synchronization between workgroups possible
__kernel
void POISSON_2_relaxation_kernel(	__global REAL *P,
									__global REAL *RHS,
									__global int *FLAG,
									int imax,
									int jmax,
									REAL rdx2,
									REAL rdy2,
									REAL beta_2,
									REAL omg)
{
	imax = imax + 2;
	jmax = jmax + 2;

	unsigned int gid = get_global_id(0);
	
	int i, j;

	//if (gid == 0) {
		for (i = 1; i <= imax-2; i++) {
			for (j = 1; j <= jmax-2; j++) {
				if ((FLAG[i*jmax + j] & C_F) && (FLAG[i*jmax + j] < 0x0100)) {
					P[i*jmax + j] = (1.-omg)*P[i*jmax + j]
						- beta_2 * 
							((P[(i+1)*jmax + j] + P[(i-1)*jmax + j]) * rdx2 +
							 (P[i*jmax + j+1] + P[i*jmax + j-1]) * rdy2
						- RHS[i*jmax + j]);
				}
			}
		}
	//}
}

//__kernel
//void POISSON_2_relaxation_kernel(	__global REAL *P,
//									__global REAL *RHS,
//									__global int *FLAG,
//									int imax,
//									int jmax,
//									REAL delx,
//									REAL dely,
//									REAL omg,
//									__local REAL *shared)
//{
//
//	int xmax = imax + 2;
//	int ymax = jmax + 2;
//
//	int x = get_global_id(0);
//	int y = get_global_id(1);
//	int gid = y * xmax + x;
//
//	int tx = get_local_id(0);
//	int ty = get_local_id(1);
//
//	REAL rdx2, rdy2;
//	REAL beta_2;
//
//	rdx2 = 1./delx/delx;
//	rdy2 = 1./dely/dely;
//	beta_2 = -omg/(2.0*(rdx2+rdy2));
//	
//	__local REAL a[16][16][3];
//
//	int num_groups_y = get_num_groups(1);
//	int local_size_y = get_local_size(1);
//	int group_id_y = get_group_id(1);
//
//	int current_y, current_group = 2;
//	if ((x > 0 && x < xmax - 1) && (y > 0 && y < ymax - 1)) {
//		if ((FLAG[gid] & C_F) && (FLAG[gid] < 0x0100)) {
//			//for (current_group = 0; current_group < num_groups_y; current_group++) {
//				if (group_id_y == current_group) {
//					for (current_y = 0; current_y < local_size_y; current_y++) {
//						if (current_y == ty) {
//							P[gid] = P[(y-1)*xmax + x] + 1 + P[(y+1)*xmax + x];
//						}
//						barrier(CLK_GLOBAL_MEM_FENCE);
//						
//					}
//					barrier(CLK_GLOBAL_MEM_FENCE);
//					if (current_y == local_size_y)
//						current_group++;
//				} 
//				//mem_fence(CLK_GLOBAL_MEM_FENCE);
//			//}
//		}
//	}
//
//	//		int current_y, current_group;
//	//if ((x > 0 && x < xmax - 1) && (y > 0 && y < ymax - 1)) {
//	//	if ((FLAG[gid] & C_F) && (FLAG[gid] < 0x0100)) {
//	//		for (current_group = 0; current_group < get_num_groups(1); current_group++) {
//	//			if (current_group == get_group_id(1)) {
//	//				for (current_y = 0; current_y < get_local_size(1); current_y++) {
//	//				//if ((tx > 0 && tx < (get_local_size(0) - 1)) && (ty > 0 && ty < (get_local_size(1) - 1))) {
//	//					if (current_y == ty && ) {
//
//	//						//a[ty][tx][0] = P[(y-1)*xmax + x];
//	//						//a[ty][tx][1] = 1;
//	//						//a[ty][tx][2] = P[(y+1)*xmax + x];
//	//		
//	//						//barrier(CLK_LOCAL_MEM_FENCE);
//
//	//						//P[gid] = a[ty][tx][1] + a[ty][tx][2]+ a[ty][tx][0];
//	//						P[gid] = P[(y-1)*xmax + x] + 1 + P[(y+1)*xmax + x];
//
//	//						//P[gid] = (1.-omg)*a[ty][tx][1] - 
//	//						//	beta_2*((a[ty][tx][2]+a[ty][tx][0])*rdx2 +
//	//						//	(a[ty][tx+1][1]+a[ty][tx-1][1])*rdy2 - RHS[gid]);
//	//						//P[gid] = 1;
//	//						//barrier(CLK_GLOBAL_MEM_FENCE);
//	//					}
//	//					barrier(CLK_GLOBAL_MEM_FENCE);
//	//			}
//	//			barrier(CLK_GLOBAL_MEM_FENCE);
//	//				//}
//	//				//barrier(CLK_LOCAL_MEM_FENCE);
//	//		}
//	//	}
//	//}
//
//
//	//if (y==0) {
//	//if (/*(x > 0 && x < xmax - 1) && */(/*y > 0 &&*/ y < ymax - 1)) {	
//		//if ((FLAG[gid] & C_F) && (FLAG[gid] < 0x0100)) {
//			//for (current_y = 0; current_y < get_global_size(1); current_y++) {
//
//				//if (current_y == ty) {
//					//a[ty][tx][0] = P[(y-1)*xmax + x];
//					//a[ty][tx][1] = P[y*xmax + x];
//					//a[ty][tx][2] = P[(y+1)*xmax + x];
//
//					//barrier(CLK_LOCAL_MEM_FENCE);
//
//					//P[gid] = gid /*a[current_y][tx][1] + a[current_y][tx][2] +a[current_y][tx][0]*/;
//
//					//P[gid] = (1.-omg)*a[current_y][tx][1] - 
//					//	beta_2*((a[current_y][tx][2]+a[current_y][tx][0])*rdx2 +
//					//	(a[current_y][tx+1][1]+a[current_y][tx-1][1])*rdy2 - RHS[gid]);
//
//					//barrier(CLK_LOCAL_MEM_FENCE);
//				//}
//				//barrier(CLK_LOCAL_MEM_FENCE);
//			//}
//		//}
//	//} else {
//	//	a[ty][tx][0] = 0.0;
//	//	a[ty][tx][1] = P[y*xmax + x];
//	//	a[ty][tx][2] = P[(y+1)*xmax + x];
//	//}
//
//		//P[i*jmax + j] = (1.-omg)*P[i*jmax + j] - 
//	//beta_2*((P[(i+1)*jmax + j]+P[(i-1)*jmax + j])*rdx2 +
//	//(P[i*jmax + j+1]+P[i*jmax + j-1])*rdy2 - RHS[i*jmax + j]);
//
//	//if ((i > 0 && i < imax - 1) && (j > 0 && j < jmax - 1)) {	
//	//	if ((FLAG[i*jmax + j] & C_F) && (FLAG[i*jmax + j] < 0x0100)) {
//	//		for (group_id = 0; group_id < get_global_size(0); group_id += 16) {
//	//			liid = group_id*get_local_size(0)+ti;
//	//			ljid = group_id*get_local_size(1)+tj;
//
//	//			a[ti][tj][0] = P[(liid-1)*jmax + ljid];
//	//			a[ti][tj][1] = P[liid*jmax + ljid];
//	//			a[ti][tj][2] = P[(liid+1)*jmax + ljid];
//
//	//			barrier(CLK_LOCAL_MEM_FENCE);
//
//	//			//P[liid*jmax + ljid] = (1.-omg)*a[ti][tj][1] - 
//	//			//	beta_2*((a[ti][tj][2]+a[ti][tj][0])*rdx2 +
//	//			//	(a[ti][tj+1][0]+a[liid][tj-1][0])*rdy2 - RHS[liid*jmax + ljid]);
//	//	
//	//			barrier(CLK_LOCAL_MEM_FENCE);
//	//		}
//	//	}
//	//}
//}

//__kernel
//void POISSON_2_comp_res_kernel(	__global REAL *P,
//								__global REAL *RHS,
//								__global int *FLAG,
//								int imax,
//								int jmax,
//								REAL rdx2,
//								REAL rdy2,
//								__global REAL *res)
//{
//	// shared reduction kernel 0 - naive
//	imax = imax + 2;
//	jmax = jmax + 2;
//
//	unsigned int gid = get_global_id(0);
//
//	int i, j;
//	REAL add;
//
//	if (gid == 0) {
//		*res = 0.0;
//
//		for (i = 1; i <= imax-2; i++) {
//			for (j = 1; j <= jmax-2; j++) {
//				if ((FLAG[i*jmax + j] & C_F) && (FLAG[i*jmax + j] < 0x0100)) {
//					add =	(P[(i+1)*jmax + j]-2*P[i*jmax + j]+P[(i-1)*jmax + j])*rdx2
//						+ (P[i*jmax + j+1]-2*P[i*jmax + j]+P[i*jmax + j-1])*rdy2
//						- RHS[i*jmax + j];
//
//					*res += add * add;
//
//				}
//			}
//		}
//	}
//
//}


__kernel
void POISSON_2_comp_res_kernel(	__global REAL *P,
								__global REAL *RHS,
								__global int *FLAG,
								int imax,
								int jmax,
								REAL rdx2,
								REAL rdy2,
								__global REAL *res_result,
								__local REAL *scratch)
{
	imax = imax + 2;
	jmax = jmax + 2;

	unsigned int gid = get_global_id(0);
	unsigned int tid = get_local_id(0);
	
	// slow on GPU
	//int i = gid % imax;
	//int j = gid / jmax;

	int i = gid & (imax - 1);
	int j = gid >> (int)log2((float)jmax);

	REAL add;

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

	// shared reduction kernel 1
	//for (unsigned int s = 1; s < get_local_size(0); s *= 2) {
	//	if ((tid % (2*s)) == 0) {
	//		scratch[tid] += scratch[tid + s];
	//	}
	//	barrier(CLK_LOCAL_MEM_FENCE);
	//}

	// shared reduction kernel 2
	//for (int s = 1; s < get_local_size(0); s *= 2) {
	//	int index = 2 * s * tid;
	//	if (index < get_local_size(0)) {
	//		scratch[index] += scratch[index + s];
	//	}
	//	barrier(CLK_LOCAL_MEM_FENCE);
	//}

	// shared reduction kernel 3
	for (int s = get_local_size(0)/2; s>0; s >>= 1) {
		if (tid < s) {
			scratch[tid] += scratch[tid + s];
		}
		barrier(CLK_LOCAL_MEM_FENCE);
	}

	if (tid == 0) {
		res_result[get_group_id(0)] = scratch[0];
	}

}

//TODO bad results
//__kernel
//void POISSON_2_comp_res_kernel(	__global REAL *P,
//								__global REAL *RHS,
//								__global int *FLAG,
//								int imax,
//								int jmax,
//								REAL rdx2,
//								REAL rdy2,
//								__global REAL *res_result,
//								__local REAL *scratch)
//{
//	// shared reduction kernel 4
//	imax = imax + 2;
//	jmax = jmax + 2;
//
//	unsigned int local_size = get_local_size(0);
//	unsigned int tid = get_local_id(0);
//	unsigned int gid = get_group_id(0)*(local_size*2) + tid;
//	
//	int i = gid & (imax - 1);
//	int j = gid >> (int)log2((float)jmax);
//
//	REAL add1, add2;
//
//	if ((i > 0 && i < imax - 1) && (j > 0 && j < jmax - 1)) {
//		/* only fluid cells */
//		if ((FLAG[i*jmax + j] & C_F) && (FLAG[i*jmax + j] < 0x0100)) {
//			add1 =	(P[(i+1)*jmax + j]-2*P[i*jmax + j]+P[(i-1)*jmax + j])*rdx2
//					+ (P[i*jmax + j+1]-2*P[i*jmax + j]+P[i*jmax + j-1])*rdy2
//					- RHS[i*jmax + j];
//			add2 = (P[(i+1)*jmax + j + (local_size)]
//						- 2*P[i*jmax + j + local_size]
//						+ P[(i-1)*jmax + j + (local_size)])*rdx2
//					+ (P[i*jmax + j+1 + local_size]
//						- 2*P[i*jmax + j + local_size]
//						+ P[i*jmax + j-1 + local_size])*rdy2
//					- RHS[i*jmax + j + local_size];
//						
//			scratch[tid] = add1*add1 + add2*add2;
//		} else {
//		scratch[tid] = 0.0;
//		}
//	} else {
//		scratch[tid] = 0.0;
//	}
//
//	barrier(CLK_LOCAL_MEM_FENCE);
//	
//	for (int s = local_size/2; s>0; s >>= 1) {
//		if (tid < s) {
//			scratch[tid] += scratch[tid + s];
//		}
//		barrier(CLK_LOCAL_MEM_FENCE);
//	}
//
//	if (tid == 0) {
//		res_result[get_group_id(0)] = scratch[0];
//	}
//}

//TODO bad results
//__kernel
//void POISSON_2_comp_res_kernel(	__global REAL *P,
//								__global REAL *RHS,
//								__global int *FLAG,
//								int imax,
//								int jmax,
//								REAL rdx2,
//								REAL rdy2,
//								__global REAL *res_result,
//								__local REAL *scratch)
//{
//	// shared reduction kernel 5
//	imax = imax + 2; // imax == get_global_size(0)
//	jmax = jmax + 2; // jmax == get_global_size(1)
//
//	unsigned int local_size = get_local_size(0);
//	unsigned int tid = get_local_id(0);
//	unsigned int gid = get_group_id(0)*(local_size*2) + tid;
//	
//	int i = gid & (imax - 1);
//	int j = gid >> (int)log2((float)jmax);
//
//	REAL add1, add2;
//
//	if ((i > 0 && i < imax - 1) && (j > 0 && j < jmax - 1)) {
//		/* only fluid cells */
//		if ((FLAG[i*jmax + j] & C_F) && (FLAG[i*jmax + j] < 0x0100)) {
//			add1 =	(P[(i+1)*jmax + j]-2*P[i*jmax + j]+P[(i-1)*jmax + j])*rdx2
//					+ (P[i*jmax + j+1]-2*P[i*jmax + j]+P[i*jmax + j-1])*rdy2
//					- RHS[i*jmax + j];
//			add2 = (P[(i+1)*jmax + j + local_size]-2*P[i*jmax + j + local_size]+P[(i-1)*jmax + j + local_size])*rdx2
//					+ (P[i*jmax + j+1 + local_size]-2*P[i*jmax + j + local_size]+P[i*jmax + j-1 + local_size])*rdy2
//					- RHS[i*jmax + j + local_size];
//						
//			scratch[tid] = add1*add1 + add2*add2;
//		} else {
//		scratch[tid] = 0.0;
//		}
//	} else {
//		scratch[tid] = 0.0;
//	}
//
//	barrier(CLK_LOCAL_MEM_FENCE);
//	
//	for (int s = local_size/2; s>32; s >>= 1) {
//		if (tid < s) {
//			scratch[tid] += scratch[tid + s];
//		}
//		barrier(CLK_LOCAL_MEM_FENCE);
//	}
//
//	if (tid < 32) {
//		scratch[tid] += scratch[tid + 32];
//		scratch[tid] += scratch[tid + 16];
//		scratch[tid] += scratch[tid + 8];
//		scratch[tid] += scratch[tid + 4];
//		scratch[tid] += scratch[tid + 2];
//		scratch[tid] += scratch[tid + 1];
//	}
//
//
//	if (tid == 0) {
//		res_result[get_group_id(0)] = scratch[0];
//	}
//}

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

	int i = get_global_id(1);
	int j = get_global_id(0);

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

	int i = get_global_id(1);
	int j = get_global_id(0);
	
	if ((i >= 0 && i <= imax - 1) && (j >= 0 && j <= jmax - 1)) { // TODO important
		// western and eastern boundary
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

		// northern and southern boundary
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

	int i = get_global_id(1);
	int j = get_global_id(0);

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

	int i = get_global_id(1);
	int j = get_global_id(0);

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