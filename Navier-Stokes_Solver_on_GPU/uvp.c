#include <stdio.h>
#include <math.h>
#include "datadef.h"
#include "init.h"

/*---------------------------------------------------------------*/
/* Computation of new temperature                                */
/*---------------------------------------------------------------*/
void COMP_TEMP(REAL **U,REAL **V,REAL **TEMP,int **FLAG,
	       int imax,int jmax,REAL delt,REAL delx,REAL dely,
	       REAL gamma,REAL Re,REAL Pr)
{
 int  i,j;
 REAL LAPLT, DUTDX, DVTDY,indelx2,indely2;
 REAL **T2;

 T2 = RMATRIX(0,imax+1,0,jmax+1);
 indelx2 = 1./delx/delx;
 indely2 = 1./dely/dely;

 for(i=1;i<=imax;i++)
    for(j=1;j<=jmax;j++)
       if( (FLAG[i][j] & C_F) && (FLAG[i][j] < C_E) )
	 {
	  LAPLT = (TEMP[i+1][j]-2.0*TEMP[i][j]+TEMP[i-1][j])*indelx2 +
		  (TEMP[i][j+1]-2.0*TEMP[i][j]+TEMP[i][j-1])*indely2;
	  DUTDX = ( (U[i][j]*0.5*(TEMP[i][j]+TEMP[i+1][j]) -
				U[i-1][j]*0.5*(TEMP[i-1][j]+TEMP[i][j])) +
		     gamma*(fabs(U[i][j])*0.5*(TEMP[i][j]-TEMP[i+1][j]) -
				fabs(U[i-1][j])*0.5*(TEMP[i-1][j]-TEMP[i][j]))
		  )/delx;
	  DVTDY = ( (V[i][j]*0.5*(TEMP[i][j]+TEMP[i][j+1]) -
				V[i][j-1]*0.5*(TEMP[i][j-1]+TEMP[i][j])) +
		     gamma*(fabs(V[i][j])*0.5*(TEMP[i][j]-TEMP[i][j+1]) -
				fabs(V[i][j-1])*0.5*(TEMP[i][j-1]-TEMP[i][j]))
		  )/dely;
	  T2[i][j] = TEMP[i][j]+delt*(LAPLT/Re/Pr - DUTDX - DVTDY);
	 }

 for(i=1;i<=imax;i++)
    for(j=1;j<=jmax;j++)
       if( (FLAG[i][j] & C_F) && (FLAG[i][j] < C_E) )
	  TEMP[i][j] = T2[i][j];

 FREE_RMATRIX(T2,0,imax+1,0,jmax+1);
}

/*----------------------------------------------------------------*/
/* Computation of tentative velocity field (F,G)                  */
/*----------------------------------------------------------------*/
void COMP_FG(REAL **U,REAL **V,REAL **TEMP,REAL **F,REAL **G,int **FLAG,
	     int imax,int jmax,REAL delt,REAL delx,REAL dely,
	     REAL GX,REAL GY,REAL gamma,REAL Re,REAL beta)
{
 int  i,j;
 REAL DU2DX,DUVDY,DUVDX,DV2DY,LAPLU,LAPLV;

 for (i=1;i<=imax-1;i++)
    for (j=1;j<=jmax;j++)
      {
       /* only if both adjacent cells are fluid cells */
       if( ((FLAG[i][j] & C_F) && (FLAG[i][j] < C_E)) &&
           ((FLAG[i+1][j] & C_F) && (FLAG[i+1][j] < C_E)) )
         {
          DU2DX = ((U[i][j]+U[i+1][j])*(U[i][j]+U[i+1][j])+
	                 gamma*fabs(U[i][j]+U[i+1][j])*(U[i][j]-U[i+1][j])-
	           (U[i-1][j]+U[i][j])*(U[i-1][j]+U[i][j])-
	   	      gamma*fabs(U[i-1][j]+U[i][j])*(U[i-1][j]-U[i][j]))
                  /(4.0*delx);
          DUVDY = ((V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1])+
                         gamma*fabs(V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1])-
	           (V[i][j-1]+V[i+1][j-1])*(U[i][j-1]+U[i][j])-
	                 gamma*fabs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j]))
                  /(4.0*dely);
          LAPLU = (U[i+1][j]-2.0*U[i][j]+U[i-1][j])/delx/delx+
	          (U[i][j+1]-2.0*U[i][j]+U[i][j-1])/dely/dely;
   
          F[i][j] = U[i][j]+delt*(LAPLU/Re-DU2DX-DUVDY+GX)
		           -delt*beta*GX*(TEMP[i][j]+TEMP[i+1][j])/2;
         }
       else
          F[i][j] = U[i][j];
      }

 for (i=1;i<=imax;i++)
    for (j=1;j<=jmax-1;j++)
      {
       /* only if both adjacent cells are fluid cells */
       if( ((FLAG[i][j] & C_F) && (FLAG[i][j] < C_E)) &&
           ((FLAG[i][j+1] & C_F) && (FLAG[i][j+1] < C_E)) )
         {
          DUVDX = ((U[i][j]+U[i][j+1])*(V[i][j]+V[i+1][j])+
	   	      gamma*fabs(U[i][j]+U[i][j+1])*(V[i][j]-V[i+1][j])-
	   	(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]+V[i][j])-
	   	      gamma*fabs(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]-V[i][j]))
	          /(4.0*delx);
          DV2DY = ((V[i][j]+V[i][j+1])*(V[i][j]+V[i][j+1])+
	   	      gamma*fabs(V[i][j]+V[i][j+1])*(V[i][j]-V[i][j+1])-
	           (V[i][j-1]+V[i][j])*(V[i][j-1]+V[i][j])-
	   	      gamma*fabs(V[i][j-1]+V[i][j])*(V[i][j-1]-V[i][j]))
	          /(4.0*dely);

          LAPLV = (V[i+1][j]-2.0*V[i][j]+V[i-1][j])/delx/delx+
	          (V[i][j+1]-2.0*V[i][j]+V[i][j-1])/dely/dely;

          G[i][j] = V[i][j]+delt*(LAPLV/Re-DUVDX-DV2DY+GY)
		           -delt*beta*GY*(TEMP[i][j]+TEMP[i][j+1])/2;;		      
         }
       else
          G[i][j] = V[i][j];
      }
 /* F und G at external boundary */
 /*------------------------------*/ 
 for (j=1;j<=jmax;j++)
   {
    F[0][j]    = U[0][j];
    F[imax][j] = U[imax][j];
   }
 for (i=1;i<=imax;i++)
   {
    G[i][0]    = V[i][0];
    G[i][jmax] = V[i][jmax];
   }
}


/*-------------------------------------------------------------*/
/* Computation of the right hand side of the pressure equation */
/*-------------------------------------------------------------*/
void COMP_RHS(REAL **F,REAL **G,REAL **RHS,int **FLAG,int imax,int jmax,
              REAL delt,REAL delx,REAL dely)
{
 int i,j;

 for (i=1;i<=imax;i++)
    for (j=1;j<=jmax;j++)
       if ((FLAG[i][j] & C_F) && (FLAG[i][j] < 0x0100))  
	 /* only for fluid and non-surface cells */
          RHS[i][j] = ((F[i][j]-F[i-1][j])/delx+(G[i][j]-G[i][j-1])/dely)/delt;
}


/*-------------------------------------------------------------*/
/* SOR iteration for the poisson equation for the pressure     */
/*-------------------------------------------------------------*/
int POISSON(REAL **P,REAL **RHS,int **FLAG,
            int imax,int jmax,REAL delx,REAL dely,
            REAL eps,int itermax,REAL omg,REAL *res,int ifull,int p_bound)
{
 int i,j,iter;
 REAL rdx2,rdy2;
 REAL add,beta_2,beta_mod;
 REAL p0 = 0.0;

 rdx2 = 1./delx/delx;
 rdy2 = 1./dely/dely;
 beta_2 = -omg/(2.0*(rdx2+rdy2));

 for (i=1;i<=imax;i++)
   for (j=1;j<=jmax;j++)
     if (FLAG[i][j] & C_F)
       p0 += P[i][j]*P[i][j];

 p0 = sqrt(p0/ifull);
 if (p0 < 0.0001)
   p0 = 1.0;

                                            /* SOR-iteration */
                                            /*---------------*/
 for (iter=1;iter<=itermax;iter++)
   {
    if (p_bound == 1)
                        /* modify the equation at the boundary */
                        /*-------------------------------------*/
      {
                        /* relaxation for fluid cells */
                        /*----------------------------*/
       for (i=1;i<=imax;i+=1)
          for (j=1;j<=jmax;j+=1)   
                        /* five point star for interior fluid cells */
             if (FLAG[i][j] == 0x001f) 
                P[i][j] = (1.-omg)*P[i][j] - 
                          beta_2*((P[i+1][j]+P[i-1][j])*rdx2 +
                                  (P[i][j+1]+P[i][j-1])*rdy2 - RHS[i][j]);
                        /* modified star near boundary */
             else if ((FLAG[i][j] & C_F) && (FLAG[i][j] < 0x0100)){ 
                beta_mod = -omg/((eps_E+eps_W)*rdx2+(eps_N+eps_S)*rdy2);
                P[i][j] = (1.-omg)*P[i][j] -
                          beta_mod*( (eps_E*P[i+1][j]+eps_W*P[i-1][j])*rdx2 +
                                     (eps_N*P[i][j+1]+eps_S*P[i][j-1])*rdy2 -
                                     RHS[i][j]);
	      }
                         /* computation of residual */
                         /*-------------------------*/
       *res = 0.0;
       for (i=1;i<=imax;i++)
          for (j=1;j<=jmax;j++)
             if ((FLAG[i][j] & C_F) && (FLAG[i][j] < 0x0100))   
                         /* only fluid cells */
                         /*------------------*/
                {
                 add =  (eps_E*(P[i+1][j]-P[i][j]) - 
                         eps_W*(P[i][j]-P[i-1][j])) * rdx2  +
                        (eps_N*(P[i][j+1]-P[i][j]) -
                         eps_S*(P[i][j]-P[i][j-1])) * rdy2  -  RHS[i][j];
                 *res += add*add;
 	        }

        *res = sqrt((*res)/ifull)/p0;
                         /* convergence? */
                         /*--------------*/

        if (*res<eps)
           return iter;
       }

    else if (p_bound == 2)
      {
                         /* copy values at external boundary */
                         /*----------------------------------*/
	for (i=1;i<=imax;i+=1)
	  {
	    P[i][0]      = P[i][1];
	    P[i][jmax+1] = P[i][jmax];
	  }
	for (j=1;j<=jmax;j+=1)
	  {
	    P[0][j]      = P[1][j];
	    P[imax+1][j] = P[imax][j];
	  }
	/* and at interior boundary cells */
	/*--------------------------------*/
	for (i=1;i<=imax;i+=1)
          for (j=1;j<=jmax;j+=1)
	    if (FLAG[i][j] >=B_N && FLAG[i][j] <=B_SO) 
	      switch (FLAG[i][j])
		{
		case B_N:{  P[i][j] = P[i][j+1];                 break;}
		case B_O:{  P[i][j] = P[i+1][j];                 break;}
		case B_S:{  P[i][j] = P[i][j-1];                 break;} 
		case B_W:{  P[i][j] = P[i-1][j];                 break;}
		case B_NO:{ P[i][j] = 0.5*(P[i][j+1]+P[i+1][j]); break;}
		case B_SO:{ P[i][j] = 0.5*(P[i][j-1]+P[i+1][j]); break;}
		case B_SW:{ P[i][j] = 0.5*(P[i][j-1]+P[i-1][j]); break;}
		case B_NW:{ P[i][j] = 0.5*(P[i][j+1]+P[i-1][j]); break;}
		default:                                         break;
		}
	
	/* relaxation for fluid cells */
	/*----------------------------*/
	for (i=1;i<=imax;i+=1)
          for (j=1;j<=jmax;j+=1)
	    if ((FLAG[i][j] & C_F) && (FLAG[i][j] < 0x0100))
	      P[i][j] = (1.-omg)*P[i][j] - 
		beta_2*((P[i+1][j]+P[i-1][j])*rdx2 +
			(P[i][j+1]+P[i][j-1])*rdy2 - RHS[i][j]);
	
	/* computation of residual */
	/*-------------------------*/
	*res = 0.0;
	for (i=1;i<=imax;i++)
          for (j=1;j<=jmax;j++)
	    if ((FLAG[i][j] & C_F) && (FLAG[i][j] < 0x0100))   
	      /* only fluid cells */
	      /*------------------*/
	      {
		add =  (P[i+1][j]-2*P[i][j]+P[i-1][j])*rdx2+
		  (P[i][j+1]-2*P[i][j]+P[i][j-1])*rdy2-RHS[i][j];
		*res += add*add;
	      }
	
        *res = sqrt((*res)/ifull)/p0;
	/* convergence? */
	/*--------------*/
	
        if (*res<eps)
	  return iter;
      }
   }
 return iter;
}


/*---------------------------------------------------------------*/
/* Computation of new velocity values                            */
/*---------------------------------------------------------------*/
void ADAP_UV (REAL **U,REAL **V,REAL **F,REAL **G,REAL **P,int **FLAG,
              int imax,int jmax,REAL delt,REAL delx,REAL dely)
{
 int i,j;

 for (i=1;i<=imax-1;i++)
    for (j=1;j<=jmax;j++)
       /* only if both adjacent cells are fluid cells */
       if( ((FLAG[i][j] & C_F) && (FLAG[i][j] < C_E)) &&
           ((FLAG[i+1][j] & C_F) && (FLAG[i+1][j] < C_E)) )
          U[i][j] = F[i][j]-(P[i+1][j]-P[i][j])*delt/delx;

 for (i=1;i<=imax;i++)
    for (j=1;j<=jmax-1;j++)
       /* only if both adjacent cells are fluid cells */
       if( ((FLAG[i][j] & C_F) && (FLAG[i][j] < C_E)) &&
           ((FLAG[i][j+1] & C_F) && (FLAG[i][j+1] < C_E)) )
          V[i][j] = G[i][j]-(P[i][j+1]-P[i][j])*delt/dely;
}


/*------------------------------------------------------------*/
/* Computation of adaptive time stepsize satisfying           */
/* the CFL stability criteria                                 */
/* and set the flag "write" if some data has to be written    */
/* into a file.                                               */
/*------------------------------------------------------------*/
void COMP_delt(REAL *delt, REAL t, int imax, int jmax, REAL delx, REAL dely,
               REAL **U, REAL **V, REAL Re, REAL Pr, REAL tau, int *write,
               REAL del_trace, REAL del_inj, REAL del_streak, REAL del_vec){
  int i, j;
  REAL umax, vmax, deltu, deltv, deltRePr; 
  REAL t_trace, t_inj, t_streak, t_vec, t_neu; 

 /* delt satisfying CFL conditions */
 /*--------------------------------*/
  if(tau >= 1.0e-10){ /* else no time stepsize control */
    umax = 1.0e-10; vmax = 1.0e-10; 
    for(i=0; i<=imax+1; i++) for(j=1; j<=jmax+1; j++)
      if(fabs(U[i][j]) > umax)
        umax = fabs(U[i][j]);

    for(i=1; i<=imax+1; i++) for(j=0; j<=jmax+1; j++)
      if(fabs(V[i][j]) > vmax)
        vmax = fabs(V[i][j]);

    deltu = delx/umax; deltv = dely/vmax; 
    if(Pr < 1)	deltRePr = 1/(1/(delx*delx)+1/(dely*dely))*Re*Pr/2.;
    else	deltRePr = 1/(1/(delx*delx)+1/(dely*dely))*Re/2.;

    if(deltu<deltv) 
      if(deltu<deltRePr) *delt = deltu;
      else 	       *delt = deltRePr;
    else
      if(deltv<deltRePr) *delt = deltv;
      else 	       *delt = deltRePr;
    *delt = tau*(*delt); /* multiply by safety factor */
  }

 /* look if some data has to be written to a file in the next time step */ 
 /*---------------------------------------------------------------------*/
  *write = 0;
  t_neu = t + (*delt);
  t_trace = t_inj = t_streak = t_vec = t_neu + 1.0e+10;

  if( (int)(t/del_trace)!=(int)(t_neu/del_trace) ){
    t_trace = (int)(t_neu/del_trace) * del_trace;
    *write += 1;
  }
  if( (int)(t/del_inj)!=(int)(t_neu/del_inj) ){
    t_inj = (int)(t_neu/del_inj) * del_inj;
    *write += 2;
  }
  if( (int)(t/del_streak)!=(int)(t_neu/del_streak) ){
    t_streak = (int)(t_neu/del_streak) * del_streak;
    *write += 4;
  }
  if( (int)(t/del_vec)!=(int)(t_neu/del_vec) ){
    t_vec = (int)(t_neu/del_vec) * del_vec;
    *write += 8;
  }
}
