#include <stdio.h>
#include <stdlib.h>
#include "datadef.h"
#include "visual.h"

/*---------------------------------------------------------------------------*/
/* Writing U,V,P,PSI, and ZETA into "vecfile" for visualization              */
/*---------------------------------------------------------------------------*/
void OUTPUTVEC_bin(REAL **U,REAL **V,REAL **P,REAL **TEMP,
		   REAL **PSI,REAL **ZETA,REAL **HEAT,int **FLAG,
                   REAL xlength,REAL ylength,int imax,int jmax,
                   char* vecfile)
{
 int i,j;
 float temp;
 FILE *fp;
 fp = fopen(vecfile, "wb");

 temp = xlength;
 fwrite(&temp, sizeof(float), 1, fp);
 temp = ylength;
 fwrite(&temp, sizeof(float), 1, fp);
 temp = imax;
 fwrite(&temp, sizeof(float), 1, fp);
 temp = jmax;
 fwrite(&temp, sizeof(float), 1, fp);

 for(j=1;j<=jmax;j+=1)
  for(i=1;i<=imax;i+=1){
   if( (FLAG[i][j] & C_F) && (FLAG[i][j] < C_E) ) 
      temp = (U[i][j]+U[i-1][j])/2.0;
   else			 temp = 0.0;
   fwrite(&temp, sizeof(float), 1, fp);
  }
 for(j=1;j<=jmax;j+=1)
  for(i=1;i<=imax;i+=1){
   if( (FLAG[i][j] & C_F) && (FLAG[i][j] < C_E) ) 
      temp = (V[i][j]+V[i][j-1])/2.0;
   else			 temp = 0.0;
   fwrite(&temp, sizeof(float), 1, fp);
  }
 for(j=1;j<=jmax;j+=1)
  for(i=1;i<=imax;i+=1){
   if( (FLAG[i][j] & C_F) && (FLAG[i][j] < C_E) ) 
      temp = P[i][j];
   else			 temp = 0.0;
   fwrite(&temp, sizeof(float), 1, fp);
  }
 for(j=1;j<=jmax;j+=1)
  for(i=1;i<=imax;i+=1){
   if( (FLAG[i][j] & C_F) && (FLAG[i][j] < C_E) ) 
        temp = TEMP[i][j];
   else			 temp = -0.5;

   temp = TEMP[i][j];
   fwrite(&temp, sizeof(float), 1, fp);
  }
 for(j=1;j<=jmax-1;j+=1)
  for(i=1;i<=imax-1;i+=1){
    temp = ZETA[i][j];
    fwrite(&temp, sizeof(float), 1, fp);
  }
 for(j=0;j<=jmax;j+=1)
  for(i=0;i<=imax;i+=1){
    temp = PSI[i][j];
    fwrite(&temp, sizeof(float), 1, fp);
  }
 for(j=0;j<=jmax;j+=1)
  for(i=0;i<=imax;i+=1){
    temp = HEAT[i][j];
    fwrite(&temp, sizeof(float), 1, fp);
  }

 fclose(fp);
}

/*-----------------------------------------------------*/
/* Computation of stream function and vorticity        */
/*-----------------------------------------------------*/
void COMPPSIZETA(REAL **U,REAL **V,REAL **PSI,REAL **ZETA,int **FLAG,
                int imax,int jmax,REAL delx,REAL dely)
{
 int i,j;

 /* Computation of the vorticity zeta at the upper right corner     */
 /* of cell (i,j) (only if the corner is surrounded by fluid cells) */
 /*-----------------------------------------------------------------*/
 for (i=1;i<=imax-1;i++)
    for (j=1;j<=jmax-1;j++)
        if( ((FLAG[i][j] & C_F) && (FLAG[i][j] < C_E))  &&
            ((FLAG[i+1][j] & C_F) && (FLAG[i+1][j] < C_E))  &&
            ((FLAG[i][j+1] & C_F) && (FLAG[i][j+1] < C_E))  &&
            ((FLAG[i+1][j+1] & C_F) && (FLAG[i+1][j+1] < C_E)) )
          ZETA[i][j] = (U[i][j+1]-U[i][j])/dely - (V[i+1][j]-V[i][j])/delx;
       else
          ZETA[i][j] = 0.0;

 /* Computation of the stream function at the upper right corner    */
 /* of cell (i,j) (only if bother lower cells are fluid cells)      */
 /*-----------------------------------------------------------------*/
 for (i=0;i<=imax;i++)
   {
    PSI[i][0] = 0.0;
    for(j=1;j<=jmax;j++)
        if( ((FLAG[i][j] & C_F) && (FLAG[i][j] < C_E))  ||
            ((FLAG[i+1][j] & C_F) && (FLAG[i+1][j] < C_E)) )
          PSI[i][j] = PSI[i][j-1] + U[i][j]*dely;
       else
          PSI[i][j] = PSI[i][j-1];
  } 
}

/*-----------------------------------------------------*/
/* Computation of the heat function                    */
/*-----------------------------------------------------*/
void COMP_HEAT(REAL **U,REAL **V,REAL **TEMP,REAL **HEAT,int **FLAG,
               REAL Re,REAL Pr,int imax,int jmax,REAL delx,REAL dely)
{
 int i,j;

 /* Computation at the upper right corner of cell (i,j) */
 /*-----------------------------------------------------*/
 for (i=0;i<=imax;i++)
   {
    HEAT[i][0] = 0.0;
    for(j=1;j<=jmax;j++)
       if( ((FLAG[i][j] & C_F) && (FLAG[i][j] < C_E))  ||
           ((FLAG[i+1][j] & C_F) && (FLAG[i+1][j] < C_E)) )
          HEAT[i][j] = HEAT[i][j-1] + 
                    dely*(U[i][j]*0.5*(1.0+TEMP[i+1][j]+TEMP[i][j])*Re*Pr-
                          (TEMP[i+1][j]-TEMP[i][j])/delx );
   }
}

/*-----------------------------------------------------------------*/
/* Allocate memory for a particle and set the coordinates to (x,y) */
/*-----------------------------------------------------------------*/
struct particle *PARTALLOC(REAL x, REAL y){
  struct particle *part;

  if((part=(struct particle *)malloc(sizeof(struct particle))) == NULL)
    {	
     printf("no memory\n");
     exit(0);
    }
  part->x = x; part->y = y;
  part->next = NULL;
  return( part );
}

/*------------------------------------------------------------------*/
/* Allocate memory for a particleline and set coordinates where     */
/* particles are injected.                                          */
/*------------------------------------------------------------------*/
struct particleline *SET_PARTICLES(int N,REAL pos1x,REAL pos1y,
				   	 REAL pos2x,REAL pos2y){
   int i;
   REAL hx,hy;
   struct particleline *Partlines;

   if((Partlines=(struct particleline *)
                 malloc((unsigned)(N) * sizeof(struct particleline))) == NULL)
     {	
      printf("no memory\n");
      exit(0);
     }
   Partlines -= 1;  

   if (N>=2)
     {
      hx  = (pos2x-pos1x)/(N-1);
      hy  = (pos2y-pos1y)/(N-1);
      for(i=1; i<=N; i++){
         Partlines[i].length = 0;
         Partlines[i].Particles = 
             PARTALLOC(pos1x+hx*(i-1),pos1y+hy*(i-1));
         Partlines[i].Particles->next = 
             PARTALLOC(pos1x+hx*(i-1),pos1y+hy*(i-1));
         Partlines[i].length++;
       }
     }
   return(Partlines);
}/*End SET_PARTICLES*/

/*-------------------------------------------------------------*/
/* Moving particles                                            */
/*-------------------------------------------------------------*/
void ADVANCE_PARTICLES(int imax,int jmax,REAL delx,REAL dely,REAL delt,
                       REAL **U,REAL **V,int **FLAG,
                       int N,struct particleline *Partlines)
{
  int i, j, k;
  REAL x, y, x1, y1, x2, y2, u, v; 
  struct particle *part,*help;

  for(k=1;k<=N;k++){
     for(part=Partlines[k].Particles; part->next != NULL; part=part->next){
        /* first element is only a dummy element           */
        /*-------------------------------------------------*/
        x = part->next->x; y = part->next->y;

        /* Computation of new x-coordinates by discretizing dx/dt=u */
        /*----------------------------------------------------------*/
        i = (int)(x/delx)+1;  j = (int)((y+0.5*dely)/dely)+1;

        x1 = (i-1)*delx;    y1 = ((j-1)-0.5)*dely;
        x2 = i*delx;        y2 = (j-0.5)*dely;

        /* bilinear interpolation */
        /*------------------------*/
        u = ((x2-x)*(y2-y)*U[i-1][j-1] +
  	     (x-x1)*(y2-y)*U[i][j-1]   +
	     (x2-x)*(y-y1)*U[i-1][j]   +
	     (x-x1)*(y-y1)*U[i][j])/delx/dely;

        /* Computation of new y-coordinates by discretizing dy/dt=v */
        /*----------------------------------------------------------*/
        i = (int)((x+0.5*delx)/delx)+1; j = (int)(y/dely)+1;

        x1 = ((i-1)-0.5)*delx;    y1 = (j-1)*dely;
        x2 = (i-0.5)*delx;        y2 = j*dely;

        /* bilinear interpolation */
        /*------------------------*/
        v = ((x2-x)*(y2-y)*V[i-1][j-1] +
	     (x-x1)*(y2-y)*V[i][j-1]   +
	     (x2-x)*(y-y1)*V[i-1][j]   +
	     (x-x1)*(y-y1)*V[i][j])/delx/dely;

        x += delt*u;   y += delt*v; 

        /*-------------------------------------*/
        /* determine new cell for the particle */
        /*-------------------------------------*/
        i = (int)(x/delx)+1;   j = (int)(y/dely)+1;

        /* if particle left the fluid domain, delete it */
        /*----------------------------------------------*/
        if (x>=imax*delx || y>=jmax*dely || x<=0 || y<=0){
          help = part->next->next;
          free(part->next);
          part->next = help;
          Partlines[k].length--;
	  if (help == NULL)
	    break; 
        }
        else{
          /*-------------------------------------*/
          /* special treatment if particle would */
	  /* be in an inner obstacle cell        */
          /*-------------------------------------*/
          if (FLAG[i][j] < C_F)
             ADVANCE_AT_BOUND(i,j,&x,&y,u,v,U,V,FLAG,delx,dely,delt);
            
          part->next->x = x; part->next->y = y;
         }
      }
   }
}/*End ADVANCE_PARTICLES*/


/*-------------------------------------------------------------------------*/
/* Computation of new particle location of a particle near a no-slip wall, */
/* guaranteeing, that the new position is not in the obstacle cell.        */
/* Here a modified interpolation algorithm is applied, using the fact that */
/* at no-skip walls, the velocity is not only given at the midpoint of the */
/* edge but on the whole edge                                              */
/*-------------------------------------------------------------------------*/
void ADVANCE_AT_BOUND(int i,int j,REAL *x,REAL *y,REAL u,REAL v,
                      REAL **U,REAL **V,int **FLAG,
                      REAL delx,REAL dely,REAL delt)
{
 int ialt,jalt;
 REAL xalt,yalt;
 REAL ul,ur,vo,vu;     
 REAL x1,x2,y1,y2;

 /* get old particle position */
 xalt = (*x)-delt*u;          yalt = (*y)-delt*v;
 ialt = (int)(xalt/delx)+1;   jalt = (int)(yalt/dely)+1; 

 if (i != ialt){      /* compute new x */
   if (FLAG[ialt+1][jalt] < C_F)
      ur = 0.0;      else{
      if (yalt>= (jalt-0.5)*dely)
         if (FLAG[ialt+1][jalt+1] < C_F){
            y2 = jalt *dely;
            ur = U[ialt][jalt]*(y2-yalt)*2.0/dely;
           }
         else{  
            y1 = (jalt-0.5)*dely;
            y2 = (jalt+0.5)*dely;
            ur = (U[ialt][jalt]*(y2-yalt)+U[ialt][jalt+1]*(yalt-y1))/dely;
	   }
      else      
         if (FLAG[ialt+1][jalt-1] < C_F){
            y1 = (jalt-1.0)*dely;
            ur = U[ialt][jalt]*(yalt-y1)*2.0/dely;
	   }
         else{ 
            y1 = (jalt-1.5)*dely;
            y2 = (jalt-0.5)*dely;
            ur = (U[ialt][jalt-1]*(y2-yalt)+U[ialt][jalt]*(yalt-y1))/dely;  
	   }
     }
   if (FLAG[ialt-1][jalt] < C_F)
      ul = 0.0;  
   else{
      if (yalt>= (jalt-0.5)*dely) 
         if (FLAG[ialt-1][jalt+1] < C_F){
            y2 = jalt *dely;
            ul = U[ialt-1][jalt]*(y2-yalt)*2.0/dely;
           }
         else{   
            y1 = (jalt-0.5)*dely;
            y2 = (jalt+0.5)*dely;
            ul = (U[ialt-1][jalt]*(y2-yalt)+U[ialt-1][jalt+1]*(yalt-y1))/dely;
	   }
      else       
         if (FLAG[ialt-1][jalt-1] < C_F){
            y1 = (jalt-1.0)*dely;
            ul = U[ialt-1][jalt]*(yalt-y1)*2.0/dely;
	   }
         else{ 
            y1 = (jalt-1.5)*dely;
            y2 = (jalt-0.5)*dely;
            ul = (U[ialt-1][jalt-1]*(y2-yalt)+U[ialt-1][jalt]*(yalt-y1))/dely;
	   }
     }
   u = (ul*(ialt*delx-xalt)+ur*(xalt-(ialt-1)*delx))/delx;
   (*x) = xalt+u*delt;
  }      /* end new x */

 if (j != jalt){        /* copute new y */
   if (FLAG[ialt][jalt+1] < C_F)
      vo = 0.0;   
   else{
      if (xalt>= (ialt-0.5)*delx) 
         if (FLAG[ialt+1][jalt+1] < C_F){
            x2 = ialt*delx;
            vo = V[ialt][jalt]*(x2-xalt)*2.0/delx;
           }
         else{  
            x1 = (ialt-0.5)*delx;
            x2 = (ialt+0.5)*delx;
            vo = (V[ialt][jalt]*(x2-xalt)+V[ialt+1][jalt]*(xalt-x1))/delx;
	   }
      else      
         if (FLAG[ialt-1][jalt+1] < C_F){
            x1 = (ialt-1.0)*delx;
            vo = V[ialt][jalt]*(xalt-x1)*2.0/delx;
	   }
         else{   
            x1 = (ialt-1.5)*delx;
            x2 = (ialt-0.5)*delx;
            vo = (V[ialt-1][jalt]*(x2-xalt)+V[ialt][jalt]*(xalt-x1))/delx;
	   }
     }
   if (FLAG[ialt][jalt-1] < C_F)
      vu = 0.0;  
   else{
      if (xalt>= (ialt-0.5)*delx) 
         if (FLAG[ialt+1][jalt-1] < C_F){
            x2 = ialt*delx;
            vu = V[ialt][jalt-1]*(x2-xalt)*2.0/delx;
           }
         else{   
            x1 = (ialt-0.5)*delx;
            x2 = (ialt+0.5)*delx;
            vu = (V[ialt][jalt-1]*(x2-xalt)+V[ialt+1][jalt-1]*(xalt-x1))/delx;
	   }
      else       
         if (FLAG[ialt-1][jalt-1] < C_F){
            x1 = (ialt-1.0)*delx;
            vu = V[ialt][jalt-1]*(xalt-x1)*2.0/delx;
	   }
         else{ 
            x1 = (ialt-1.5)*delx;
            x2 = (ialt-0.5)*delx;
            vu = (V[ialt-1][jalt-1]*(x2-xalt)+V[ialt][jalt-1]*(xalt-x1))/delx;
	   }
     }
   v = (vu*(jalt*dely-yalt)+vo*(yalt-(jalt-1)*dely))/dely;
   (*y) = yalt+v*delt;
  }    /* end new y */
}/* End ADVANCE_AT_BOUND */


/*--------------------------------------------------------------------*/
/* Injection of new particles for streaklines                         */
/*--------------------------------------------------------------------*/
void INJECT_PARTICLES(int N, struct particleline *Partlines){
  int i;
  struct particle *part;

  for(i=1; i<=N; i++){
    part = PARTALLOC(Partlines[i].Particles->x,Partlines[i].Particles->y);
    part->next = Partlines[i].Particles->next;
    Partlines[i].Particles->next = part;
    Partlines[i].length++;       
  }
}/*End INJECT_PARTICLES*/

/*--------------------------------------------------------------*/
/* Append particle positions to file "partfile" in ascii format */
/*--------------------------------------------------------------*/
void WRITE_PARTICLES(char *partfile, int N, struct particleline *Partlines){
  int i;
  FILE *fp;
  struct particle *part;

  fp = fopen(partfile,"ab");

  for(i=1; i<=N; i++){
    fprintf(fp,"%d\n",Partlines[i].length);
    for(part=Partlines[i].Particles; part->next != NULL; part=part->next)
      fprintf(fp,"%3.3f %3.3f\n", part->next->x, part->next->y);
  }

  fclose(fp);
}/*End  WRITE_PARTICLES*/


/*---------------------------------------------------------------------*/
/* Append particle positions to file "partfile" in binary format       */
/*---------------------------------------------------------------------*/
void WRITE_PARTICLES_bin(char *partfile, int N, struct particleline *Partlines)
{
 int i;
 FILE *fp;
 float temp, temp2[2];
 struct particle *part;

 fp = fopen(partfile, "ab");
 for(i=1; i<=N; i++){
    temp=Partlines[i].length;
    fwrite(&temp, sizeof(float), 1, fp);
    part=Partlines[i].Particles;
    for(; part->next != NULL; part=part->next){
      temp2[0]=part->next->x;
      temp2[1]=part->next->y;
      fwrite(temp2, sizeof(float), 2, fp);
    }
 }
 fclose(fp);
}

/*----------------------------------------------------------------------*/
/* Move particle positions and append them to a file if wanted          */
/*----------------------------------------------------------------------*/
void PARTICLE_TRACING(char* outputfile,REAL t,int imax,int jmax,
		      REAL delx,REAL dely,REAL delt,
                      REAL **U,REAL **V,int **FLAG,
                      int N, struct particleline *Partlines, int write)
{
  FILE *fp;

  if(t == 0){
    fp = fopen(outputfile, "wb");
    fprintf(fp,"%d\n%d\n%f\n%f\n%d\n", imax, jmax, delx, dely, N);
    fclose(fp);
    WRITE_PARTICLES(outputfile,N,Partlines);
  }

  ADVANCE_PARTICLES(imax,jmax,delx,dely,delt,U,V,FLAG,N,Partlines);

  if(write & 1)
    WRITE_PARTICLES(outputfile,N,Partlines);
    
}/*End PARTICLE_TRACING*/

/*---------------------------------------------------------------------*/
/* Move particles for streaklines, inject and write particle positions */
/*---------------------------------------------------------------------*/
void STREAKLINES(char* streakfile,int write,
                 int imax, int jmax, REAL delx, REAL dely,REAL delt,REAL t, 
                 REAL **U,REAL **V,int **FLAG,
                 int N, struct particleline *Partlines)
{
  FILE *fp;

  if(t==0){
    fp = fopen(streakfile, "wb");
    fprintf(fp, "%d\n", imax);
    fprintf(fp, "%d\n", jmax);
    fprintf(fp, "%g\n", delx);
    fprintf(fp, "%g\n", dely);
    fprintf(fp, "%d\n", N);
    fclose(fp);
    WRITE_PARTICLES_bin(streakfile,N,Partlines);
   }

  ADVANCE_PARTICLES(imax,jmax,delx,dely,delt,U,V,FLAG,N,Partlines);

  if(write & 2)  INJECT_PARTICLES(N,Partlines);
  if(write & 4)  WRITE_PARTICLES_bin(streakfile,N,Partlines);

}/*End STREAKLINES*/
