#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include "datadef.h"

/*------------------------------------------------------------------*/
/* Reading the Input parameters from "Inputfile"                    */
/*------------------------------------------------------------------*/
int READ_PARAMETER(char *Inputfile,char *problem,
		   REAL *xlength,REAL *ylength,int *imax,int *jmax,
		   REAL *delx,REAL *dely,
		   REAL *t_end,REAL *delt, REAL *tau,
		   REAL *del_trace,REAL *del_inj,REAL *del_streak,REAL *del_vec,
                   char *vecfile,char *tracefile,char *streakfile,
                   char *infile,char *outfile,
                   int *N,REAL *pos1x,REAL *pos1y,REAL *pos2x, REAL *pos2y,
		   int *itermax,REAL *eps,REAL *omg,REAL *gamma,int *p_bound,
		   REAL *Re,REAL *Pr,REAL *beta,REAL *GX,REAL *GY,
		   REAL *UI,REAL *VI,REAL *TI,
		   int *wW,int *wE,int *wN,int *wS){
  char c;
  FILE *fp;

  if ((fp = fopen(Inputfile, "r")) == NULL)
    {
     printf("Error while opening Inputfile %s\n",Inputfile);
     return(1);
    }

 /* Reading the type of the problem and checking if defined or not */
 /*----------------------------------------------------------------*/
  fscanf(fp, "%s", problem); 
     for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c)); /* find end of line */
  if( strcmp(problem, "convection") && strcmp(problem, "rayleigh") &&
      strcmp(problem, "fluidtrap") &&
      strcmp(problem, "dcavity") && strcmp(problem, "backstep") &&
      strcmp(problem, "plate") && strcmp(problem, "circle") &&
      strcmp(problem, "dam") && strcmp(problem, "drop") &&
      strcmp(problem, "molding") && strcmp(problem, "wave"))
  {
    printf("Problem %s not defined!\n", problem);
    printf("Choose dcavity\n");
    printf("       convection\n");
    printf("       rayleigh\n");
    printf("       fluidtrap\n");
    printf("       backstep\n");
    printf("       plate\n");
    printf("       circle\n");
    printf("       drop\n");
    printf("       dam\n");
    printf("       molding\n");
    printf("       wave\n");
    return(1);
  }

 /* reading "Inputfile" line for line */
 /*-----------------------------------*/
  fscanf(fp, "%s", infile); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, "%s", outfile); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
 
  fscanf(fp, INREAL, xlength); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, ylength); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, "%d", imax); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, "%d", jmax);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));

  (*delx) = (*xlength)/(*imax); (*dely) = (*ylength)/(*jmax);

  fscanf(fp, INREAL, t_end); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, delt); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, tau); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));

  fscanf(fp, INREAL, del_trace); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, del_inj); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, del_streak); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, del_vec); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));

  fscanf(fp, "%s", vecfile); 
     for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, "%s", tracefile); 
     for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, "%s", streakfile); 
     for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));

  fscanf(fp, "%d", N); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, pos1x); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, pos1y); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, pos2x); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, pos2y); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));

  fscanf(fp, "%d", itermax); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, eps); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, omg); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, gamma); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, "%d", p_bound);
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));

  fscanf(fp, INREAL, Re); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, Pr); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, beta); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, GX); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, GY); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, UI); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, VI); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, INREAL, TI); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));

  fscanf(fp, "%d", wW); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, "%d", wE); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, "%d", wN); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));
  fscanf(fp, "%d", wS); 
        for(fscanf(fp,"%c",&c);c!='\n';fscanf(fp,"%c",&c));

 /* Printing problem parameters */
 /*-----------------------------*/
  printf("\nProblem: %s\n",problem);
  printf("\nxlength=%.3e  ylength=%.3e  imax=%d  jmax=%d\n",
					(*xlength),(*ylength),(*imax),(*jmax));
  printf("delx=%.3e dely=%.3e\n", (*delx),(*dely));
  printf("t_end=%.3e delt=%.3e tau=%.3e\n", (*t_end),(*delt),(*tau));
  printf("del_trace=%.3e del_inj=%.3e del_streak=%.3e del_vec=%.3e\n",
                            (*del_trace),(*del_inj),(*del_streak),(*del_vec));
  printf("vecfile: %s tracefile: %s streakfile: %s\n",
                            vecfile,tracefile,streakfile);
  printf("infile: %s outfile: %s\n",infile,outfile);
  printf("N=%d pos1x=%.3e pos1y=%.3e pos2x=%.3e pos2y=%.3e\n",
                            (*N),(*pos1x),(*pos1y),(*pos2x),(*pos2y));
  printf("itermax=%d eps=%.3e omg=%.3e gamma=%.3e, p_bound=%d\n",
			   (*itermax),(*eps),(*omg),(*gamma),(*p_bound));
  printf("Re=%.3e Pr=%.3e beta=%.3e GX=%.3e GY=%.3e\n", 
					  (*Re),(*Pr),(*beta),(*GX),(*GY));
  printf("UI=%.3e VI=%.3e TI=%.3e\n", (*UI),(*VI),(*TI));
  printf("wW=%d wE=%d wN=%d wS=%d\n\n", (*wW),(*wE),(*wN),(*wS));

 /* Closing "Inputfile" */
  fclose(fp);  

  if ((*p_bound !=1) && (*p_bound != 2)){
    printf("p_bound must be 1 or 2!\n");
    return(1);
  }    
  if ((*wW > 4)||(*wW < 1)){
    printf("wW must be 1,2,3, or 4\n");
    return(1);
  }
  if ((*wE > 4)||(*wE < 1)){
    printf("wE must be 1,2,3, or 4\n");
    return(1);
  }
  if ((*wN > 4)||(*wN < 1)){
    printf("wN must be 1,2,3, or 4\n");
    return(1);
  }
  if ((*wS > 4)||(*wS < 1)){
    printf("wS must be 1,2,3, or 4\n");
    return(1);
  }
  if (((*wW==4) && (*wE!=4)) || (*wW!=4) && (*wE==4)){
    printf("Periodic boundary conditions need wW=wE=4\n");
    return(1);
  }
  if (((*wS==4) && (*wN!=4)) || (*wS!=4) && (*wN==4)){
    printf("Periodic boundary conditions need wS=wN=4\n");
    return(1);
  }

  return(0);
}

/*---------------------------------------------------------------*/
/* Setting the initial values for U,V,P, and TEMP                 */
/*---------------------------------------------------------------*/
void INIT_UVP(char *problem,
              REAL **U,REAL **V,REAL **P,REAL **TEMP,int imax,int jmax,
              REAL UI,REAL VI,REAL TI)
{
  int i,j;
 /* loop through all cells */
 /*------------------------*/
  for(i=0;i<=imax+1;i++)
    for(j=0;j<=jmax+1;j++)
      {
	U[i][j] = UI;
	V[i][j] = VI;
	P[i][j] = 0.;
	TEMP[i][j] = TI;
      }
  /* Set U=0.0 in the lower half for the flow past a backward facing step */
  /*----------------------------------------------------------------------*/
  if(strcmp(problem, "backstep")==0)
     for(i=0;i<=imax+1;i++)
        for(j=0;j<=jmax/2;j++)
           U[i][j] = 0.0;
}

/*----------------------------------------------------------------------*/
/* Initializing the integer array FLAG, dependent of the problem type   */
/*----------------------------------------------------------------------*/
void INIT_FLAG(char *problem,int **FLAG,int imax,int jmax,REAL delx,REAL dely,
               int *ibound)
{
  int i,j;
  int low,up;
  REAL mx,my,x,y,rad1;

                   /* boundary strip to C_B */
                   /*-----------------------*/
  for(i=0;i<=imax+1;i++)
    {
     FLAG[i][0]      = C_B;
     FLAG[i][jmax+1] = C_B;
    }
  for(j=1;j<=jmax;j++)
    {
     FLAG[0][j]      = C_B;
     FLAG[imax+1][j] = C_B;
    }
                   /* all inner cells fluid cells */
                   /*-----------------------------*/
  for(i=1;i<=imax;i++)
     for(j=1;j<=jmax;j++)
        FLAG[i][j] = C_F;

                   /* problem dependent obstacle cells in the interior */
                   /*--------------------------------------------------*/
  if(strcmp(problem,"fluidtrap")==0)
    {
     for(i=9*imax/22+1;i<=13*imax/22;i++)
       {
	for(j=1;j<=4*jmax/11;j++)
           FLAG[i][j] = C_B;
        for(j=8*jmax/11+1;j<=jmax;j++)
           FLAG[i][j] = C_B;          
       }
    }

  if(strcmp(problem,"plate")==0)
    {
                   /* flow past an inclined plate */
                   /*---------------------------------------*/
     low = 2*jmax/5;          /* lower and upper bound of the plate */
     up  = 3*jmax/5;
     FLAG[low][low]     = C_B;
     FLAG[low][low+1]   = C_B;
     FLAG[up][up-1]     = C_B;
     FLAG[up][up]       = C_B;
     for (i=low+1;i<=up-1;i++)
        for (j=i-1;j<=i+1;j++)
           FLAG[i][j] = C_B;
    }

  if(strcmp(problem,"backstep")==0 || strcmp(problem,"wave")==0)
    {
                    /* flow past a backward facing step */
                    /*----------------------------------*/
     for (i=1;i<=jmax;i++)
        for (j=1;j<=jmax/2;j++)
           FLAG[i][j] = C_B;
    }

  if(strcmp(problem,"circle")==0)
    {
                   /* flow past a cylinder/circle */
                   /*-----------------------------*/  
     mx = 20.0/41.0*jmax*dely;
     my = mx;
     rad1 = 5.0/41.0*jmax*dely;
     for (i=1;i<=imax;i++)
        for (j=1;j<=jmax;j++)   
          {
           x = (i-0.5)*delx;
           y = (j-0.5)*dely;
           if ((x-mx)*(x-mx)+(y-my)*(y-my) <= rad1*rad1)
              FLAG[i][j] = C_B;
	  }
    }

  if(strcmp(problem,"molding")==0)
    {
                   /* circular obstacle */
                   /*-------------------*/  
     mx = jmax*dely/2;
     my = jmax*dely/2;
     rad1 = jmax*dely/6;
     for (i=1;i<=imax;i++)
        for (j=1;j<=jmax;j++)   
          {
           x = (i-0.5)*delx;
           y = (j-0.5)*dely;
           if ((x-mx)*(x-mx)+(y-my)*(y-my) <= rad1*rad1)
              FLAG[i][j] = C_B;
	  }
    }

                     /* Printing the geometry of the fluid domain */
                     /*-------------------------------------------*/
  printf ("\nGeometry of the fluid domain:\n\n");
  for(j=jmax+1;j>=0;j--)
    {
     for(i=0;i<=imax+1;i++)
        if (!(FLAG[i][j] & C_F))
           printf("**");
        else      
           printf("  ");                                    
     printf ("\n");
    }
  printf ("\n");
  printf ("\n");
                    /* FLAGs for boundary cells */
                    /*--------------------------*/
  (*ibound) = 0;
  for(i=1;i<=imax;i++)
     for(j=1;j<=jmax;j++){
        if (!(FLAG[i][j] & C_F))
           (*ibound)++;
        FLAG[i][j] += ((FLAG[i-1][j] & C_F)*B_W + (FLAG[i+1][j] & C_F)*B_O +
                      (FLAG[i][j-1] & C_F)*B_S + (FLAG[i][j+1] & C_F)*B_N)/C_F;
        switch (FLAG[i][j]){
           case 0x0003:
           case 0x0007:
           case 0x000b:
           case 0x000c:
           case 0x000d:
           case 0x000e:
           case 0x000f:{           
                     printf("Illegal obstacle cell [%d][%d]\n",i,j);
                     exit(0);
                    }  
	 }
      }
}


/*-------------------------------------------------------------------------*/
/* Writing the values of U,V,P,TEMP,FLAG into a file for subsequent        */ 
/* calculations                                                            */
/*-------------------------------------------------------------------------*/
void WRITE_bin(REAL **U,REAL **V,REAL **P,REAL **TEMP,int **FLAG,
		      int imax,int jmax,char* file)
{
 int i;
 FILE *fp;

 fp = fopen(file, "wb"); 

 fwrite(&imax, sizeof(int), 1, fp);
 fwrite(&jmax, sizeof(int), 1, fp);

 for(i=0;i<=imax+1;i+=1)
   fwrite(U[i], sizeof(REAL), jmax+2, fp);
 for(i=0;i<=imax+1;i+=1)
   fwrite(V[i], sizeof(REAL), jmax+2, fp);
 for(i=0;i<=imax+1;i+=1)
   fwrite(P[i], sizeof(REAL), jmax+2, fp);
 for(i=0;i<=imax+1;i+=1)
   fwrite(TEMP[i], sizeof(REAL), jmax+2, fp);
 for(i=0;i<=imax+1;i+=1)
   fwrite(FLAG[i], sizeof(int), jmax+2, fp);
 fclose(fp);
}

/*-------------------------------------------------------------------------*/
/* Reading initial values from a file                                      */
/*-------------------------------------------------------------------------*/
int READ_bin(REAL **U,REAL **V,REAL **P,REAL **TEMP,int **FLAG,
		 int imax,int jmax,char* file)
{
 int i,j;
 FILE *fp;

 if(strcmp(file, "none") == 0) return(-1);

 if( (fp = fopen(file,"rb")) == NULL){
   printf("Error in READ_bin: File %s cannot be opened!\n", file);
   return(1);
 }

 fread(&i, sizeof(int), 1, fp);
 fread(&j, sizeof(int), 1, fp);

 if(i!=imax || j!=jmax){
    printf("ATTENTION: imax and jmax have wrong values in %s\n",file);
    printf("imax = %d  jmax = %d\n", i, j);
    return(2);
 }

 for(i=0;i<=imax+1;i+=1)
   fread(U[i], sizeof(REAL), jmax+2, fp);
 for(i=0;i<=imax+1;i+=1)
   fread(V[i], sizeof(REAL), jmax+2, fp);
 for(i=0;i<=imax+1;i+=1)
   fread(P[i], sizeof(REAL), jmax+2, fp);
 for(i=0;i<=imax+1;i+=1)
   fread(TEMP[i], sizeof(REAL), jmax+2, fp);
 for(i=0;i<=imax+1;i+=1)
   fread(FLAG[i], sizeof(int), jmax+2, fp);
 fclose(fp);

 return(0);
}


/*-------------------------------------------------------------------------*/
/* RMATRIX allocates memory for a [nrl,nrh]x[ncl,nch]-array of REAL-type   */
/*-------------------------------------------------------------------------*/
REAL **RMATRIX(int nrl,int nrh,int ncl,int nch)
{
  int i;
  REAL **m;
  if((m = (REAL**) malloc((unsigned) (nrh-nrl+1)*sizeof(double*))) == NULL)
     {	
      printf("no memory\n");
      exit(0);
     }
  m -= nrl;
  for(i=nrl;i<=nrh;i++)
    {
     if((m[i] = (REAL*) malloc((unsigned) (nch-ncl+1)*sizeof(double)))==NULL)
       {	
        printf("no memory\n");
        exit(0);
       }
     m[i] -= ncl;
    }
  return m;
} 

/*-------------------------------------------------------------------------*/
/* FREE_RMATRIX frees the memory of an array allocated with RMATRIX        */
/*-------------------------------------------------------------------------*/
void FREE_RMATRIX(REAL** m,int nrl,int nrh,int ncl,int nch)
{
  int i;
  for (i=nrh;i>=nrl;i--) free((void*) (m[i]+ncl));
  free((char*) (m+nrl));
}

/*-------------------------------------------------------------------------*/
/* IMATRIX allocates memory for a [nrl,nrh]x[ncl,nch]-array of integer-type*/
/*-------------------------------------------------------------------------*/
int **IMATRIX(int nrl,int nrh,int ncl,int nch)
{
  int i;
  int **m;
  if((m = (int**) malloc((unsigned) (nrh-nrl+1)*sizeof(int*))) == NULL)
     {	
      printf("no memory\n");
      exit(0);
     }
  m -= nrl;
  for(i=nrl;i<=nrh;i++)
    {
     if((m[i] = (int*) malloc((unsigned) (nch-ncl+1)*sizeof(int)))==NULL)
       {	
        printf("no memory\n");
        exit(0);
       }
     m[i] -= ncl;
    }
  return m;
} 

/*-------------------------------------------------------------------------*/
/* FREE_IMATRIX frees the memory of an array allocated with IMATRIX        */
/*-------------------------------------------------------------------------*/
void FREE_IMATRIX(int** m,int nrl,int nrh,int ncl,int nch)
{
  int i;
  for (i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
  free((char*) (m+nrl));
}
