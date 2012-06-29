void COMP_TEMP(REAL **U,REAL **V,REAL **TEMP,int **FLAG,
	       int imax,int jmax,REAL delt,REAL delx,REAL dely,
	       REAL alpha,REAL Re,REAL Pr);

void COMP_FG(REAL **U,REAL **V,REAL **TEMP,REAL **F,REAL **G,int **FLAG,
	     int imax,int jmax,REAL delt,REAL delx,REAL dely,
	     REAL GX,REAL GY,REAL alpha,REAL Re,REAL beta);

void COMP_RHS(REAL **F,REAL **G,REAL **RHS,int **FLAG,int imax,int jmax,
              REAL delt,REAL delx,REAL dely);

int POISSON(REAL **P,REAL **RHS,int **FLAG,
            int imax,int jmax,REAL delx,REAL dely,
            REAL eps,int itermax,REAL omg,REAL *res,int ifull,int p_bound);

void ADAP_UV(REAL **U,REAL **V,REAL **F,REAL **G,REAL **P,int **FLAG,
             int imax,int jmax,REAL delt,REAL delx,REAL dely);

void COMP_delt(REAL *delt, REAL T, int imax, int jmax, REAL delx, REAL dely,
               REAL **U, REAL **V, REAL Re, REAL Pr, REAL tau, int *write,
               REAL del_trace, REAL del_inj, REAL del_streak, REAL del_vec);
