struct particleline *INIT_PARTICLES(int *N,int imax,int jmax,
                                    REAL delx,REAL dely,
                                    int ppc,char *problem,REAL **U,REAL **V);

void SET_PART(struct particleline *Partline,REAL x,REAL y);

void MARK_CELLS(int **FLAG,int imax,int jmax,REAL delx,REAL dely,
                int *ifull,int *isurf,
                int N, struct particleline *Particlelines);

void SET_UVP_SURFACE(REAL **U,REAL **V,REAL **P,int **FLAG,REAL GX,REAL GY,
                     int imax,int jmax,REAL Re,REAL delx,REAL dely,REAL delt);
