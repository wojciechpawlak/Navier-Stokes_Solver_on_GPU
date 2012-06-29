void OUTPUTVEC_bin(REAL **U,REAL **V,REAL **P,REAL **TEMP,
		   REAL **PSI,REAL **ZETA,REAL **HEAT,int **FLAG,
                   REAL,REAL,int im,int jm,char* outputfile);

void COMPPSIZETA(REAL **U,REAL **V,REAL **PSI,REAL **ZETA,int **FLAG,
                int imax,int jmax,REAL delx,REAL dely);

void COMP_HEAT(REAL **U,REAL **V,REAL **TEMP,REAL **HEAT,int **FLAG,
               REAL Re,REAL Pr,int imax,int jmax,REAL delx,REAL dely);

void ADVANCE_PARTICLES(int ibar,int jbar,REAL delx,REAL dely,REAL delt,
                       REAL **U,REAL **V,int **FLAG,
                       int N,struct particleline *Partlines);

void ADVANCE_AT_BOUND(int i,int j,REAL *x,REAL *y,REAL u,REAL v,
                      REAL **U,REAL **V,int **FLAG,
                      REAL delx,REAL dely,REAL delt);

void INJECT_PARTICLES(int N, struct particleline *Partlines);

void WRITE_PARTICLES(char *outputfile,int N,struct particleline *Partlines);

void WRITE_PARTICLES_bin(char *partfile,int N,struct particleline *Partlines);

void PARTICLE_TRACING(char* outputfile,REAL t,int imax,int jmax,
		      REAL delx,REAL dely,REAL delt,
                      REAL **U,REAL **V,int **FLAG,
                      int N, struct particleline *Partlines, int write);

void STREAKLINES(char* streakfile,int write,
                 int imax, int jmax, REAL delx, REAL dely,REAL delt,REAL t, 
                 REAL **U,REAL **V,int **FLAG,
                 int N, struct particleline *Partlines);

struct particle *PARTALLOC(REAL x, REAL y);

struct particleline *SET_PARTICLES(int N,REAL pos1x,REAL pos1y,
					 REAL pos2x,REAL pos2y);
