#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <CL/cl.h>

#include "alloc.h"
#include "utils.h"
#include "utils_opencl.h"
#include "datadef.h"
#include "init.h"
#include "boundary.h"
#include "uvp.h"
#include "visual.h"
#include "surface.h"

// Workgroup sizes
#define BLOCK_SIZE_2D_1 16
#define BLOCK_SIZE_2D_2 16
#define BLOCK_SIZE_1D 256

// Problem options
#define DCAV // used as benchmark problem (simplest)

// Parts of simulations
#define GPU
#define CPU
#define VERIFY
//#define PRINT
//#define LOG
//#define VISUAL

// Parts of GPU execution
#define ON_GPU_P0
#define ON_CPU_1_RELAX
#define ON_CPU_1_RES
#define ON_GPU_2_COPY
#define ON_CPU_2_RELAX
//#define ON_GPU_2_RES // commented - no residual computation

// Time calculation
#define TIME_WITH_MEM
//#define TIME_WITHOUT_MEM

// Name of source file with OpenCL kernels
const char *sourceFile = "kernels.cl";

int main(int argc, char *argv[])
{
#ifdef LOG
	
	FILE* log_file_cpu;
	log_file_cpu = fopen("logCPU.txt", "w");
	
	FILE* log_file_gpu;
	log_file_gpu = fopen("logGPU.txt", "w");
#endif
	
	/*
	 * Navier-Stokes Solver Initialization (both for CPU and GPU implementations)
	 */
	char problem[30];
	char infile[30], outfile[30];
	REAL xlength, ylength;
	int  imax, jmax;
	REAL delx, dely;
	REAL t_end, delt, tau;
	REAL del_trace, del_inj, del_streak, del_vec;
	char vecfile[30], tracefile[30], streakfile[30];
	int N;
	REAL pos1x, pos2x, pos1y, pos2y;
	int itermax;
	REAL eps, omg, gamma;
	int  p_bound;
	REAL Re, Pr, beta, GX, GY, UI, VI, TI;
	int wW, wE, wN, wS;
	int itersor=0, write;

	REAL t;
	REAL res;
	REAL **U, **V, **P, **PSI, **ZETA, **RHS, **F, **G, **TEMP, **HEAT;
	int  **FLAG;
	int  ppc, ifull=0, isurf=0, ibound=0;
	struct particleline *Particlelines;
	int init_case, cycle;

	/* READ the parameters of the problem.                */
	/* Stop if problem type or inputfile are not defined  */       
	/*----------------------------------------------------*/
	if (READ_PARAMETER(argv[1],problem,
			&xlength, &ylength, &imax, &jmax, &delx, &dely,
			&t_end, &delt, &tau, 
			&del_trace, &del_inj, &del_streak, &del_vec,
			vecfile,tracefile,streakfile,
			infile, outfile,
			&N, &pos1x, &pos1y, &pos2x, &pos2y,
			&itermax,&eps,&omg,&gamma,&p_bound,
			&Re, &Pr, &beta, &GX, &GY, &UI, &VI, &TI,
			&wW, &wE, &wN, &wS) != 0 ) {
		
		return(1); 
	}

	/* Allocate memory for the arrays */
	/*--------------------------------*/
	U    = RMATRIX(0, imax+1, 0, jmax+1);
	V    = RMATRIX(0, imax+1, 0, jmax+1);
	F    = RMATRIX(0, imax+1, 0, jmax+1);
	G    = RMATRIX(0, imax+1, 0, jmax+1);
	P    = RMATRIX(0, imax+1, 0, jmax+1);
	TEMP = RMATRIX(0, imax+1, 0, jmax+1);
	PSI  = RMATRIX(0, imax,	  0, jmax);
	ZETA = RMATRIX(1, imax-1, 1, jmax-1);
	HEAT = RMATRIX(0, imax,   0, jmax);
	RHS  = RMATRIX(0, imax+1, 0, jmax+1); 
	FLAG = IMATRIX(0, imax+1, 0, jmax+1);
	ppc  = 4;                             

	/* Read initial values from file "infile" */
	/*----------------------------------------*/
	init_case = READ_bin(U, V, P, TEMP, FLAG, imax, jmax, infile); 

	
	if (init_case > 0) {
		/* Error while reading "infile" */
		return(1); 
	}
	

	if (init_case < 0) {            
		/* Set initial values if        */
		/* "infile" is not specified    */
		INIT_UVP(problem, U, V, P, TEMP, imax, jmax, UI, VI, TI); 
		INIT_FLAG(problem, FLAG, imax, jmax, delx, dely, &ibound);
	}
#ifdef VISUAL
	/* Initialize particles for streaklines or particle tracing */
	/*----------------------------------------------------------*/
	//if (strcmp(streakfile, "none") || strcmp(tracefile, "none")) {
	//	Particlelines = SET_PARTICLES(N, pos1x, pos1y, pos2x, pos2y);
	//}

	/* Initialize particles for free boundary problems */      
	/*-------------------------------------------------*/
	//if (!strcmp(problem, "drop") || !strcmp(problem, "dam")) {
	//	Particlelines = INIT_PARTICLES(&N, imax, jmax, delx, dely,
	//		ppc, problem, U, V);
	//}
#endif

	/* Set initial values for boundary conditions				*/
	/* and specific boundary conditions, depending on "problem" */      
	/*----------------------------------------------------------*/
	SETBCOND(U, V, P, TEMP, FLAG, imax, jmax, wW, wE, wN, wS);
	SETSPECBCOND(problem, U, V, P, TEMP, imax, jmax, UI, VI);

#ifdef GPU

	/*
	 * OpenCL Initialization
	 */

	cl_int status;

	//-----------------------------------------------------
	// Discover and initialize the platforms
	//-----------------------------------------------------
	cl_uint numPlatforms = 0;
	cl_platform_id *platforms = NULL;

	// Retrieve the number of platforms
	status = clGetPlatformIDs(0, NULL, &numPlatforms);
	if (status != CL_SUCCESS) {
		printf("clGetPlatformIDs failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	// Make sure some platforms were found 
	if (numPlatforms == 0) {
		printf("No platforms detected.\n");
		exit(-1);
	}

	// Allocate enough space for each platform
	platforms = (cl_platform_id*) malloc(numPlatforms * sizeof(cl_platform_id));
	if (platforms == NULL) {
		perror("malloc failed");
		exit(-1);
	}

	// Fill in platforms
	status = clGetPlatformIDs(numPlatforms, platforms, NULL);
	if (status != CL_SUCCESS) {
		printf("clGetPlatformIDs failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	// Print basic information about each platform
	//printPlatforms(platforms, numPlatforms, status); 	// uncommment for platform information

	//-----------------------------------------------------
	// Discover and initialize the devices
	//----------------------------------------------------- 

	cl_uint numDevices = 0;
	cl_device_id *devices = NULL;

	// Retrive the number of devices present
	status = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_GPU, //TODO CL_DEVICE_TYPE_ALL for all types
		0, NULL, &numDevices);						
	if (status != CL_SUCCESS) {
		printf("clGetDeviceIDs failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	// Make sure some devices were found
	if (numDevices == 0) {
		printf("No devices detected.\n");
		exit(-1);
	}

	// Allocate enough space for each device
	devices = (cl_device_id*)malloc(numDevices * sizeof(cl_device_id));
	if(devices == NULL) {
		perror("malloc failed");
		exit(-1);
	}

	// Fill in devices
	status = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_GPU, // CL_DEVICE_TYPE_ALL
		numDevices,	devices, NULL);
	if(status != CL_SUCCESS) {
		printf("clGetDeviceIDs failed: %s\n", cluErrorString(status));
		exit(-1);
	} 

	// Print basic information about each device
	// printDevices(devices, numDevices, status); 	// uncomment to get device information
	// getDeviceInfo(devices, numDevices); 	// uncomment for detailed device information

	//-----------------------------------------------------
	// Create a context
	//----------------------------------------------------- 

	cl_context context = NULL;

	// Create a context and associate it with the devices
	context = clCreateContext(NULL, numDevices, devices, NULL, NULL, &status);
	if (status != CL_SUCCESS || context == NULL) {
		printf("clCreateContext failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	//-----------------------------------------------------
	// Create a command queue
	//----------------------------------------------------- 

	cl_command_queue cmdQueue = NULL;

	// Create a command queue and associate it with the device to execute on
	cmdQueue = clCreateCommandQueue(context, devices[0], 0, &status);
	if (status != CL_SUCCESS || cmdQueue == NULL) {
		printf("clCreateCommandQueue failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	//-----------------------------------------------------
	// Create device buffers
	//----------------------------------------------------- 

	// Computer the number of elements in grids used
	const int elements = (imax+2) * (jmax+2); // size of grid + ghost cells on boundaries

	// Compute the size of the data 
	size_t datasize = sizeof(REAL)*elements;
	size_t datasize_int = sizeof(int)*elements;

	// Check limitations of first available device
	char buf[100];
	printf("Device selected:\n");
	printf("Device %u: \n", 0);
	status = clGetDeviceInfo(devices[0], CL_DEVICE_VENDOR,
		sizeof(buf), buf, NULL);
	printf("\tDevice: %s\n", buf);
	status = clGetDeviceInfo(devices[0], CL_DEVICE_NAME,
		sizeof(buf), buf, NULL);
	printf("\tName: %s\n", buf);

	size_t buf2 = 0;
	status = clGetDeviceInfo(devices[0], CL_DEVICE_MAX_WORK_GROUP_SIZE,
		sizeof(size_t), &buf2, NULL);
	printf("\tMax Work Group Size: %u threads\n\n", buf2);
	if (status != CL_SUCCESS) {
		printf("clGetDeviceInfo failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	int MAX_THREADS_PR_WORKGROUP = buf2;
	int THREADS_PR_WORKGROUP = argc > 2 ? nextPow2(atoi(argv[2])) : MAX_THREADS_PR_WORKGROUP;

	printf("Threads per work group = %d.\n\n", THREADS_PR_WORKGROUP);

	// Compute the number of workgroups for 1D and 2D worksizes
	int NUM_WORKGROUPS_2D = (getGlobalSize(BLOCK_SIZE_2D_1, imax) * getGlobalSize(BLOCK_SIZE_2D_2, jmax))
		/ (BLOCK_SIZE_2D_1*BLOCK_SIZE_2D_2);
	int NUM_WORKGROUPS_1D = (getGlobalSize(BLOCK_SIZE_2D_1, imax) * getGlobalSize(BLOCK_SIZE_2D_2, jmax))
		/ BLOCK_SIZE_1D;

	printf("Work groups allocated for 2D block size (%d, %d , 1) = %d\n", BLOCK_SIZE_2D_1, BLOCK_SIZE_2D_2, NUM_WORKGROUPS_2D);
	printf("Work groups allocated for 1D block size (%d, 1, 1) = %d\n\n", BLOCK_SIZE_1D, NUM_WORKGROUPS_1D);

	//////////////////////////////////////////////////////////////////////////

	//-----------------------------------------------------
	// Allocate memory for host buffers for GPU execution
	//-----------------------------------------------------
	int *FLAG_h;
	REAL *U_h, *V_h, *TEMP_h, *F_h, *G_h, *RHS_h, *P_h, *p0_result_h;
#ifdef ON_GPU_2_RES 
	REAL *res_result_h;
#endif

	FLAG_h	= (int *) malloc(datasize_int);
	
	U_h		= (REAL *) malloc(datasize);
	V_h		= (REAL *) malloc(datasize);
	TEMP_h	= (REAL *) malloc(datasize);
	F_h		= (REAL *) malloc(datasize);
	G_h		= (REAL *) malloc(datasize);
	RHS_h	= (REAL *) malloc(datasize);
	P_h		= (REAL *) malloc(datasize);
	p0_result_h = (REAL *) malloc(NUM_WORKGROUPS_1D*sizeof(REAL));

#ifdef ON_GPU_2_RES
	res_result_h = (REAL *) malloc(NUM_WORKGROUPS_1D*sizeof(REAL));
#endif

	//TODO initialize them separately
	// Copy initial contents for host buffers from original arrays
	copy_array_int_2d_to_1d(FLAG,	FLAG_h, imax+2, jmax+2);

	copy_array_real_2d_to_1d(U,		U_h,	imax+2, jmax+2);
	copy_array_real_2d_to_1d(V,		V_h,	imax+2, jmax+2);
	copy_array_real_2d_to_1d(TEMP,	TEMP_h, imax+2, jmax+2);
	copy_array_real_2d_to_1d(P,		P_h,	imax+2, jmax+2);

	// does not need to copy values to F_h, G_h, RHS_h, because they are raw

	REAL **PSI_h, **ZETA_h, **HEAT_h;

	PSI_h	= RMATRIX(0, imax,	 0, jmax);
	ZETA_h	= RMATRIX(1, imax-1, 1, jmax-1);
	HEAT_h	= RMATRIX(0, imax,   0, jmax);

	copy_array_real(PSI,	PSI_h,	imax+1,	jmax+1);
	copy_array_real(ZETA,	ZETA_h, imax-1, jmax-1);
	copy_array_real(HEAT,	HEAT_h, imax+1,	jmax+1);

	//////////////////////////////////////////////////////////////////////////

	//-----------------------------------------------------
	// Allocate memory for data on device
	//-----------------------------------------------------

	// Input and Output arrays on the device
	cl_mem U_d;
	cl_mem V_d;
	cl_mem FLAG_d;
	cl_mem TEMP_d;
	cl_mem TEMP_new_d;
	cl_mem F_d;
	cl_mem G_d;
	cl_mem RHS_d;
	cl_mem P_d;
	cl_mem p0_result_d;

#ifdef ON_GPU_2_RES 
	cl_mem res_result_d;
#endif

#ifdef ON_GPU_2_RES_NAIVE
	cl_mem res_d;
#endif

	// Create a buffer object (U_d)
	U_d = clCreateBuffer(context, CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,
		datasize, U_h, &status);
	if (status != CL_SUCCESS || U_d == NULL) {
		printf("clCreateBuffer failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	// Create a buffer object (V_d)
	V_d = clCreateBuffer(context, CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,
		datasize, V_h, &status);
	if (status != CL_SUCCESS || V_d == NULL) {
		printf("clCreateBuffer failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	// Create a buffer object (FLAG_d)
	FLAG_d = clCreateBuffer(context, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR,
		datasize_int, FLAG_h, &status);
	if (status != CL_SUCCESS || FLAG_d == NULL) {
		printf("clCreateBuffer failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	// Create a buffer object (TEMP_d)
	TEMP_d = clCreateBuffer(context, CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR,
		datasize, TEMP_h, &status);
	if (status != CL_SUCCESS || TEMP_d == NULL) {
		printf("clCreateBuffer failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	// Create a buffer object (TEMP_new_d)
	TEMP_new_d = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
		datasize, NULL, &status);
	if (status != CL_SUCCESS || TEMP_new_d == NULL) {
		printf("clCreateBuffer failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	// Create a buffer object (F_d)
	F_d = clCreateBuffer(context, CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,
		datasize, F_h, &status);
	if (status != CL_SUCCESS || F_d == NULL) {
		printf("clCreateBuffer failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	// Create a buffer object (G_d)
	G_d = clCreateBuffer(context, CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,
		datasize, G_h, &status);
	if (status != CL_SUCCESS || G_d == NULL) {
		printf("clCreateBuffer failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	// Create a buffer object (RHS_d)
	RHS_d = clCreateBuffer(context, CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,
		datasize, RHS_h, &status);
	if (status != CL_SUCCESS || RHS_d == NULL) {
		printf("clCreateBuffer failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	// Create a buffer object (P_d)
	P_d = clCreateBuffer(context, CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,
		datasize, P_h, &status);
	if (status != CL_SUCCESS || P_d == NULL) {
		printf("clCreateBuffer failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	// Create a buffer object for vector of partial results
	// for initial pressure value computation
	p0_result_d = clCreateBuffer(context, CL_MEM_WRITE_ONLY|CL_MEM_COPY_HOST_PTR,
		NUM_WORKGROUPS_1D*sizeof(REAL), p0_result_h, &status);
	if (status != CL_SUCCESS || p0_result_d == NULL) {
		printf("clCreateBuffer failed: %s\n", cluErrorString(status));
		exit(-1);
	}

#ifdef ON_GPU_2_RES_NAIVE
	// Create a buffer object for scalar value for residual computation
	res_d = clCreateBuffer(context, CL_MEM_WRITE_ONLY|CL_MEM_COPY_HOST_PTR,
		sizeof(REAL), &res, &status);
	if (status != CL_SUCCESS || res_d == NULL) {
		printf("clCreateBuffer failed: %s\n", cluErrorString(status));
		exit(-1);
	}
#endif

#ifdef ON_GPU_RES 
	// Create a buffer object for vector of partial results
	// for residual computation
	res_result_d = clCreateBuffer(context, CL_MEM_WRITE_ONLY|CL_MEM_COPY_HOST_PTR,
		NUM_WORKGROUPS_1D*sizeof(REAL), res_result_h, &status);
	if (status != CL_SUCCESS || res_result_d == NULL) {
		printf("clCreateBuffer failed: %s\n", cluErrorString(status));
		exit(-1);
	}
#endif

	//////////////////////////////////////////////////////////////////////////

	//-----------------------------------------------------
	// Create and compile the program
	//----------------------------------------------------- 

	cl_program program;

	char *programSource;

	// Read in the source code of the program
	programSource = READ_kernelSource(sourceFile);

	// Create a program
	program = clCreateProgramWithSource(context, 1, (const char**)&programSource, 
		NULL, &status);
	if (status != CL_SUCCESS) {
		printf("clCreateProgramWithSource failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	// Build (compile & link) the program for the devices
	status = clBuildProgram(program, numDevices, devices, NULL, NULL, NULL);

	// Print build errors if any
	if (status != CL_SUCCESS) {
		printf("Program failed to build.\n");
		cl_build_status buildStatus;
		
		for (unsigned int i = 0; i < numDevices; i++) {
			clGetProgramBuildInfo(program, devices[i], CL_PROGRAM_BUILD_STATUS,
				sizeof(cl_build_status), &buildStatus, NULL);
			if (buildStatus == CL_SUCCESS) {
				continue;
			}

			char *buildLog;
			size_t buildLogSize;
			clGetProgramBuildInfo(program, devices[i], CL_PROGRAM_BUILD_LOG,
				0, NULL, &buildLogSize);
			buildLog = (char*)malloc(buildLogSize);
			if (buildLog == NULL) {
				perror("malloc");
				exit(-1);
			}

			clGetProgramBuildInfo(program, devices[i], CL_PROGRAM_BUILD_LOG,
				buildLogSize, buildLog, NULL);
			buildLog[buildLogSize-1] = '\0';

			printf("Device %u Build Log:\n%s\n", i, buildLog);   
			free(buildLog);
		}

		exit(0);
	}
	else {
		printf("No kernel build errors\n\n");
	}

	//////////////////////////////////////////////////////////////////////////
	
	/*
	 * GPU Execution
	 */

	// Create timers
	double timer_gpu;
	clock_t start_gpu, stop_gpu;

	//-----------------------------------------------------
	// Create the kernel
	//----------------------------------------------------- 

	cl_kernel TEMP_kernel = NULL;
	cl_kernel FG_kernel = NULL;
	cl_kernel RHS_kernel = NULL;
	cl_kernel POISSON_p0_kernel = NULL;
	cl_kernel POISSON_1_relaxation_kernel = NULL;
	cl_kernel POISSON_1_comp_res_kernel = NULL;
	cl_kernel POISSON_2_copy_boundary_kernel = NULL;
	cl_kernel POISSON_2_relaxation_kernel = NULL;
	cl_kernel POISSON_2_comp_res_kernel = NULL;
	cl_kernel ADAP_UV_kernel = NULL;
	cl_kernel SETBCOND_outer_kernel = NULL;
	cl_kernel SETBCOND_inner_kernel = NULL;
	cl_kernel SETSPECBCOND_kernel = NULL;

	// Create a kernel for computation of temperature TEMP
	TEMP_kernel = clCreateKernel(program, "COMP_TEMP_kernel", &status);
	if (status != CL_SUCCESS) {
		printf("clCreateKernel failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	// Create a kernel for computation of velocity vectors F, G
	FG_kernel = clCreateKernel(program, "COMP_FG_kernel", &status);
	if (status != CL_SUCCESS) {
		printf("clCreateKernel failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	// Create a kernel for computation of right-hand side for pressure equation RHS
	RHS_kernel = clCreateKernel(program, "COMP_RHS_kernel", &status);
	if (status != CL_SUCCESS) {
		printf("clCreateKernel failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	// Create a kernel for computation of initial pressure values for pressure P
	POISSON_p0_kernel = clCreateKernel(program, "POISSON_p0_kernel", &status);
	if (status != CL_SUCCESS) {
		printf("clCreateKernel failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	// Create a kernel for relaxation for pressure P (method 2)
	POISSON_1_relaxation_kernel = clCreateKernel(program, "POISSON_1_relaxation_kernel", &status);
	if (status != CL_SUCCESS) {
		printf("clCreateKernel failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	// Create a kernel for computation of norm of residual for pressure P (method 1)
	POISSON_1_comp_res_kernel = clCreateKernel(program, "POISSON_1_comp_res_kernel", &status);
	if (status != CL_SUCCESS) {
		printf("clCreateKernel failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	// Create a kernel for copying values a boundary cells for pressure P (method 2)
	POISSON_2_copy_boundary_kernel = clCreateKernel(program, "POISSON_2_copy_boundary_kernel", &status);
	if (status != CL_SUCCESS) {
		printf("clCreateKernel failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	// Create a kernel for relaxation for pressure P (method 2)
	POISSON_2_relaxation_kernel = clCreateKernel(program, "POISSON_2_relaxation_kernel", &status);
	if (status != CL_SUCCESS) {
		printf("clCreateKernel failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	// Create a kernel for computation of norm of residual for pressure P (method 2)
	POISSON_2_comp_res_kernel = clCreateKernel(program, "POISSON_2_comp_res_kernel", &status);
	if (status != CL_SUCCESS) {
		printf("clCreateKernel failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	// Create a kernel for computation of new velocity filed U, V
	ADAP_UV_kernel = clCreateKernel(program, "ADAP_UV_kernel", &status);
	if (status != CL_SUCCESS) {
		printf("clCreateKernel failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	// Create a kernel for computation of new velocity filed U, V
	SETBCOND_outer_kernel = clCreateKernel(program, "SETBCOND_outer_kernel", &status);
	if (status != CL_SUCCESS) {
		printf("clCreateKernel failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	// Create a kernel for computation of new velocity filed U, V
	SETBCOND_inner_kernel = clCreateKernel(program, "SETBCOND_inner_kernel", &status);
	if (status != CL_SUCCESS) {
		printf("clCreateKernel failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	// Create a kernel for computation of new velocity filed U, V
	SETSPECBCOND_kernel = clCreateKernel(program, "SETSPECBCOND_kernel", &status);
	if (status != CL_SUCCESS) {
		printf("clCreateKernel failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	//-----------------------------------------------------
	// Configure the work-item structure
	//----------------------------------------------------- 

	// Define an index space (global work size) of work items for execution.
	// A workgroup size (local work size) is not required, but can be used.
	size_t localWorkSize[] = {BLOCK_SIZE_2D_1, BLOCK_SIZE_2D_2};

	size_t globalWorkSize[] = {getGlobalSize(BLOCK_SIZE_2D_1, imax), getGlobalSize(BLOCK_SIZE_2D_2, jmax)};

	//size_t localWorkSize1d[] = {BLOCK_SIZE * BLOCK_SIZE};
	size_t localWorkSize1d[] = {BLOCK_SIZE_1D};

	size_t globalWorkSize1d[] = {getGlobalSize(BLOCK_SIZE_1D, imax) * getGlobalSize(BLOCK_SIZE_1D, jmax)};

	printf("Local Work Size 2D: %d %d\n", localWorkSize[0], localWorkSize[1]);
	printf("Global Work Size 2D: %d %d\n", globalWorkSize[0], globalWorkSize[1]);
	printf("Local Work Size 1D: %d\n", localWorkSize1d[0]);
	printf("Global Work Size 1D: %d\n\n", globalWorkSize1d[0]);

	//////////////////////////////////////////////////////////////////////////

	printf("GPU execution - start\n");
	start_gpu = clock();

	//-----------------------------------------------------
	// Write host data to device buffers
	//----------------------------------------------------- 

	status = clEnqueueWriteBuffer(cmdQueue, FLAG_d, CL_TRUE, 0, 
		datasize_int, FLAG_h, 0, NULL, NULL);

	status |= clEnqueueWriteBuffer(cmdQueue, U_d, CL_TRUE, 0, 
		datasize, U_h, 0, NULL, NULL);

	status |= clEnqueueWriteBuffer(cmdQueue, V_d, CL_TRUE, 0, 
		datasize, V_h, 0, NULL, NULL);

	status |= clEnqueueWriteBuffer(cmdQueue, TEMP_d, CL_TRUE, 0, 
		datasize, TEMP_h, 0, NULL, NULL);

	status |= clEnqueueWriteBuffer(cmdQueue, TEMP_new_d, CL_TRUE, 0, 
		datasize, TEMP_h, 0, NULL, NULL);

	status |= clEnqueueWriteBuffer(cmdQueue, F_d, CL_TRUE, 0, 
		datasize, F_h, 0, NULL, NULL);

	status |= clEnqueueWriteBuffer(cmdQueue, G_d, CL_TRUE, 0, 
		datasize, G_h, 0, NULL, NULL);

	status |= clEnqueueWriteBuffer(cmdQueue, RHS_d, CL_TRUE, 0, 
		datasize, RHS_h, 0, NULL, NULL);

	status |= clEnqueueWriteBuffer(cmdQueue, P_d, CL_TRUE, 0, 
		datasize, P_h, 0, NULL, NULL);
	
	status |= clEnqueueWriteBuffer(cmdQueue, p0_result_d, CL_TRUE, 0, 
		NUM_WORKGROUPS_1D*sizeof(REAL), p0_result_h, 0, NULL, NULL);

#ifdef ON_GPU_2_RES
	status |= clEnqueueWriteBuffer(cmdQueue, res_result_d, CL_TRUE, 0, 
		NUM_WORKGROUPS_1D*sizeof(REAL), res_result_h, 0, NULL, NULL);
#endif

#ifdef ON_GPU_2_RES_NAIVE
	status |= clEnqueueWriteBuffer(cmdQueue, res_result_d, CL_TRUE, 0, 
		NUM_WORKGROUPS_1D*sizeof(REAL), res_result_h, 0, NULL, NULL);
#endif

	if (status != CL_SUCCESS) {
		printf("clEnqueueWriteBuffer failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	//////////////////////////////////////////////////////////////////////////

	/*
	 * Main time loop
	 */
	for (t=0.0, cycle=0; t < t_end; t+=delt, cycle++) {
		COMP_delt(&delt, t, imax, jmax, delx, dely, U, V, Re, Pr, tau, &write,
			del_trace, del_inj, del_streak, del_vec);    


#ifndef DCAV
		//TODO make it work for other problems than dcav
		/* Determine fluid cells for free boundary problems */
		/* and set boundary values at free surface          */
		/*--------------------------------------------------*/
		//if (!strcmp(problem, "drop") || !strcmp(problem, "dam") ||
		//	!strcmp(problem, "molding") || !strcmp(problem, "wave")) {
		//	//TODO change FLAG
		//	MARK_CELLS(FLAG, imax, jmax, delx, dely, &ifull, &isurf,
		//		N, Particlelines);
		//	//TODO change U,V,P
		//	SET_UVP_SURFACE(U, V, P, FLAG, GX, GY,
		//		imax, jmax, Re, delx, dely, delt);
		//} else {
		//	ifull = imax*jmax-ibound;
		//}
#else
			ifull = imax*jmax-ibound;
#endif

		//---------------------------------------------------------------
		//  Set arguments for kernels and enqueue them for execution
		//---------------------------------------------------------------

		cl_event event;

		/* Compute new temperature */
		/*-------------------------*/

		// Associate the input and output buffers with the TEMP_kernel 
		status = clSetKernelArg(TEMP_kernel, 0, sizeof(cl_mem), &U_d);
		status |= clSetKernelArg(TEMP_kernel, 1, sizeof(cl_mem), &V_d);
		status |= clSetKernelArg(TEMP_kernel, 2, sizeof(cl_mem), &TEMP_d);
		status |= clSetKernelArg(TEMP_kernel, 3, sizeof(cl_mem), &TEMP_new_d);
		status |= clSetKernelArg(TEMP_kernel, 4, sizeof(cl_mem), &FLAG_d);
		status |= clSetKernelArg(TEMP_kernel, 5, sizeof(int), &imax);
		status |= clSetKernelArg(TEMP_kernel, 6, sizeof(int), &jmax);
		status |= clSetKernelArg(TEMP_kernel, 7, sizeof(REAL), &delt);
		status |= clSetKernelArg(TEMP_kernel, 8, sizeof(REAL), &delx);
		status |= clSetKernelArg(TEMP_kernel, 9, sizeof(REAL), &dely);
		status |= clSetKernelArg(TEMP_kernel, 10, sizeof(REAL), &gamma);
		status |= clSetKernelArg(TEMP_kernel, 11, sizeof(REAL), &Re);
		status |= clSetKernelArg(TEMP_kernel, 12, sizeof(REAL), &Pr);
		if (status != CL_SUCCESS) {
			printf("clSetKernelArg failed: %s\n", cluErrorString(status));
			exit(-1);
		}

		// Execute the TEMP_kernel kernel
		status = clEnqueueNDRangeKernel(cmdQueue, TEMP_kernel, 2, NULL, globalWorkSize, 
			localWorkSize, 0, NULL, &event);
		if(status != CL_SUCCESS) {
			printf("clEnqueueNDRangeKernel failed: %s\n", cluErrorString(status));
			exit(-1);
		}
		//clWaitForEvents(1, &event);

		/* Compute tentative velocity field (F, G) */
		/*----------------------------------------*/

		// Associate the input and output buffers with the FG_kernel 
		status = clSetKernelArg(FG_kernel, 0, sizeof(cl_mem), &U_d);
		status |= clSetKernelArg(FG_kernel, 1, sizeof(cl_mem), &V_d);
		status |= clSetKernelArg(FG_kernel, 2, sizeof(cl_mem), &TEMP_new_d);
		status |= clSetKernelArg(FG_kernel, 3, sizeof(cl_mem), &F_d);
		status |= clSetKernelArg(FG_kernel, 4, sizeof(cl_mem), &G_d);
		status |= clSetKernelArg(FG_kernel, 5, sizeof(cl_mem), &FLAG_d);
		status |= clSetKernelArg(FG_kernel, 6, sizeof(int), &imax);
		status |= clSetKernelArg(FG_kernel, 7, sizeof(int), &jmax);
		status |= clSetKernelArg(FG_kernel, 8, sizeof(REAL), &delt);
		status |= clSetKernelArg(FG_kernel, 9, sizeof(REAL), &delx);
		status |= clSetKernelArg(FG_kernel, 10, sizeof(REAL), &dely);
		status |= clSetKernelArg(FG_kernel, 11, sizeof(REAL), &GX);
		status |= clSetKernelArg(FG_kernel, 12, sizeof(REAL), &GY);
		status |= clSetKernelArg(FG_kernel, 13, sizeof(REAL), &gamma);
		status |= clSetKernelArg(FG_kernel, 14, sizeof(REAL), &Re);
		status |= clSetKernelArg(FG_kernel, 15, sizeof(REAL), &beta);
		if (status != CL_SUCCESS) {
			printf("clSetKernelArg failed: %s\n", cluErrorString(status));
			exit(-1);
		}

		// Execute the FG_kernel
		status = clEnqueueNDRangeKernel(cmdQueue, FG_kernel, 2, NULL, globalWorkSize, 
			localWorkSize, 0, NULL, &event);
		if(status != CL_SUCCESS) {
			printf("clEnqueueNDRangeKernel failed: %s\n", cluErrorString(status));
			exit(-1);
		}
		//clWaitForEvents(1, &event);

		/* Compute right hand side for pressure equation */
		/*-----------------------------------------------*/

		// Associate the input and output buffers with the RHS_kernel 
		status = clSetKernelArg(RHS_kernel, 0, sizeof(cl_mem), &F_d);
		status |= clSetKernelArg(RHS_kernel, 1, sizeof(cl_mem), &G_d);
		status |= clSetKernelArg(RHS_kernel, 2, sizeof(cl_mem), &RHS_d);
		status |= clSetKernelArg(RHS_kernel, 3, sizeof(cl_mem), &FLAG_d);
		status |= clSetKernelArg(RHS_kernel, 4, sizeof(int), &imax);
		status |= clSetKernelArg(RHS_kernel, 5, sizeof(int), &jmax);
		status |= clSetKernelArg(RHS_kernel, 6, sizeof(REAL), &delt);
		status |= clSetKernelArg(RHS_kernel, 7, sizeof(REAL), &delx);
		status |= clSetKernelArg(RHS_kernel, 8, sizeof(REAL), &dely);
		if (status != CL_SUCCESS) {
			printf("clSetKernelArg failed: %s\n", cluErrorString(status));
			exit(-1);
		}
		
		// Execute the RHS_kernel kernel
		status = clEnqueueNDRangeKernel(cmdQueue, RHS_kernel, 2, NULL, globalWorkSize, 
			localWorkSize, 0, NULL, &event);
		if(status != CL_SUCCESS) {
			printf("clEnqueueNDRangeKernel failed: %s\n", cluErrorString(status));
			exit(-1);
		}
		//clWaitForEvents(1, &event);

		//TODO because of CPU versions of Relaxation and Comp Res kernels
		status = clEnqueueReadBuffer(cmdQueue, RHS_d, CL_TRUE, 0,
			datasize, RHS_h, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			printf("clEnqueueReadBuffer failed: %s\n", cluErrorString(status));
			exit(-1);
		}

		/* Solve the pressure equation by successive over relaxation */
		/*-----------------------------------------------------------*/

		int iter = 0;

		if (ifull > 0) {
			int iimax = imax + 2;
			int jjmax = jmax + 2;

			REAL p0 = 0.0;

#ifdef ON_CPU_P0
			for (int i = 1; i <= iimax-2; i++) {
				for (int j = 1; j <= jjmax-2; j++) {
					if (FLAG_h[i*jjmax + j] & C_F) {
						p0 += P_h[i*jjmax + j]*P_h[i*jjmax + j];
					}
				}
			}	
#endif

#ifdef ON_GPU_P0
			// Associate the input and output buffers with the POISSON_p0_kernel
			status = clSetKernelArg(POISSON_p0_kernel, 0, sizeof(cl_mem), &P_d);
			status |= clSetKernelArg(POISSON_p0_kernel, 1, sizeof(cl_mem), &FLAG_d);
			status |= clSetKernelArg(POISSON_p0_kernel, 2, sizeof(int), &imax);
			status |= clSetKernelArg(POISSON_p0_kernel, 3, sizeof(int), &jmax);
			status |= clSetKernelArg(POISSON_p0_kernel, 4, sizeof(cl_mem), &p0_result_d);
			status |= clSetKernelArg(POISSON_p0_kernel, 5, sizeof(REAL)*localWorkSize1d[0], NULL);
			if (status != CL_SUCCESS) {
				printf("clSetKernelArg failed: %s\n", cluErrorString(status));
				exit(-1);
			}

			// Execute the POISSON_p0_kernel kernel
			status = clEnqueueNDRangeKernel(cmdQueue, POISSON_p0_kernel, 1,
				NULL, globalWorkSize1d, localWorkSize1d, 0, NULL, &event);
			if(status != CL_SUCCESS) {
				printf("clEnqueueNDRangeKernel failed: %s\n", cluErrorString(status));
				exit(-1);
			}
			//clWaitForEvents(1, &event);

			// Read the partial results from reduction kernel
			status = clEnqueueReadBuffer(cmdQueue, p0_result_d, CL_TRUE, 0,
				NUM_WORKGROUPS_1D*sizeof(REAL), p0_result_h, 0, NULL, NULL);
			if (status != CL_SUCCESS) {
				printf("clEnqueueReadBuffer failed: %s\n", cluErrorString(status));
				exit(-1);
			}

			// Compute the initial pressure value p0
			for (int group_id = 0; group_id < NUM_WORKGROUPS_1D; group_id++)
			{
				p0 += p0_result_h[group_id];
			}
#endif

			p0 = sqrt(p0/ifull);
			if (p0 < 0.0001)
				p0 = 1.0;

			REAL rdx2, rdy2;
			rdx2 = 1./delx/delx;
			rdy2 = 1./dely/dely;

			REAL beta_2;

			beta_2 = -omg/(2.0*(rdx2+rdy2));

			/* SOR-iteration */
			/*---------------*/
			for (iter = 1; iter <= itermax; iter++) {	
				
				res = 0.0;

				if (p_bound == 1) {

#ifdef ON_CPU_1_RELAX
					REAL beta_mod;
					
					for (int i = 1; i <= iimax-1; i++) {
						for (int j = 1; j <= jjmax-1; j++) {
							/* five point star for interior fluid cells */
							if (FLAG_h[i*jmax + j] == 0x001f) {
								P_h[i*jmax + j] = (1.-omg)*P_h[i*jmax + j] - 
									beta_2*((P_h[(i+1)*jmax + j]+P_h[(i-1)*jmax + j])*rdx2 +
									(P_h[i*jmax + j+1]+P_h[i*jmax + j-1])*rdy2 - RHS_h[i*jmax + j]);
							}
							/* modified star near boundary */
							else if ((FLAG_h[i*jmax + j] & C_F) && (FLAG_h[i*jmax + j] < 0x0100)) { 
								beta_mod = -omg/((eps_E+eps_W)*rdx2+(eps_N+eps_S)*rdy2);
								
								P_h[i*jmax + j] = (1.-omg)*P_h[i*jmax + j] -
									beta_mod*( (eps_E*P_h[(i+1)*jmax + j]+eps_W*P_h[(i-1)*jmax + j])*rdx2 +
									(eps_N*P_h[i*jmax + j+1]+eps_S*P_h[i*jmax + j-1])*rdy2 - RHS_h[i*jmax + j]);
							}
						}
					}
#endif

#ifdef ON_GPU_1_RELAX
					//TODO time dependencies
					/* relaxation for fluid cells */
					/*----------------------------*/
					// Associate the input and output buffers with the POISSON_1_relaxation_kernel 
					//status = clSetKernelArg(POISSON_1_relaxation_kernel, 0, sizeof(cl_mem), &P_d);
					//status |= clSetKernelArg(POISSON_1_relaxation_kernel, 1, sizeof(cl_mem), &RHS_d);
					//status |= clSetKernelArg(POISSON_1_relaxation_kernel, 2, sizeof(cl_mem), &FLAG_d);
					//status |= clSetKernelArg(POISSON_1_relaxation_kernel, 3, sizeof(int), &imax);
					//status |= clSetKernelArg(POISSON_1_relaxation_kernel, 4, sizeof(int), &jmax);
					//status |= clSetKernelArg(POISSON_1_relaxation_kernel, 5, sizeof(REAL), &delx);
					//status |= clSetKernelArg(POISSON_1_relaxation_kernel, 6, sizeof(REAL), &dely);
					//status |= clSetKernelArg(POISSON_1_relaxation_kernel, 7, sizeof(REAL), &omg);
					//if (status != CL_SUCCESS) {
					//	printf("clSetKernelArg failed: %s\n", cluErrorString(status));
					//	exit(-1);
					//}
					//
					//// Execute the POISSON_1_relaxation_kernel
					//status = clEnqueueNDRangeKernel(cmdQueue, POISSON_1_relaxation_kernel, 2,
					//	NULL, globalWorkSize, localWorkSize, 0, NULL, &event);
					//if(status != CL_SUCCESS) {
					//	printf("clEnqueueNDRangeKernel failed: %s\n", cluErrorString(status));
					//	exit(-1);
					//}

					////clWaitForEvents(1, &event);
#endif

#ifdef ON_CPU_1_RES
					/* computation of residual */
					/*-------------------------*/
					//REAL add;

					//res = 0.0;

					//for (int i = 1; i <= iimax-2; i++) {
					//	for (int j = 1; j <= jjmax-2; j++) {
					//		/* only fluid cells */
					//		if ((FLAG_h[i*jjmax + j] & C_F) && (FLAG_h[i*jjmax + j] < 0x0100)) {
					//			add = (	eps_E_h*(P_h[(i+1)*jjmax + j]-P_h[i*jjmax + j]) - 
					//					eps_W_h*(P_h[i*jjmax + j]-P_h[(i-1)*jjmax + j])) * rdx2
					//					+ (	eps_N_h*(P_h[i*jjmax + j+1]-P_h[i*jjmax + j]) -
					//						eps_S_h*(P_h[i*jjmax + j]-P_h[i*jjmax + j-1])) * rdy2
					//					- RHS_h[i*jjmax + j];

					//			res += add * add;	
					//		}
					//	}
					//}

					//res = sqrt(res/ifull)/p0;

					///* convergence? */
					//if (res < eps) {
					//	break;
					//}

#endif

#ifdef ON_GPU_1_RES
					//TODO because of dependencies above
					//status = clEnqueueWriteBuffer(cmdQueue, P_d, CL_FALSE, 0, 
					//	datasize, P_h, 0, NULL, NULL);
					//if (status != CL_SUCCESS) {
					//	printf("clEnqueueWriteBuffer failed: %s\n", cluErrorString(status));
					//	exit(-1);
					//}

					//// Associate the input and output buffers with the POISSON_1_comp_res_kernel 
					//status = clSetKernelArg(POISSON_1_comp_res_kernel, 0, sizeof(cl_mem), &P_d);
					//status |= clSetKernelArg(POISSON_1_comp_res_kernel, 1, sizeof(cl_mem), &RHS_d);
					//status |= clSetKernelArg(POISSON_1_comp_res_kernel, 2, sizeof(cl_mem), &FLAG_d);
					//status |= clSetKernelArg(POISSON_1_comp_res_kernel, 3, sizeof(int), &imax);
					//status |= clSetKernelArg(POISSON_1_comp_res_kernel, 4, sizeof(int), &jmax);
					//status |= clSetKernelArg(POISSON_1_comp_res_kernel, 5, sizeof(REAL), &delx);
					//status |= clSetKernelArg(POISSON_1_comp_res_kernel, 6, sizeof(REAL), &dely);
					//status |= clSetKernelArg(POISSON_1_comp_res_kernel, 7, sizeof(REAL), &res);
					//if (status != CL_SUCCESS) {
					//	printf("clSetKernelArg failed: %s\n", cluErrorString(status));
					//	exit(-1);
					//}

					//// Execute the POISSON_1_comp_res_kernel
					//status = clEnqueueNDRangeKernel(cmdQueue, POISSON_1_comp_res_kernel, 2,
					//	NULL, globalWorkSize, localWorkSize, 0, NULL, &event);
					//if(status != CL_SUCCESS) {
					//	printf("clEnqueueNDRangeKernel failed: %s\n", cluErrorString(status));
					//	exit(-1);
					//}

					////clWaitForEvents(1, &event);

					//res = sqrt(res/ifull)/p0;

					///* convergence? */
					//if (res < eps) {
					//	break;
					//}
#endif

				} else if (p_bound == 2) {

#ifdef ON_CPU_2_COPY
					//TODO because of time dependencies above
					//status = clEnqueueReadBuffer(cmdQueue, P_d, CL_TRUE, 0,
					//	datasize, P_h, 0, NULL, NULL);
					//if (status != CL_SUCCESS) {
					//	printf("clEnqueueReadBuffer failed: %s\n", cluErrorString(status));
					//	exit(-1);
					//}
					
					
					/* copy values at external boundary */
					/*----------------------------------*/
					for (int i = 1; i <= iimax-2; i++) {
						P_h[i*jjmax] = P_h[i*jjmax + 1];
						P_h[i*jjmax + (jjmax-1)] = P_h[i*jjmax + (jjmax-2)];
					}
					for (int j = 1; j <= jjmax-2; j++) {
						P_h[j] = P_h[jjmax + j];
						P_h[(iimax-1)*jjmax + j] = P_h[(iimax-2)*jjmax + j];
					}
					/* and at interior boundary cells */
					/*--------------------------------*/
					for (int i = 1; i <= iimax-2; i++) {
						for (int j = 1; j <= jjmax-2; j++) {
							if (FLAG_h[i*jjmax + j] >= B_N && FLAG_h[i*jjmax + j] <= B_SO) {
								switch (FLAG_h[i*jjmax + j]) {
								case B_N:	{ P_h[i*jjmax + j] = P_h[i*jjmax + j+1];							break; }
								case B_O:	{ P_h[i*jjmax + j] = P_h[(i+1)*jjmax + j];							break; }
								case B_S:	{ P_h[i*jjmax + j] = P_h[i*jjmax + j-1];							break; } 
								case B_W:	{ P_h[i*jjmax + j] = P_h[(i-1)*jjmax + j];							break; }
								case B_NO:	{ P_h[i*jjmax + j] = 0.5*(P_h[i*jjmax + j+1]+P_h[(i+1)*jjmax + j]);	break; }
								case B_SO:	{ P_h[i*jjmax + j] = 0.5*(P_h[i*jjmax + j-1]+P_h[(i+1)*jjmax + j]);	break; }
								case B_SW:	{ P_h[i*jjmax + j] = 0.5*(P_h[i*jjmax + j-1]+P_h[(i-1)*jjmax + j]);	break; }
								case B_NW:	{ P_h[i*jjmax + j] = 0.5*(P_h[i*jjmax + j+1]+P_h[(i-1)*jjmax + j]);	break; }
								default: break;
								}
							}
						}
					}	

					status = clEnqueueWriteBuffer(cmdQueue, P_d, CL_TRUE, 0, 
						datasize, P_h, 0, NULL, NULL);
					if (status != CL_SUCCESS) {
						printf("clEnqueueWriteBuffer failed: %s\n", cluErrorString(status));
						exit(-1);
					}
#endif

#ifdef ON_GPU_2_COPY
					/* copy values at external and interior boundary cells	*/
					/*------------------------------------------------------*/
					
					// Associate the input and output buffers with the POISSON_p0_kernel
					status = clSetKernelArg(POISSON_2_copy_boundary_kernel, 0, sizeof(cl_mem), &P_d);
					status |= clSetKernelArg(POISSON_2_copy_boundary_kernel, 1, sizeof(cl_mem), &FLAG_d);
					status |= clSetKernelArg(POISSON_2_copy_boundary_kernel, 2, sizeof(int), &imax);
					status |= clSetKernelArg(POISSON_2_copy_boundary_kernel, 3, sizeof(int), &jmax);
					if (status != CL_SUCCESS) {
						printf("clSetKernelArg failed: %s\n", cluErrorString(status));
						exit(-1);
					}
					
					// Execute the POISSON_2_copy_boundary_kernel
					status = clEnqueueNDRangeKernel(cmdQueue, POISSON_2_copy_boundary_kernel, 2,
						NULL, globalWorkSize, localWorkSize, 0, NULL, &event);
					if(status != CL_SUCCESS) {
						printf("clEnqueueNDRangeKernel failed: %s\n", cluErrorString(status));
						exit(-1);
					}
					//clWaitForEvents(1, &event);
#endif

#ifdef ON_CPU_2_RELAX					
					status |= clEnqueueReadBuffer(cmdQueue, P_d, CL_TRUE, 0,
						datasize, P_h, 0, NULL, NULL);
					if (status != CL_SUCCESS) {
						printf("clEnqueueReadBuffer failed: %s\n", cluErrorString(status));
						exit(-1);
					}

					//if (cycle == 1 && iter == 1)
					//	print_1darray_to_file(P_h, imax + 2, jmax + 2, "P_h_before.txt");

					/* relaxation for fluid cells */
					/*----------------------------*/
						
					for (int i = 1; i <= iimax-2; i++) {
						for (int j = 1; j <= jjmax-2; j++) {
							if ((FLAG_h[i*jjmax + j] & C_F) && (FLAG_h[i*jjmax + j] < 0x0100)) {
								P_h[i*jjmax + j] = (1.-omg)*P_h[i*jjmax + j] - 
									beta_2*((P_h[(i+1)*jjmax + j]+P_h[(i-1)*jjmax + j])*rdx2 +
									(P_h[i*jjmax + j+1]+P_h[i*jjmax + j-1])*rdy2 - RHS_h[i*jjmax + j]);
							}
						}
					}

					//if (cycle == 1 && iter == 1)
					//	print_1darray_to_file(P_h, imax + 2, jmax + 2, "P_h_after.txt");

					status = clEnqueueWriteBuffer(cmdQueue, P_d, CL_TRUE, 0, 
						datasize, P_h, 0, NULL, NULL);
					if (status != CL_SUCCESS) {
						printf("clEnqueueWriteBuffer failed: %s\n", cluErrorString(status));
						exit(-1);
					}

#endif

#ifdef ON_GPU_2_RELAX		
					// Associate the input and output buffers with the POISSON_2_relaxation_kernel 
					status = clSetKernelArg(POISSON_2_relaxation_kernel, 0, sizeof(cl_mem), &P_d);
					status |= clSetKernelArg(POISSON_2_relaxation_kernel, 1, sizeof(cl_mem), &RHS_d);
					status |= clSetKernelArg(POISSON_2_relaxation_kernel, 2, sizeof(cl_mem), &FLAG_d);
					status |= clSetKernelArg(POISSON_2_relaxation_kernel, 3, sizeof(int), &imax);
					status |= clSetKernelArg(POISSON_2_relaxation_kernel, 4, sizeof(int), &jmax);
					status |= clSetKernelArg(POISSON_2_relaxation_kernel, 5, sizeof(REAL), &rdx2);
					status |= clSetKernelArg(POISSON_2_relaxation_kernel, 6, sizeof(REAL), &rdy2);
					status |= clSetKernelArg(POISSON_2_relaxation_kernel, 7, sizeof(REAL), &beta_2);
					status |= clSetKernelArg(POISSON_2_relaxation_kernel, 8, sizeof(REAL), &omg);
					if (status != CL_SUCCESS) {
						printf("clSetKernelArg failed: %s\n", cluErrorString(status));
						exit(-1);
					}

					// Execute the POISSON_2_relaxation_kernel
					status = clEnqueueNDRangeKernel(cmdQueue, POISSON_2_relaxation_kernel, 2,
						NULL, globalWorkSize, localWorkSize, 0, NULL, &event);
					if(status != CL_SUCCESS) {
						printf("clEnqueueNDRangeKernel failed: %s\n", cluErrorString(status));
						exit(-1);
					}
					//clWaitForEvents(1, &event);

					//status |= clEnqueueReadBuffer(cmdQueue, P_d, CL_TRUE, 0,
					//	datasize, P_h, 0, NULL, NULL);
					//if (status != CL_SUCCESS) {
					//	printf("clEnqueueReadBuffer failed: %s\n", cluErrorString(status));
					//	exit(-1);
					//}
					//print_1darray_to_file(P_h, imax+2, jmax+2, "P_h_after.txt");
#endif

#ifdef ON_CPU_2_RES
					REAL add;

					res = 0.0;

					for (int i = 1; i <= iimax-2; i++) {
						for (int j = 1; j <= jjmax-2; j++) {
							if ((FLAG_h[i*jjmax + j] & C_F) && (FLAG_h[i*jjmax + j] < 0x0100)) {
								add =	(P_h[(i+1)*jjmax + j]-2*P_h[i*jjmax + j]+P_h[(i-1)*jjmax + j])*rdx2
										+ (P_h[i*jjmax + j+1]-2*P_h[i*jjmax + j]+P_h[i*jjmax + j-1])*rdy2
										- RHS_h[i*jjmax + j];

								res += add * add;
								
							}
						}
					}	

					res = sqrt(res/ifull)/p0;

					// convergence?
					if (res < eps) {
						break;
					}
#endif

#ifdef ON_GPU_2_RES_NAIVE
					// Associate the input and output buffers with the POISSON_2_comp_res_naive_kernel 
					status = clSetKernelArg(POISSON_2_comp_res_kernel, 0, sizeof(cl_mem), &P_d);
					status |= clSetKernelArg(POISSON_2_comp_res_kernel, 1, sizeof(cl_mem), &RHS_d);
					status |= clSetKernelArg(POISSON_2_comp_res_kernel, 2, sizeof(cl_mem), &FLAG_d);
					status |= clSetKernelArg(POISSON_2_comp_res_kernel, 3, sizeof(int), &imax);
					status |= clSetKernelArg(POISSON_2_comp_res_kernel, 4, sizeof(int), &jmax);
					status |= clSetKernelArg(POISSON_2_comp_res_kernel, 5, sizeof(REAL), &rdx2);
					status |= clSetKernelArg(POISSON_2_comp_res_kernel, 6, sizeof(REAL), &rdy2);
					status |= clSetKernelArg(POISSON_2_comp_res_kernel, 7, sizeof(cl_mem), &res_d);
					if (status != CL_SUCCESS) {
						printf("clSetKernelArg failed: %s\n", cluErrorString(status));
						exit(-1);
					}

					// Execute the POISSON_2_comp_res_kernel
					status = clEnqueueNDRangeKernel(cmdQueue, POISSON_2_comp_res_kernel, 1,
						NULL, globalWorkSize1d, localWorkSize1d, 0, NULL, &event);
					if(status != CL_SUCCESS) {
						printf("clEnqueueNDRangeKernel failed: %s\n", cluErrorString(status));
						exit(-1);
					}

					// Read the final result from reduction kernel
					status = clEnqueueReadBuffer(cmdQueue, res_d, CL_TRUE, 0,
						sizeof(REAL), &res, 0, NULL, NULL);
					if (status != CL_SUCCESS) {
						printf("clEnqueueReadBuffer failed: %s\n", cluErrorString(status));
						exit(-1);
					}

					res = sqrt(res/ifull)/p0;

					// convergence?
					if (res < eps) {
						break;
					}
#endif

#ifdef ON_GPU_2_RES
					/* computation of residual */
					/*-------------------------*/

					// Associate the input and output buffers with the POISSON_2_comp_res_kernel 
					status = clSetKernelArg(POISSON_2_comp_res_kernel, 0, sizeof(cl_mem), &P_d);
					status |= clSetKernelArg(POISSON_2_comp_res_kernel, 1, sizeof(cl_mem), &RHS_d);
					status |= clSetKernelArg(POISSON_2_comp_res_kernel, 2, sizeof(cl_mem), &FLAG_d);
					status |= clSetKernelArg(POISSON_2_comp_res_kernel, 3, sizeof(int), &imax);
					status |= clSetKernelArg(POISSON_2_comp_res_kernel, 4, sizeof(int), &jmax);
					status |= clSetKernelArg(POISSON_2_comp_res_kernel, 5, sizeof(REAL), &rdx2);
					status |= clSetKernelArg(POISSON_2_comp_res_kernel, 6, sizeof(REAL), &rdy2);
					status |= clSetKernelArg(POISSON_2_comp_res_kernel, 7, sizeof(cl_mem), &res_result_d);
					status |= clSetKernelArg(POISSON_2_comp_res_kernel, 8, sizeof(REAL)*localWorkSize1d[0], NULL);
					if (status != CL_SUCCESS) {
						printf("clSetKernelArg failed: %s\n", cluErrorString(status));
						exit(-1);
					}

					// Execute the POISSON_2_comp_res_kernel
					status = clEnqueueNDRangeKernel(cmdQueue, POISSON_2_comp_res_kernel, 1,
						NULL, globalWorkSize1d, localWorkSize1d, 0, NULL, &event);
					if(status != CL_SUCCESS) {
						printf("clEnqueueNDRangeKernel failed: %s\n", cluErrorString(status));
						exit(-1);
					}
					//clWaitForEvents(1, &event);

					// Read the partial results from reduction kernel
					status = clEnqueueReadBuffer(cmdQueue, res_result_d, CL_TRUE, 0,
						NUM_WORKGROUPS_1D*sizeof(REAL), res_result_h, 0, NULL, NULL);
					if (status != CL_SUCCESS) {
						printf("clEnqueueReadBuffer failed: %s\n", cluErrorString(status));
						exit(-1);
					}

					res = 0.0;

					// Compute residual from partial results
					for (int group_id = 0; group_id < NUM_WORKGROUPS_1D; group_id++)
					{
						res += res_result_h[group_id];
					}

					res = sqrt(res/ifull)/p0;

					// convergence?
					if (res < eps) {
						break;
					}
#endif					
				}
			}
		}

		// End of SOR iteration

		printf("t_end= %1.5g, t= %1.3e, delt= %1.1e, iterations %3d, res: %e, F-cells: %d, S-cells: %d, B-cells: %d\n",
			t_end, t+delt, delt, iter, res, ifull, isurf, ibound);  

		/* Compute the new velocity field */
		/*--------------------------------*/

		// Associate the input and output buffers with the ADAP_UV_kernel 
		status = clSetKernelArg(ADAP_UV_kernel, 0, sizeof(cl_mem), &U_d);
		status |= clSetKernelArg(ADAP_UV_kernel, 1, sizeof(cl_mem), &V_d);
		status |= clSetKernelArg(ADAP_UV_kernel, 2, sizeof(cl_mem), &F_d);
		status |= clSetKernelArg(ADAP_UV_kernel, 3, sizeof(cl_mem), &G_d);
		status |= clSetKernelArg(ADAP_UV_kernel, 4, sizeof(cl_mem), &P_d);
		status |= clSetKernelArg(ADAP_UV_kernel, 5, sizeof(cl_mem), &FLAG_d);
		status |= clSetKernelArg(ADAP_UV_kernel, 6, sizeof(int), &imax);
		status |= clSetKernelArg(ADAP_UV_kernel, 7, sizeof(int), &jmax);
		status |= clSetKernelArg(ADAP_UV_kernel, 8, sizeof(REAL), &delt);
		status |= clSetKernelArg(ADAP_UV_kernel, 9, sizeof(REAL), &delx);
		status |= clSetKernelArg(ADAP_UV_kernel, 10, sizeof(REAL), &dely);
		if (status != CL_SUCCESS) {
			printf("clSetKernelArg failed: %s\n", cluErrorString(status));
			exit(-1);
		}

		// Execute the ADAP_UV_kernel
		status = clEnqueueNDRangeKernel(cmdQueue, ADAP_UV_kernel, 2,
			NULL, globalWorkSize, localWorkSize, 0, NULL, &event);
		if(status != CL_SUCCESS) {
			printf("clEnqueueNDRangeKernel failed: %s\n", cluErrorString(status));
			exit(-1);
		}
		//clWaitForEvents(1, &event);

		/* Set boundary conditions */
		/*-------------------------*/

		// Associate the input and output buffers with the SETBCOND_outer_kernel
		status = clSetKernelArg(SETBCOND_outer_kernel, 0, sizeof(cl_mem), &U_d);
		status |= clSetKernelArg(SETBCOND_outer_kernel, 1, sizeof(cl_mem), &V_d);
		status |= clSetKernelArg(SETBCOND_outer_kernel, 2, sizeof(cl_mem), &P_d);
		status |= clSetKernelArg(SETBCOND_outer_kernel, 3, sizeof(cl_mem), &TEMP_new_d);
		status |= clSetKernelArg(SETBCOND_outer_kernel, 4, sizeof(int), &imax);
		status |= clSetKernelArg(SETBCOND_outer_kernel, 5, sizeof(int), &jmax);
		status |= clSetKernelArg(SETBCOND_outer_kernel, 6, sizeof(int), &wW);
		status |= clSetKernelArg(SETBCOND_outer_kernel, 7, sizeof(int), &wE);
		status |= clSetKernelArg(SETBCOND_outer_kernel, 8, sizeof(int), &wN);
		status |= clSetKernelArg(SETBCOND_outer_kernel, 9, sizeof(int), &wS);
		if (status != CL_SUCCESS) {
			printf("clSetKernelArg failed: %s\n", cluErrorString(status));
			exit(-1);
		}

		// Execute the SETBCOND_outer_kernel
		status = clEnqueueNDRangeKernel(cmdQueue, SETBCOND_outer_kernel, 2,
			NULL, globalWorkSize, localWorkSize, 0, NULL, &event);
		if(status != CL_SUCCESS) {
			printf("clEnqueueNDRangeKernel failed: %s\n", cluErrorString(status));
			exit(-1);
		}
		//clWaitForEvents(1, &event);

		// Associate the input and output buffers with the SETBCOND_inner_kernel 
		status = clSetKernelArg(SETBCOND_inner_kernel, 0, sizeof(cl_mem), &U_d);
		status |= clSetKernelArg(SETBCOND_inner_kernel, 1, sizeof(cl_mem), &V_d);
		status |= clSetKernelArg(SETBCOND_inner_kernel, 2, sizeof(cl_mem), &TEMP_d);
		status |= clSetKernelArg(SETBCOND_inner_kernel, 3, sizeof(cl_mem), &FLAG_d);
		status |= clSetKernelArg(SETBCOND_inner_kernel, 4, sizeof(int), &imax);
		status |= clSetKernelArg(SETBCOND_inner_kernel, 5, sizeof(int), &jmax);
		if (status != CL_SUCCESS) {
			printf("clSetKernelArg failed: %s\n", cluErrorString(status));
			exit(-1);
		}

		// Execute the SETBCOND_inner_kernel
		status = clEnqueueNDRangeKernel(cmdQueue, SETBCOND_inner_kernel, 2,
			NULL, globalWorkSize, localWorkSize, 0, NULL, &event);
		if(status != CL_SUCCESS) {
			printf("clEnqueueNDRangeKernel failed: %s\n", cluErrorString(status));
			exit(-1);
		}
		//clWaitForEvents(1, &event);

		/* Set special boundary conditions */
		/* Overwrite preset default values */
		/*---------------------------------*/

		// Associate the input and output buffers with the SETSPECBCOND_kernel 
		//status = clSetKernelArg(SETSPECBCOND_kernel, 0, sizeof(cl_mem), &problem);
		status |= clSetKernelArg(SETSPECBCOND_kernel, 0, sizeof(cl_mem), &U_d);
		status |= clSetKernelArg(SETSPECBCOND_kernel, 1, sizeof(cl_mem), &V_d);
		status |= clSetKernelArg(SETSPECBCOND_kernel, 2, sizeof(cl_mem), &P_d);
		status |= clSetKernelArg(SETSPECBCOND_kernel, 3, sizeof(cl_mem), &TEMP_d);
		status |= clSetKernelArg(SETSPECBCOND_kernel, 4, sizeof(int), &imax);
		status |= clSetKernelArg(SETSPECBCOND_kernel, 5, sizeof(int), &jmax);
		status |= clSetKernelArg(SETSPECBCOND_kernel, 6, sizeof(REAL), &UI);
		status |= clSetKernelArg(SETSPECBCOND_kernel, 7, sizeof(REAL), &VI);
		if (status != CL_SUCCESS) {
			printf("clSetKernelArg failed: %s\n", cluErrorString(status));
			exit(-1);
		}

		// Execute the SETSPECBCOND_kernel
		status = clEnqueueNDRangeKernel(cmdQueue, SETSPECBCOND_kernel, 2,
			NULL, globalWorkSize, localWorkSize, 0, NULL, &event);
		if(status != CL_SUCCESS) {
			printf("clEnqueueNDRangeKernel failed: %s\n", cluErrorString(status));
			exit(-1);
		}
		//clWaitForEvents(1, &event);

#ifndef DCAV
		//TODO other problems than dcav
		//if (!strcmp(problem, "drop") || !strcmp(problem, "dam") ||
		//	!strcmp(problem, "molding") || !strcmp(problem, "wave")) {
		//	SET_UVP_SURFACE(U_h, V_h, P_h, FLAG_h, GX, GY, imax, jmax, Re, delx, dely, delt);
		//}
#endif

#ifdef VISUAL
		//TODO Data Visualisation
		/* Write data for visualization */
		/*------------------------------*/
		//if ((write & 8) && strcmp(vecfile, "none")) {     
		//	COMPPSIZETA(U_h, V_h, PSI_h, ZETA_h, FLAG_h, imax, jmax, delx, dely);
		//	COMP_HEAT(U_h, V_h, TEMP_h, HEAT_h, FLAG_h, Re, Pr, imax, jmax, delx, dely);
		//	OUTPUTVEC_txt(U_h, V_h, P_h, TEMP_h, PSI_h, ZETA_h, HEAT_h, FLAG_h, xlength, ylength,
		//		imax, jmax, vecfile);
		//}
		//if ((write & 8) && strcmp(outfile, "none")) {
		//	WRITE_txt(U_h, V_h, P_h, TEMP_h, FLAG_h, imax, jmax, outfile);
		//}
		//if (strcmp(tracefile, "none")) {
		//	PARTICLE_TRACING(tracefile, t, imax, jmax, delx, dely, delt,
		//		U_h, V_h, FLAG_h, N, Particlelines, write);
		//}
		//if (strcmp(streakfile, "none")) {
		//	STREAKLINES(streakfile, write, imax, jmax, delx, dely, delt, t,
		//		U_h, V_h, FLAG_h, N, Particlelines);
		//}
#endif
	}  

	/*
	 * End of main time loop
	 */

	//////////////////////////////////////////////////////////////////////////
	
	//-----------------------------------------------------
	// Read the output buffers back to the host
	//----------------------------------------------------- 

	// Read the OpenCL output buffer (U_d) to the host output array (U_h)
	status = clEnqueueReadBuffer(cmdQueue, U_d, CL_TRUE, 0, 
		datasize, U_h, 0, NULL, NULL);
	// Read the OpenCL output buffer (V_d) to the host output array (V_h)
	status |= clEnqueueReadBuffer(cmdQueue, V_d, CL_TRUE, 0, 
		datasize, V_h, 0, NULL, NULL);
	// Read the OpenCL output buffer (TEMP_new_d) to the host output array (TEMP_h)
	status |= clEnqueueReadBuffer(cmdQueue, TEMP_new_d, CL_TRUE, 0, 
		datasize, TEMP_h, 0, NULL, NULL);
	// Read the OpenCL output buffer (F_d) to the host output array (F_h)
	status |= clEnqueueReadBuffer(cmdQueue, F_d, CL_TRUE, 0, 
		datasize, F_h, 0, NULL, NULL);
	// Read the OpenCL output buffer (G_d) to the host output array (G_h)
	status |= clEnqueueReadBuffer(cmdQueue, G_d, CL_TRUE, 0, 
		datasize, G_h, 0, NULL, NULL);
	// Read the OpenCL output buffer (RHS_d) to the host output array (RHS_h)
	status |= clEnqueueReadBuffer(cmdQueue, RHS_d, CL_TRUE, 0,
		datasize, RHS_h, 0, NULL, NULL);
	// Read the OpenCL output buffer (P_d) to the host output array (P_h)
	status |= clEnqueueReadBuffer(cmdQueue, P_d, CL_TRUE, 0,
		datasize, P_h, 0, NULL, NULL);
	if (status != CL_SUCCESS) {
		printf("clEnqueueReadBuffer failed: %s\n", cluErrorString(status));
		exit(-1);
	}

	stop_gpu = clock();
	printf("GPU execution - finish\n\n");

	timer_gpu = (double) (stop_gpu - start_gpu) / CLOCKS_PER_SEC;

	//////////////////////////////////////////////////////////////////////////

#ifdef VISUAL
	//TODO input for GPU
	//if (strcmp(vecfile,"none"))
	//{     
	//	COMPPSIZETA(U, V, PSI, ZETA, FLAG, imax, jmax, delx, dely);
	//	COMP_HEAT(U, V, TEMP, HEAT, FLAG, Re, Pr, imax, jmax, delx, dely);
	//	OUTPUTVEC_txt(U, V, P, TEMP, PSI, ZETA, HEAT, FLAG, xlength, ylength,
	//		imax, jmax, vecfile);
	//}
	//if (strcmp(outfile,"none")) {
	//	WRITE_txt(U, V, P, TEMP, FLAG, imax, jmax, outfile);
	//}
#endif

#ifdef PRINT
	print_1darray_to_file(U_h, imax+2, jmax+2, "U_h.txt");
	print_1darray_to_file(V_h, imax+2, jmax+2, "V_h.txt");
	print_1darray_to_file(TEMP_h, imax+2, jmax+2, "TEMP_h.txt");
	print_1darray_to_file(F_h, imax+2, jmax+2, "F_h.txt");
	print_1darray_to_file(G_h, imax+2, jmax+2, "G_h.txt");
	print_1darray_to_file(RHS_h, imax+2, jmax+2, "RHS_h.txt");
	print_1darray_to_file(P_h, imax+2, jmax+2, "P_h.txt");
	//print_1darray_to_file(PSI_h, imax+1, jmax+1, "PSI_h.txt");
	//print_1darray_to_file(ZETA_h, imax, jmax, "ZETA_h.txt");
	//print_1darray_to_file(HEAT_h, imax+1, jmax+1, "HEAT_h.txt");
	print_1darray_int_to_file(FLAG_h, imax+2, jmax+2, "FLAG_h.txt");
#endif

#endif

	//////////////////////////////////////////////////////////////////////////

#ifdef CPU

	/*
	 * CPU Execution
	 */

	// Create timers
	double timer_cpu;
	clock_t start_cpu, stop_cpu;

	printf("CPU execution - start\n");
	start_cpu = clock();

	res = 0.0;

	/*
	 * Main time loop
	 */
	for (t=0.0, cycle=0; t < t_end; t+=delt, cycle++) {
		COMP_delt(&delt, t, imax, jmax, delx, dely, U, V, Re, Pr, tau, &write,
			del_trace, del_inj, del_streak, del_vec);    


		//TODO other problems than dcav
		/* Determine fluid cells for free boundary problems */
		/* and set boundary values at free surface          */
		/*--------------------------------------------------*/
		//if (!strcmp(problem, "drop") || !strcmp(problem, "dam") ||
		//	!strcmp(problem, "molding") || !strcmp(problem, "wave")) {
		//	MARK_CELLS(FLAG, imax, jmax, delx, dely, &ifull, &isurf,
		//		N, Particlelines);
		//	SET_UVP_SURFACE(U, V, P, FLAG, GX, GY,
		//		imax, jmax, Re, delx, dely, delt);
		//} else {
			ifull = imax*jmax-ibound;
		//}

		/* Compute new temperature */
		/*-------------------------*/
		COMP_TEMP(U, V, TEMP, FLAG, imax, jmax,
			delt, delx, dely, gamma, Re, Pr);

		/* Compute tentative velocity field (F, G) */
		/*----------------------------------------*/
		COMP_FG(U, V, TEMP, F, G, FLAG, imax, jmax,
			delt, delx, dely, GX, GY, gamma, Re, beta);

		/* Compute right hand side for pressure equation */
		/*-----------------------------------------------*/
		COMP_RHS(F, G, RHS, FLAG, imax, jmax, delt, delx, dely);

		/* Solve the pressure equation by successive over relaxation */
		/*-----------------------------------------------------------*/
		if (ifull > 0) {
			itersor = POISSON(P, RHS, FLAG, imax, jmax, delx, dely,
			eps, itermax, omg, &res, ifull, p_bound);
		}

		printf("t_end= %1.5g, t= %1.3e, delt= %1.1e, iterations %3d, res: %e, F-cells: %d, S-cells: %d, B-cells: %d\n",
			t_end, t+delt, delt, itersor, res, ifull, isurf, ibound);  

		/* Compute the new velocity field */
		/*--------------------------------*/
		ADAP_UV(U, V, F, G, P, FLAG, imax, jmax, delt, delx, dely);
		
		/* Set boundary conditions */
		/*-------------------------*/
		SETBCOND(U, V, P, TEMP, FLAG, imax, jmax, wW, wE, wN, wS);
		/* Set special boundary conditions */
		/* Overwrite preset default values */
		/*---------------------------------*/
		SETSPECBCOND(problem, U, V, P, TEMP, imax, jmax, UI, VI);

		//TODO other problems than dcav
		//	if (!strcmp(problem, "drop") || !strcmp(problem, "dam") ||
		//		!strcmp(problem, "molding") || !strcmp(problem, "wave")) {
		//		SET_UVP_SURFACE(U, V, P, FLAG, GX, GY, imax, jmax, Re, delx, dely, delt);
		//	}

#ifdef VISUAL

		//TODO Data Visualisation
		//	/* Write data for visualization */
		//	/*------------------------------*/
		//	if ((write & 8) && strcmp(vecfile, "none")) {     
		//		COMPPSIZETA(U, V, PSI, ZETA, FLAG, imax, jmax, delx, dely);
		//		COMP_HEAT(U, V, TEMP, HEAT, FLAG, Re, Pr, imax, jmax, delx, dely);
		//		OUTPUTVEC_txt(U, V, P, TEMP, PSI, ZETA, HEAT, FLAG, xlength, ylength,
		//			imax, jmax, vecfile);
		//	}
		//	if ((write & 8) && strcmp(outfile, "none")) {
		//		WRITE_txt(U, V, P, TEMP, FLAG, imax, jmax, outfile);
		//	}
		//	if (strcmp(tracefile, "none")) {
		//		PARTICLE_TRACING(tracefile, t, imax, jmax, delx, dely, delt,
		//			U, V, FLAG,	N, Particlelines, write);
		//	}
		//	if (strcmp(streakfile, "none")) {
		//		STREAKLINES(streakfile, write, imax, jmax, delx, dely, delt, t,
		//			U, V, FLAG, N, Particlelines);
		//	}

#endif
	}           
	
	/*
	 * End of main time loop
	 */

	stop_cpu = clock();
	printf("CPU execution - finish\n\n");

	timer_cpu = (double) (stop_cpu - start_cpu) / CLOCKS_PER_SEC;

	//////////////////////////////////////////////////////////////////////////

#ifdef VISUAL

	//if (strcmp(vecfile, "none"))
	//{     
	//	COMPPSIZETA(U, V, PSI, ZETA, FLAG, imax, jmax, delx, dely);
	//	COMP_HEAT(U, V, TEMP, HEAT, FLAG, Re, Pr, imax, jmax, delx, dely);
	//	OUTPUTVEC_txt(U, V, P, TEMP, PSI, ZETA, HEAT, FLAG, xlength, ylength,
	//		imax, jmax, vecfile);
	//}
	//if (strcmp(outfile, "none")) {
	//	WRITE_txt(U, V, P, TEMP, FLAG, imax, jmax, outfile);
	//}
#endif

	#ifdef PRINT

	print_array_to_file(U, imax+2, jmax+2, "U.txt");
	print_array_to_file(V, imax+2, jmax+2, "V.txt");
	print_array_to_file(TEMP, imax+2, jmax+2, "TEMP.txt");
	print_array_to_file(F, imax+2, jmax+2, "F.txt");
	print_array_to_file(G, imax+2, jmax+2, "G.txt");
	print_array_to_file(RHS, imax+2, jmax+2, "RHS.txt");
	print_array_to_file(P, imax+2, jmax+2, "P.txt");
	//print_array_to_file(PSI, imax+1, jmax+1, "PSI.txt");
	//print_array_to_file(ZETA, imax, jmax, "ZETA.txt");
	//print_array_to_file(HEAT, imax+1, jmax+1, "HEAT.txt");
	print_array_int_to_file(FLAG, imax+2, jmax+2, "FLAG.txt");

	#endif

#endif

	/*--------------------------------------------------------------------------------*/

#ifdef VERIFY

	/*
	 * Verification
	 */

	printf("U: ");
	if (compare_array(U, U_h, imax+2, jmax+2)) {
		printf("PASSED\n");
	} else {
		printf("FAILED\n");
	}

	printf("V: ");
	if (compare_array(V, V_h, imax+2, jmax+2)) {
		printf("PASSED\n");
	} else {
		printf("FAILED\n");
	}

	printf("TEMP: ");
	if (compare_array(TEMP, TEMP_h, imax+2, jmax+2)) {
		printf("PASSED\n");
	} else {
		printf("FAILED\n");
	}

	printf("F: ");
	if (compare_array(F, F_h, imax+2, jmax+2)) {
		printf("PASSED\n");
	} else {
		printf("FAILED\n");
	}

	printf("G: ");
	if (compare_array(G, G_h, imax+2, jmax+2)) {
		printf("PASSED\n");
	} else {
		printf("FAILED\n");
	}

	printf("RHS: ");
	if (compare_array(RHS, RHS_h, imax+2, jmax+2)) {
		printf("PASSED\n");
	} else {
		printf("FAILED\n");
	}

	printf("P: ");
	if (compare_array(P, P_h, imax+2, jmax+2)) {
		printf("PASSED\n");
	} else {
		printf("FAILED\n");
	}

	// Print timings
	printf("\n");
	printf("  GPU time    : %.4f (ms)\n", timer_gpu);
	printf("  CPU time    : %.4f (ms)\n", timer_cpu);
	printf("  Speedup %.2fx\n\n", timer_cpu/timer_gpu);
	
#endif

	//////////////////////////////////////////////////////////////////////////

	//-----------------------------------------------------
	// Release resources
	//----------------------------------------------------- 

#ifdef GPU

	// Free OpenCL resources
	clReleaseMemObject(FLAG_d);
	clReleaseMemObject(U_d);
	clReleaseMemObject(V_d);
	clReleaseMemObject(TEMP_d);
	clReleaseMemObject(TEMP_new_d);
	clReleaseMemObject(F_d);
	clReleaseMemObject(G_d);
	clReleaseMemObject(P_d);
	clReleaseMemObject(p0_result_d);

#ifdef ON_GPU_2_RES
	clReleaseMemObject(res_result_d);
#endif

#ifdef ON_GPU_2_RES_NAIVE
	clReleaseMemObject(res_d);
#endif
	clReleaseKernel(FG_kernel);
	clReleaseKernel(TEMP_kernel);
	clReleaseProgram(program);
	clReleaseCommandQueue(cmdQueue);
	clReleaseContext(context);


	// Free host resources
	free(FLAG_h);
	free(U_h);
	free(V_h);
	free(TEMP_h);
	free(F_h);
	free(G_h);
	free(RHS_h);
	free(P_h);
	free(p0_result_h);
#ifdef ON_GPU_2_RES
	free(res_result_h);
#endif

	FREE_RMATRIX(PSI_h,		0, imax,	0, jmax);
	FREE_RMATRIX(ZETA_h,	1, imax-1,	1, jmax-1);
	FREE_RMATRIX(HEAT_h,	0, imax,	0, jmax);

	free(platforms);
	free(devices);

#endif

	FREE_RMATRIX(U,		0, imax+1,	0, jmax+1);
	FREE_RMATRIX(V,		0, imax+1,	0, jmax+1);
	FREE_RMATRIX(F,		0, imax+1,	0, jmax+1);
	FREE_RMATRIX(G,		0, imax+1,	0, jmax+1);
	FREE_RMATRIX(P,		0, imax+1,	0, jmax+1);
	FREE_RMATRIX(TEMP,	0, imax+1,	0, jmax+1);
	FREE_RMATRIX(PSI,	0, imax,	0, jmax);
	FREE_RMATRIX(ZETA,	1, imax-1,	1, jmax-1);
	FREE_RMATRIX(HEAT,	0, imax,	0, jmax);
	FREE_RMATRIX(RHS,	0, imax+1,	0, jmax+1);
	FREE_IMATRIX(FLAG,	0, imax+1,	0, jmax+1);

#ifdef LOG
	fclose(log_file_cpu);
	fclose(log_file_gpu);
#endif

	printf("End of program\n");
	return(0);
 }
