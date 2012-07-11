#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

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

#define BLOCK_SIZE 16

#define GPU
#define CPU
#define VERIFY
#define PRINT

int main(int argc, char *argv[])
{
	/*
	 * Navier-Stokes Solver Initialization
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
	// STEP 1: Discover and initialize the platforms
	//-----------------------------------------------------
	cl_uint numPlatforms = 0;
	cl_platform_id *platforms = NULL;

	// Retrieve the number of platforms
	status = clGetPlatformIDs(0, NULL, &numPlatforms);
	if (status != CL_SUCCESS) {
		printf("clGetPlatformIDs failed\n");
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
		printf("clGetPlatformIDs failed\n");
		exit(-1);
	}

	// Print basic information about each platform
	printPlatforms(platforms, numPlatforms, status);

	//-----------------------------------------------------
	// STEP 2: Discover and initialize the devices
	//----------------------------------------------------- 

	cl_uint numDevices = 0;
	cl_device_id *devices = NULL;

	// Retrive the number of devices present
	status = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_GPU, // CL_DEVICE_TYPE_ALL
		0, NULL, &numDevices);						
	if (status != CL_SUCCESS) {
		printf("clGetDeviceIDs failed\n");
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
		printf("clGetDeviceIDs failed\n");
		exit(-1);
	} 

	// Print basic information about each device
	printDevices(devices, numDevices, status);

	getDeviceInfo(devices, numDevices);

	//-----------------------------------------------------
	// STEP 3: Create a context
	//----------------------------------------------------- 

	cl_context context = NULL;

	// Create a context and associate it with the devices
	context = clCreateContext(NULL, numDevices, devices, NULL, NULL, &status);
	if (status != CL_SUCCESS || context == NULL) {
		printf("clCreateContext failed\n");
		exit(-1);
	}

	//-----------------------------------------------------
	// STEP 4: Create a command queue
	//----------------------------------------------------- 

	cl_command_queue cmdQueue = NULL;

	// Create a command queue and associate it with the device you 
	// want to execute on
	cmdQueue = clCreateCommandQueue(context, devices[0], 0, &status);
	if (status != CL_SUCCESS || cmdQueue == NULL) {
		printf("clCreateCommandQueue failed\n");
		exit(-1);
	}

	//-----------------------------------------------------
	// STEP 5: Create device buffers
	//----------------------------------------------------- 

	// Elements in each array - U, V, F, G
	const int elements = (imax+2) * (jmax+2);

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
		printf("clGetDeviceInfo failed\n");
		exit(-1);
	}

	int MAX_THREADS_PR_WORKGROUP = buf2;
	int THREADS_PR_WORKGROUP = argc > 2 ? atoi(argv[2]) : MAX_THREADS_PR_WORKGROUP;

	printf("Threads per work group = %d. \n", THREADS_PR_WORKGROUP);

	int WORKGROUPS = (elements + THREADS_PR_WORKGROUP)/THREADS_PR_WORKGROUP; 
	printf("Work groups allocated = %d\n\n", WORKGROUPS);

	// Input and Output arrays on the device
	cl_mem U_d;
	cl_mem V_d;
	cl_mem FLAG_d;
	cl_mem TEMP_d;
	cl_mem TEMP_new_d;
	cl_mem F_d;
	cl_mem G_d;
	cl_mem RHS_d;

	// Allocate memory for data on device
	
	// Create a buffer object (U_d)
	U_d = clCreateBuffer(context, CL_MEM_READ_ONLY,
		datasize, NULL, &status);
	if (status != CL_SUCCESS || U_d == NULL) {
		printf("clCreateBuffer failed\n");
		exit(-1);
	}

	// Create a buffer object (V_d)
	V_d = clCreateBuffer(context, CL_MEM_READ_ONLY,
		datasize, NULL, &status);
	if (status != CL_SUCCESS || V_d == NULL) {
		printf("clCreateBuffer failed\n");
		exit(-1);
	}

	// Create a buffer object (FLAG_d)
	FLAG_d = clCreateBuffer(context, CL_MEM_READ_ONLY,
		datasize_int, NULL, &status);
	if (status != CL_SUCCESS || FLAG_d == NULL) {
		printf("clCreateBuffer failed\n");
		exit(-1);
	}

	// Create a buffer object (TEMP_d)
	TEMP_d = clCreateBuffer(context, CL_MEM_READ_ONLY,
		datasize, NULL, &status);
	if (status != CL_SUCCESS || TEMP_d == NULL) {
		printf("clCreateBuffer failed\n");
		exit(-1);
	}

	// Create a buffer object (TEMP_new_d)
	TEMP_new_d = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
		datasize, NULL, &status);
	if (status != CL_SUCCESS || TEMP_new_d == NULL) {
		printf("clCreateBuffer failed\n");
		exit(-1);
	}

	// Create a buffer object (F_d)
	F_d = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
		datasize, NULL, &status);
	if (status != CL_SUCCESS || F_d == NULL) {
		printf("clCreateBuffer failed\n");
		exit(-1);
	}

	// Create a buffer object (G_d)
	G_d = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
		datasize, NULL, &status);
	if (status != CL_SUCCESS || G_d == NULL) {
		printf("clCreateBuffer failed\n");
		exit(-1);
	}

	// Create a buffer object (RHS_d)
	RHS_d = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
		datasize, NULL, &status);
	if (status != CL_SUCCESS || RHS_d == NULL) {
		printf("clCreateBuffer failed\n");
		exit(-1);
	}

	int *FLAG_h;
	REAL *U_h, *V_h, *TEMP_h, *F_h, *G_h,  *RHS_h;

	REAL **P_h, **PSI_h, **ZETA_h, **HEAT_h;

	FLAG_h	= (int *) malloc(datasize_int);
	
	U_h		= (REAL *) malloc(datasize);
	V_h		= (REAL *) malloc(datasize);
	TEMP_h	= (REAL *) malloc(datasize);
	F_h		= (REAL *) malloc(datasize);
	G_h		= (REAL *) malloc(datasize);
	RHS_h	= (REAL *) malloc(datasize);

	P_h		= RMATRIX(0, imax+1, 0, jmax+1);
	PSI_h	= RMATRIX(0, imax,	 0, jmax);
	ZETA_h	= RMATRIX(1, imax-1, 1, jmax-1);
	HEAT_h	= RMATRIX(0, imax,   0, jmax);

	copy_array_int_2d_to_1d(FLAG,	FLAG_h, imax+2, jmax+2);

	copy_array_real_2d_to_1d(U,		U_h,	imax+2, jmax+2);
	copy_array_real_2d_to_1d(V,		V_h,	imax+2, jmax+2);
	copy_array_real_2d_to_1d(TEMP,	TEMP_h, imax+2, jmax+2);

	copy_array_real(P,		P_h,	imax+2, jmax+2);
	copy_array_real(PSI,	PSI_h,	imax+1,	jmax+1);
	copy_array_real(ZETA,	ZETA_h, imax-1, jmax-1);
	copy_array_real(HEAT,	HEAT_h, imax+1,	jmax+1);

	//-----------------------------------------------------
	// STEP 7: Create and compile the program
	//----------------------------------------------------- 

	cl_program program;

	char *programSource;
	const char *sourceFile = "kernels.cl";

	// Read in the source code of the program
	programSource = READ_kernelSource(sourceFile);

	// Create a program.
	// The 'source' string is the code from the FG_kernel.cl file.
	program = clCreateProgramWithSource(context, 1, (const char**)&programSource, 
		NULL, &status);
	if (status != CL_SUCCESS) {
		printf("clCreateProgramWithSource failed\n");
		exit(-1);
	}

	// Build (compile & link) the program for the devices.
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
		printf("No build errors\n");
	}


	
	/*
	 * GPU Execution
	 */

	// Create timers
	double timer_gpu;
	clock_t start_gpu, stop_gpu;

	//-----------------------------------------------------
	// STEP 8: Create the kernel
	//----------------------------------------------------- 

	cl_kernel TEMP_kernel = NULL;
	cl_kernel FG_kernel = NULL;
	cl_kernel RHS_kernel = NULL;

	// Create a kernel for computation of temperature TEMP
	TEMP_kernel = clCreateKernel(program, "COMP_TEMP_kernel", &status);
	if (status != CL_SUCCESS) {
		printf("clCreateKernel failed\n");
		exit(-1);
	}

	// Create a kernel for computation of velocity vectors F, G
	FG_kernel = clCreateKernel(program, "COMP_FG_kernel", &status);
	if (status != CL_SUCCESS) {
		printf("clCreateKernel failed\n");
		exit(-1);
	}

	// Create a kernel for computation of right-hand side for pressure equation RHS
	RHS_kernel = clCreateKernel(program, "COMP_RHS_kernel", &status);
	if (status != CL_SUCCESS) {
		printf("clCreateKernel failed\n");
		exit(-1);
	}



	//-----------------------------------------------------
	// STEP 10: Configure the work-item structure
	//----------------------------------------------------- 

	// Define an index space (global work size) of work items for execution.
	// A workgroup size (local work size) is not required, but can be used.
	size_t localWorkSize[2];
	localWorkSize[0] = BLOCK_SIZE;
	localWorkSize[1] = BLOCK_SIZE;

	size_t globalWorkSize[2];
	//globalWorkSize[0] = WORKGROUPS*THREADS_PR_WORKGROUP;
	globalWorkSize[0] = THREADS_PR_WORKGROUP;
	globalWorkSize[1] = THREADS_PR_WORKGROUP;

	start_gpu = clock();

	/*
	 * Main time loop
	 */
	t=0.0;
	//for (t=0.0, cycle=0; t < t_end; t+=delt, cycle++) {
		COMP_delt(&delt, t, imax, jmax, delx, dely, U, V, Re, Pr, tau, &write,
			del_trace, del_inj, del_streak, del_vec);    

		//TODO other problems than dcav
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

		//-----------------------------------------------------
		// STEP 9: Set arguments for kernels
		//----------------------------------------------------- 

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
			printf("clSetKernelArg failed\n");
			exit(-1);
		}


		//-----------------------------------------------------
		// STEP 6: Write host data to device buffers
		//----------------------------------------------------- 

		status = clEnqueueWriteBuffer(cmdQueue, FLAG_d, CL_FALSE, 0, 
			datasize_int, FLAG_h, 0, NULL, NULL);

		status |= clEnqueueWriteBuffer(cmdQueue, U_d, CL_FALSE, 0, 
			datasize, U_h, 0, NULL, NULL);

		status |= clEnqueueWriteBuffer(cmdQueue, V_d, CL_FALSE, 0, 
			datasize, V_h, 0, NULL, NULL);

		status |= clEnqueueWriteBuffer(cmdQueue, TEMP_d, CL_FALSE, 0, 
			datasize, TEMP_h, 0, NULL, NULL);

		status |= clEnqueueWriteBuffer(cmdQueue, TEMP_new_d, CL_FALSE, 0, 
			datasize, TEMP_h, 0, NULL, NULL);

		status |= clEnqueueWriteBuffer(cmdQueue, F_d, CL_FALSE, 0, 
			datasize, F_h, 0, NULL, NULL);

		status |= clEnqueueWriteBuffer(cmdQueue, G_d, CL_FALSE, 0, 
			datasize, G_h, 0, NULL, NULL);

		status |= clEnqueueWriteBuffer(cmdQueue, RHS_d, CL_FALSE, 0, 
			datasize, RHS_h, 0, NULL, NULL);

		if (status != CL_SUCCESS) {
			printf("clEnqueueWriteBuffer failed\n");
			exit(-1);
		}

		/* Compute new temperature */
		/*-------------------------*/
		cl_event event;

		// Execute the kernel
		status = clEnqueueNDRangeKernel(cmdQueue, TEMP_kernel, 2, NULL, globalWorkSize, 
			localWorkSize, 0, NULL, &event);
		if(status != CL_SUCCESS) {
			printf("clEnqueueNDRangeKernel failed: %s\n", cluErrorString(status));
			exit(-1);
		}

		clWaitForEvents(1, &event);

		// Read the OpenCL output buffer (TEMP_new_d) to the host output array (TEMP_h)
		//status = clEnqueueReadBuffer(cmdQueue, TEMP_new_d, CL_TRUE, 0,
		//	datasize, TEMP_h, 0, NULL, NULL);
		//if (status != CL_SUCCESS) {
		//	printf("clEnqueueReadBuffer failed\n");
		//	exit(-1);
		//}

		/* Compute tentative velocity field (F, G) */
		/*----------------------------------------*/
		//-----------------------------------------------------
		// STEP 11: Enqueue the kernel for execution
		//-----------------------------------------------------

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
			printf("clSetKernelArg failed\n");
			exit(-1);
		}
		
		//status = clEnqueueWriteBuffer(cmdQueue, TEMP_d, CL_FALSE, 0, 
		//	datasize, TEMP_h, 0, NULL, NULL);
		//if (status != CL_SUCCESS) {
		//	printf("clEnqueueWriteBuffer failed\n");
		//	exit(-1);
		//}

		// Execute the kernel
		status = clEnqueueNDRangeKernel(cmdQueue, FG_kernel, 2, NULL, globalWorkSize, 
			localWorkSize, 0, NULL, &event);
		if(status != CL_SUCCESS) {
			printf("clEnqueueNDRangeKernel failed: %s\n", cluErrorString(status));
			exit(-1);
		}

		clWaitForEvents(1, &event);

		//-----------------------------------------------------
		// STEP 12: Read the output buffer back to the host
		//----------------------------------------------------- 

		// Read the OpenCL output buffer (F_d) to the host output array (F_h)
		//status = clEnqueueReadBuffer(cmdQueue, F_d, CL_TRUE, 0,
		//	datasize, F_h, 0, NULL, NULL);

		// Read the OpenCL output buffer (G_d) to the host output array (G_h)
		//status |= clEnqueueReadBuffer(cmdQueue, G_d, CL_TRUE, 0,
		//	datasize, G_h, 0, NULL, NULL);

		//if (status != CL_SUCCESS) {
		//	printf("clEnqueueReadBuffer failed\n");
		//	exit(-1);
		//}

		/* Compute right hand side for pressure equation */
		/*-----------------------------------------------*/
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
			printf("clSetKernelArg failed\n");
			exit(-1);
		}
		
		// Execute the kernel
		status = clEnqueueNDRangeKernel(cmdQueue, RHS_kernel, 2, NULL, globalWorkSize, 
			localWorkSize, 0, NULL, &event);
		if(status != CL_SUCCESS) {
			printf("clEnqueueNDRangeKernel failed: %s\n", cluErrorString(status));
			exit(-1);
		}

		clWaitForEvents(1, &event);

		// Read the OpenCL output buffer (F_d) to the host output array (F_h)
		status = clEnqueueReadBuffer(cmdQueue, TEMP_new_d, CL_TRUE, 0, 
			datasize, TEMP_h, 0, NULL, NULL);
		status = clEnqueueReadBuffer(cmdQueue, F_d, CL_TRUE, 0, 
			datasize, F_h, 0, NULL, NULL);
		status |= clEnqueueReadBuffer(cmdQueue, G_d, CL_TRUE, 0, 
			datasize, G_h, 0, NULL, NULL);
		status |= clEnqueueReadBuffer(cmdQueue, RHS_d, CL_TRUE, 0,
			datasize, RHS_h, 0, NULL, NULL);
		if (status != CL_SUCCESS) {
			printf("clEnqueueReadBuffer failed\n");
			exit(-1);
		}

		/* Solve the pressure equation by successive over relaxation */
		/*-----------------------------------------------------------*/
		//if (ifull > 0) {
		//	itersor = POISSON(P_h, RHS_h, FLAG_h, imax, jmax, delx, dely,
		//		eps, itermax, omg, &res, ifull, p_bound);
		//}

		//printf("t_end= %1.5g, t= %1.3e, delt= %1.1e, iterations %3d, res: %e, F-cells: %d, S-cells: %d, B-cells: %d\n",
		//	t_end, t+delt, delt, itersor, res, ifull, isurf, ibound);  

		/* Compute the new velocity field */
		/*--------------------------------*/
		//ADAP_UV(U_h, V_h, F_h, G_h, P_h, FLAG_h, imax, jmax, delt, delx, dely);

		/* Set boundary conditions */
		/*-------------------------*/
		//SETBCOND(U_h, V_h, P_h, TEMP_h, FLAG_h, imax, jmax, wW, wE, wN, wS);
		/* Set special boundary conditions */
		/* Overwrite preset default values */
		/*---------------------------------*/
		//SETSPECBCOND(problem, U_h, V_h, P_h, TEMP_h, imax, jmax, UI, VI);

		//if (!strcmp(problem, "drop") || !strcmp(problem, "dam") ||
		//	!strcmp(problem, "molding") || !strcmp(problem, "wave")) {
		//	SET_UVP_SURFACE(U_h, V_h, P_h, FLAG_h, GX, GY, imax, jmax, Re, delx, dely, delt);
		//}

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
	//}           
	
	/*
	 * End of main time loop
	 */

	stop_gpu = clock();

	timer_gpu = (double) (stop_gpu - start_gpu) / CLOCKS_PER_SEC;

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

	printf("  GPU time    : %.4f (ms)\n", timer_gpu);

	#ifdef PRINT

	print_1darray_to_file(TEMP_h, imax+2, jmax+2, "TEMP_h.txt");
	print_1darray_to_file(F_h, imax+2, jmax+2, "F_h.txt");
	print_1darray_to_file(G_h, imax+2, jmax+2, "G_h.txt");
	print_1darray_to_file(RHS_h, imax+2, jmax+2, "RHS_h.txt");

	#endif

	#endif

	#ifdef CPU

	/*
	 * CPU Execution
	 */

	// Create timers
	double timer_cpu;
	clock_t start_cpu, stop_cpu;

	start_cpu = clock();

	/*
	 * Main time loop
	 */

	t=0.0;
	//for (t=0.0, cycle=0; t < t_end; t+=delt, cycle++) {
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
		//	ifull = imax*jmax-ibound;
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

	//	/* Solve the pressure equation by successive over relaxation */
	//	/*-----------------------------------------------------------*/
	//	if (ifull > 0) {
	//		itersor = POISSON(P, RHS, FLAG, imax, jmax, delx, dely,
	//		eps, itermax, omg, &res, ifull, p_bound);
	//	}

	//	printf("t_end= %1.5g, t= %1.3e, delt= %1.1e, iterations %3d, res: %e, F-cells: %d, S-cells: %d, B-cells: %d\n",
	//		t_end, t+delt, delt, itersor, res, ifull, isurf, ibound);  

	//	/* Compute the new velocity field */
	//	/*--------------------------------*/
	//	ADAP_UV(U, V, F, G, P, FLAG, imax, jmax, delt, delx, dely);

	//	/* Set boundary conditions */
	//	/*-------------------------*/
	//	SETBCOND(U, V, P, TEMP, FLAG, imax, jmax, wW, wE, wN, wS);
	//	/* Set special boundary conditions */
	//	/* Overwrite preset default values */
	//	/*---------------------------------*/
	//	SETSPECBCOND(problem, U, V, P, TEMP, imax, jmax, UI, VI);

	//	if (!strcmp(problem, "drop") || !strcmp(problem, "dam") ||
	//		!strcmp(problem, "molding") || !strcmp(problem, "wave")) {
	//		SET_UVP_SURFACE(U, V, P, FLAG, GX, GY, imax, jmax, Re, delx, dely, delt);
	//	}

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
	//}           
	
	/*
	 * End of main time loop
	 */

	stop_cpu = clock();

	timer_cpu = (double) (stop_cpu - start_cpu) / CLOCKS_PER_SEC;

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

	printf("  CPU time    : %.4f (ms)\n",timer_cpu);

	#ifdef PRINT

	//print_array_to_file(U, imax+2, jmax+2, "U.txt");
	//print_array_to_file(V, imax+2, jmax+2, "V.txt");
	print_array_to_file(F, imax+2, jmax+2, "F.txt");
	print_array_to_file(G, imax+2, jmax+2, "G.txt");
	//print_array_to_file(P, imax+2, jmax+2, "P.txt");
	print_array_to_file(TEMP, imax+2, jmax+2, "TEMP.txt");
	//print_array_to_file(PSI, imax+1, jmax+1, "PSI.txt");
	//print_array_to_file(ZETA, imax, jmax, "ZETA.txt");
	//print_array_to_file(HEAT, imax+1, jmax+1, "HEAT.txt");
	print_array_to_file(RHS, imax+2, jmax+2, "RHS.txt");
	//print_array_int_to_file(FLAG, imax+2, jmax+2, "FLAG.txt");

	#endif

	#endif

	#ifdef VERIFY

	/*
	 * Verification
	 */

	if (compare_array(TEMP, TEMP_h, imax+2, jmax+2)) {
		printf("PASSED\n");
	} else {
		printf("FAILED\n");
	}

	if (compare_array(F, F_h, imax+2, jmax+2)) {
		printf("PASSED\n");
	} else {
		printf("FAILED\n");
	}

	if (compare_array(G, G_h, imax+2, jmax+2)) {
		printf("PASSED\n");
	} else {
		printf("FAILED\n");
	}

	if (compare_array(RHS, RHS_h, imax+2, jmax+2)) {
		printf("PASSED\n");
	} else {
		printf("FAILED\n");
	}

	// Print timings
	printf("  Speedup %.2fx\n\n", timer_cpu/timer_gpu);
	
	#endif

	//-----------------------------------------------------
	// STEP 13: Release OpenCL resources
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

	FREE_RMATRIX(P_h,		0, imax+1,	0, jmax+1);
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

	printf("End of program\n");
	return(0);
}
