#include <stdio.h>
#include <stdlib.h>

#include <CL\cl.h>

const char *cluErrorString(cl_int err) {

	switch(err) {

	case CL_SUCCESS: return "CL_SUCCESS";
	case CL_DEVICE_NOT_FOUND: return "CL_DEVICE_NOT_FOUND";
	case CL_DEVICE_NOT_AVAILABLE: return "CL_DEVICE_NOT_AVAILABLE";
	case CL_COMPILER_NOT_AVAILABLE: return
										"CL_COMPILER_NOT_AVAILABLE";
	case CL_MEM_OBJECT_ALLOCATION_FAILURE: return
											   "CL_MEM_OBJECT_ALLOCATION_FAILURE";
	case CL_OUT_OF_RESOURCES: return "CL_OUT_OF_RESOURCES";
	case CL_OUT_OF_HOST_MEMORY: return "CL_OUT_OF_HOST_MEMORY";
	case CL_PROFILING_INFO_NOT_AVAILABLE: return
											  "CL_PROFILING_INFO_NOT_AVAILABLE";
	case CL_MEM_COPY_OVERLAP: return "CL_MEM_COPY_OVERLAP";
	case CL_IMAGE_FORMAT_MISMATCH: return "CL_IMAGE_FORMAT_MISMATCH";
	case CL_IMAGE_FORMAT_NOT_SUPPORTED: return
											"CL_IMAGE_FORMAT_NOT_SUPPORTED";
	case CL_BUILD_PROGRAM_FAILURE: return "CL_BUILD_PROGRAM_FAILURE";
	case CL_MAP_FAILURE: return "CL_MAP_FAILURE";
	case CL_INVALID_VALUE: return "CL_INVALID_VALUE";
	case CL_INVALID_DEVICE_TYPE: return "CL_INVALID_DEVICE_TYPE";
	case CL_INVALID_PLATFORM: return "CL_INVALID_PLATFORM";
	case CL_INVALID_DEVICE: return "CL_INVALID_DEVICE";
	case CL_INVALID_CONTEXT: return "CL_INVALID_CONTEXT";
	case CL_INVALID_QUEUE_PROPERTIES: return "CL_INVALID_QUEUE_PROPERTIES";
	case CL_INVALID_COMMAND_QUEUE: return "CL_INVALID_COMMAND_QUEUE";
	case CL_INVALID_HOST_PTR: return "CL_INVALID_HOST_PTR";
	case CL_INVALID_MEM_OBJECT: return "CL_INVALID_MEM_OBJECT";
	case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR: return
												 "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
	case CL_INVALID_IMAGE_SIZE: return "CL_INVALID_IMAGE_SIZE";
	case CL_INVALID_SAMPLER: return "CL_INVALID_SAMPLER";
	case CL_INVALID_BINARY: return "CL_INVALID_BINARY";
	case CL_INVALID_BUILD_OPTIONS: return "CL_INVALID_BUILD_OPTIONS";
	case CL_INVALID_PROGRAM: return "CL_INVALID_PROGRAM";
	case CL_INVALID_PROGRAM_EXECUTABLE: return
											"CL_INVALID_PROGRAM_EXECUTABLE";
	case CL_INVALID_KERNEL_NAME: return "CL_INVALID_KERNEL_NAME";
	case CL_INVALID_KERNEL_DEFINITION: return "CL_INVALID_KERNEL_DEFINITION";
	case CL_INVALID_KERNEL: return "CL_INVALID_KERNEL";
	case CL_INVALID_ARG_INDEX: return "CL_INVALID_ARG_INDEX";
	case CL_INVALID_ARG_VALUE: return "CL_INVALID_ARG_VALUE";
	case CL_INVALID_ARG_SIZE: return "CL_INVALID_ARG_SIZE";
	case CL_INVALID_KERNEL_ARGS: return "CL_INVALID_KERNEL_ARGS";
	case CL_INVALID_WORK_DIMENSION: return "CL_INVALID_WORK_DIMENSION";
	case CL_INVALID_WORK_GROUP_SIZE: return "CL_INVALID_WORK_GROUP_SIZE";
	case CL_INVALID_WORK_ITEM_SIZE: return "CL_INVALID_WORK_ITEM_SIZE";
	case CL_INVALID_GLOBAL_OFFSET: return "CL_INVALID_GLOBAL_OFFSET";
	case CL_INVALID_EVENT_WAIT_LIST: return "CL_INVALID_EVENT_WAIT_LIST";
	case CL_INVALID_EVENT: return "CL_INVALID_EVENT";
	case CL_INVALID_OPERATION: return "CL_INVALID_OPERATION";
	case CL_INVALID_GL_OBJECT: return "CL_INVALID_GL_OBJECT";
	case CL_INVALID_BUFFER_SIZE: return "CL_INVALID_BUFFER_SIZE";
	case CL_INVALID_MIP_LEVEL: return "CL_INVALID_MIP_LEVEL";

	default: return "UNKNOWN CL ERROR CODE";
	}
}

void printPlatforms(cl_platform_id* platforms, cl_uint numPlatforms, cl_int status) 
{
	printf("%u platforms detected\n", numPlatforms);
	for (unsigned int i = 0; i < numPlatforms; i++) {
		char buf[100];
		printf("Platform %u: \n", i);
		status = clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR,
			sizeof(buf), buf, NULL);
		printf("\tVendor: %s\n", buf);
		status |= clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME,
			sizeof(buf), buf, NULL);
		printf("\tName: %s\n", buf);

		if (status != CL_SUCCESS) {
			printf("clGetPlatformInfo failed\n");
			exit(-1);
		}
	}
	printf("\n");
}

void printDevices(cl_device_id * devices, cl_uint numDevices, cl_int status) 
{
	printf("%u devices detected\n", numDevices);
	for (unsigned int i = 0; i < numDevices; i++) {
		char buf[100];
		printf("Device %u: \n", i);
		status = clGetDeviceInfo(devices[i], CL_DEVICE_VENDOR,
			sizeof(buf), buf, NULL);
		printf("\tDevice: %s\n", buf);
		status |= clGetDeviceInfo(devices[i], CL_DEVICE_NAME,
			sizeof(buf), buf, NULL);
		printf("\tName: %s\n", buf);

		if (status != CL_SUCCESS) {
			printf("clGetDeviceInfo failed\n");
			exit(-1);
		}
	}
	printf("\n");
}

void getDeviceInfo(cl_device_id* ID, cl_uint deviceCount)
{
	cl_int err; size_t size;

	for (int d = 0; d < deviceCount; d++){
		// Get device name and print to screen
		err = clGetDeviceInfo(ID[d], CL_DEVICE_VENDOR, 0, NULL, &size);
		char* vendor = (char*) malloc(sizeof(char) * size);
		err = clGetDeviceInfo(ID[d], CL_DEVICE_VENDOR, size, vendor, NULL);

		err = clGetDeviceInfo(ID[d], CL_DEVICE_NAME, 0, NULL, &size);
		char* name = (char*) malloc(sizeof(char) * size);
		err = clGetDeviceInfo(ID[d], CL_DEVICE_NAME, size, name, NULL);

		printf("Compute capabilities of %s, %s\n\n" ,name, vendor);

		err = clGetDeviceInfo(ID[d], CL_DEVICE_VERSION, 0, NULL, &size);
		char* version = (char*) malloc(sizeof(char) * size);
		err = clGetDeviceInfo(ID[d], CL_DEVICE_VERSION, size, version, NULL);
		printf("\t OpenCL Device Support \t \t %s\n", version);

		err = clGetDeviceInfo(ID[d], CL_DRIVER_VERSION, 0, NULL, &size);
		char* driver = (char*) malloc(sizeof(char) * size);
		err = clGetDeviceInfo(ID[d], CL_DRIVER_VERSION, size, driver, NULL);
		printf("\t OpenCL Software Driver \t %s\n", driver);

		err = clGetDeviceInfo(ID[d], CL_DEVICE_PROFILE, 0, NULL, &size);
		char* profile = (char*) malloc(sizeof(char) * size);
		err = clGetDeviceInfo(ID[d], CL_DEVICE_PROFILE, size, profile, NULL);
		printf("\t OpenCL Profile Supported \t %s\n\n", profile);

		cl_uint maxComputeUnits;
		err = clGetDeviceInfo(ID[d], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint),
			&maxComputeUnits, NULL);
		printf("\t Compute Units Available \t %d\n", maxComputeUnits);

		size_t maxWorkGroupItems;
		err = clGetDeviceInfo(ID[d], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t),
			&maxWorkGroupItems, NULL);
		printf("\t Max Work-items per Work-group \t %d\n", maxWorkGroupItems);

		cl_ulong localMem;
		err = clGetDeviceInfo(ID[d], CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), 
			&localMem, NULL);
		printf("\t Compute Unit Local Memory \t %dKb\n\n", localMem/1024);

		cl_ulong globalMem;
		err = clGetDeviceInfo(ID[d], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), 
			&globalMem, NULL);
		printf("\t Global Memory Size \t \t %dMB\n", globalMem/1048576);

		cl_bool ECC;
		err = clGetDeviceInfo(ID[d], CL_DEVICE_ERROR_CORRECTION_SUPPORT, 
			sizeof(cl_bool), &ECC, NULL);
		switch (ECC) {
		case CL_TRUE:
			printf("\t Error Correction Support \t Yes\n");
			break;
		case CL_FALSE:
			printf("\t Error Correction Support \t No\n");
			break;
		}

		cl_ulong globalCache;
		err = clGetDeviceInfo(ID[d], CL_DEVICE_GLOBAL_MEM_CACHE_SIZE,
			sizeof(cl_ulong), &globalCache, NULL);
		printf("\t Global Memory Cache \t \t %dKB\n", globalCache/1025);

		cl_uint globalCacheLine;
		err = clGetDeviceInfo(ID[d], CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE, 
			sizeof(cl_uint), &globalCacheLine, NULL);
		printf("\t Global Memory Cache line \t %dBytes\n\n", globalCacheLine);
	}
}

//int getUserDefPlatform(cl_platform_id** ID, char** Name, int* P)
//{
//	// Get number of platforms
//	cl_uint platformCount = 0;
//	cl_int err = clGetPlatformIDs(0, NULL, &platformCount);
//	if (err != CL_SUCCESS) {
//		printf("Failed to get number of platforms available\n");
//		return 1;
//	}
//
//	// Allocate and get platforms
//	*ID = (cl_platform_id*) malloc(platformCount * sizeof(cl_platform_id));
//	err = clGetPlatformIDs(platformCount, *ID, NULL);
//	if (err != CL_SUCCESS) {
//		printf("Failed to get platform ids available\n");
//		return 1;
//	}
//
//	// Iterate through platforms and display the name of each
//	printf("Platforms available:\n");
//	for (int p = 0; p < platformCount; p++) {
//		// Get info length
//		size_t required;
//		err = clGetPlatformInfo(* ID[p], CL_PLATFORM_NAME, 0, NULL, &required);
//		if (err != CL_SUCCESS) {
//			printf("Failed to get length of info string\n");
//			return 1;
//		}
//
//		// Get info
//		*Name = (char*) malloc(required * sizeof(char));
//		err = clGetPlatformInfo(* ID[p], CL_PLATFORM_NAME, required, *Name, NULL);
//		if (err != CL_SUCCESS) {
//			printf("Failed to get info string\n");
//			return 1;
//		}
//
//		// Print name
//		printf("[%d] \t %s \n" , p+1, *Name);
//
//		// Free name
//		free(*Name);
//	}
//
//	// Ask user for which platform to use
//	printf("Choose a platform:");
//	scanf("%d" ,P);
//	*P = *P - 1;
//
//	// Get info length
//	size_t required;
//	err = clGetPlatformInfo(*ID[*P], CL_PLATFORM_NAME, 0, NULL, &required);
//	if (err != CL_SUCCESS) {
//		printf("Failed to get length of info string\n");
//		return 1;
//	}
//
//	// Get info
//	*Name = (char*) malloc(required * sizeof(char));
//	err = clGetPlatformInfo(*ID[*P], CL_PLATFORM_NAME, required, *Name, NULL);
//	if (err != CL_SUCCESS) {
//		printf("Failed to get info string\n");
//		return 1;
//	}
//
//	return 0;
//}
//
//int getUserDefDevice(cl_platform_id* platformIds, int p, cl_device_id** ID ,
//	char** Name, int* D, cl_uint* deviceCount)
//{
//	// Get number of devices
//	cl_int err;
//	// cl_uint deviceCount = *Count;
//	err = clGetDeviceIDs(platformIds[p], CL_DEVICE_TYPE_ALL, 0, NULL,
//		deviceCount);
//	if (err != CL_SUCCESS) {
//		printf("Failed to get number of devices in platform\n");
//		return 1;
//	}
//
//	// Get device ids
//	*ID = (cl_device_id*) malloc(*deviceCount * sizeof(cl_device_id));
//	err = clGetDeviceIDs (platformIds[p], CL_DEVICE_TYPE_ALL, *deviceCount, *ID,
//		NULL);
//	if (err != CL_SUCCESS) {
//		printf("Failed to get device ids in platform\n");
//		return 1;
//	}
//	cl_device_id * id = *ID;
//
//	// Iterate through devices and display the name of each
//	printf(" \ nDevices available on platform %d:\n", p+1);
//	for (int d = 0; d < *deviceCount; d++) {
//		// Get info length
//		size_t required;
//		err = clGetDeviceInfo(id[d], CL_DEVICE_NAME, 0, NULL, &required);
//		if (err != CL_SUCCESS) {
//			printf("Failed to get l eng th of info string\n");
//			return 1;
//		}
//
//		// Get info
//		*Name = (char*)malloc(r e qui r ed *sizeof(char));
//		err = clGetDeviceInfo(id [ d ], CL_DEVICE_NAME, r equi r ed, *Name, NULL);
//		if (err != CL_SUCCESS) {
//			printf("Failed to get length of info string\n");
//			return 1;
//		}
//
//		// Print name
//		printf("[%d] \t %s \n" , d+1, *Name);
//
//		// Free name
//		free(*Name);
//	}
//
//	// Ask the user to choose a device
//	printf(" Choose a device : ");
//	scanf("%d" ,D);
//	*D = *D - 1;
//
//	// Get info length
//	size_t required;
//	err = clGetDeviceInfo(id[*D], CL_DEVICE_NAME, 0, NULL, &required);
//	if (err != CL_SUCCESS) {
//		printf("Failed to get length of info string\n");
//		return 1;
//	}
//
//	// Get info
//	*Name = (char*) malloc(required * sizeof(char));
//	err = clGetDeviceInfo(id[*D], CL_DEVICE_NAME, required, *Name, NULL);
//	if (err != CL_SUCCESS) {
//		printf("Failed to get length of info string\n");
//		return 1;
//	}
//
//	return 0;
//}
//


//
//int buildProgramFromSources(cl_context context, cl_uint deviceCount,
//	cl_device_id *deviceIds, const char* filname, cl_program* program) {
//		// Open file
//		FILE* file = fopen(filename, "r");
//
//		// Find size and go back to start
//		fseek(file, 0, SEEK_END);
//		size_t length = ftell(file);
//		fseek(file, 0, SEEK_SET);
//
//		// Allocate and read content
//		char* source = (char*)malloc(length * sizeof(char));
//		fread(source, length, 1, file);
//
//		// Close file
//		fclose(file);
//
//		// Create program
//		cl_int err;
//		*program = clCreateProgramWithSource(context, 1, &source, &length, &err);
//		if (err != CL_SUCCESS) {
//			printf("Failed to create program in context.\n");
//			return 1;
//		}
//
//		// Build sources
//		err = clBuildProgram(*program, deviceCount, deviceIds, NULL, NULL, NULL);
//		if (err != CL_SUCCESS) {
//			printf(" Failed to build program for devices in platform.\n");
//			return 1;
//		}
//		// Free content
//		free(source);
//		return 0;
//}

size_t getGlobalWorkSize()
{



}