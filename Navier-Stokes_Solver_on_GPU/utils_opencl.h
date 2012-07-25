#pragma once
#ifndef UTILS_OPENCL_H_
#define UTILS_OPENCL_H_

const char *cluErrorString(cl_int err);

void printPlatforms(cl_platform_id* platforms, cl_uint numPlatforms, cl_int status);

void printDevices(cl_device_id * devices, cl_uint numDevices, cl_int status);

void getDeviceInfo(cl_device_id* ID, cl_uint deviceCount);

unsigned int nextPow2(unsigned int x);

size_t getGlobalSize(int group_size, int global_size) ;

#endif /* UTILS_OPENCL_H_ */