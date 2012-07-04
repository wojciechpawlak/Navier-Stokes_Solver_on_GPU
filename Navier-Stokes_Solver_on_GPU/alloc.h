#pragma once
#ifndef ALLOC_H_
#define ALLOC_H_

#include "datadef.h"

REAL **RMATRIX(int nrl,int nrh,int ncl,int nch);

void FREE_RMATRIX(REAL** m,int nrl,int nrh,int ncl,int nch);

int **IMATRIX(int nrl,int nrh,int ncl,int nch);

void FREE_IMATRIX(int** m,int nrl,int nrh,int ncl,int nch);

#endif /* ALLOC_H_ */