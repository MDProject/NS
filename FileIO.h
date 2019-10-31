#ifndef DEFINE_FILEIO
#define DEFINE_FILEIO
#include "Define.h"

void WriteCompactArray(char* fpath, int* I, int* J, Scalar* K, int Nx, int Nz);

// parameters order:	L,H,Nx,Nz
void WriteParameter(char* fpath);

void WriteSlipLength(char* fpath, Scalar* ls1, Scalar* ls2);

#endif
