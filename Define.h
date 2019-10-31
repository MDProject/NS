#ifndef DEFINE_DEFINE
#define DEFINE_DEFINE
#include <iostream>
#include <stdio.h>

typedef double Scalar;

extern Scalar R;
extern Scalar L;
extern Scalar H;
extern int Nx;
extern int Nz;

void assemble_ls_array(Scalar** ls1, Scalar** ls2, int Nx);


#endif