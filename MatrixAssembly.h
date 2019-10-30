#ifndef DEFINE_BHASSEMBLY
#define DEFINE_BHASSEMBLY
#include "DifferenceOperator.h"

/*
	* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	Nz-1
	^
	|
	2
	* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	* * * * * * * * * * * * * * (1 --> Nx-1) * * * * * * * * * * * * * 
	i and j increase along the direction [-L --> L](1~Nx-1) and [-H --> H](2~Nz-1)

	shrink above inner grid points from 2D to 1D sequence according to the j-principal order [(j-2)(Nx-1)+i] 
*/

void assemble_biharmonic_matrix(Scalar** C, Scalar* ls1, Scalar* ls2, int** I, int** J, Scalar** K, int Nx, int Nz, Scalar DX, Scalar DZ);

void free_biharmonic_matrix(int** I, int** J, Scalar** K, int Nx, int Nz);

#endif
