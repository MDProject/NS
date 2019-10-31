#include "Define.h"

Scalar Ls1(Scalar x) {
	return 1;
}

Scalar Ls2(Scalar x) {
	return 1;
}

void assemble_ls_array(Scalar** ls1, Scalar** ls2, int Nx) {
	Scalar DX = 2.*L / (Nx - 1);
	(*ls1) = (Scalar*)calloc(Nx, sizeof(Scalar));
	(*ls2) = (Scalar*)calloc(Nx, sizeof(Scalar));
	for (int i = 1; i <= Nx - 1; i++) {
		(*ls1)[i] = Ls1(-L + (i - 1)*DX);
		(*ls2)[i] = Ls2(-L + (i - 1)*DX);
	}
}


Scalar R = 0.4;
Scalar L = 10;
Scalar H = 5;
int Nx = 201;
int Nz = 101;