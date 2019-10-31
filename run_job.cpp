#include "DifferenceOperator.h"
#include "FileIO.h"
#include "MatrixAssembly.h"
#include "Define.h"

int main() {
	Scalar DX = 2.*L / (Nx - 1);
	Scalar DZ = 2.*H / (Nz - 1);

	Scalar* ls1, *ls2;
	assemble_ls_array(&ls1, &ls2, Nx);
	Scalar** C, ** L;
	biharmonic_4th_center(&C, DX, DZ);
	laplacian_4th_center(&L, DX, DZ);

	int* I, *J;
	Scalar* K;

	assemble_biharmonic_matrix_2nd(C, ls1, ls2, &I, &J, &K, Nx, Nz);
	char fpath_biharmonic[] = "D:\\WinSCP\\PDE\\Matrix\\biharmonic.txt";
	WriteCompactArray(fpath_biharmonic, I, J, K, Nx, Nz);
	free_biharmonic_matrix(C, I, J, K);
	
	assemble_laplacian_matrix_2nd(L, ls1, ls2, &I, &J, &K, Nx, Nz);
	char fpath_laplacian[] = "D:\\WinSCP\\PDE\\Matrix\\laplacian.txt";
	WriteCompactArray(fpath_laplacian, I, J, K, Nx, Nz);
	free_laplacian_matrix(L, I, J, K);

	char fpath_setting[] = "D:\\WinSCP\\PDE\\Matrix\\setting.txt";
	WriteParameter(fpath_setting);

	char fpath_sliplen[] = "D:\\WinSCP\\PDE\\Matrix\\ls.txt";
	WriteSlipLength(fpath_sliplen, ls1, ls2);
	system("pause");
}