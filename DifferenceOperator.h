#ifndef DEFINE_DIFFOPERATOR
#define DEFINE_DIFFOPERATOR
#include "Define.h"

/* C[i][j]/L[i][j]
		*	*	*	*	*	(left --> right: j increase)
		*	*	*	*	*	(up --> down:	i increase)
		*	*	*	*	*
		*	*	*	*	*
		*	*	*	*	*
*/
void biharmonic_4th_center(Scalar*** C, Scalar DZ, Scalar DX);
void free_biharmonic_4th_center(Scalar** C);


void laplacian_4th_center(Scalar*** L, Scalar DZ, Scalar DX);
void free_laplacian_4th_center(Scalar** L);

#endif