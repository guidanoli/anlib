#include "sistlinear.h"


double absoluto(double x)
{
	if(x >= 0.0)
		return x;
	else
		return x*(-1);
}

void fatoracao (int n, double** A, int* p)
{
	int i, j, k, piv, tempI;
	double f, tempD;

	/* Inicializa o vetor de permutação */
	for( i = 0 ; i < n ; i++ )
	{
		p[i] = i;
	}

	/* Eliminaçao de Gauss */
	for( i = 0 ; i < n-1 ; i++ )
	{

		/* Pivotamento */
		piv = i;
		
		for( k = i + 1 ; k < n ; k++ )
		{
			if(absoluto(A[k][i]) > absoluto(A[piv][i]))
				piv = k;
		}

		if(piv != i)
			for( k = i ; k < n ; k++ )
			{
				tempD = A[i][k];
				A[i][k] = A[piv][k];
				A[piv][k] = tempD;
			}

		tempI = p[i];
		p[i] = p[piv];
		p[piv] = tempI;

		for( k = i + 1 ; k < n ; k++ )
		{
			f = A[k][i]/A[i][i];

			/* Construção da meia-matriz L */
			A[k][i] = f;

			for( j = i+1 ; j < n ; j++ )
			{

				/* Construção da meia-matriz U */
				A[k][j] = A[k][j] - f*(A[i][j]);
			}
		}
	}


}