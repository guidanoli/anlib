/*	Interface de matriz.h 
	Programado por Guilherme Dantas
	15/8/18
*/

/* Cria um vetor de dimensão n */
double* vetcria (int n);

/* Libera o espaço na memória ocupado por um vetor v */
void vetlibera (double* v);

/* Retorna o produto escalar entre dois vetores v e w de dimensão n */
double escalar (int n, double* v, double* w);

/* Retorna a norma-2 de um vetor v de dimensão n */
double norma2 (int n, double* v);

/* Retorna uma matriz de dimensões m por n */
double** matcria (int m, int n);

/* Libera o espaço na memória ocupado por uma matriz A com m linhas */
void matlibera (int m, double** A);

/* Transpõe uma matriz A m por n e armazena na matriz T n por m */
void transposta (int m, int n, double** A, double** T);

/* Multiplica uma matriz A m por n por um vetor v de dimensão n e armazena o produto no vetor w */
void multmv (int m, int n, double** A, double* v, double* w);

/* Multiplica uma matriz A (m por n) por uma matriz B (n por q) e armazena o produto na matriz C (m por q) */
void multmm (int m, int n, int q, double** A, double** B, double** C);