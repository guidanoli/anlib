#include <stdio.h>

/* Retorna o n�mero de passos para chegar a uma boa aproxima��o da raiz da fun��o f, e armazena em r a raiz aproximada */
int bissecao (double a, double b, int p, double (*f) (double x), double *r);

int IQI (double x0, double x1, double x2, int p, double (*f) (double x), double* r);