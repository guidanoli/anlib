#include <stdio.h>

/* Retorna o número de passos para chegar a uma boa aproximação da raiz da função f, e armazena em r a raiz aproximada */
int bissecao (double a, double b, int p, double (*f) (double x), double *r);

int IQI (double x0, double x1, double x2, int p, double (*f) (double x), double* r);