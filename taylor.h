/* Interface do módulo taylor */

#define M_PI 3.14159265358979323846

/* Aproxima cos(x) para o domínio recomendado de 0 a PI centrado em 0 por uma série de Taylor */
double tcos (double x);

/* Calcula o resíduo da aproximação de cos(x) para o domínio recomendado de 0 a PI centrado em 0 por uma série de Taylor */
double tcos_erro (double x);

/* Aproxima sqrt(x) para o domínio recomendado de 1 a 2 centrado em 1 por uma série de Taylor */
double tsqrt (double x);

/* Calcula o resíduo da aproximação sqrt(x) para o domínio recomendado de 1 a 2 centrado em 1 por uma série de Taylor */
double tsqrt_erro (double x);

/* Averigua o menor valor, retornando 1 se x >= 0 e 0 caso contrário */
int menor (double x, double y);

/* Retorna o valor absoluto de x */
double fabs (double x);