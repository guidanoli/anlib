/* Interface do m�dulo taylor */

#define M_PI 3.14159265358979323846

/* Aproxima cos(x) para o dom�nio recomendado de 0 a PI centrado em 0 por uma s�rie de Taylor */
double tcos (double x);

/* Calcula o res�duo da aproxima��o de cos(x) para o dom�nio recomendado de 0 a PI centrado em 0 por uma s�rie de Taylor */
double tcos_erro (double x);

/* Aproxima sqrt(x) para o dom�nio recomendado de 1 a 2 centrado em 1 por uma s�rie de Taylor */
double tsqrt (double x);

/* Calcula o res�duo da aproxima��o sqrt(x) para o dom�nio recomendado de 1 a 2 centrado em 1 por uma s�rie de Taylor */
double tsqrt_erro (double x);

/* Averigua o menor valor, retornando 1 se x >= 0 e 0 caso contr�rio */
int menor (double x, double y);

/* Retorna o valor absoluto de x */
double fabs (double x);