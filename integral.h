/* Defini��o do m�dulo integral
   Programado por Guilherme Dantas
   10/10/2018 */

#include <stdio.h>

/* Recebe o ponteiro para fun��o a ser derivada, o ponto
   em que a derivada ser� calculada, e a janela h.
   Retorna a derivada da fun��o pelo m�todo de segunda ordem.
   Se o ponteiro pra fun��o for nulo ou h for zero, retorna 0.0 */
double derivada (double (*f) (double x), double x, double h);

/* Recebe o ponteiro para fun��o, o ponteiro para a fun��o que �
   a derivada anal�tica da primeira fun��o, e o ponto em que a derivada
   ser� calculada.
   Retorna o h para o qual o erro da derivada num�rica � menor.
   Se o ponteiro para alguma das fun��es for nulo, retorna 0.0 */
double h_otimo (double (*f) (double x), double (*fl) (double x), double x);

/* Recebe o ponteiro para fun��o, os limites inferior e superior
   e a quantidade de passos na integra��o.
   Retorna a integra��o da fun��o f no intervalo [a,b] com n passos
   calculada pelo m�todo de Simpson.
   Se o ponteiro for nulo, retorna 0.0 */
double simpson (double (*f) (double), double a, double b, int n);

/* Recebe o ponteiro para fun��o, os limites inferior e superior
   e a quantidade de passos na integra��o.
   Retorna a integra��o da fun��o f no intervalo [a,b] com n passos
   calculada pelo m�todo do Ponto M�dio.
   Se o ponteiro for nulo, retorna 0.0 */
double pontomedio (double (*f) (double), double a, double b, int n);