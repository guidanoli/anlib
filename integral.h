/* Definição do módulo integral
   Programado por Guilherme Dantas
   10/10/2018 */

#include <stdio.h>

/* Recebe o ponteiro para função a ser derivada, o ponto
   em que a derivada será calculada, e a janela h.
   Retorna a derivada da função pelo método de segunda ordem.
   Se o ponteiro pra função for nulo ou h for zero, retorna 0.0 */
double derivada (double (*f) (double x), double x, double h);

/* Recebe o ponteiro para função, o ponteiro para a função que é
   a derivada analítica da primeira função, e o ponto em que a derivada
   será calculada.
   Retorna o h para o qual o erro da derivada numérica é menor.
   Se o ponteiro para alguma das funções for nulo, retorna 0.0 */
double h_otimo (double (*f) (double x), double (*fl) (double x), double x);

/* Recebe o ponteiro para função, os limites inferior e superior
   e a quantidade de passos na integração.
   Retorna a integração da função f no intervalo [a,b] com n passos
   calculada pelo método de Simpson.
   Se o ponteiro for nulo, retorna 0.0 */
double simpson (double (*f) (double), double a, double b, int n);

/* Recebe o ponteiro para função, os limites inferior e superior
   e a quantidade de passos na integração.
   Retorna a integração da função f no intervalo [a,b] com n passos
   calculada pelo método do Ponto Médio.
   Se o ponteiro for nulo, retorna 0.0 */
double pontomedio (double (*f) (double), double a, double b, int n);