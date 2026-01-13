/* mymath.h */
#ifndef MYMATH_H
#define MYMATH_H

#include <complex.h>
#include <stddef.h>

/* 
 * DFT - Discrete Fourier Transform (implementación directa)
 * Parámetros:
 *   x - puntero a la señal de entrada compleja
 *   X - puntero al resultado de la transformada compleja
 *   N - tamaño de la señal (número de muestras)
 */
void dft(const double complex *x, double complex *X, size_t N);

/* 
 * vec_mat - Multiplicación vector-matriz
 * y = x * A
 * donde x es vector fila de tamaño m
 *       A es matriz de m x n
 *       y es vector fila resultante de tamaño n
 * Parámetros:
 *   x - puntero al vector de entrada (tamaño m)
 *   A - puntero a la matriz (tamaño m x n, almacenada por filas)
 *   y - puntero al vector resultante (tamaño n)
 *   m - número de filas de A / tamaño de x
 *   n - número de columnas de A / tamaño de y
 */
void vec_mat(const double *x, const double *A, double *y, size_t m, size_t n);

#endif /* MYMATH_H */