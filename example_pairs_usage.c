/*
 * Ejemplo de uso de la función read_pairs_file()
 * 
 * Compilar junto con pairs_reader.c:
 * gcc -std=c11 -Wall -I./include src/pairs_reader.c example_pairs_usage.c -o example_pairs
 * 
 * Uso:
 * ./example_pairs 99_pairs.txt
 */

#include <stdio.h>
#include <stdlib.h>
#include "pairs_reader.h"

int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("Uso: %s <archivo_de_pares>\n", argv[0]);
        printf("Ejemplo: %s 99_pairs.txt\n", argv[0]);
        return 1;
    }

    size_t rows, cols;
    
    // Leer el archivo de pares
    int **pairs = read_pairs_file(argv[1], &rows, &cols);
    
    if (!pairs) {
        printf("ERROR: no se pudo leer el archivo\n");
        return 1;
    }

    // Imprimir el contenido
    print_pairs((const int**)pairs, rows, cols);
    
    // Ejemplo: acceder a elementos específicos
    printf("\n=== Acceso a elementos específicos ===\n");
    if (rows > 0 && cols > 0) {
        printf("Primer par: [%d, %d]\n", pairs[0][0], pairs[0][1]);
    }
    if (rows > 0 && cols > 1) {
        printf("Primer par (completo): ");
        for (size_t j = 0; j < cols; j++) {
            printf("%d ", pairs[0][j]);
        }
        printf("\n");
    }
    
    // Limpiar memoria
    free_pairs(pairs, rows);
    
    printf("\nMemoria liberada exitosamente.\n");
    return 0;
}
