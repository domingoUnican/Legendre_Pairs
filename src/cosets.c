/* cosets_v2.c - Adaptación usando cyclotomic_cosets.c */
#include "cyclotomic_cosets.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

/*
CosetList* all_cosets(size_t N, size_t k) {
    // Simplemente usar cyclotomic_cosets directamente
    CosetList *result = (CosetList*)malloc(sizeof(CosetList));
    if (!result) {
        fprintf(stderr, "ERROR: sin memoria para resultado\n");
        return NULL;
    }
    
    *result = cyclotomic_cosets(k, N);
    
    if (result->len == 0) {
        fprintf(stderr, "ERROR: no se pudieron calcular cosets\n");
        free(result);
        return NULL;
    }
    
    return result;
}
*/

int main() {
    size_t N = 7;
    size_t k = 2;
    
    
    // Calcular todos los cosets
    //CosetList cosets = cyclotomic_cosets(k, N);
    
    // Imprimir los resultados
    //cosetlist_print(&cosets);
    
    // Liberar memoria
    //free_cosetlist(&cosets);
    printf("2");
    getchar();
    return 0;
}