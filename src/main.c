
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cyclotomic_cosets.h"
#include "mymath.h"
#include <complex.h>
#include <math.h>


static void print_bits(const uint8_t *v, size_t n, const char *name) {
    printf("%s = [", name);
    for (size_t i = 0; i < n; ++i) {
        printf("%u", (unsigned)v[i]);
        if (i + 1 < n) printf(", ");
    }
    printf("]\n");
}

static void print_vector(const uint8_t *v, size_t n) {
    printf("[");
    for (size_t i = 0; i < n; i++) {
        printf("%u", (unsigned)v[i]);
        if (i + 1 < n) printf(", ");
    }
    printf("]");
}


/* Versión modificada de BinaryCombinations que genera vectores completos */
typedef struct {
    uint8_t **vectors;        // Array de vectores de longitud N
    uint8_t **combinations;   // Array de combinaciones de cosets (para referencia)
    size_t num_vectors;       // Número total de vectores (2^num_cosets)
    size_t vector_length;     // Longitud de cada vector (N)
    size_t num_cosets;        // Número de cosets
} CosetVectors;

/* Función auxiliar: busca en qué coset está un elemento */
static int find_element_in_cosets(const CosetList *cl, size_t element) {
    for (size_t i = 0; i < cl->len; i++) {
        const Coset *c = &cl->data[i];
        for (size_t j = 0; j < c->len; j++) {
            if (c->data[j] == element) {
                return i;  // Devuelve el índice del coset
            }
        }
    }
    return -1;  // No encontrado (no debería pasar)
}


/* Genera vector completo para una combinación dada de cosets */
uint8_t* generate_vector_for_combination(const CosetList *cl, 
                                        const uint8_t *combination, 
                                        size_t N) {
    uint8_t *vector = (uint8_t*)calloc(N, sizeof(uint8_t));
    if (!vector) {
        fprintf(stderr, "ERROR: sin memoria para vector\n");
        return NULL;
    }
    
    // Para cada elemento 0..N-1, asignamos el valor según su coset
    for (size_t elem = 0; elem < N; elem++) {
        int coset_idx = find_element_in_cosets(cl, elem);
        if (coset_idx >= 0) {
            vector[elem] = combination[coset_idx];
        }
        // Si no se encuentra (no debería pasar), queda en 0
    }
    
    return vector;
}



/* Genera todos los vectores basados en combinaciones de cosets */
CosetVectors* generate_coset_vectors(const CosetList *cl, size_t N) {
    if (!cl || cl->len == 0 || N == 0) {
        fprintf(stderr, "ERROR: lista de cosets inválida\n");
        return NULL;
    }
    
    if (cl->len > 30) {  // Límite razonable
        fprintf(stderr, "ERROR: demasiados cosets (%zu)\n", cl->len);
        return NULL;
    }
    
    CosetVectors *result = (CosetVectors*)malloc(sizeof(CosetVectors));
    if (!result) {
        fprintf(stderr, "ERROR: sin memoria para CosetVectors\n");
        return NULL;
    }
    
    result->num_cosets = cl->len;
    result->num_vectors = 1ULL << cl->len;  // 2^num_cosets
    result->vector_length = N;
    
    // Reservar memoria para combinaciones
    result->combinations = (uint8_t**)malloc(result->num_vectors * sizeof(uint8_t*));
    result->vectors = (uint8_t**)malloc(result->num_vectors * sizeof(uint8_t*));
    
    if (!result->combinations || !result->vectors) {
        fprintf(stderr, "ERROR: sin memoria para arrays\n");
        if (result->combinations) free(result->combinations);
        if (result->vectors) free(result->vectors);
        free(result);
        return NULL;
    }
    
    // Inicializar todo a NULL
    memset(result->combinations, 0, result->num_vectors * sizeof(uint8_t*));
    memset(result->vectors, 0, result->num_vectors * sizeof(uint8_t*));
    
    // Generar cada combinación y su vector correspondiente
    for (size_t i = 0; i < result->num_vectors; i++) {
        // 1. Generar combinación binaria para los cosets
        result->combinations[i] = (uint8_t*)malloc(cl->len * sizeof(uint8_t));
        if (!result->combinations[i]) {
            fprintf(stderr, "ERROR: sin memoria para combination %zu\n", i);
            goto error_cleanup;
        }
        
        for (size_t j = 0; j < cl->len; j++) {
            result->combinations[i][j] = (i >> j) & 1;
        }
        
        // 2. Generar vector completo para esta combinación
        result->vectors[i] = generate_vector_for_combination(cl, result->combinations[i], N);
        if (!result->vectors[i]) {
            fprintf(stderr, "ERROR: sin memoria para vector %zu\n", i);
            goto error_cleanup;
        }
    }
    
    return result;
    
error_cleanup:
    // Limpiar memoria parcialmente asignada
    for (size_t j = 0; j < result->num_vectors; j++) {
        if (result->combinations && result->combinations[j]) {
            free(result->combinations[j]);
        }
        if (result->vectors && result->vectors[j]) {
            free(result->vectors[j]);
        }
    }
    if (result->combinations) free(result->combinations);
    if (result->vectors) free(result->vectors);
    free(result);
    return NULL;
}

/* Imprime todas las combinaciones y sus vectores */
void print_coset_vectors(const CosetVectors *cv, const CosetList *cl) {
    if (!cv || !cl) {
        printf("(datos nulos)\n");
        return;
    }
    
    printf("=== Vectores basados en cosets ===\n");
    printf("Cosets: %zu\n", cl->len);
    printf("Vectores: %zu (2^%zu)\n\n", cv->num_vectors, cl->len);
    
    // Imprimir mapeo coset → elementos
    printf("Mapeo cosets:\n ");
    for (size_t i = 0; i < cl->len; i++) {
        const Coset *c = &cl->data[i];
        printf(" C%zu = { ", i);
        for (size_t j = 0; j < c->len; j++) {
            printf("%zu", c->data[j]);
            if (j + 1 < c->len) printf(", ");
        }
        printf(" }\n ");
    }
    printf("\n");
    
    
    for (size_t i = 0; i < cv->num_vectors; i++) {
        printf("Vector %zu: ", i);
        printf("[");
        for (size_t j = 0; j < cv->num_cosets; j++) {
            printf("%u", cv->combinations[i][j]);
            if (j + 1 < cv->num_cosets) printf(", ");
        }
        printf("] -> ");
        print_vector(cv->vectors[i], cv->vector_length);
        printf("\n");
    }
    printf("\n");
}


void print_first_k_vectors(const CosetVectors *cv, size_t k) {
    if (!cv) {
        printf("(datos nulos)\n");
        return;
    }
    
    size_t limit = (k < cv->num_vectors) ? k : cv->num_vectors;
    
    printf("=== Primeros %zu vectores ===\n", limit);
    for (size_t i = 0; i < limit; i++) {
        printf("Vector %zu: ", i);
        print_vector(cv->vectors[i], cv->vector_length);
        printf("\n");
    }
}


void free_coset_vectors(CosetVectors *cv) {
    if (!cv) return;
    
    if (cv->combinations) {
        for (size_t i = 0; i < cv->num_vectors; i++) {
            if (cv->combinations[i]) free(cv->combinations[i]);
        }
        free(cv->combinations);
    }
    
    if (cv->vectors) {
        for (size_t i = 0; i < cv->num_vectors; i++) {
            if (cv->vectors[i]) free(cv->vectors[i]);
        }
        free(cv->vectors);
    }
    
    free(cv);
}

int save_vectors(const CosetVectors *cv, const char *filename) {
    if (!cv || !filename) {
        fprintf(stderr, "ERROR: parámetros inválidos\n");
        return 0;
    }
    
    FILE *file = fopen(filename, "w");
    if (!file) {
        fprintf(stderr, "ERROR: no se pudo abrir archivo '%s'\n", filename);
        return 0;
    }
    
    for (size_t i = 0; i < cv->num_vectors; i++) {
        for (size_t j = 0; j < cv->vector_length; j++) {
            fprintf(file, "%u", cv->vectors[i][j]);
        }
        fprintf(file, "\n");
    }
    
    fclose(file);
    printf("Vectores guardados\n");
    return 1;
}

/* Convierte vector binario a complejo (-1/+1) */
static void binary_to_complex(const uint8_t *binary, double complex *complex_arr, size_t N) {
    for (size_t i = 0; i < N; i++) {
        complex_arr[i] = (binary[i] == 0) ? 1.0 + 0.0*I : 0 + 0.0*I;
    }
}

/* Versión 2: Guarda módulo DFT redondeado Y segundo fichero con transformación */
int save_dft_both(const CosetVectors *cv, const char *filename1,  const char *filename2, size_t N) {
    if (!cv || !filename1 || !filename2) {
        fprintf(stderr, "ERROR: parámetros inválidos\n");
        return 0;
    }
    
    // Calcular valor constante (N-3)/2
    double constant = (N+1)/2;
    
    FILE *file1 = fopen(filename1, "w");
    FILE *file2 = fopen(filename2, "w");
    
    if (!file1 || !file2) {
        fprintf(stderr, "ERROR: no se pudo abrir archivos\n");
        if (file1) fclose(file1);
        if (file2) fclose(file2);
        return 0;
    }
    
    double complex *time_domain = (double complex*)malloc(N * sizeof(double complex));
    double complex *freq_domain = (double complex*)malloc(N * sizeof(double complex));
    
    if (!time_domain || !freq_domain) {
        fprintf(stderr, "ERROR: sin memoria para DFT\n");
        fclose(file1); fclose(file2);
        if (time_domain) free(time_domain);
        if (freq_domain) free(freq_domain);
        return 0;
    }
    
    //printf("Procesando %zu vectores...\n", cv->num_vectors);
    //printf("Valor constante: (N-3)/2 = (%.0f-3)/2 = %.1f\n", (double)N, constant);
    
    for (size_t i = 0; i < cv->num_vectors; i++) {
        
        binary_to_complex(cv->vectors[i], time_domain, N);
        dft(time_domain, freq_domain, N);
        
        
        // Escribir en primer fich: módulo redondeado
        printf("********\n");
        for (size_t j = 1; j < N; j++) {

            double modulo = pow(cabs(freq_domain[j]),2);
            double rounded = rint(modulo);

            printf("%f : %f : %f : %f\n",creal(time_domain[j]), creal(freq_domain[j]), modulo, rounded);
            fprintf(file1, "%.0f", rounded);
            if (j + 1 < N) fprintf(file1, " ");
        }
        fprintf(file1, "\n");
        
        // Escribir en segundo fich: (N-3)/2 - módulo redondeado
        for (size_t j = 1; j < N; j++) {
            double modulo = pow(cabs(freq_domain[j]),2);
            double rounded = rint(modulo);
            double transformed = constant - rounded;
            
            fprintf(file2, "%.0f", rint(transformed));
            
            if (j + 1 < N) fprintf(file2, " ");
        }
        fprintf(file2, "\n");
    }
    
    free(time_domain);
    free(freq_domain);
    fclose(file1);
    fclose(file2);

    printf("Modulo DFT redondeado guardado\n");
    printf("(|(N-3)/2 - modulo)| redondeado guardado\n\n");
    
    return 1;
}

/* Versión simple: muestra solo las líneas solicitadas */
void show_specific_lines(const char *combinations_path, size_t line1, size_t line2) {
    FILE *file = fopen(combinations_path, "r");
    if (!file) {
        printf("ERROR: No se pudo abrir '%s'\n", combinations_path);
        return;
    }
    
    char *lines[2] = {NULL, NULL};
    size_t current_line = 0;
    char buffer[4096];
    
    // Buscar las dos líneas
    while (fgets(buffer, sizeof(buffer), file)) {
       //printf("Leyendo indice %zu: %s", current_line, buffer);
        if (current_line == line1 || current_line == line2) {
            // Asignar a la posición correcta (0 para line1, 1 para line2)
            int index = (current_line == line1) ? 0 : 1;
            
            // Eliminar newline y asignar memoria
            size_t len = strlen(buffer);
            if (len > 0 && buffer[len-1] == '\n') buffer[len-1] = '\0';
            
            lines[index] = strdup(buffer);
            if (!lines[index]) {
                printf("ERROR: Sin memoria\n");
                fclose(file);
                return;
            }
        }
        
        // Si ya encontramos ambas, salir3
        if (lines[0] && lines[1]) break;
        current_line++;
    }
    
    // Mostrar resultados
    
    if (lines[0]) {
        printf("  - Linea %zu del primer fichero: \t%s\n", line1, lines[0]);
        free(lines[0]);
    }

    if (lines[1]) {
        printf("  - Linea %zu del segundo fichero: \t%s\n", line2, lines[1]);
        free(lines[1]);
    }
    
    fclose(file);
}


void find_matches_files(const char *file1_path, const char *file2_path, const char *combinations_path) {
    
    FILE *file1 = fopen(file1_path, "r");
    FILE *file2 = fopen(file2_path, "r");
    
    if (!file1 || !file2) {
        fprintf(stderr, "ERROR: No se pudo abrir los archivos\n");
        if (file1) fclose(file1);
        if (file2) fclose(file2);
        return;
    }
    
    printf("=== BUSCANDO LPS ===\n");
    
    // Primero, leer todo el archivo 2 en un array de hashes simples
    
    char line1[4096];
    char line2[4096];
    size_t line1_num = 0;
    int found_any = 0;
    
    // Contar líneas en archivo 2 para estimar
    size_t total_lines_2 = 0;
    rewind(file2);
    while (fgets(line2, sizeof(line2), file2)) {
        total_lines_2++;
    }
    rewind(file2);
    
    // Para cada línea en archivo 1
    while (fgets(line1, sizeof(line1), file1)) {
        
        size_t len1 = strlen(line1);
        if (len1 > 0 && line1[len1-1] == '\n') line1[len1-1] = '\0';
        rewind(file2);
        size_t line2_num = 0;
        
        while (fgets(line2, sizeof(line2), file2)) {
            
            size_t len2 = strlen(line2);
            if (len2 > 0 && line2[len2-1] == '\n') line2[len2-1] = '\0';
            
            if (strcmp(line1, line2) == 0) {
                printf("LP:\n");
                printf("  - Linea %zu en primer archivo y Linea %zu en segundo archivo: %s\n", line1_num,line2_num, line1);
                
                show_specific_lines(combinations_path, line1_num, line2_num);
                
                found_any = 1;
                // No break para encontrar todas
            }
            line2_num++;
        }
        
        line1_num++;
        
    }
    
    fclose(file1);
    fclose(file2);
}



int main(void) {
    size_t N = 45;
    size_t k = 7;
    printf("N=%d k=%d\n\n",N,k);

    // Calcular cosets
    CosetList cl = cyclotomic_cosets(k, N);
    //cosetlist_print(&cl);

    // Generar vectores
    CosetVectors *cv = generate_coset_vectors(&cl, N);
    
    if (cv) {
        print_coset_vectors(cv, &cl);
        
        // Guardar cosets
        save_vectors(cv,"combinations.txt");
        // Guardar ficheros de dft y cte-cft
        save_dft_both(cv, "dft.txt","cte-dft.txt", N);

        // Buscar LPS en los ficheros
        find_matches_files("dft.txt","cte-dft.txt", "combinations.txt");
        
        free_coset_vectors(cv);
    }
    
    free_cosetlist(&cl);
    print("===FIN===")
    getchar();
    return 0;
}