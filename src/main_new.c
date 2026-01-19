
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

void show_specific_lines(const char *combinations_path, size_t line1, size_t line2, FILE *lp_file) {
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
        if (current_line == line1 || current_line == line2) {
            int index = (current_line == line1) ? 0 : 1;
            
            size_t len = strlen(buffer);
            if (len > 0 && buffer[len-1] == '\n') buffer[len-1] = '\0';
            
            lines[index] = strdup(buffer);
        }
        
        if (lines[0] && lines[1]) break;
        current_line++;
    }
    
    // Mostrar resultados y guardar en fichero
    if (lines[0]) {
        printf("  - Linea %zu del primer fichero: \t%s\n", line1, lines[0]);
    }

    if (lines[1]) {
        printf("  - Linea %zu del segundo fichero: \t%s\n", line2, lines[1]);
    }

    // Escritura en lp.txt con el formato par1\npar2\n\n
    if (lines[0] && lines[1] && lp_file) {
        fprintf(lp_file, "%s\n", lines[0]);
        fprintf(lp_file, "%s\n\n", lines[1]);
    }
    
    // Limpieza de memoria
    if (lines[0]) free(lines[0]);
    if (lines[1]) free(lines[1]);
    
    fclose(file);
}

typedef struct {
    char *data;
    size_t original_index;
} DFTLine;

void find_matches_files(const char *file1_path, const char *file2_path, const char *combinations_path) {
    FILE *f1 = fopen(file1_path, "r");
    FILE *f2 = fopen(file2_path, "r");
    FILE *lp_file = fopen("lp.txt", "w");
    
    if (!f1 || !f2 || !lp_file) {
        if (f1) fclose(f1); if (f2) fclose(f2);
        return;
    }

    // 1. Cargar fichero 2 en memoria para evitar accesos a disco
    size_t capacity = 1000;
    size_t total_lines2 = 0;
    DFTLine *lines2 = malloc(capacity * sizeof(DFTLine));
    char buffer[4096];

    while (fgets(buffer, sizeof(buffer), f2)) {
        if (total_lines2 >= capacity) {
            capacity *= 2;
            lines2 = realloc(lines2, capacity * sizeof(DFTLine));
        }
        size_t len = strlen(buffer);
        if (len > 0 && buffer[len-1] == '\n') buffer[len-1] = '\0';
        lines2[total_lines2].data = strdup(buffer);
        lines2[total_lines2].original_index = total_lines2;
        total_lines2++;
    }

    printf("=== BUSCANDO LPS ===\n");
    char line1[4096];
    size_t line1_num = 0;
    while (fgets(line1, sizeof(line1), f1)) {
        size_t len = strlen(line1);
        if (len > 0 && line1[len-1] == '\n') line1[len-1] = '\0';

        // Solo buscamos desde line1_num en adelante para evitar duplicados (i,j) y (j,i)
        // Y para no comparar un vector consigo mismo si los ficheros fueran iguales
        for (size_t j = 0; j < total_lines2; j++) {
            if (strcmp(line1, lines2[j].data) == 0) {
                // Si el problema exige que el par sea (A, B) tal que i < j para no duplicar:
                if (line1_num < lines2[j].original_index) { 
                    show_specific_lines(combinations_path, line1_num, lines2[j].original_index, lp_file);
                }
            }
        }
        line1_num++;
    }

    // Limpieza
    for (size_t i = 0; i < total_lines2; i++) free(lines2[i].data);
    free(lines2);
    fclose(f1); fclose(f2); fclose(lp_file);
}
void generar_opciones_comprimidas(int N, int C) {
    if (N % C != 0) {
        printf("Error: N=%d no es divisible por el factor de compresion C=%d\n", N, C);
        return;
    }

    int L = N / C; // Longitud del vector resultante
    printf("N=%d, C=%d -> Longitud comprimida L=%d\n", N, C, L);
    printf("Cada posicion puede variar entre 0 y %d\n", C);
    printf("--- Generando combinaciones ---\n");

    // número total de combinaciones (C + 1)^L
    long long total_combinaciones = 1;
    for (int i = 0; i < L; i++) total_combinaciones *= (C + 1);

    printf("Total de secuencias posibles: %lld\n\n", total_combinaciones);


    int *secuencia = (int *)calloc(L, sizeof(int));

    for (long long i = 0; i < total_combinaciones; i++) {
        // Imprimir secuencia actual
        printf("[");
        for (int j = 0; j < L; j++) {
            printf("%d", secuencia[j]);
            if (j < L - 1) printf(", ");
        }
        printf("]\n");

        for (int j = L - 1; j >= 0; j--) {
            if (secuencia[j] < C) {
                secuencia[j]++;
                break;
            } else {
                secuencia[j] = 0;
            }
        }
    }

    free(secuencia);
}
void analizar_cosets_comprimidos(size_t N, size_t k, size_t C) {
    if (N % C != 0) {
        printf("Error: N=%zu no es divisible por el factor C=%zu\n", N, C);
        return;
    }

    size_t L = N / C; // Nueva longitud del vector

    //Calcular cosets para la longitud reducida L
    CosetList cl_comprimida = cyclotomic_cosets(k, L);

    printf("\n=== Cosets del Espacio Comprimido (Modulo %zu) ===\n", L);
    for (size_t i = 0; i < cl_comprimida.len; i++) {
        printf("C_comp[%zu] = { ", i);
        for (size_t j = 0; j < cl_comprimida.data[i].len; j++) {
            printf("%zu", cl_comprimida.data[i].data[j]);
            if (j + 1 < cl_comprimida.data[i].len) printf(", ");
        }
        printf(" }\n");
    }

    printf("\nCada uno de estos %zu cosets puede tomar valores en el rango [0, %zu]\n", 
            cl_comprimida.len, C);
    
    
    free_cosetlist(&cl_comprimida);
}


int main(void) {
    size_t N = 10;
    size_t k = 2;
    size_t C = 5;
    printf("N=%d k=%d C=%d\n\n",N,k,C);
    generar_opciones_comprimidas(N, C);
    analizar_cosets_comprimidos(N,k,C);
    /*
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
    printf("===FIN===");
    */
    getchar();
    return 0;
}