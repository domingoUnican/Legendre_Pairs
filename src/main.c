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

    /*
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
    */
}


uint8_t** PSD(const CosetVectors * cv, size_t N) {
    double psd_total = 0.0;
    uint8_t resultado[cv->num_vectors][N];
    complex double dft_temp[N];

    for (size_t i = 0; i < cv->num_vectors; i++) {
        dft(dft_temp, (double complex *)cv->vectors[i], N);
        for (size_t j = 0; j < N; j++) {
            
            resultado[i][j] = rint(creal(dft_temp[j]) * creal(dft_temp[j]) + cimag(dft_temp[j]) * cimag(dft_temp[j]));
        }
    }

    return resultado;   
    
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

static void legendre_sequence(int p,  uint8_t *sequence)
{
    sequence[0] = 5;
    for (int n = 1; n < p; n++) {
        uint8_t legendre_symbol = 3;
        for (int k = 1; k <= (p - 1) / 2; k++) {
            if ((n % p) == (k * k % p)) {
                legendre_symbol = 6;
                break;
            }
        }
        sequence[n] = legendre_symbol;
    }
}



bool is_compression(int N, int p, uint8_t *sequence, uint8_t *compressed)
{
    int temp[p];
    for (int i = 0; i < p; i++) {
        temp[i] = compressed[i];
    }
    for (int i = 0; i < p; i++)
    {
        for (int j = 0; p * j < N; j++)
        {
            temp[i] -= sequence[(i + p * j) % N];
        }
        if (temp[i] != 0)
        {
            return false;
        }
    }
    return true;

}

bool is_less_than_compression(int N, int p,uint8_t *sequence, uint8_t *compressed, uint8_t *bound_sequence)
{
    uint8_t temp;
    uint8_t temp_bound;
    for (int i = 0; i < p; i++)
    {
        temp = compressed[i];
        temp_bound = 0;
        for (int j = 0; p * j < N; j++)
        {
            temp -= sequence[(i + p * j) % N];
            temp_bound += bound_sequence[(i + p * j) % N];
        }
        if (temp < 0 )
        {
            return false;
        }
        if (false)
        {
            printf("Bound exceeded at index %d: temp = %u, bound = %u\n", i, temp, temp_bound);
            print_bits(sequence, N, "Sequence0");
            print_bits(bound_sequence, N, "Sequence1");
            for (int j = 0; p * j < N; j++) {
                printf("Index %d: sequence contribution = %u, bound contribution = %u\n", (i + p * j) % N, sequence[(i + p * j) % N], bound_sequence[(i + p * j) % N]);
            }
            print_bits(compressed, p, "Compressed0");

        }
    }
    return true;

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
        //printf("********\n");
        for (size_t j = 1; j < N; j++) {

            double modulo = pow(cabs(freq_domain[j]),2);
            double rounded = rint(modulo);

            //printf("%f : %f : %f : %f\n",creal(time_domain[j]), creal(freq_domain[j]), modulo, rounded);
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
    printf("(|cte - modulo)| redondeado guardado\n\n");
    
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

    // Escritura en lp.txt con el formato par1\npar2\n\n
    if (lines[0] && lines[1] && lp_file) {
        printf(" - LP ENCONTRADO\n");
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

void process_and_filter_vectors(const CosetList *cl, size_t N) {
    if (!cl) return;

    size_t num_cosets = cl->len;
    size_t num_vectors = 1ULL << num_cosets;
    
    // Umbral de la condición: (N + 1) / 2
    double threshold = ((double)N + 1.0) / 2.0;
    double constant = ((double)N + 1.0) / 2.0; 
    size_t matches_found = 0;

    // Abrir archivos para guardar solo lo que cumpla la condición
    FILE *f_comb = fopen("combinations.txt", "a");
    FILE *f_psd  = fopen("dft.txt", "a");
    FILE *f_cte  = fopen("cte-dft.txt", "a");

    if (!f_comb || !f_psd || !f_cte) {
        printf("ERROR: No se pudieron abrir los archivos.\n");
        return;
    }

    // Buffers de trabajo (RAM constante)
    double complex *time_domain = malloc(N * sizeof(double complex));
    double complex *freq_domain = malloc(N * sizeof(double complex));
    uint8_t *combination = malloc(num_cosets * sizeof(uint8_t));
    double *current_psd = malloc(N * sizeof(double));

    printf("=== INICIANDO PROCESO (N=%zu) ===\n", N);
    printf("Condicion: Max(PSD) < %.2f\n", threshold);

    for (size_t i = 0; i < num_vectors; i++) {
        // 1. Generar combinación binaria de los cosets
        for (size_t j = 0; j < num_cosets; j++) {
            combination[j] = (i >> j) & 1;
        }

        // 2. Mapear cosets al vector de bits de longitud N
        uint8_t *vector_bits = generate_vector_for_combination(cl, combination, N);
        
        // 3. Transformada de Fourier
        binary_to_complex(vector_bits, time_domain, N);
        dft(time_domain, freq_domain, N);

        // 4. Calcular PSD y evaluar el máximo
        double max_psd = -1.0;
        for (size_t j = 1; j < N; j++) { // Empezamos en 1 para omitir componente DC
            current_psd[j] = rint(pow(cabs(freq_domain[j]), 2));
            if (current_psd[j] > max_psd) {
                max_psd = current_psd[j];
            }
        }

        // ==========================================================
        // SI SE CUMPLE: Max(PSD) < (N + 1) / 2
        // ==========================================================
        if (max_psd < threshold) {
            matches_found++;

            // Guardar en combinations.txt
            for (size_t j = 0; j < N; j++) fprintf(f_comb, "%u", vector_bits[j]);
            fprintf(f_comb, "\n");

            // Guardar en dft.txt (PSD) y cte-dft.txt (N+1/2 - PSD)
            for (size_t j = 1; j < N; j++) {
                fprintf(f_psd, "%.0f%s", current_psd[j], (j + 1 < N) ? " " : "");
                
                double transformed = rint(constant - current_psd[j]);
                fprintf(f_cte, "%.0f%s", transformed, (j + 1 < N) ? " " : "");
            }
            fprintf(f_psd, "\n");
            fprintf(f_cte, "\n");
        }

        // Liberar el vector de esta iteración para mantener la RAM limpia
        free(vector_bits);
    }

    // Cierre y limpieza
    fclose(f_comb); fclose(f_psd); fclose(f_cte);
    free(time_domain); free(freq_domain); free(combination); free(current_psd);

    printf("\n=== PROCESO FINALIZADO ===\n");
    printf("Vectores que cumplen la condicion: %zu\n", matches_found);
    printf("Resultados guardados en los ficheros .txt\n");
}

int gcd(int a, int b)
{
    int temp;
    while (b != 0)
    {
        temp = a % b;

        a = b;
        b = temp;
    }
    return a;
}

typedef struct {
    FILE *f_comb;
    FILE *f_psd;
    FILE *f_cte;
    size_t matches_found;
    double threshold;
    double constant;
    double complex **dft_matrix;
    int **psd_matrix;
    size_t N;
    size_t num_cosets;
    double complex *current_dft;
    int *current_psd;
    int *bound_psd;
    uint8_t *compression_a;
    uint8_t *compression_b;
    int p;
    int coset_idx;
    CosetList *cl;
    uint8_t *current_combination;
} DFSContext;

bool check_bound(const DFSContext *ctx) {
    uint8_t *vector_bits = generate_vector_for_combination(ctx->cl, ctx->current_combination, ctx->N);
    uint8_t *temp0 = malloc(ctx->cl->len * sizeof(uint8_t));
    for (size_t j = ctx->coset_idx; j < ctx->cl->len; j++) {
        temp0[j] = 1;
    }
    uint8_t *temp1 = generate_vector_for_combination(ctx->cl, temp0, ctx->N);
    if  ( !is_less_than_compression(ctx->N, ctx->p, vector_bits, ctx->compression_a, temp1))
    {
        if (!is_less_than_compression(ctx->N, ctx->p, vector_bits, ctx->compression_b, temp1)) {
            free(vector_bits);
            free(temp1);
            free(temp0);
            return false;
        }
    }
    free(vector_bits);
    free(temp1);
    free(temp0);
    /*if ( !is_less_compression_a && !is_less_compression_b) {
        return false;
    }*/
    int max_diff = -1;
    for (size_t j = 1; j < ctx->N; j++) {
        int bound_remaining_j = 0;
        size_t start_idx =  ctx->coset_idx + 1;
        for (size_t i = start_idx; i < ctx->num_cosets; i++) {
            bound_remaining_j += ctx->psd_matrix[i][j];
        }
        int diff = ctx->current_psd[j] - bound_remaining_j;
        if (diff > max_diff) {
            max_diff = diff;
        }
    }
    return max_diff <= (int)ctx->threshold;
}

bool is_valid_combination(const DFSContext *ctx) {
    uint8_t *vector_bits = generate_vector_for_combination(ctx->cl, ctx->current_combination, ctx->N);
    uint8_t temp[ctx->N];
    for (size_t j = 0; j < ctx->N; j++) {
        temp[j] = vector_bits[j];
    }
    free(vector_bits);
    bool is_compressed = is_compression(ctx->N, ctx->p, temp, ctx->compression_a);
    if (is_compressed ||  is_compression(ctx->N, ctx->p, temp, ctx->compression_b)){
        int max_psd = -1;
        for (size_t j = 1; j < ctx->N; j++) {
            if (ctx->current_psd[j] > max_psd) {
                max_psd = ctx->current_psd[j];
            }
        }
        return max_psd <= (int)ctx->threshold;
    }
    return false;
}

/* DFS recursivo para explorar combinaciones de cosets */
void dfs_explore_combinations(
    DFSContext *ctx
) {
    // Caso base: hemos asignado todos los cosets
    if (ctx->coset_idx == ctx->cl->len) {
        
        
        // PODA: Si cumple la condición, guardar
        if (is_valid_combination(ctx)) {
            ctx->matches_found++;
            /*for (size_t j = 0; j < cl->len; j++) {
                printf("%u", current_combination[j]);
            //}
            printf(" -> Max PSD: %d\n", max_psd);
            printf("Vector bits: ");*/
            uint8_t *vector_bits = generate_vector_for_combination(ctx->cl, ctx->current_combination, ctx->N);            
            /*for (size_t j = 0; j < ctx->N; j++) {
                printf("%u", vector_bits[j]);
            //}
            printf("\n");*/
            // Guardar vector binario
            for (size_t j = 0; j < ctx->N; j++) {
                fprintf(ctx->f_comb, "%u", vector_bits[j]);
            }
            fprintf(ctx->f_comb, "\n");
            
            // Guardar PSD y transformación
            for (size_t j = 1; 2 * j < ctx->N; j++) {
                int psd_val = (int)rint(pow(cabs(ctx->current_dft[j]), 2));
                fprintf(ctx->f_psd, "%d%s", psd_val, (j + 1 < ctx->N) ? " " : "");
                
                int transformed = (int)rint(ctx->constant - psd_val);
                fprintf(ctx->f_cte, "%d%s", transformed, (j + 1 < ctx->N) ? " " : "");
            }
            fprintf(ctx->f_psd, "\n");
            fprintf(ctx->f_cte, "\n");
            
            free(vector_bits);
        }
        return;
    }
    
    // Rama 1: No incluir el coset actual (valor 0)
    
    if (check_bound(ctx)) {
        /*         printf("Entra por 0 en el %zu:\n", coset_idx);
        for (size_t prueba = 0; prueba < ctx->num_cosets; prueba++) {
            printf("%u", current_combination[prueba]);
        }
        printf("\n"); */
        ctx->current_combination[ctx->coset_idx] = 0;
        ctx->coset_idx++;
        dfs_explore_combinations(ctx);
        ctx->coset_idx--;
    }
    
    // Rama 2: Incluir el coset actual (valor 1)
    for (size_t j = 1; j < ctx->N; j++) {
        ctx->current_dft[j] = ctx->current_dft[j] + ctx->dft_matrix[ctx->coset_idx][j];
        ctx->current_psd[j] = (int)rint(pow(cabs(ctx->current_dft[j]), 2));
    }
    
    if (check_bound(ctx)) {
        //printf("Entra por 1 en el coset %zu.\n", coset_idx);
        ctx->current_combination[ctx->coset_idx] = 1;     
        /*for (size_t prueba = 0; prueba < ctx->num_cosets; prueba++) {
            printf("%u", current_combination[prueba]);
        
        printf("\n");*/
        ctx->coset_idx++;
        dfs_explore_combinations(ctx);
        ctx->coset_idx--;
    }
    
    // Retroceso: deshacer cambios
    for (size_t j = 1; j < ctx->N; j++) {
        ctx->current_dft[j] = ctx->current_dft[j] - ctx->dft_matrix[ctx->coset_idx][j];
        ctx->current_psd[j] = (int)rint(pow(cabs(ctx->current_dft[j]), 2));
    }
}

/* Versión DFS del procesamiento */
void process_and_filter_vectors_dfs(CosetList *cl, size_t N, int p) {
    if (!cl || cl->len == 0) return;

    double threshold = ((double)N + 1.0) / 2.0;
    double constant = ((double)N + 1.0) / 2.0;
    
    // Asignar matrices dinámicamente (OPCIÓN 1: Recomendada)
    double complex **dft_matrix = malloc(cl->len * sizeof(double complex *));
    int **psd_matrix = malloc(cl->len * sizeof(int *));
    
    for (size_t i = 0; i < cl->len; i++) {
        dft_matrix[i] = malloc(N * sizeof(double complex));
        psd_matrix[i] = malloc(N * sizeof(int));
    }
    
    double complex *time_domain = malloc(N * sizeof(double complex));
    double complex *freq_domain = malloc(N * sizeof(double complex));
    int *current_psd = malloc(N * sizeof(int));
    int *bound_psd = malloc(N * sizeof(int));
    uint8_t *compression_a = malloc(p * sizeof(uint8_t));
    uint8_t *compression_b = malloc(p * sizeof(uint8_t));
    double complex *current_dft = malloc(N * sizeof(double complex));
    uint8_t *combination = malloc(cl->len * sizeof(uint8_t));
    legendre_sequence(p, compression_a);
    for (int i = 0; i < p; i++) {
        if (compression_a[i] == 5) {
            compression_b[i] = 5;
        } else if (compression_a[i] == 3) {
            compression_b[i] = 6;
        } else if (compression_a[i] == 6) {
            compression_b[i] = 3;
        }
    }
    // Inicializar
    for (size_t i = 0; i < cl->len; i++) {
        combination[i] = 0;
    }
    for (size_t j = 0; j < N; j++) {
        current_dft[j] = 0.0 + 0.0*I;
        current_psd[j] = 0;
        bound_psd[j] = 0;
    }
    //printf("El valor de cl->len es: %zu\n", cl->len);
    // Calcular DFT y PSD para cada coset
    for (size_t i = 0; i < cl->len; i++) {
        combination[i] = 1;
        uint8_t *vector_bits = generate_vector_for_combination(cl, combination, N);
        for (size_t j = 0; j < N; j++) {
            printf("%u", vector_bits[j]);
        }
        //printf("\n Vector generado para el coset %zu\n", i);
        binary_to_complex(vector_bits, time_domain, N);
        dft(time_domain, freq_domain, N);
        
        for (size_t j = 0; j < N; j++) {
            dft_matrix[i][j] = freq_domain[j];
            psd_matrix[i][j] = (int)rint(pow(cabs(freq_domain[j]), 2));
            bound_psd[j] += psd_matrix[i][j];
        }
        free(vector_bits);
        combination[i] = 0;
    }
    
    free(time_domain);
    free(freq_domain);
    
    // Abrir archivos
    FILE *f_comb = fopen("combinations.txt", "a");
    FILE *f_psd = fopen("dft.txt", "a");
    FILE *f_cte = fopen("cte-dft.txt", "a");

    if (!f_comb || !f_psd || !f_cte) {
        printf("ERROR: No se pudieron abrir los archivos.\n");
        if (f_comb) fclose(f_comb);
        if (f_psd) fclose(f_psd);
        if (f_cte) fclose(f_cte);
        
        // Limpiar memoria antes de salir
        for (size_t i = 0; i < cl->len; i++) {
            free(dft_matrix[i]);
            free(psd_matrix[i]);
        }
        free(dft_matrix);
        free(psd_matrix);
        free(current_dft);
        free(current_psd);
        free(bound_psd);
        free(combination);
        free(compression_a);
        free(compression_b);
        return;
    }

    // Preparar contexto
    DFSContext ctx = {
        .f_comb = f_comb,
        .f_psd = f_psd,
        .f_cte = f_cte,
        .matches_found = 0,
        .threshold = threshold,
        .constant = constant,
        .N = N,
        .p = p,
        .num_cosets = cl->len,
        .dft_matrix = dft_matrix,
        .psd_matrix = psd_matrix,
        .current_dft = current_dft,
        .current_psd = current_psd,
        .bound_psd = bound_psd,
        .compression_a = compression_a,
        .compression_b = compression_b,
        .cl = cl, 
        .current_combination = combination,
        .coset_idx = 0
    };

    printf("=== EXPLORACIÓN DFS (N=%zu, cosets=%zu) ===\n", N, cl->len);
    printf("Condicion: Max(PSD) < %.2f\n\n", threshold);
    // Iniciar DFS desde coset 0
    dfs_explore_combinations(&ctx);
    
    // Limpieza completa de memoria
    for (size_t i = 0; i < cl->len; i++) {
        free(dft_matrix[i]);
        free(psd_matrix[i]);
    }
    free(dft_matrix);
    free(psd_matrix);
    free(current_dft);
    free(current_psd);
    free(bound_psd);
    free(combination);
    
    fclose(f_comb);
    fclose(f_psd);
    fclose(f_cte);

    printf("\n=== EXPLORACIÓN FINALIZADA ===\n");
    printf("Vectores encontrados: %zu\n\n", ctx.matches_found);
}



int main(void) {
    int p = 5;
    size_t N = (size_t)(9*p);
    size_t k;
    int tamanos[N];
    for (size_t i = 0; i < N; i++) {
        tamanos[i] = -1;
    }
    for (size_t k = 1; k<N; k++) {
        if (gcd(N, k) == 1) {
            CosetList cl = cyclotomic_cosets(k, N);
            tamanos[k] = cl.len;
            free_cosetlist(&cl);
        }
    }
    for (size_t minimo_size = 1; minimo_size <= N; minimo_size++) {
        for (int k=2; k<N; k++){
            if (tamanos[k]==minimo_size){
                printf("k=%d tiene %d cosets\n",k,tamanos[k]);
                printf("N=%d k=%d\n\n",N,k);
                CosetList cl = cyclotomic_cosets(k, N);
                printf("Número de cosets: %zu\n", cl.len);
                // Usar versión DFS en lugar de la iterativa
                process_and_filter_vectors_dfs(&cl, N, p);
                free_cosetlist(&cl);
            }
        }
    }
    printf("===FIN===");
    getchar();
    return 0;
}
