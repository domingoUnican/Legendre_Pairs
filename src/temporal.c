#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cyclotomic_cosets.h"
#include "mymath.h"
#include <complex.h>
#include <math.h>


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


void legendre_sequence(int p, int q,  int *sequence, int flag)
{
    sequence[0] = (q*q + 1)>> 1; // esto es dividir por 2, pero usando bit shift para enteros 
    for (int n = 1; n < p; n++) {
        int legendre_symbol = flag;
        for (int k = 1; k <= (p - 1) / 2; k++) {
            if ((n % p) == (k * k % p)) {
                legendre_symbol =  - legendre_symbol ; // cambiar el signo
                break;
            }
        }
        sequence[n] = (q*q + legendre_symbol * q) >> 1;
    }
}

int* generate_vector_for_combination(const CosetList *cl, 
                                        const int *combination, 
                                        size_t N) {
    int *vector = (int*)calloc(N,sizeof(int));
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

static void binary_to_complex(const int *binary, double complex *complex_arr, size_t N) {
    for (size_t i = 0; i < N; i++) {
        complex_arr[i] = (binary[i] == 0) ? 1.0 : 0.0;
    }
}


bool is_less_than_compression(int N, int p, const int *sequence, const int *compressed, const int *bound_sequence)
{
    //return true;
    int *temp = calloc(p, sizeof(int));
    int *temp_bound = calloc(p, sizeof(int));
    for (int i = 0; i < p; i++)
    {
        for (int j = 0; p * j < N; j++)
        {
            temp[i] += sequence[(i + p * j) % N];
            temp_bound[i] += bound_sequence[(i + p * j) % N];
        }
        if (temp[i] > compressed[i] || temp_bound[i] + temp[i] < compressed[i]) 
        {
            /*bool is_lesser = true;
            for (int k = 0; k < 45; k++) {
                is_lesser = is_lesser && (prueba[k] <= bound_sequence[k]);
            }
            if (is_lesser)
            {
                printf("La secuencia no es menor que la compresión, pero sí menor que el bound.\n");
                for (int k = 0; k < 45; k++) {
                    printf("%d ", sequence[k]);
                }
                printf("\n");
                for (int k = 0; k < 45; k++) {
                    printf("%d ", prueba[k]);
                }
                printf("\n");
            }
            printf("\nBound exceeded at index %d: temp = %d, comprosed=%d, bound = %d\n", i, temp[i], compressed[i], temp_bound[i]);
            /*printf("Sequence: ");
            for (int k = 0; k < N; k++) {
                    printf("%d ", sequence[k]);
                }
            printf("\nBound Sequence: ");
            for (int k = 0; k < N; k++) {
                    printf("%d ", bound_sequence[k]);
                }
            printf("\n");*/
            return false;
        }
    }
    return true;

}

typedef struct {
    FILE *f_comb;
    FILE *f_psd;
    FILE *f_cte;
    size_t matches_found;
    double threshold;
    //double constant; constant is the same as threshold. 
    double complex **dft_matrix;
    int **psd_matrix;
    size_t N;
    size_t num_cosets;
    double complex *current_dft;
    int *current_psd;
    int *bound_psd;
    int *compression_a;
    int *compression_b;
    int p;
    int q;
    int coset_idx;
    CosetList *cl;
    int *current_combination;
} DFSContext;

bool is_compression(int N, int p, int *sequence, int *compressed)
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

bool check_bound(const DFSContext *ctx) {
    int *vector_bits = generate_vector_for_combination(ctx->cl, ctx->current_combination, ctx->N);
    int *temp0 = malloc(ctx->cl->len * sizeof(int));
    for (size_t j = ctx->coset_idx; j < ctx->cl->len; j++) {
        temp0[j] = 1;
    }
    int *temp1 = generate_vector_for_combination(ctx->cl, temp0, ctx->N);
    if  ( !is_less_than_compression(ctx->N, ctx->p, vector_bits, ctx->compression_a, temp1))
    {
        if (!is_less_than_compression(ctx->N, ctx->p, vector_bits, ctx->compression_b, temp1)) {
            return false;
            free(vector_bits);
            free(temp1);
            free(temp0);
            return false;
        }
    }
    free(vector_bits);
    free(temp1);
    free(temp0);
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
    int *vector_bits = generate_vector_for_combination(ctx->cl, ctx->current_combination, ctx->N);
    int *temp= malloc(ctx->N * sizeof(int));
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
        free(temp);
        return max_psd <= (int)ctx->threshold;
    }
    free(temp);
    return false;
}

void dfs_explore_combinations( DFSContext *ctx) 
{
    // Caso base: hemos asignado todos los cosets
    if (ctx->coset_idx == ctx->cl->len) {
        
        
        // PODA: Si cumple la condición, guardar
        if (is_valid_combination(ctx)) {
            ctx->matches_found++;
            int *vector_bits = generate_vector_for_combination(ctx->cl, ctx->current_combination, ctx->N);            
            for (size_t j = 0; j < ctx->N; j++) {
                fprintf(ctx->f_comb, "%u", vector_bits[j]);
            }
            fprintf(ctx->f_comb, "\n");
            
            // Guardar PSD y transformación
            for (size_t j = 1; 2 * j < ctx->N; j++) {
                int psd_val = (int)rint(pow(cabs(ctx->current_dft[j]), 2));
                fprintf(ctx->f_psd, "%d%s", psd_val, (j + 1 < ctx->N) ? " " : "");
                
                int transformed = (int)rint(ctx->threshold - psd_val);
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
    ctx->current_combination[ctx->coset_idx] = 1; 
    if (check_bound(ctx)) {
        //printf("Entra por 1 en el coset %zu.\n", coset_idx);
            
        /*for (size_t prueba = 0; prueba < ctx->num_cosets; prueba++) {
            printf("%u", current_combination[prueba]);
        
        printf("\n");*/
        ctx->coset_idx++;
        dfs_explore_combinations(ctx);
        ctx->coset_idx--;
    }
    ctx->current_combination[ctx->coset_idx] = 0; 
    // Retroceso: deshacer cambios
    for (size_t j = 1; j < ctx->N; j++) {
        ctx->current_dft[j] = ctx->current_dft[j] - ctx->dft_matrix[ctx->coset_idx][j];
        ctx->current_psd[j] = (int)rint(pow(cabs(ctx->current_dft[j]), 2));
    }
}


/* Versión DFS del procesamiento */
void process_and_filter_vectors_dfs(CosetList *cl, size_t N, int p, int q) {
    if (!cl || cl->len == 0) return;

    double threshold = ((double)N + 1.0) / 2.0;
    
    // Asignar matrices dinámicamente (OPCIÓN 1: Recomendada)
    double complex **dft_matrix = malloc(cl->len * sizeof(double complex *));
    int **psd_matrix = malloc(cl->len * sizeof(int *));
    
    for (size_t i = 0; i < cl->len; i++) {
        dft_matrix[i] = malloc(N * sizeof(double complex));
        psd_matrix[i] = malloc(N * sizeof(int));
    }
    
    double complex *time_domain = malloc(N * sizeof(double complex));
    double complex *freq_domain = malloc(N * sizeof(double complex));
    int *current_psd = calloc(N, sizeof(int)); // inicializar a 0
    int *bound_psd = calloc(N, sizeof(int)); // inicializar a 0
    int *compression_a = malloc(p * sizeof(int)); 
    int *compression_b = malloc(p * sizeof(int));
    double complex *current_dft = calloc(N, sizeof(double complex)); // inicializar a 0
    int *combination = calloc(cl->len, sizeof(int)); // inicializar a 0
    legendre_sequence(p, q, compression_a, 1);
    legendre_sequence(p, q, compression_b, -1);
    //printf("El valor de cl->len es: %zu\n", cl->len);
    // Calcular DFT y PSD para cada coset
    for (size_t i = 0; i < cl->len; i++) {
        combination[i] = 1;
        int *vector_bits = generate_vector_for_combination(cl, combination, N);
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
        .N = N,
        .p = p,
        .q = q,
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
    int p = 11;
    int q = 3;
    size_t N = (size_t)(p*q*q);
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
    for (int minimo_size = 2; minimo_size < N; minimo_size++) {
        for (int k=1; k< N; k++){
            if (tamanos[k]==minimo_size){
                printf("k=%d tiene %d cosets\n",k,tamanos[k]);
                printf("N=%d k=%d\n\n",N,k);
                CosetList cl = cyclotomic_cosets(k, N);
                // Usar versión DFS en lugar de la iterativa
                process_and_filter_vectors_dfs(&cl, N, p, q);
                free_cosetlist(&cl);
            }
        }
    }
    printf("===FIN===");
    getchar();
    return 0;
}