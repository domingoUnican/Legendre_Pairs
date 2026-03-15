#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cyclotomic_cosets.h"
#include "mymath.h"
#include <complex.h>
#include <math.h>
#include "pairs_reader.h"
#define ALFABET_SIZE 2
// New version of gcd that works for negative numbers as well
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
            free(temp);
            free(temp_bound);
            return false;
        }
    }
    free(temp);
    free(temp_bound);
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
    int **candidate_pairs;
    size_t num_candidate_pairs;
    size_t dimension_candidate_pairs;
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
bool is_less_than_candidates(const DFSContext *ctx, int *sequence, int *bound_sequence)
{
    int compresion_count =(int) ( ctx->N / ctx->dimension_candidate_pairs );
    for (size_t i=0; i<ctx->dimension_candidate_pairs; i++) {
        for (size_t j=1; j<compresion_count; j++) {
            sequence[i] += sequence[i + j * ctx->dimension_candidate_pairs];
            bound_sequence[i] += bound_sequence[i + j * ctx->dimension_candidate_pairs];
        }
    }
    bool matches_a = false;
    for (size_t pos = 0; (pos < ctx->num_candidate_pairs) && !matches_a; pos++) {
        matches_a = true;
        for (size_t j = 0; (j < ctx->dimension_candidate_pairs) && matches_a; j++) {
            if (sequence[j] > ctx->candidate_pairs[pos][j] || bound_sequence[j] + sequence[j] < ctx->candidate_pairs[pos][j]) {
                matches_a = false;
            }
        }
    }
    return matches_a; // No se encontró ningún par que coincida
}


bool is_compression_of_candidates(const DFSContext *ctx, int *sequence)
{
    int compresion_count =(int) ( ctx->N / ctx->dimension_candidate_pairs );
    for (size_t i=0; i<ctx->dimension_candidate_pairs; i++) {
        for (size_t j=1; j<compresion_count; j++) {
            sequence[i] += sequence[i + j * ctx->dimension_candidate_pairs];
        }
    }
    bool matches_a = false;
    for (size_t pos = 0; (pos < ctx->num_candidate_pairs) && !matches_a; pos++) {
        bool matches_a = true;
        for (size_t j = 0; (j < ctx->dimension_candidate_pairs) && matches_a; j++) {
            if (sequence[j] != ctx->candidate_pairs[pos][j]) {
                matches_a = false;
            }
        }
    }
    return matches_a; // No se encontró ningún par que coincida
}

bool check_bound(const DFSContext *ctx) {
    int *vector_bits = generate_vector_for_combination(ctx->cl, ctx->current_combination, ctx->N);
    int *temp0 = malloc(ctx->cl->len * sizeof(int));
    for (size_t j = ctx->coset_idx; j < ctx->cl->len; j++) {
        temp0[j] = 1;
    }
    int *temp1 = generate_vector_for_combination(ctx->cl, temp0, ctx->N);
    bool result = is_less_than_candidates(ctx, vector_bits, temp1);
    
    /* This is for using the compression check*/
    /*
    if  ( !is_less_than_compression(ctx->N, ctx->p, vector_bits, ctx->compression_a, temp1))
    {
        if (!is_less_than_compression(ctx->N, ctx->p, vector_bits, ctx->compression_b, temp1)) {
            result = false;
        }
    }
    */
    free(vector_bits);
    free(temp1);
    free(temp0);
    
    if (!result) {
        return false;
    }
    int max_diff = -1;
    for (size_t j = 1; j < ctx->N; j++) {
        int bound_remaining_j = 0;
        size_t start_idx =  ctx->coset_idx ;
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
    
    bool result = false;
    /* This does not uses the candidates, now we are going to use them*/
    //bool is_compressed = is_compression(ctx->N, ctx->p, temp, ctx->compression_a);
    //if (is_compressed ||  is_compression(ctx->N, ctx->p, temp, ctx->compression_b)){
    if (is_compression_of_candidates(ctx,temp))
    {
        int max_psd = -1;
        for (size_t j = 1; j < ctx->N; j++) {
            if (ctx->current_psd[j] > max_psd) {
                max_psd = ctx->current_psd[j];
            }
        }
        result = (max_psd <= (int)ctx->threshold);
    }
    
    free(temp);
    return result;
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
    for (int alfabet_val = 0; alfabet_val < ALFABET_SIZE; alfabet_val++) {
        ctx->current_combination[ctx->coset_idx] = alfabet_val;
        for (size_t j = 1; j < ctx->N; j++) {
           ctx->current_dft[j] = ctx->current_dft[j] + alfabet_val * ctx->dft_matrix[ctx->coset_idx][j];
            ctx->current_psd[j] = (int)rint(pow(cabs(ctx->current_dft[j]), 2));
        }
        ctx->coset_idx++;
        if (check_bound(ctx)) {
            dfs_explore_combinations(ctx);
        }
        ctx->coset_idx--;
        ctx->current_combination[ctx->coset_idx] = 0;
        for (size_t j = 1; j < ctx->N; j++) {
            ctx->current_dft[j] = ctx->current_dft[j] - alfabet_val * ctx->dft_matrix[ctx->coset_idx][j];
            ctx->current_psd[j] = (int)rint(pow(cabs(ctx->current_dft[j]), 2));
        }
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
    size_t rows, cols;
    int **pairs = read_pairs_file("single_sequences.txt", &rows, &cols);
    if (cols == 33)
    {
        printf("Archivo de pares leído correctamente: %zu filas, %zu columnas\n", rows, cols);
    }
    if (!f_comb || !f_psd || !f_cte || !pairs) {
        printf("ERROR: No se pudieron abrir los archivos o leer los pares.\n");
        if (f_comb) fclose(f_comb);
        if (f_psd) fclose(f_psd);
        if (f_cte) fclose(f_cte);
        if (pairs) free_pairs(pairs, rows);
        
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
    //printf("Nueva secuencia de pares:\n");
    int compresion_count =(int) ( N / cols ); 
    for (size_t i=0; i<rows; i++) {
        for (size_t j=0; j<cols; j++) {
            pairs[i][j] = (compresion_count + pairs[i][j]) >> 1;
            //printf("%d ", pairs[i][j]);
        }
        //printf("\n");
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
        .coset_idx = 0,
        .candidate_pairs = pairs,
        .num_candidate_pairs = rows,
        .dimension_candidate_pairs = cols
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
    free(compression_a);
    free(compression_b);
    if (pairs) free_pairs(pairs, rows);
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
    for (int minimo_size = 2; minimo_size <= N; minimo_size++) {
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