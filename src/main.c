
#include <stdio.h>
#include <stdlib.h>
#include "cyclotomic_cosets.h"

static void print_bits(const uint8_t *v, size_t n, const char *name) {
    printf("%s = [", name);
    for (size_t i = 0; i < n; ++i) {
        printf("%u", (unsigned)v[i]);
        if (i + 1 < n) printf(", ");
    }
    printf("]\n");
}

int main(void) {
    size_t N = 15, k = 2;
    CosetList cl = cyclotomic_cosets(k, N);
    printf("=== Cosets para N=%zu, k=%zu ===\n", N, k);
    printf("Número de cosets: %zu\n", cl.len);
    cosetlist_print(&cl);

    /* Tomamos el primer coset de la lista y preguntamos en qué cosets está contenido. */
    if (cl.len > 0) {
        const Coset *c = &cl.data[0];
        uint8_t *vec = coset_membership_vector(c, &cl);
        if (vec) {
            print_bits(vec, cl.len, "membership(c, cl)");
            free(vec);
        }
    }

    free_cosetlist(&cl);
    return 0;
}
