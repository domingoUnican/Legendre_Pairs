#ifndef CYCLOTOMIC_COSETS_H
#define CYCLOTOMIC_COSETS_H

#include <stdint.h>
#include <stdbool.h>

/* --------- Estructuras de datos --------- */

typedef struct {
    size_t *data;       /* elementos del coset */
    size_t len;         /* longitud actual */
    size_t cap;         /* capacidad reservada */
} Coset;

typedef struct {
    Coset *data;        /* array de cosets */
    size_t len;         /* número de cosets */
    size_t cap;         /* capacidad reservada */
} CosetList;

/* --------- Funciones públicas principales --------- */

/**
 * Calcula los cosets ciclotómicos módulo N bajo multiplicación por k.
 * @param k factor de multiplicación
 * @param N módulo
 * @return lista de cosets
 */
CosetList cyclotomic_cosets(size_t k, size_t N);

/**
 * Libera toda la memoria asociada a una CosetList.
 * @param cl puntero a la lista a liberar
 */
void free_cosetlist(CosetList *cl);

/**
 * Imprime los cosets en formato legible.
 * @param cl puntero a la lista de cosets
 */
void cosetlist_print(const CosetList *cl);

/* --------- Funciones de pertenencia --------- */

/**
 * Comprueba si el coset a está contenido en el coset b.
 * @param a primer coset
 * @param b segundo coset
 * @return true si a ⊆ b, false en caso contrario
 */
bool coset_is_subset_of(const Coset *a, const Coset *b);

/**
 * Calcula un vector de pertenencia que indica en qué cosets de la lista
 * está contenido el coset c.
 * @param c coset a verificar
 * @param cl lista de cosets
 * @return vector de 0/1 (debe liberarse con free), NULL en error
 */
uint8_t* coset_membership_vector(const Coset *c, const CosetList *cl);

#endif /* CYCLOTOMIC_COSETS_H */