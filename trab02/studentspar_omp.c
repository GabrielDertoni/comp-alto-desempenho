#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include <omp.h>

#include <stdbool.h>
#include <assert.h>

#include "next_permutation.h"

// Generic swap
#define SWAP(type, a, b)      \
    do {                      \
        type *__ptr_a = a;    \
        type *__ptr_b = b;    \
        type tmp = *__ptr_a;  \
        *__ptr_a = *__ptr_b;  \
        *__ptr_b = tmp;       \
    } while(0)

#define DBG(arr, n)                          \
    do {                                     \
        printf("[ %d", (arr)[0]);              \
        for (int __i = 1; __i < (n); __i++)    \
            printf(", %d", (arr)[__i]);        \
        printf(" ]\n");                      \
    } while(0)

typedef struct {
    int cost;
    int n;
    int steps[];
} search_result_t;

static inline search_result_t *min_result(search_result_t *omp_in, search_result_t *omp_out) {
    if (omp_in->cost < omp_out->cost)
        memcpy(omp_out, omp_in, sizeof(search_result_t) + omp_in->n * sizeof(int));
    free(omp_in);
    return omp_out;
}

#pragma omp declare reduction(min_result : search_result_t* :         \
        omp_out = min_result(omp_in, omp_out)) \
        initializer (omp_priv=(omp_orig))

static inline int_fast64_t factorial(int n) {
    int_fast64_t ret = 1;
    for (int_fast64_t i = 2; i <= n; i++)
        ret *= i;
    return ret;
}

// Get the nth lexicographical ordering of `values`, where `index = 0` will permute `values` to be
// the smallest lexicographical permutation and `index = n! - 1` will be the biggest. Thus, `index`
// must be in the range [0, n!). Assumes `values` is sorted in non-decreasing order. O(n^2).
//
// source: https://stackoverflow.com/questions/7918806/finding-n-th-permutation-without-computing-others
void nth_permutation(int *values, int n, int_fast64_t index) {
    int_fast64_t fact = factorial(n);

    for (int i = 0; i < n; i++) {
        fact /= n - i;
        int_fast64_t idx = index / fact;
        index %= fact;
        int val = values[idx];

        memmove(values + 1, values, idx * sizeof(int));
        *(values++) = val;
    }
}

static inline int path_cost(const int *restrict adj_mat, const int *restrict path, int n) {
    int cost = 0;
    for (int i = 1; i <= n; i++)
        cost += adj_mat[path[i-1] * n + path[i]];
    return cost;
}

search_result_t *min_from_permutations(const int *adj_mat, int n, int_fast64_t perm_start, int_fast64_t n_perm) {
    search_result_t *best = (search_result_t*)malloc(sizeof(search_result_t) + (n + 1) * sizeof(int));
    best->cost = INT_MAX;
    best->n = n + 1;

    search_result_t *curr = (search_result_t*)malloc(sizeof(search_result_t) + (n + 1) * sizeof(int));
    for (int i = 0; i < n; i++)
        curr->steps[i] = i;
    curr->steps[n] = 0;
    curr->n = n + 1;

    nth_permutation(curr->steps + 1, n - 1, perm_start);

    for (int_fast64_t perm = perm_start; perm < perm_start + n_perm; perm++) {
        curr->cost = path_cost(adj_mat, curr->steps, n);
        if (curr->cost < best->cost) {
            memcpy(best->steps, curr->steps, (n + 1) * sizeof(int));
            best->cost = curr->cost;
        }

        if (perm + 1 < perm_start + n_perm)
            next_permutation_opt(curr->steps + 1, n - 1);
    }

    free(curr);
    return best;
}

void print_time(double time_in_secs) {
    if (time_in_secs < 1e-12) printf("<too small to measure>");
    static char *time_scale_names[] = { "s", "ms", "us", "ns", "ps" };
    double scaled = time_in_secs;
    int scale = 0;
    while ((int)scaled == 0) {
        scale++;
        scaled *= 1000.0;
    }
    printf("%.02lf%s", scaled, time_scale_names[scale]);
}

int main(int argc, char *argv[]) {
    if (argc < 2) return 1;
    srand(42);

    int n = atoi(argv[1]);
    int *adj_mat = (int*)malloc(n * n * sizeof(int));
    for (int i = 0; i < n * n; i++)
        adj_mat[i] = random() % 1000;

    for (int i = 0; i < n; i++)
        DBG(((int*)&adj_mat[i * n]), n);

    double time = omp_get_wtime();
    int_fast64_t fact = factorial(n - 1);
    int t = omp_get_max_threads();
    printf("Computando %ld caminhos possíveis, utilizando %d threads\n", fact, t);

    search_result_t *min = (search_result_t*)malloc(sizeof(search_result_t) + (n + 1) * sizeof(int));
    min->cost = INT_MAX;
    min->n = n + 1;

    #pragma omp parallel reduction(min_result:min)
    {
        int self = omp_get_thread_num();
        int_fast64_t start_perm = self * (fact / t);
        int_fast64_t n_perm;
        if (self == t - 1)
            n_perm = fact - start_perm;
        else
            n_perm = fact / t;

        min = min_from_permutations(adj_mat, n, start_perm, n_perm);
    }

    time = omp_get_wtime() - time;
    printf("Took: ");
    print_time(time);
    printf("\n");
    printf("Path cost: %d\n", min->cost);
    printf("Caminho: ");
    DBG(min->steps, n + 1);

    free(min);
    free(adj_mat);
    return 0;
}
