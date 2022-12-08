#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include <mpi.h>

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

typedef struct {
    int cost;
    int steps[];
} search_result_t;

search_result_t *min_from_permutations(const int *adj_mat, int n, int_fast64_t perm_start, int_fast64_t n_perm) {
    search_result_t *best = (search_result_t*)malloc(sizeof(search_result_t) + (n + 1) * sizeof(int));
    best->cost = INT_MAX;

    search_result_t *curr = (search_result_t*)malloc(sizeof(search_result_t) + (n + 1) * sizeof(int));
    for (int i = 0; i < n; i++)
        curr->steps[i] = i;
    curr->steps[n] = 0;

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

void min_reduce_mpi(void *in, void *inout, int *len, MPI_Datatype *_dptr) {
    search_result_t *in_as_result = (search_result_t*)in;
    search_result_t *inout_as_result = (search_result_t*)inout;
    if (in_as_result->cost < inout_as_result->cost)
        memcpy(inout, in, *len);
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    if (argc < 2) return 1;
    int n = atoi(argv[1]);
    srand(42);

    int nproc, myrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int *adj_mat = (int*)malloc(n * n * sizeof(int));

    if (myrank == 0) {
        for (int i = 0; i < n * n; i++)
            adj_mat[i] = random() % 1000;

        int_fast64_t fact = factorial(n - 1);
        printf("Computando %ld caminhos possíveis, utilizando %d processos\n", fact, nproc);
    }

    MPI_Bcast(adj_mat, n * n, MPI_INT, 0, MPI_COMM_WORLD);

    double time;
    if (myrank == 0)
        time = MPI_Wtime();

    int_fast64_t total_perm = factorial(n - 1);
    int_fast64_t perm_per_proc = total_perm / nproc;
    int_fast64_t perm_start = myrank * perm_per_proc;

    // The last process may get a bit more permutations to compute.
    if (myrank == nproc - 1)
        perm_per_proc = total_perm - perm_start;

    search_result_t *local_min = min_from_permutations(adj_mat, n, perm_start, perm_per_proc);
    search_result_t *global_min;
    if (myrank == 0)
        global_min = (search_result_t*)malloc(sizeof(search_result_t) + (n + 1) * sizeof(int));

    MPI_Op min_reduce_op;
    MPI_Op_create(min_reduce_mpi, true, &min_reduce_op);
    MPI_Reduce(local_min, global_min, sizeof(search_result_t) + (n + 1) * sizeof(int),
               MPI_BYTE, min_reduce_op, 0, MPI_COMM_WORLD);

    if (myrank == 0) {
        time = MPI_Wtime() - time;

        printf("Took: ");
        print_time(time);
        printf("\n");
        printf("Path cost: %d\n", global_min->cost);
        printf("Caminho: [ %d", global_min->steps[0]);
        for (int i = 1; i <= n; i++)
            printf(", %d", global_min->steps[i]);
        printf(" ]\n");

        free(global_min);
    }

    free(local_min);

    MPI_Finalize();

    return 0;
}
