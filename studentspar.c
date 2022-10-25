#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#include <unistd.h>

#define N_GRADES 101
#define OPT_SUMS_SZ 128

typedef int pref_sum_t[OPT_SUMS_SZ];

typedef struct {
    double val;
    size_t index;
} indexed_val_t;

#pragma omp declare reduction(argmax : indexed_val_t :         \
        omp_out = omp_out.val > omp_in.val ? omp_out : omp_in) \
        initializer (omp_priv=(omp_orig))

static inline void merge_sums_128(int *restrict sums, const int *restrict other) {
    #pragma omp simd
    for (int i = 0; i < OPT_SUMS_SZ; i++)
        sums[i] += other[i];
}

void fill_random_vector(int_fast8_t *mat, int n) {
    for (int i = 0; i < n; i++)
        mat[i] = rand() % N_GRADES;
}

static inline void pref_sum(int *sums) {
    for (int i = 1; i < OPT_SUMS_SZ; i++)
        sums[i] += sums[i-1];
}

static inline void counting_sort_accum_freq(const int_fast8_t *restrict mat, int *restrict sums, size_t mat_len) {
    for (int i = 0; i < mat_len; i++)
        sums[mat[i]]++;

    pref_sum(sums);
}

const int* lower_bound_128(const int *first, int value) {
    first += 64 * (first[63] < value);
    first += 32 * (first[31] < value);
    first += 16 * (first[15] < value);
    first +=  8 * (first[ 7] < value);
    first +=  4 * (first[ 3] < value);
    first +=  2 * (first[ 1] < value);
    first +=  1 * (first[ 0] < value);
    first +=  1 * (first[ 0] < value);
    return first;
}

// Source: https://en.cppreference.com/w/cpp/algorithm/upper_bound
const int* upper_bound(const int* first, size_t n, int value) {
    const int* it;
    size_t step;

    while (n > 0) {
        it = first;
        step = n / 2;
        it += step;
        if (value >= *it) {
            first = ++it;
            n -= step + 1;
        }
        else
            n = step;
    }
    return first;
}

static inline double median_from_sums_128(const pref_sum_t sums) {
    int count = sums[127];
    if (count % 2 == 0) {
        const int median_up = lower_bound_128(sums, count/2) - sums;
        const int median_down = lower_bound_128(sums, (count/2)+1) - sums;
        return (double)(median_down + median_up) / 2.0;
    } else {
        const int median = lower_bound_128(sums, count/2) - sums;
        return (double)median;
    }
}

static inline void compute_statistics_from_sums(const int *restrict sums, int_fast8_t *restrict min,
                                                int_fast8_t *restrict max, double *restrict median,
                                                double *restrict mean, double *restrict stdev) {
    const int n_occur = sums[N_GRADES-1];
    int64_t total = 0;
    int64_t total_sq = 0;
    for (size_t i = 1; i < N_GRADES; i++) {
        int count = sums[i] - sums[i-1];
        total += count * i;
        total_sq += count * i * i;
    }

    *min = upper_bound(sums, N_GRADES, 0) - sums;
    *max = lower_bound_128(sums, n_occur) - sums;
    *median = median_from_sums_128(sums);
    const double m = (double)total / (double)n_occur;
    *mean = m;
    *stdev = sqrt((double)total_sq / (double)n_occur - m * m);
}

void compute_all_statistics(const int_fast8_t *mat, int r, int c, int a,
                            // City
                            int_fast8_t *restrict min_city, int_fast8_t *restrict max_city,
                            double *restrict median_city, double *restrict mean_city,
                            double *restrict stdev_city,
                            // Region
                            int_fast8_t *restrict min_reg, int_fast8_t *restrict max_reg,
                            double *restrict median_reg, double *restrict mean_reg,
                            double *restrict stdev_reg,
                            // Brasil
                            int_fast8_t *restrict min_total, int_fast8_t *restrict max_total,
                            double *restrict median_total, double *restrict mean_total,
                            double *restrict stdev_total, int *restrict best_reg,
                            int *restrict best_city_reg, int *restrict best_city) {
    const size_t ngrades_per_region = c * a;

    pref_sum_t sums_total;
    __builtin_memset(sums_total, 0, sizeof(sums_total));

    indexed_val_t reg_argmax  = { .index = -1, .val = -1 };
    indexed_val_t city_argmax = { .index = -1, .val = -1 };

    int t = omp_get_max_threads();

    #pragma omp parallel for             \
        reduction(argmax:reg_argmax)     \
        reduction(argmax:city_argmax)    \
        reduction(+:sums_total[:128])    \
        schedule(dynamic)                \
        if(t <= r)
    for (int reg = 0; reg < r; reg++) {
        pref_sum_t sums_reg;
        __builtin_memset(sums_reg, 0, sizeof(sums_reg));

        #pragma omp parallel for             \
            reduction(argmax:city_argmax)    \
            reduction(+:sums_reg[:128])      \
            schedule(dynamic)
        for (int city = 0; city < c; city++) {
            const int i = reg * ngrades_per_region + city * a;
            const int j = reg * c + city;

            pref_sum_t sums;
            __builtin_memset(sums, 0, sizeof(sums));
            counting_sort_accum_freq(mat + i, sums, a);

            compute_statistics_from_sums(sums, &min_city[j], &max_city[j],
                                         &median_city[j], &mean_city[j], &stdev_city[j]);

            merge_sums_128(sums_reg, sums);
            if (city_argmax.val < mean_city[j]) {
                city_argmax.index = j;
                city_argmax.val = mean_city[j];
            }
        }

        compute_statistics_from_sums(sums_reg, &min_reg[reg], &max_reg[reg],
                                     &median_reg[reg], &mean_reg[reg], &stdev_reg[reg]);

        merge_sums_128(sums_total, sums_reg);
        if (reg_argmax.val < mean_reg[reg]) {
            reg_argmax.index = reg;
            reg_argmax.val = mean_reg[reg];
        }
    }

    compute_statistics_from_sums(sums_total, min_total, max_total,
                                 median_total, mean_total, stdev_total);

    *best_reg = reg_argmax.index;
    *best_city_reg = city_argmax.index / c;
    *best_city = city_argmax.index % c;
}

int main(int argc, char *argv[]) {
    // Leitura dos dados de entrada
    size_t r, c, a;
    int seed;
#ifndef PERF
    if (!scanf("%zu %zu %zu %d", &r, &c, &a, &seed)) return 1;
#else
    if (argc < 5) return 1;
    r = atoi(argv[1]);
    c = atoi(argv[2]);
    a = atoi(argv[3]);
    seed = atoi(argv[4]);
#endif

    const size_t n = r * c * a;
    const size_t ncity = r * c;
    const size_t ngrades_per_region = c * a;
    int_fast8_t* mat = (int_fast8_t*)malloc(n * sizeof(int_fast8_t));
    // City
    int_fast8_t* min_city = (int_fast8_t*)malloc(ncity * sizeof(int_fast8_t));
    int_fast8_t* max_city = (int_fast8_t*)malloc(ncity * sizeof(int_fast8_t));
    double* median_city = (double*)malloc(ncity * sizeof(double));
    double* mean_city = (double*)malloc(ncity * sizeof(double));
    double* stdev_city = (double*)malloc(ncity * sizeof(double));

    // Region
    int_fast8_t* min_reg = (int_fast8_t*)malloc(r * sizeof(int_fast8_t));
    int_fast8_t* max_reg = (int_fast8_t*)malloc(r * sizeof(int_fast8_t));
    double* median_reg = (double*)malloc(r * sizeof(double));
    double* mean_reg = (double*)malloc(r * sizeof(double));
    double* stdev_reg = (double*)malloc(r * sizeof(double));

    // Brasil
    int best_reg, best_city_reg, best_city;
    int_fast8_t min_total, max_total;
    double median_total, mean_total, stdev_total;

    srand(seed);
    fill_random_vector(mat, n);

#ifdef DEBUG
    for (int reg = 0; reg < r; reg++) {
        printf("Região %d\n", reg);
        for (int city = 0; city < c; city++) {
            for (int grade = 0; grade < a; grade++) {
                const int i = reg * ngrades_per_region + city * a + grade;
                printf("%d ", mat[i]);
            }
            printf("\n");
        }
        printf("\n");
    }
#endif

    double start_time = omp_get_wtime();
    compute_all_statistics(mat, r, c, a, min_city, max_city, median_city, mean_city, stdev_city,
                                min_reg, max_reg, median_reg, mean_reg, stdev_reg,
                                &min_total, &max_total, &median_total, &mean_total, &stdev_total, &best_reg,
                                &best_city_reg, &best_city);
    double time_taken = omp_get_wtime() - start_time;

#ifndef PERF
    for (int reg = 0; reg < r; reg++) {
        for (int city = 0; city < c; city++) {
            int i = reg * c + city;
            printf(
                "Reg %d - Cid %d: menor: %d, maior: %d, mediana: %.02lf, média: %.02lf e DP: %.02lf\n",
                reg, city, min_city[i], max_city[i], median_city[i], mean_city[i], stdev_city[i]
            );
        }
        printf("\n");
    }

    for (int reg = 0; reg < r; reg++) {
        printf(
            "Reg %d: menor: %d, maior: %d, mediana: %.02lf, média: %.02lf e DP: %.02lf\n",
            reg, min_reg[reg], max_reg[reg], median_reg[reg], mean_reg[reg], stdev_reg[reg]
        );
    }
    printf("\n");

    printf(
        "Brasil: menor: %d, maior: %d, mediana: %.02lf, média: %.02lf e DP: %.02lf\n",
        min_total, max_total, median_total, mean_total, stdev_total
    );
    printf("\n");

    printf("Melhor região: Região %d\n", best_reg);
    printf("Melhor cidade: Região %d, Cidade: %d\n", best_city_reg, best_city);
    printf("\n");
#else
    // Use variables to prevent them beeing optimized out

    // City
    *(volatile int_fast8_t*)min_city;
    *(volatile int_fast8_t*)max_city;
    *(volatile double*)median_city;
    *(volatile double*)mean_city;
    *(volatile double*)stdev_city;

    // Region
    *(volatile int_fast8_t*)min_reg;
    *(volatile int_fast8_t*)max_reg;
    *(volatile double*)median_reg;
    *(volatile double*)mean_reg;
    *(volatile double*)stdev_reg,

    // Brasil
    *(volatile int_fast8_t*)mat;
    *(volatile int_fast8_t*)&min_total;
    *(volatile int_fast8_t*)&max_total;
    *(volatile double*)&median_total;
    *(volatile double*)&mean_total;
    *(volatile double*)&stdev_total;
    *(volatile int*)&best_reg;
    *(volatile int*)&best_city_reg;
    *(volatile int*)&best_city;
#endif

    printf("Tempo de resposta sem considerar E/S, em segundos: %.03lfs\n", time_taken);

    free(mat);
    free(min_city);
    free(max_city);
    free(median_city);
    free(mean_city);
    free(stdev_city);
    free(min_reg);
    free(max_reg);
    free(median_reg);
    free(mean_reg);
    free(stdev_reg);

    return 0;
}
