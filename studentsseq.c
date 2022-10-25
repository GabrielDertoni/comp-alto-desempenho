#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <omp.h>

// Número de notas diferentes
#define N_GRADES 101

// Número ótimo para ser colocado num vetor de frequencias de notas. Esse número precisa ser no
// mínimo igual a `N_GRADES`. No caso, o valor de 128 nos permite otimizar algumas funções
// auxiliares.
#define OPT_SUMS_SZ 128

// Tipo utilizado para variáveis que armazenam somas de prefixo e/out frequencias de notas.
typedef int pref_sum_t[OPT_SUMS_SZ];

// Preenche um vetor de valores aleatórios.
void fill_random_vector(int_fast8_t *mat, int n) {
    for (int i = 0; i < n; i++)
        mat[i] = rand() % N_GRADES;
}

// Faz uma soma de prefixo no vetor `sums`, inplace.
static inline void pref_sum(pref_sum_t sums) {
    for (int i = 1; i < OPT_SUMS_SZ; i++)
        sums[i] += sums[i-1];
}

// Conta o número de ocorrências e armazena no vetor de frequencias `freq`.
static inline void count_freq(const int_fast8_t *restrict mat, pref_sum_t freq, size_t mat_len) {
    for (int i = 0; i < mat_len; i++)
        freq[mat[i]]++;
}

// Retorna um ponteiro que aponta dentro de `sums` ou no elemento seguinte. O ponteiro retornado
// aponta para o primeiro elemento maior ou igual a `value` que se econtra no vetor ordenado
// `sums`. Caso todos os elementos de `sums` sejam menores que `value`, retorna `&sums[128]`.
const int* lower_bound_128(const pref_sum_t sums, int value) {
    // Essa função pôde ser otimizada agressivamente por sabermos o tamanho do vetor `sums`. Dessa
    // forma, podemos fazer unroll de todos os loops e remover branches. Como 128 é uma potência de
    // 2, o algoritmo (de busca binária) sempre levará exatamente o mesmo número de iterações.
    sums += 64 * (sums[63] < value);
    sums += 32 * (sums[31] < value);
    sums += 16 * (sums[15] < value);
    sums +=  8 * (sums[ 7] < value);
    sums +=  4 * (sums[ 3] < value);
    sums +=  2 * (sums[ 1] < value);
    sums +=  1 * (sums[ 0] < value);
    sums +=  1 * (sums[ 0] < value);
    return sums;
}

// Retorna um ponteiro que aponta dentro de `sums` ou no elemento seguinte. O ponteiro retornado
// aponta para o primeiro elemento maior que `value` que se econtra no vetor ordenado `sums`. Caso
// todos os elementos de `sums` sejam menores ou iguais a `value`, retorna `&sums[128]`.
const int* upper_bound_128(const pref_sum_t sums, int value) {
    sums += 64 * (sums[63] <= value);
    sums += 32 * (sums[31] <= value);
    sums += 16 * (sums[15] <= value);
    sums +=  8 * (sums[ 7] <= value);
    sums +=  4 * (sums[ 3] <= value);
    sums +=  2 * (sums[ 1] <= value);
    sums +=  1 * (sums[ 0] <= value);
    sums +=  1 * (sums[ 0] <= value);
    return sums;
}

// Obtém a mediana a partir de uma soma de prefixo das frequencias das notas.
static inline double median_from_sums_128(const pref_sum_t sums) {
    int count = sums[OPT_SUMS_SZ-1];
    if (count % 2 == 0) {
        const int median_up = lower_bound_128(sums, count/2) - sums;
        const int median_down = lower_bound_128(sums, (count/2)+1) - sums;
        return (double)(median_down + median_up) / 2.0;
    } else {
        const int median = lower_bound_128(sums, count/2) - sums;
        return (double)median;
    }
}

// Computa estatísticas de mínimo, máximo, mediana, média e desvio padrão a partir de uma soma de
// prefixo das frequencias das notas.
static inline void compute_statistics_from_sums(const pref_sum_t sums, int_fast8_t *restrict min,
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

    *min = upper_bound_128(sums, 0) - sums;
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

    double best_reg_mean = 0, best_city_mean = 0;

    pref_sum_t sums_total;
    __builtin_memset(sums_total, 0, sizeof(sums_total));
    for (size_t reg = 0; reg < r; reg++) {
        pref_sum_t sums_reg;
        __builtin_memset(sums_reg, 0, sizeof(sums_reg));

        for (size_t city = 0; city < c; city++) {
            const size_t i = reg * ngrades_per_region + city * a;
            const size_t j = reg * c + city;

            pref_sum_t sums;
            __builtin_memset(sums, 0, sizeof(sums));
            count_freq(mat + i, sums, a);
            pref_sum(sums);

            compute_statistics_from_sums(sums, &min_city[j], &max_city[j],
                                         &median_city[j], &mean_city[j], &stdev_city[j]);

            if (best_city_mean < mean_city[j]) {
                best_city_mean = mean_city[j];
                *best_city_reg = reg;
                *best_city = city;
            }

            for (int i = 0; i < OPT_SUMS_SZ; i++)
                sums_reg[i] += sums[i];
        }

        compute_statistics_from_sums(sums_reg, &min_reg[reg], &max_reg[reg],
                                     &median_reg[reg], &mean_reg[reg], &stdev_reg[reg]);

        if (best_reg_mean < mean_reg[reg]) {
            best_reg_mean = mean_reg[reg];
            *best_reg = reg;
        }

        for (int i = 0; i < OPT_SUMS_SZ; i++)
            sums_total[i] += sums_reg[i];
    }

    compute_statistics_from_sums(sums_total, min_total, max_total,
                                 median_total, mean_total, stdev_total);
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
