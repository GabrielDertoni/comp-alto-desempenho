#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>


#define MAX_GRADE 100
#define SUMS_SIZE 128 // Optmized for some algorithms, must be greater than MAX_GRADE + 1


typedef int64_t Total;
typedef int_fast8_t Grade;

typedef int Sums[SUMS_SIZE];

typedef struct Sizes {
    size_t grades;
    size_t cities;
    size_t regions;
    size_t students;
} Sizes;

typedef struct Best {
    size_t index;
    double average;
} Best;

typedef struct Stats {
    Best best;
    Grade* restrict min;
    Grade* restrict max;
    double* restrict stdev;
    double* restrict median;
    double* restrict average;
} Stats;

typedef struct Country {
    Grade min;
    Grade max;
    double stdev;
    double median;
    double average;
} Country;


#define zero_sums(sums) __builtin_memset(sums, 0, SUMS_SIZE * sizeof(int))


static inline void free_stats(Stats* stats);
static inline Stats allocate_stats(size_t size);
static inline void merge_sums(const Sums merged, Sums output);
static inline Grade* random_grades(unsigned seed, size_t number_grades);
static inline void calculate_sums(const Grade* grades, size_t size, Sums output);
static inline void print_results(double time_taken, const Sizes* sizes, const Stats* restrict cities, const Stats* restrict regions, const Country* country);
static inline void calculate_stats(const Sums sums, Grade* restrict min, Grade* restrict max, double* restrict median, double* restrict average, double*  restrict stdev);


static inline void compute(const Grade* grades, const Sizes* sizes, Stats* restrict cities, Stats* restrict regions, Country* country) {
    Sums country_sums;
    zero_sums(country_sums);

    Best best_city = (Best){ .index = 0, .average = 0.0 };
    Best best_region  = (Best){ .index = 0, .average = 0.0 };

    for (size_t i = 0; i < sizes->regions; i++) {
        Sums region_sums;
        zero_sums(region_sums);

        for (size_t j = 0; j < sizes->cities; j++) {
            const size_t city = i * sizes->cities + j;

            Sums city_sums;
            zero_sums(city_sums);

            calculate_sums(grades + i * sizes->cities * sizes->students + j * sizes->students, sizes->students, city_sums);
            calculate_stats(city_sums, &cities->min[city], &cities->max[city], &cities->median[city], &cities->average[city], &cities->stdev[city]);

            merge_sums(city_sums, region_sums);

            if (cities->average[city] > best_city.average) {
                best_city = (Best){ .index = city, .average = cities->average[city] };
            }
        }

        calculate_stats(region_sums, &regions->min[i], &regions->max[i], &regions->median[i], &regions->average[i], &regions->stdev[i]);
        merge_sums(region_sums, country_sums);

        if (regions->average[i] > best_region.average) {
            best_region = (Best){ .index = i, .average = regions->average[i] };
        }
    }

    calculate_stats(country_sums, &country->min, &country->max, &country->median, &country->average, &country->stdev);

    cities->best = best_city;
    regions->best = best_region;
}

int main() {
    Sizes sizes;
    unsigned seed;

    if (scanf("%zu %zu %zu %u", &sizes.regions, &sizes.cities, &sizes.students, &seed) != 4) {
        fputs("Failed to read input", stderr);
        exit(EXIT_FAILURE);
    }

    sizes.grades = sizes.regions * sizes.cities * sizes.students;
    Grade* grades = random_grades(seed, sizes.grades);

    Country country;
    Stats regions = allocate_stats(sizes.regions);
    Stats cities = allocate_stats(sizes.regions * sizes.cities);

    double start_time = omp_get_wtime();
    compute(grades, &sizes, &cities, &regions, &country);

    print_results(omp_get_wtime() - start_time, &sizes, &cities, &regions, &country);

    free(grades);
    free_stats(&cities);
    free_stats(&regions);

    return EXIT_SUCCESS;
}


// Calculate

const int* lower_bound(const int* begin, int value) { // Only works with size 125
    begin += 64 * (begin[63] < value);
    begin += 32 * (begin[31] < value);
    begin += 16 * (begin[15] < value);
    begin +=  8 * (begin[ 7] < value);
    begin +=  4 * (begin[ 3] < value);
    begin +=  2 * (begin[ 1] < value);
    begin +=  1 * (begin[ 0] < value);
    begin +=  1 * (begin[ 0] < value);

    return begin;
}

const int* upper_bound(const int* begin, int value) { // Only works with size 125
    begin += 64 * (begin[63] <= value);
    begin += 32 * (begin[31] <= value);
    begin += 16 * (begin[15] <= value);
    begin +=  8 * (begin[ 7] <= value);
    begin +=  4 * (begin[ 3] <= value);
    begin +=  2 * (begin[ 1] <= value);
    begin +=  1 * (begin[ 0] <= value);
    begin +=  1 * (begin[ 0] <= value);

    return begin;
}

static inline size_t distance(const int* begin, const int* it) {
    return it - begin;
}

static inline void calculate_totals(const Sums sums, Total* total, Total* total_sq) {
    *total = *total_sq = 0;

    for (size_t i = 1; i < SUMS_SIZE; i++) {
        size_t count = sums[i] - sums[i - 1];

        *total += count * i;
        *total_sq += count * i * i;
    }
}

static inline double calculate_median(const Sums sums, size_t size) {
    size_t half_size = size / 2;

    if (size % 2 == 0) { // Even
        const size_t up = distance(sums, lower_bound(sums, half_size));
        const size_t down = distance(sums, lower_bound(sums, half_size + 1));

        return (up + down) / 2.0;
    }

    return distance(sums, lower_bound(sums, half_size));
}

static inline void calculate_stats(const Sums sums, Grade* restrict min, Grade* restrict max, double* restrict median, double* restrict average, double* restrict stdev) {
    size_t size = sums[MAX_GRADE];

    Total total, total_sq;
    calculate_totals(sums, &total, &total_sq);

    *min = distance(sums, upper_bound(sums, 0));
    *max = distance(sums, lower_bound(sums, size));

    *median = calculate_median(sums, size);

    *average = (double)total / size;
    *stdev = sqrt((double)total_sq / size - *average * *average);
}


// Sums

static inline void calculate_sums(const Grade* grades, size_t size, Sums output) {
    for (size_t i = 0; i < size; i++) { // Count
        output[grades[i]]++;
    }

    for (size_t i = 1; i < SUMS_SIZE; i++) { // Prefix Sum
        output[i] += output[i - 1];
    }
}

static inline void merge_sums(const Sums merged, Sums output) {
    for (size_t i = 0; i < SUMS_SIZE; i++) {
        output[i] += merged[i];
    }
}


// Grades

static inline Grade* random_grades(unsigned seed, size_t number_grades) {
    Grade* grades = malloc(number_grades * sizeof(Grade));
    srand(seed);

    if (grades == NULL) {
        fputs("Failed to allocate memory", stderr);
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < number_grades; i++) {
        grades[i] = rand() % (MAX_GRADE + 1);
    }

    return grades;
}


// Stats

static inline Stats allocate_stats(size_t size) {
    Stats stats;

    if (
        (stats.min     = malloc(size * sizeof(Grade)))  == NULL ||
        (stats.max     = malloc(size * sizeof(Grade)))  == NULL ||
        (stats.stdev   = malloc(size * sizeof(double))) == NULL ||
        (stats.median  = malloc(size * sizeof(double))) == NULL ||
        (stats.average = malloc(size * sizeof(double))) == NULL
    ) {
        fputs("Failed to allocate memory", stderr);
        exit(EXIT_FAILURE);
    }

    return stats;
}

static inline void free_stats(Stats* stats) {
    free(stats->min);
    free(stats->max);
    free(stats->stdev);
    free(stats->median);
    free(stats->average);
}


// Results

static inline void print_results(double time_taken, const Sizes* sizes, const Stats* restrict cities, const Stats* restrict regions, const Country* country) {
    for (size_t i = 0; i < sizes->regions; i++) {
        for (size_t j = 0; j < sizes->cities; j++) {
            const size_t city = i * sizes->cities + j;

            printf(
                "Reg %zu - Cid %zu: menor: %d, maior: %d, mediana: %.02lf, média: %.02lf e DP: %.02lf\n",
                i, j, cities->min[city], cities->max[city], cities->median[city], cities->average[city], cities->stdev[city]
            );
        }
        puts("");
    }

    for (size_t i = 0; i < sizes->regions; i++) {
        printf(
            "Reg %zu: menor: %d, maior: %d, mediana: %.02lf, média: %.02lf e DP: %.02lf\n",
            i, regions->min[i], regions->max[i], regions->median[i], regions->average[i], regions->stdev[i]
        );
    }
    puts("");

    printf(
        "Brasil: menor: %d, maior: %d, mediana: %.02lf, média: %.02lf e DP: %.02lf\n",
        country->min, country->max, country->median, country->average, country->stdev
    );
    puts("");

    printf("Melhor região: Região %zu\n", regions->best.index);
    printf("Melhor cidade: Região %zu, Cidade: %zu\n", cities->best.index / sizes->cities, cities->best.index % sizes->cities);
    puts("");

    printf("Tempo de resposta sem considerar E/S, em segundos: %.03lfs\n", time_taken);
}
