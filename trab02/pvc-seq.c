#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>


#define MIN_WEIGHT 1
#define MAX_WEIGHT 999


typedef int_fast8_t node_t;
typedef int_fast16_t edge_t;

typedef struct Minimum {
    edge_t cost;
    node_t* path;
} minimum_t;


static inline node_t* generate_path(size_t nodes);
static inline edge_t* generate_graph(size_t nodes);
static inline minimum_t allocate_minimum(size_t nodes);
static inline void next_path(node_t* path, size_t nodes);
static inline void print_path(node_t* path, size_t nodes);
static inline int_fast64_t factorial(int_fast64_t number);
static inline edge_t path_cost(node_t* path, edge_t* graph, size_t nodes);


int main(int argc, char** argv) {
    if (argc <= 1) {
        fputs("Usage:\n  pvc-seq <nodes> [<seed>]\n", stderr);
        exit(EXIT_FAILURE);
    }

    size_t nodes;
    sscanf(argv[1], "%zu", &nodes);

    if (nodes <= 1) {
        fputs("Invalid number of nodes\n", stderr);
        exit(EXIT_FAILURE);
    }

    if (argc >= 3) {
        int seed = atoi(argv[2]);
        srand(seed);
    }


    node_t* path = generate_path(nodes);
    edge_t* graph = generate_graph(nodes);
    minimum_t minimum = allocate_minimum(nodes);

    clock_t time = clock();

    for (int_fast64_t i = 0; i < factorial(nodes - 1); i++) {
        edge_t cost = path_cost(path, graph, nodes);

        if (cost < minimum.cost) {
            minimum.cost = cost;
            memcpy(minimum.path, path, (nodes + 1) * sizeof(node_t));
        }

        next_path(path, nodes);
    }

    time = clock() - time;


    printf("Response Time: %.02lfs\n", time / (double)CLOCKS_PER_SEC);
    printf("Toal Cost: %" PRIdFAST16 "\n\n", minimum.cost);
    print_path(minimum.path, nodes);

    free(path);
    free(graph);
    free(minimum.path);

    return EXIT_SUCCESS;
}


// Utils

static inline int_fast64_t factorial(int_fast64_t number) {
    if (number > 20) {
        fputs("Overflow detected in factorial\n", stderr);
        exit(EXIT_FAILURE);
    }

    int_fast64_t factorial = 1;

    for (int_fast64_t i = 2; i <= number; i++) {
        factorial *= i;
    }

    return factorial;
}

static inline void swap(node_t* first, node_t* second) {
    node_t value = *first;
    *first = *second;
    *second = value;
}

static inline void reverse(node_t* array, size_t size) {
    for (size_t i = 0; i < size / 2; i++) {
        swap(&array[i], &array[size - i - 1]);
    }
}

static inline node_t* upper_bound(node_t value, node_t* array, size_t size) {
    node_t* last = array + size - 1;

    while (size > 0) {
        size_t step = size / 2;

        if (value >= last[-step]) {
            last -= step + 1;
            size -= step + 1;
        }
        else {
            size = step;
        }
    }

    return last;
}


// Graph

static inline edge_t generate_edge() {
    return MIN_WEIGHT + rand() % (MAX_WEIGHT - MIN_WEIGHT + 1);
}

static inline edge_t* generate_graph(size_t nodes) {
    edge_t* graph = malloc(nodes * nodes * sizeof(edge_t));

    if (graph == NULL) {
        fputs("Failed to allocate memory\n", stderr);
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < nodes; i++) {
        for (size_t j = i + 1; j < nodes; j++) {
            graph[i * nodes + j] = generate_edge();
            graph[j * nodes + i] = generate_edge();
        }

        graph[i * nodes + i] = 0;
    }

    return graph;
}


// Minumum

static inline minimum_t allocate_minimum(size_t nodes) {
    minimum_t minimum;
    minimum.cost = nodes * MAX_WEIGHT;

    if ((minimum.path = malloc((nodes + 1) * sizeof(node_t))) == NULL) {
        fputs("Failed to allocate memory\n", stderr);
        exit(EXIT_FAILURE);
    }

    return minimum;
}


// Path

static inline node_t* generate_path(size_t nodes) {
    node_t* path = malloc((nodes + 1) * sizeof(node_t));

    if (path == NULL) {
        fputs("Failed to allocate memory\n", stderr);
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < nodes; i++) {
        path[i] = i;
    }
    path[nodes] = 0; // Loop back

    return path;
}

static inline void next_path(node_t* path, size_t nodes) {
    node_t* permutation = &path[1];
    size_t size = nodes - 1;
    size_t pivot;

    for (pivot = size - 1; pivot > 0; pivot--) {
        if (permutation[pivot - 1] < permutation[pivot]) {
            break;
        }
    }

    if (pivot-- == 0) { // Pivot not found
        reverse(permutation, size);
        return;
    }

    swap(&permutation[pivot], upper_bound(permutation[pivot], &permutation[pivot + 1], size - pivot - 1));
    reverse(&permutation[pivot + 1], size - pivot - 1);
}

static inline edge_t path_cost(node_t* path, edge_t* graph, size_t nodes) {
    edge_t cost = 0;

    for (size_t i = 0; i < nodes; i++) {
        cost += graph[path[i] * nodes + path[i + 1]];
    }

    return cost;
}

static inline void print_path(node_t* path, size_t nodes) {
    for (size_t i = 0; i <= nodes; i++) {
        printf("%" PRIdFAST8 " ", path[i]);
    }

    putchar('\n');
}
