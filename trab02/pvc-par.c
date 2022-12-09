#include <mpi.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <inttypes.h>


#define MIN_WEIGHT 1
#define MAX_WEIGHT 999


typedef int node_t;
typedef int edge_t;

typedef struct Minimum {
    edge_t cost;
    node_t path[];
} minimum_t;


static inline edge_t* allocate_graph(size_t nodes);
static inline minimum_t* allocate_minimum(size_t nodes);
static inline void next_path(node_t* path, size_t nodes);
static inline void print_path(node_t* path, size_t nodes);
static inline int_fast64_t factorial(int_fast64_t number);
static inline void generate_graph(edge_t* graph, size_t nodes);
static inline node_t* generate_path(size_t nodes, int_fast64_t index);
static inline edge_t path_cost(node_t* path, edge_t* graph, size_t nodes);
void reduce_minimum(void* in, void* inout, int* length, MPI_Datatype* data_type);


int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &processes);


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


    double time = 0;
    minimum_t* result = NULL;
    edge_t* graph = allocate_graph(nodes);
    minimum_t* minimum = allocate_minimum(nodes);


    if (rank == 0) { // Process 0 generates the graph and broadcasts it
        if (argc >= 3) {
            int seed = atoi(argv[2]);
            srand(seed);
        }

        result = allocate_minimum(nodes);
        generate_graph(graph, nodes);

        time = MPI_Wtime();
    }

    MPI_Bcast(graph, nodes * nodes, MPI_INT, 0, MPI_COMM_WORLD);


    int_fast64_t total_paths = factorial(nodes - 1);
    int_fast64_t block_size = total_paths / processes;
    int_fast64_t start_path = rank * block_size;

    if (rank == processes - 1) { // Last process gets the remainder of the paths
        block_size = total_paths - start_path;
    }


    node_t* path = generate_path(nodes, start_path);

    for (int_fast64_t i = 0; i < block_size; i++) {
        edge_t cost = path_cost(path, graph, nodes);

        if (cost < minimum->cost) {
            minimum->cost = cost;
            memcpy(minimum->path, path, (nodes + 1) * sizeof(node_t));
        }

        next_path(path, nodes);
    }


    MPI_Op minimum_reduction;
    MPI_Op_create(reduce_minimum, true, &minimum_reduction);
    MPI_Reduce(minimum, result, sizeof(minimum_t) + (nodes + 1) * sizeof(node_t), MPI_BYTE, minimum_reduction, 0, MPI_COMM_WORLD);


    if (rank == 0) { // Process 0 prints the results
        time = MPI_Wtime() - time;

        printf("Response Time: %.02lfs\n", time);
        printf("Toal Cost: %d\n\n", result->cost);
        print_path(result->path, nodes);
    }


    free(path);
    free(graph);
    free(result);
    free(minimum);

    MPI_Finalize();
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

// Source: https://stackoverflow.com/questions/7918806/finding-n-th-permutation-without-computing-others
static inline void get_permutation(node_t* permutation, size_t size, int_fast64_t index) {
    int_fast64_t factorial_ = factorial(size);

    for (size_t i = 0; i < size; i++) {
        factorial_ /= size - i;

        int_fast64_t j = index / factorial_;
        index %= factorial_;

        int value = permutation[j];

        memmove(&permutation[1], permutation, j * sizeof(node_t));
        *(permutation++) = value;
    }
}


// Graph

static inline edge_t* allocate_graph(size_t nodes) {
    edge_t* graph = malloc(nodes * nodes * sizeof(edge_t));

    if (graph == NULL) {
        fputs("Failed to allocate memory\n", stderr);
        exit(EXIT_FAILURE);
    }

    return graph;
}

static inline edge_t generate_edge() {
    return MIN_WEIGHT + rand() % (MAX_WEIGHT - MIN_WEIGHT + 1);
}

static inline void generate_graph(edge_t* graph, size_t nodes) {
    for (size_t i = 0; i < nodes; i++) {
        for (size_t j = i + 1; j < nodes; j++) {
            graph[i * nodes + j] = generate_edge();
            graph[j * nodes + i] = generate_edge();
        }

        graph[i * nodes + i] = 0;
    }
}


// Minumum

static inline minimum_t* allocate_minimum(size_t nodes) {
    minimum_t* minimum = malloc(sizeof(minimum_t) + (nodes + 1) * sizeof(node_t));

    if (minimum == NULL) {
        fputs("Failed to allocate memory\n", stderr);
        exit(EXIT_FAILURE);
    }

    minimum->cost = nodes * MAX_WEIGHT;
    return minimum;
}

void reduce_minimum(void* in, void* inout, int* length, MPI_Datatype* data_type) {
    (void)data_type; // Unused parameter

    if (((minimum_t*)in)->cost < ((minimum_t*)inout)->cost) {
        memcpy(inout, in, *length);
    }
}


// Path

static inline node_t* generate_path(size_t nodes, int_least64_t index) {
    node_t* path = malloc((nodes + 1) * sizeof(node_t));

    if (path == NULL) {
        fputs("Failed to allocate memory\n", stderr);
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < nodes; i++) {
        path[i] = i;
    }
    path[nodes] = 0; // Loop back

    get_permutation(&path[1], nodes - 1, index);
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
        printf("%d ", path[i]);
    }

    putchar('\n');
}
