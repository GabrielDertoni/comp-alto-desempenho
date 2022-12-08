#ifndef _NEXT_PERM_H
#define _NEXT_PERM_H

// Generic swap
#define SWAP(type, a, b)      \
    do {                      \
        type *__ptr_a = a;    \
        type *__ptr_b = b;    \
        type tmp = *__ptr_a;  \
        *__ptr_a = *__ptr_b;  \
        *__ptr_b = tmp;       \
    } while(0)

// Reverse `upper_bound`. Returns a pointer to the first value (right to left) larger than `value`.
// Assumes `values` is sorted in non-increasing order.
//
// reference: https://en.cppreference.com/w/cpp/algorithm/upper_bound
int *upper_bound_r(int *values, int size, int value) {
    int *last = values + size - 1;
    while (size > 0) {
        int step = size / 2;
        if (value >= last[-step]) {
            last -= step + 1;
            size -= step + 1;
        } else {
            size = step;
        }
    }

    return last;
}

// Reverses the range delimited by [begin, end), i.e. `begin` points to the first element of the
// array and `end` is the past-the-end pointer of the array.
void reverse(int *begin, int *end) {
    while (begin < end)
        SWAP(int, begin++, --end);
}

// Get the next lexicographical permutation, i.e. the smallest lexicographical permutation of
// `values` that is bigger than the permutation in `values` before the function call. O(n)
void next_permutation(int *values, int size) {
    if (size <= 1) return;

    // Find a non-increasing sorted slice in the end of `values`. `sorted_start` indicates the
    // start of this slice that ends in the end of `values`.
    int *sorted_start = values + size - 1;
    while (sorted_start > values && sorted_start[-1] >= sorted_start[0])
        sorted_start--;

    if (sorted_start != values) {
        // Find the value to swap with `sorted_start[-1]`. It has to be bigger than it in order to
        // increase the lexicographical order.
        int *bound = upper_bound_r(sorted_start, values + size - sorted_start, sorted_start[-1]);
        SWAP(int, bound, sorted_start - 1);
    }

    // We know the slice from `sorted_start` to the end of `values` is still sorted in
    // lexicographical non-increasing order. But now we want to make it the least possible
    // lexicographical order (non-decreasing order). Just reverse it.
    reverse(sorted_start, values + size);
}

// Get the next lexicographical permutation, i.e. the smallest lexicographical permutation of
// `values` that is bigger than the permutation in `values` before the function call. O(n)
// It is undefined behaviour to call this function whith a sorted array.
void next_permutation_opt(int *values, int size) {
    if (size <= 1) return;

    // Find a non-increasing sorted slice in the end of `values`. `sorted_start` indicates the
    // start of this slice that ends in the end of `values`.
    int *sorted_start = values + size - 1;
    while (sorted_start[-1] >= sorted_start[0])
        sorted_start--;

    // Find the value to swap with `sorted_start[-1]`. It has to be bigger than it in order to
    // increase the lexicographical order.
    int *bound = upper_bound_r(sorted_start, values + size - sorted_start, sorted_start[-1]);
    SWAP(int, bound, sorted_start - 1);

    // We know the slice from `sorted_start` to the end of `values` is still sorted in
    // lexicographical non-increasing order. But now we want to make it the least possible
    // lexicographical order (non-decreasing order). Just reverse it.
    reverse(sorted_start, values + size);
}

/*
void permutation_iter(int *values, int size, int **state) {
    if (size <= 1) return;

    int *cycles;
    if (*state == NULL) {
        *state = cycles = (int*)malloc(size * sizeof(int));
        for (int i = 0; i < size; i++)
            cycles[i] = size - i;
        return;
    }
    cycles = *state;

    for (int i = size - 1; i >= 0; i--) {
        cycles[i]--;
        int j = cycles[i];
        if (j == 0) {
            int tmp = values[i];
            memmove(values + i, values + i + 1, (size - i - 1) * sizeof(int));
            values[size - 1] = tmp;
            cycles[i] = size - i;
        } else {
            SWAP(int, &values[i], &values[size - j]);
            return;
        }
    }

    // If the we got here, we reached the very last permutation
    free(*state);
    *state = NULL;
}
*/

#undef SWAP

#endif
