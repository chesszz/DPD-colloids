#include "double_vector.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

Double_Vector *new_vec(void) {
    Double_Vector *vec = malloc(sizeof(Double_Vector));
    if (vec == NULL) {
        fprintf(stderr, "Failed to allocate memory.\n");
        exit(1);
    }

    vec->length = 0;
    vec->max_length = 16;

    double *new_items = malloc(vec->max_length * sizeof(double));
    if (new_items == NULL) {
        fprintf(stderr, "Failed to allocate memory.\n");
        exit(1);
    }

    vec->items = new_items;

    return vec;
}

void add_item(Double_Vector *vec, double item) {

    assert(vec != NULL);
    /* Check that the array is at most *barely* full (when equal). */
    assert(vec->max_length >= vec->length);

    /* If we have space in the array, then just add it at the end. */
    if (vec->max_length > vec->length) {
        vec->items[vec->length] = item;
        vec->length++;
    }

    /* Else, if we have reached the maximum value, we double the max size. */
    else if (vec->max_length == vec->length) {
        double *new_items = 
                    realloc(vec->items, 2 * vec->max_length * sizeof(double));
        if (new_items == NULL) {
            fprintf(stderr, "Failed to allocate memory.\n");
            exit(1);
        }

        /* Assign the new aray to the vector. */
        vec->items = new_items;
        vec->max_length *= 2;
        /* Check that now we have space. */
        assert(vec->max_length > vec->length);

        vec->items[vec->length] = item;
        vec->length++;
    }
    
    /* Will never reach here. */
    else {
        assert(0);
    }
}

/* Makes the vector think that it is empty, but max_length and the items are 
 * unchanged, so the memory is not reclaimed. This is to avoid having to call
 * realloc again later.
 */
void clear_vec(Double_Vector *vec) {
    vec->length = 0;
}
