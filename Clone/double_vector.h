#ifndef DOUBLE_VECTOR_H
#define DOUBLE_VECTOR_H

typedef struct {
    int length;
    int max_length;
    double *items;
} Double_Vector;

Double_Vector *new_vec(void);
void add_item(Double_Vector *vec, double item);
void clear_vec(Double_Vector *vec);

#endif /* DOUBLE_VECTOR_H */