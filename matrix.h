#ifndef _MY_MATRIX_H
#define _MY_MATRIX_H

typedef struct Matrix {
	double *data;
	// global dimensions
	int n_rows;
	int n_cols;
	// local dimensions
	int local_rows;
	int local_cols;
} Matrix;

typedef struct GridProcess {
	// global info
	int n_proc_rows; // equal to nproc
	int n_proc_cols; // always 1
	int nproc;
	// local info
	int myrow; // always 1
	int mycol; // equal to myrank
	int myrank;
} GridProcess;

typedef struct {
	double *data;
	int size;
} Vector;

Vector vector_new (int size);
void vector_free (Vector);
Vector vector_new_and_fill (int size, double (*func)(int));
void vector_print (Vector);
double vector_abs_diff (Vector U, Vector V);


void matrix_init (const int block_size); // call before all others
void matrix_exit (); // call after all others

// matrix creation
Matrix matrix_new(const int dim);
Matrix matrix_new_copy (Matrix A);

void matrix_free (Matrix A);

// matrix initialization
void matrix_fill (Matrix A, double (*func)(int, int));

// matrix fillers
double matrix_filler_identity (int m, int n);

// matrix output
void matrix_print (Matrix A);
// void matrix_print_diag (Matrix A);

// matrix operations
double matrix_abs_diff (Matrix A, Matrix B);

// matrix transformations
void matrix_transform (Matrix A, double (*func)(int, int, double));


#endif