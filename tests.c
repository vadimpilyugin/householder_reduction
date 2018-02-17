#include <stdio.h>
#include <stdlib.h>
// Libraries
#include "mpi.h"
#include "matrix.h"

#define HUGE 23
#define BIG 8
#define MEDIUM 5
#define SMALL 3
#define OK NULL
#define ERROR 2
#define FAIL "FAIL!"
#define MASTER 0

#define BLOCK_SIZE 2

int i_am_the_master;

void print_test (const char *test) {
	if (i_am_the_master) {
		if (test != NULL)
			printf ("%s...",test);
		else
			printf ("pass!\n");
	}
}

double vector_filler_natural (int row) {
	return row+1;
}
double matrix_filler_hermitian (int row, int col) {
	double a [MEDIUM][MEDIUM] = {
		{1,5,9,5,8},
		{5,6,1,4,3},
		{9,1,2,3,4},
		{5,4,3,5,6},
		{8,3,4,6,9}
	};
	return a[row][col];
}
double matrix_filler_to_print (int row, int col) {
	double a [MEDIUM][MEDIUM] = {
		{1,2,3,4,5},
		{6,7,8,9,10},
		{11,12,13,14,15},
		{16,17,18,19,20},
		{21,22,23,24,25}
	};
	return a[row][col];
}

double matrix_filler_to_multiply_a (int row, int col) {
	double a[SMALL][SMALL] = {
		{1,2,3},
		{4,5,6},
		{7,8,9}
	};
	return a[row][col];
}
double matrix_filler_to_multiply_b (int row, int col) {
	double a[SMALL][SMALL] = {
		{5,0,0},
		{0,2,0},
		{0,0,3}
	};
	return a[row][col];
}
double matrix_filler_to_multiply_c (int row, int col) {
	double a[SMALL][SMALL] = {
		{5,4,9},
		{20,10,18},
		{35,16,27}
	};
	return a[row][col];
}

double matrix_filler_to_abs_diff_a (int row, int col) {
	double a[SMALL][SMALL] = {
		{1.01,1,1},
		{1,1,1.001},
		{1,1,1}
	};
	return a[row][col];
}
double matrix_filler_to_abs_diff_b (int row, int col) {
	double b[SMALL][SMALL] = {
		{1,1,1},
		{1,1,1},
		{1,1,1}
	};
	return b[row][col];
}

double matrix_transformer_inverse (int i, int j, double elem) {
	return 1-elem;
}

double matrix_filler_inverse (int row, int col) {
	double a[MEDIUM][MEDIUM] = {
		{0,1,1,1,1},
		{1,0,1,1,1},
		{1,1,0,1,1},
		{1,1,1,0,1},
		{1,1,1,1,0}
	};
	return a[row][col];
}

double matrix_filler_print_big (int row, int col) {
	return row*BIG+col+1;
}
double matrix_filler_huge (int row, int col) {
	// от 1 до HUGE^2
	return row*HUGE+col+1;
}

void test_create_identity () {
	print_test ("test_create_identity");
	// ======= test code ==========

	Matrix A = matrix_new (MEDIUM);
	matrix_fill (A, matrix_filler_identity);
	matrix_free (A);
	
	// ======= test code ==========
	print_test (OK);
}

void test_matrix_print () {
	print_test ("test_matrix_print");
	// ======= test code ==========

	Matrix A = matrix_new (BIG);
	matrix_fill (A, matrix_filler_print_big);
	if (i_am_the_master)
		printf ("\n\n");
	matrix_print (A);
	matrix_free (A);

	// Matrix A = matrix_new (HUGE);
	// matrix_fill (A, matrix_filler_huge);
	// matrix_print (A);
	// matrix_free (A);

	// ======= test code ==========
	print_test (OK);
}

void test_matrix_transform () {
	print_test ("test_matrix_transform");
	// ======= test code ==========

	Matrix A = matrix_new (MEDIUM);
	matrix_fill (A, matrix_filler_identity);
	matrix_transform (A, matrix_transformer_inverse);
	Matrix A_true = matrix_new (MEDIUM);
	matrix_fill (A_true, matrix_filler_inverse);
	if (matrix_abs_diff (A, A_true) > 1e-5) {
		print_test (FAIL);
		exit (ERROR);
	}
	if (i_am_the_master)
		printf ("\n\n");
	matrix_print (A);
	matrix_free (A);
	matrix_free (A_true);
	
	// ======= test code ==========
	print_test (OK);
}

void test_abs_difference () {
	print_test ("test_abs_difference");
	// ======= test code ==========

	Matrix A = matrix_new (BIG);
	matrix_fill (A, matrix_filler_print_big);
	Matrix B = matrix_new_copy (A);
	
	double diff = matrix_abs_diff (A,B);
	if (diff > 1) {
		print_test (FAIL);
		exit (ERROR);
	}
	matrix_free(A);
	matrix_free(B);
	
	// ======= test code ==========
	print_test (OK);
}

void test_create_and_fill_vector () {
	print_test ("test_create_and_fill_vector");
	// ======= test code ==========

	Vector V = vector_new_and_fill (MEDIUM, vector_filler_natural);
	vector_print (V);
	vector_free (V);
	
	// ======= test code ==========
	print_test (OK);
}

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	i_am_the_master = rank == MASTER;

	matrix_init(BLOCK_SIZE);
	// =================== Tests ======================

	test_create_identity ();
	test_abs_difference ();
	test_matrix_print ();
	test_matrix_transform ();
	test_create_and_fill_vector ();

	// =================== Tests ======================
	matrix_exit();

	MPI_Finalize();
	return 0;
}
