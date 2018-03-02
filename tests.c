#include <stdio.h>
#include <stdlib.h>
// time
#include <time.h>
// rand
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
#define EPS 1e-6

int i_am_the_master;
int g_n; // для некоторых генераторов

void print_test (const char *test) {
	if (i_am_the_master) {
		if (test != NULL)
			printf ("%s...",test);
		else
			printf ("pass!\n");
	}
}
void assert (int expr, const char *msg) {
	if (!expr) {
		print_test (FAIL);
		if (msg != NULL)
			fprintf (stderr, "--- %s\n", msg);
		exit (ERROR);
	}
}
double vector_filler_natural (int row) {
	return row+1;
}
double matrix_filler_natural (int row, int col) {
	return row*g_n+col+1;
}
double matrix_filler_inverse_sum (int row, int col) {
	// row и col от 0
	return 1/(row+col+1+0.0);
}
double vector_filler_inverse_sum (int row) {
	double sum = 0;
	for (int col = 0; col < g_n; col++)
		if (col % 2 == 1)
			sum += matrix_filler_inverse_sum (row, col);
	return sum;
}
double vector_filler_inverse_sum_solution (int row) {
	return row % 2;
}
double rand_d () {
	return (rand () - RAND_MAX/2)/((double)RAND_MAX)*100;
}
double vector_filler_random (int row) {
	return rand_d();
}
double matrix_filler_random (int row, int col) {
	return rand_d();
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
double matrix_filler_slau (int row, int col) {
	double a [MEDIUM][MEDIUM] = {
		{ 45.763030, 44.269726, 31.264127, 24.050056, 82.452115 },
		{ 12.769607, 29.755301, 50.650619, 95.263367, 1.872111 },
		{ 73.776901, 45.171356, 40.387992, 31.608756, 53.029842 },
		{ 1.276430, 2.474284, 78.646368, 75.875657, 9.722007 },
		{ 82.509788, 58.273965, 55.112964, 9.725793, 56.041680 }
	};
	return a[row][col];
}
double vector_filler_slau (int row) {
	double a [MEDIUM] = {
		53.629906, 91.537775, 3.781911, 14.883654, 42.956326
	};
	return a[row];
}
double vector_m_v_mult_solution (int row) {
	double a [MEDIUM] = {
		15, 40, 65, 90, 115
	};
	return a[row];
}
double vector_filler_slau_solution (int row) {
	double a [MEDIUM] = {
		-2.295191,4.060002,0.270309,-0.137051,-0.318066
	};
	return a[row];
}
double matrix_filler_antidiagonal (int row, int col) {
	return row+col == g_n-1;
}
double vector_filler_countdown (int row) {
	return g_n - row;
}
double vector_filler_ones (int row) {
	return 1;
}
double matrix_transformer_triagonalizer (int i, int j, double elem) {
	if (i > j)
		return 0;
	else
		return elem;
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
	if (matrix_abs_diff (A, A_true) > EPS) {
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

void test_triagonalization () {
	print_test ("test_triagonalization");
	// ======= test code ==========

	if (i_am_the_master)
		printf ("\n\n");
	Matrix A = matrix_new_and_fill (MEDIUM, matrix_filler_slau);
	Vector V;
	if (i_am_the_master) {
		V = vector_new_and_fill (MEDIUM, vector_filler_slau);
	}
	
	matrix_triagonalize (A, V);
	matrix_print (A);

	Matrix A_copy = matrix_new_copy (A);
	matrix_transform (A_copy, matrix_transformer_triagonalizer);

	double diff = matrix_abs_diff (A, A_copy);
	if (i_am_the_master)
		assert (diff < EPS, "Матрица не треугольная!");

	matrix_free (A);
	if (i_am_the_master)
		vector_free (V);

	// ======= test code ==========
	print_test (OK);
}

void test_triagonalization_identity () {
	print_test ("test_triagonalization_identity");
	// ======= test code ==========

	if (i_am_the_master)
		printf ("\n\n");
	Matrix A = matrix_new_and_fill (MEDIUM, matrix_filler_identity);
	Matrix A_copy = matrix_new_copy (A);
	Vector V;
	if (i_am_the_master)
		V = vector_new (A.n_rows);
	matrix_triagonalize (A, V);

	matrix_print (A);

	double diff = matrix_abs_diff (A, A_copy);
	if (i_am_the_master)
		assert (diff < EPS, "Матрица не диагональная!");

	matrix_free (A);
	matrix_free (A_copy);
	if (i_am_the_master)
		vector_free (V);

	// ======= test code ==========
	print_test (OK);
}

void test_slau_solving () {
	print_test ("test_slau_solving");
	if (i_am_the_master)
		printf ("\n\n");
	// ======= test code ==========

	Matrix A = matrix_new_and_fill (MEDIUM, matrix_filler_slau);
	Vector V;
	Vector True_X;
	if (i_am_the_master) {
		V = vector_new_and_fill (MEDIUM, vector_filler_slau);
		True_X = vector_new_and_fill (MEDIUM, vector_filler_slau_solution);
	}
	
	matrix_triagonalize (A, V);
	Vector V_copy; // обратный ход испортит вектор V
	if (i_am_the_master)
		V_copy = vector_new_copy (V);
	matrix_print (A);
	vector_print (V);
	Vector X = gaussian_elimination (A, V);
	matrix_print (A);
	if (i_am_the_master) {
		vector_print (X);
		vector_print (V);
	}
	if (i_am_the_master) {
		assert (vector_abs_diff (X, True_X) < EPS, "Решения не совпадают!");
	}

	double diff = matrix_slau_difference (A, X, V_copy);
	if (i_am_the_master)
		printf("--- Невязка: %f\n\n", diff);
	assert (diff < EPS, "Невязка не нулевая");

	matrix_free (A);
	if (i_am_the_master) {
		vector_free (V);
		vector_free (V_copy);
		vector_free (True_X);
	}
	vector_free (X);

	// ======= test code ==========
	print_test (OK);
}

void test_slau_antidiagonal () {
	print_test ("test_slau_antidiagonal");
	// ======= test code ==========

	g_n = MEDIUM;
	Matrix A = matrix_new_and_fill (MEDIUM, matrix_filler_antidiagonal);
	Vector V;
	Vector True_X;

	if (i_am_the_master) {
		V = vector_new_and_fill (MEDIUM, vector_filler_natural);
		True_X = vector_new_and_fill (MEDIUM, vector_filler_countdown);
	}

	matrix_triagonalize (A, V);
	Vector V_copy;
	if (i_am_the_master)
		V_copy = vector_new_copy (V);
	Vector X = gaussian_elimination (A, V);
	if (i_am_the_master)
		vector_print (X);

	if (i_am_the_master) {
		assert (vector_abs_diff (X, True_X) < EPS, "Решения не совпадают!");
	}
	double diff = matrix_slau_difference (A, X, V_copy);
	assert (diff < EPS, "Невязка не нулевая");
	if (i_am_the_master)
		printf("--- Невязка: %f\n\n", diff);

	matrix_free (A);
	if (i_am_the_master) {
		vector_free (V);
		vector_free (V_copy);
		vector_free (True_X);
	}
	vector_free (X);

	// ======= test code ==========
	print_test (OK);
}

void test_slau_degenerate_case () {
	print_test ("test_slau_degenerate_case");
	// ======= test code ==========

	g_n = MEDIUM;
	if (i_am_the_master)
		printf ("\n\n");
	Matrix A = matrix_new_and_fill (MEDIUM, matrix_filler_natural);
	Vector V;
	if (i_am_the_master)
		V = vector_new_and_fill (MEDIUM, vector_filler_natural);
	
	matrix_triagonalize (A, V);
	matrix_print (A);
	Vector X = gaussian_elimination (A, V);
	assert (X.data == NULL, "Не освободили память!");
	assert (X.size == 0, "Не уменьшили число элементов!");

	matrix_free (A);
	if (i_am_the_master) {
		vector_free (V);
	}
	vector_free (X);

	// ======= test code ==========
	print_test (OK);
}

void test_slau_size (int size) {
	print_test ("test_slau_size");
	if (i_am_the_master)
		printf ("%d...", size);
	// ======= test code ==========

	Matrix A = matrix_new_and_fill (size, matrix_filler_random);
	Vector V;
	if (i_am_the_master)
		V = vector_new_and_fill (size, vector_filler_random);
	matrix_triagonalize (A, V);
	Vector V_copy;
	if (i_am_the_master)
		V_copy = vector_new_copy (V);
	Vector X = gaussian_elimination (A, V);
	if (X.data != NULL) {
		double diff = matrix_slau_difference (A,X,V_copy);
		assert (diff < EPS, "Невязка не нулевая");
		// if (i_am_the_master)
		// 	printf("--- Невязка: %f\n\n", diff);
	}

	matrix_free (A);
	if (i_am_the_master) {
		vector_free (V);
		vector_free (V_copy);
	}
	vector_free (X);

	// ======= test code ==========
	print_test (OK);
}

void test_matrix_vector_mult () {
	print_test ("test_matrix_vector_mult");
	if (i_am_the_master)
		printf ("\n\n");
	// ======= test code ==========

	g_n = MEDIUM;
	Matrix A = matrix_new_and_fill (MEDIUM, matrix_filler_natural);
	Vector V = vector_new_and_fill (MEDIUM, vector_filler_ones);
	matrix_print (A);
	Vector Res = matrix_vector_mult (A, V);
	Vector True_Res;
	if (i_am_the_master) {
		True_Res = vector_new_and_fill (MEDIUM, vector_m_v_mult_solution);
		vector_print (Res);
		assert (vector_abs_diff (Res, True_Res) < EPS, NULL);
		vector_free (V);
		vector_free (Res);
		vector_free (True_Res);
	}
	matrix_free (A);

	// ======= test code ==========
	print_test (OK);
}

void test_slau_diff () {
	print_test ("test_slau_diff");
	// if (i_am_the_master)
	// 	printf ("\n\n");
	// ======= test code ==========
	g_n = MEDIUM;
	Matrix A = matrix_new_and_fill (MEDIUM, matrix_filler_natural);
	Vector X = vector_new_and_fill (MEDIUM, vector_filler_ones);
	Vector V;
	if (i_am_the_master) {
		V = vector_new_and_fill (MEDIUM, vector_m_v_mult_solution);
	}
	double diff = matrix_slau_difference (A, X, V);
	if (i_am_the_master)
		printf ("\n\n--- Невязка: %lf\n", diff);
	assert (diff < EPS, "matrix_slau_difference failed");

	matrix_free (A);
	vector_free (X);
	if (i_am_the_master)
		vector_free (V);
	// ======= test code ==========
	print_test (OK);
}


void test_slau_from_formula (int n) {
	print_test ("test_slau_from_formula");
	// if (i_am_the_master)
	// 	printf ("\n\n");
	// ======= test code ==========

	g_n = n;
	Matrix A = matrix_new_and_fill (n, matrix_filler_inverse_sum);
	Vector V = vector_new_and_fill (n, vector_filler_inverse_sum);
	Vector True_X = vector_new_and_fill (n, vector_filler_inverse_sum_solution);
	// matrix_print (A);
	// vector_print (V);
	matrix_triagonalize (A, V);
	Vector X = gaussian_elimination (A, V);
	if (X.data != NULL) {
		double diff = vector_abs_diff (X, True_X);
		assert (diff < EPS, "Невязка не нулевая");
		// vector_print (X);
		vector_free (X);
	}

	matrix_free (A);
	vector_free (V);
	vector_free (True_X);

	// ======= test code ==========
	print_test (OK);
}


int main(int argc, char **argv) {

	MPI_Init(&argc, &argv);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	srand(time(NULL) + 3600*rank);

	i_am_the_master = rank == MASTER;

	matrix_init();
	// =================== Tests ======================

	test_create_identity ();
	test_abs_difference ();
	test_matrix_print ();
	test_matrix_transform ();
	test_create_and_fill_vector ();
	test_triagonalization ();
	test_triagonalization_identity ();
	test_slau_solving ();
	test_slau_antidiagonal ();
	test_slau_degenerate_case ();
	for (int i = 0; i < 30; i++)
		test_slau_size (i);
	test_matrix_vector_mult ();
	test_slau_diff ();
	for (int i = 0; i < 7; i++)
		test_slau_from_formula (i);

	// =================== Tests ======================
	matrix_exit();

	MPI_Finalize();
	return 0;
}
