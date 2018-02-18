// cpp header
#include "matrix.h"

// Standard headers
#include <stdio.h>
// memcpy
#include <string.h>
// errno
#include <errno.h>
// sqrt
#include <math.h>
// exit
#include <stdlib.h>

// Libs
#include "mpi.h"

// Definitions
#define ERROR 2
#define MASTER 0
#define NO_PARAM -1
#define GET_DEFAULT_CONTEXT 0
#define DEBUG 0
#define DESC_N_ELEM 9
#define USE_LOCAL_STORAGE 0
#define STR_LEN 50
#define BLOCK_SIZE 1

#define SUCCESS 0
#define ALREADY_TRIAGONAL 1
#define ALREADY_ROTATED 2
#define ZERO_COEFFICIENT 3

#define EPS 1e-7

// Macros
#define INDEX(row, col, size) (indxg2lr((row))+(size)*indxg2lc((col)))
#define sqr(x) ((x)*(x))

// Global variables start with 'g_'
static int g_block_size = 0;
static GridProcess g_grid_proc; // информация о процессе внутри грида
static int i_am_the_master;
double g_alpha; // нужна для matrix_transform_exp
static char g_str [STR_LEN];

 // INDXG2L = NB*((INDXGLOB-1)/(NB*NPROCS))+MOD(INDXGLOB-1,NB)+1

int indxg2l (int global_ind, int block_size, int n_proc_dim) {
	return block_size*(global_ind/(block_size*n_proc_dim))+global_ind%block_size;
}
// global index to local row
int indxg2lr (int global_row) {
	if (indxg2l (global_row,g_block_size,g_grid_proc.n_proc_rows) != (g_block_size*(global_row/(g_block_size*g_grid_proc.n_proc_rows))+global_row%g_block_size)) {
		perror ("Не равны");
		exit (ERROR);
	}

	return indxg2l (global_row,g_block_size,g_grid_proc.n_proc_rows);
}
// global index to local column
int indxg2lc (int global_col) {
	return indxg2l (global_col, g_block_size, g_grid_proc.n_proc_cols);
}

 // INDXG2P = MOD( ISRCPROC + (INDXGLOB - 1) / NB, NPROCS )

// global index to local process number
int indxg2p (int global_ind, int block_size, int n_proc_dim) {
	return (global_ind / block_size) % n_proc_dim;
}
int indxg2pr (int global_row) {
	return indxg2p (global_row, g_block_size, g_grid_proc.n_proc_rows);
}
int indxg2pc (int global_col) {
	return indxg2p (global_col, g_block_size, g_grid_proc.n_proc_cols);
}

// INDXL2G = NPROCS*NB*((INDXLOC-1)/NB) + MOD(INDXLOC-1,NB) + MOD(NPROCS+IPROC-ISRCPROC, NPROCS)*NB + 1

// local index to global
int indxl2g (int local_ind, int mydim, int block_size, int n_proc_dim) {
	return n_proc_dim * block_size * (local_ind / block_size) + 
		local_ind % block_size + mydim * block_size;
}
int indxl2gr (int local_row) {
	return indxl2g (local_row, g_grid_proc.myrow, g_block_size, g_grid_proc.n_proc_rows);
}
int indxl2gc (int local_col) {
	return indxl2g (local_col, g_grid_proc.mycol, g_block_size, g_grid_proc.n_proc_cols);
}

// число строк или столбцов в локальной матрице
int numroc (int dim, int block_size, int proc_ind, int n_proc_dim) {
	int nblocks = dim / block_size;
	int result = (nblocks/n_proc_dim) * block_size;
	int extra_blocks = nblocks % n_proc_dim;
	if (proc_ind < extra_blocks)
		result += block_size;
	else if (proc_ind == extra_blocks)
		result += dim % block_size;
	return result;
}

// печать строки, находящейся в памяти процесса sender, на процессе receiver
void print_g_str (int sender, int receiver) {
	int tag = 0;
	if (g_grid_proc.myrank == sender) {
		MPI_Request request;

		MPI_Isend (g_str,STR_LEN,MPI_CHAR,receiver,tag,MPI_COMM_WORLD,&request);
	}
	if (g_grid_proc.myrank == receiver) {
		MPI_Status status;

		MPI_Recv (g_str,STR_LEN,MPI_CHAR,sender,tag,MPI_COMM_WORLD,&status);
		printf ("%s",g_str);
		// fflush (stdout);
	}
}

static GridProcess new_grid (int n_proc_rows, int n_proc_cols, int myrank, int nproc) {

	GridProcess grid_proc;
	grid_proc.n_proc_rows = n_proc_rows;
	grid_proc.n_proc_cols = n_proc_cols;
	grid_proc.myrow = 0;
	grid_proc.mycol = myrank;

	grid_proc.nproc = nproc;
	grid_proc.myrank = myrank;
	return grid_proc;
}

void matrix_init () {

	int myrank, nproc;
	MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
	MPI_Comm_size (MPI_COMM_WORLD, &nproc);
	i_am_the_master = myrank == MASTER;

	g_block_size = BLOCK_SIZE;
	int n_proc_rows = 1;
	int n_proc_cols = nproc;
	// guess_process_grid_dimensions (&n_proc_rows, &n_proc_cols);
	g_grid_proc = new_grid (n_proc_rows, n_proc_cols, myrank, nproc);
	#if DEBUG
	if (i_am_the_master) {

		printf ("\n\tПроцессорная решетка\n\n");
		printf ("\t┌---------┐\n");
		printf ("\t|    |    |\n");
		printf ("\t|    |    |\n");
		printf ("\t|----|----|  %d rows\n", g_grid_proc.n_proc_rows);
		printf ("\t|    |    |\n");
		printf ("\t|    |    |\n");
		printf ("\t└---------┘\n");
		printf ("\t   %d cols\n\n", g_grid_proc.n_proc_cols);
	}

	snprintf (g_str, STR_LEN, "--- Процесс %d/%d -> /%d,%d/\n", 
		myrank,nproc,g_grid_proc.myrow,g_grid_proc.mycol);
	int proc_num;
	for (proc_num = 0; proc_num < nproc; proc_num++) {
		MPI_Barrier (MPI_COMM_WORLD);

		print_g_str (proc_num, MASTER);

		MPI_Barrier (MPI_COMM_WORLD);
	}
	if (i_am_the_master)
		printf ("\n\n");
	#endif
}

void matrix_exit () {

}

Matrix matrix_new(const int dim_) {
	Matrix A;
	int dim = dim_;

	A.n_rows = dim;
	A.n_cols = dim;

	A.local_rows = numroc (
		dim,
		g_block_size,
		g_grid_proc.myrow,
		g_grid_proc.n_proc_rows
	);
	A.local_cols = numroc (
		dim,
		g_block_size,
		g_grid_proc.mycol,
		g_grid_proc.n_proc_cols
	);

	A.data = (double *) malloc (A.local_rows*A.local_cols*sizeof (double));

	#if DEBUG
		if (i_am_the_master) {
			printf ("\n\n\tРаспределение матрицы\n\n");
			printf ("\t┌---------┐\n");
			printf ("\t|    |    |\n");
			printf ("\t|    |    |\n");
			printf ("\t|----|----|  %d rows\n", A.n_rows);
			printf ("\t|    |    |\n");
			printf ("\t|    |    |\n");
			printf ("\t└---------┘\n");
			printf ("\t   %d cols\n\n", A.n_cols);
			printf ("--- Размер блока: %dx%d\n\n", g_block_size, g_block_size);
		}

		if (i_am_the_master)
			printf ("--- Процессорная решетка: %dx%d\n\n", 
				g_grid_proc.n_proc_rows, g_grid_proc.n_proc_cols
			);

		snprintf (g_str, STR_LEN, "(%d,%d)  ", g_grid_proc.myrow, g_grid_proc.mycol);

		int grid_row, grid_col;
		for (grid_row = 0; grid_row < g_grid_proc.n_proc_rows; grid_row++) {
			for (grid_col = 0; grid_col < g_grid_proc.n_proc_cols; grid_col++) {
				MPI_Barrier (MPI_COMM_WORLD);

				if (i_am_the_master && grid_col == 0)
					printf ("  ");
				print_g_str (grid_row*g_grid_proc.n_proc_cols+grid_col, MASTER);
				if (i_am_the_master && grid_col == g_grid_proc.n_proc_cols-1)
					printf ("\n");

				MPI_Barrier (MPI_COMM_WORLD);
			}
		}
		if (i_am_the_master) {
			printf ("\n\n");
			printf ("--- Распределение матрицы по процессам: \n\n");
		}

		snprintf (g_str, STR_LEN, "{%3d,%3d} | ", A.local_rows, A.local_cols);

		for (grid_row = 0; grid_row < g_grid_proc.n_proc_rows; grid_row++) {
			for (grid_col = 0; grid_col < g_grid_proc.n_proc_cols; grid_col++) {
				MPI_Barrier (MPI_COMM_WORLD);

				if (i_am_the_master && grid_col == 0)
					printf ("  | ");
				print_g_str (grid_row * g_grid_proc.n_proc_cols + grid_col, MASTER);
				if (i_am_the_master && grid_col == g_grid_proc.n_proc_cols-1)
					printf ("\n");

				MPI_Barrier (MPI_COMM_WORLD);
			}
		}

		if (i_am_the_master)
			printf ("\n\n");

		MPI_Barrier (MPI_COMM_WORLD);

	#endif

	return A;
}

void matrix_free (Matrix A) {
	if (A.data != NULL) {
		free (A.data);
	}
}

Matrix matrix_new_from_vector (double * A_vec, int dim) {
	Matrix A = matrix_new (dim);
	int row, col;
	for (col = 0; col < A.local_cols; col++) {
		for (row = 0; row < A.local_rows; row++) {
			A.data[row + A.local_rows*col] = A_vec[(indxl2gr(row))*dim+indxl2gc(col)];
		}
	}
	return A;
}

Matrix matrix_new_diag_from_vector (double * A_vec, int dim) {
	Matrix A = matrix_new (dim);
	int row, col;
	for (col = 0; col < A.local_cols; col++) {
		for (row = 0; row < A.local_rows; row++) {
			if (indxl2gr(row) == indxl2gc(col))
				A.data[row + A.local_rows*col] = A_vec[indxl2gr(row)];
			else
				A.data[row + A.local_rows*col] = 0;
		}
	}
	return A;
}

Matrix matrix_new_copy (Matrix A) {
	if (A.n_rows != A.n_cols) {
		errno = EINVAL;
		perror ("Матрица должна быть квадратной!");
		exit (ERROR);
	}
	Matrix A_copy = matrix_new (A.n_rows);
	memcpy (A_copy.data, A.data, A.local_rows*A.local_cols*sizeof(double));
	return A_copy;
}

void matrix_fill (Matrix A, double (*func)(int, int)) {
	int row, col;
	for (col = 0; col < A.local_cols; col++) {
		for (row = 0; row < A.local_rows; row++) {
			A.data[row + A.local_rows*col] = func (indxl2gr(row), indxl2gc(col));
		}
	}
}

Matrix matrix_new_and_fill (int size, double (*func)(int, int)) {
	Matrix A = matrix_new (size);
	matrix_fill (A, func);
	return A;
}

void matrix_transform (Matrix A, double (*func)(int, int, double)) {
	int row, col;
	for (col = 0; col < A.local_cols; col++) {
		for (row = 0; row < A.local_rows; row++) {
			A.data[row + A.local_rows*col] = func (indxl2gr(row), indxl2gc(col), A.data[row + A.local_rows*col]);
		}
	}
}

double matrix_filler_identity (int m, int n) {
	if (m == n) 
		return 1;
	else
		return 0;
}

void matrix_print (Matrix A) {
	MPI_Barrier (MPI_COMM_WORLD);

	int global_row;
	int global_col;

	#if DEBUG
	if (i_am_the_master) {
		printf ("\n\n");
		printf ("\n============ Matrix %dx%d =============\n",A.n_rows,A.n_cols);
	}

	if (i_am_the_master)
		printf ("\n\n--- Распределение по процессам\n\n");

	for (global_row = 0; global_row < A.n_rows; global_row++) {
		for (global_col = 0; global_col < A.n_cols; global_col++) {
			MPI_Barrier (MPI_COMM_WORLD);

			int process_row = (global_row/g_block_size)%g_grid_proc.n_proc_rows;
			int process_col = (global_col/g_block_size)%g_grid_proc.n_proc_cols;

			if (i_am_the_master){
				printf ("%4d",process_row*g_grid_proc.n_proc_cols+process_col);
				if (global_col == A.n_cols-1)
					printf ("\n");
				fflush (stdout);
			}

			MPI_Barrier (MPI_COMM_WORLD);
		}
	}

	if (i_am_the_master)
		printf ("\n\n--- Нумерация внутри процессов\n\n");

	for (global_row = 0; global_row < A.n_rows; global_row++) {
		for (global_col = 0; global_col < A.n_cols; global_col++) {
			MPI_Barrier (MPI_COMM_WORLD);

			int process_row = (global_row/g_block_size)%g_grid_proc.n_proc_rows;
			int process_col = (global_col/g_block_size)%g_grid_proc.n_proc_cols;
			
			if (process_row == g_grid_proc.myrow && process_col == g_grid_proc.mycol) {
				{
					int local_row = indxg2lr (global_row);
					int local_col = indxg2lc (global_col);
					snprintf (g_str, STR_LEN, "%d,%d  ", local_row,local_col);
				}
				MPI_Request request;
				MPI_Isend (g_str,STR_LEN,MPI_CHAR,MASTER,0,MPI_COMM_WORLD, &request);
			}
			if (i_am_the_master) {
				int source = process_row*g_grid_proc.n_proc_cols+process_col;
				MPI_Status status;

				MPI_Recv (g_str,STR_LEN,MPI_CHAR,source,0,MPI_COMM_WORLD,&status);
				printf ("%s",g_str);
				fflush (stdout);
				if (global_col == A.n_cols-1)
					printf ("\n");
			}
			MPI_Barrier (MPI_COMM_WORLD);
		}
	}

	if (i_am_the_master)
		printf ("\n=============================\n\n");


	#endif

	for (global_row = 0; global_row < A.n_rows; global_row++) {
		for (global_col = 0; global_col < A.n_cols; global_col++) {
			MPI_Barrier (MPI_COMM_WORLD);

			int process_row = (global_row/g_block_size)%g_grid_proc.n_proc_rows;
			int process_col = (global_col/g_block_size)%g_grid_proc.n_proc_cols;
			
			if (process_row == g_grid_proc.myrow && process_col == g_grid_proc.mycol) {
				{
					int local_row = indxg2lr (global_row);
					int local_col = indxg2lc (global_col);
					double local_elem = A.data[local_col*A.local_rows+local_row];
					snprintf (g_str, STR_LEN, "%8.4lf ", local_elem);
				}
				MPI_Request request;
				MPI_Isend (g_str,STR_LEN,MPI_CHAR,MASTER,0,MPI_COMM_WORLD, &request);
			}
			if (i_am_the_master) {
				int source = process_row*g_grid_proc.n_proc_cols+process_col;
				MPI_Status status;

				MPI_Recv (g_str,STR_LEN,MPI_CHAR,source,0,MPI_COMM_WORLD,&status);
				printf ("%s",g_str);
				fflush (stdout);
				if (global_col == A.n_cols-1)
					printf ("\n");
			}


			MPI_Barrier (MPI_COMM_WORLD);
		}
	}
	#if DEBUG
	if (i_am_the_master) {
		printf ("\n=====================================\n");
		printf ("\n\n");
	}
	#endif

	MPI_Barrier (MPI_COMM_WORLD);
}

double matrix_abs_diff (Matrix A, Matrix B) {
	if (A.n_rows != B.n_rows || A.n_cols != B.n_cols) {
		errno = EINVAL;
		perror ("abs_diff: Матрицы не совпадают по размерности!");
		exit (ERROR);
	}
	if (A.local_rows != B.local_rows || A.local_cols != B.local_cols) {
		errno = EINVAL;
		perror ("abs_diff: Матрицы не совпадают по распределению!");
		exit (ERROR);
	}
	int rows = A.local_rows;
	int cols = A.local_cols;
	double diff = 0.0;
	int row;
	int col;
	// считаем разницу по локальной части
	for (row = 0; row < rows; row++) {
		for (col = 0; col < cols; col++) {
			double abs_val = abs (
				A.data[row*cols+col]-B.data[row*cols+col]
			);
			if (abs_val > diff)
				diff = abs_val;
		}
	}
	// находим максимум и раздаем его всем
	double total_diff = 0;
	MPI_Reduce(&diff, &total_diff, 1, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);
  MPI_Bcast (&total_diff, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	return total_diff;
}

Vector vector_new (int size) {
	Vector V;
	V.data = (double *) malloc (size*sizeof(double));
	if (V.data == NULL) {
		perror ("vector_new");
		fprintf (stderr, "Значение size: %d\n", size);
		exit (ERROR);
	}
	V.size = size;
	return V;
}

Vector vector_new_copy (Vector V) {
	Vector V_copy = vector_new (V.size);
	memcpy (V_copy.data, V.data, V.size*sizeof(double));
	V_copy.size = V.size;
	return V_copy;
}

void vector_free (Vector V) {
	if (V.data != NULL)
		free (V.data);
}

void vector_fill (Vector V, double (*func)(int)) {
	int row;

	for (row = 0; row < V.size; row++) {
		V.data[row] = func (row);
	}
}

Vector vector_new_and_fill (int size, double (*func)(int)) {
	Vector V = vector_new (size);
	vector_fill (V, func);
	return V;
}

void vector_print (Vector V) {
	if (i_am_the_master) {
		int row;

		printf("\n\n");
		for (row = 0; row < V.size; row++) {
			printf ("%8.4lf\n", V.data[row]);
		}
		printf("\n\n");
	}
}

double vector_abs_diff (Vector U, Vector V) {
	if (U.size != V.size) {
		errno = EINVAL;
		perror ("vector_abs_diff: размеры не совпадают!");
		exit (ERROR);
	}

	int size = V.size;
	int row;
	double diff = 0;

	for (row = 0; row < size; row++) {
		int index = row;
		double tmp = fabs (U.data[index]-V.data[index]);
		if (tmp > diff)
			diff = tmp;
	}

	return diff;
}


int calculate_x (int n, int k, Vector X, Matrix A) {
	int row;
	int column_start = INDEX(0,k,n);
	// считаем вектор x матрицы отражения
	// x = (a - ||a||*e) / ||a - ||a||*e||

	// S_k = sum for i = k+1..n (A[i,k]^2)
	double s_k = 0;
	for (row = k+1; row < n; row++)
		s_k += sqr (A.data[column_start+row]);

	// ||a|| = sqrt (S_k + A[k,k]^2)
	double norm_a = sqrt (s_k + sqr(A.data[column_start+k]));
	if (norm_a < EPS) {
		// матрица приведена к треугольному виду
		return ALREADY_TRIAGONAL;
	}

	// X[k] = A[k,k]-norm(A[k])
	X.data[k] = A.data[column_start+k] - norm_a;

	// X[k+1..n] = A[k+1..n,k]
	for (row = k+1; row < n; row++)
		X.data[row] = A.data[column_start+row];

	// ||X|| = sqrt (S_k + X[k]^2)
	double norm_x = sqrt (s_k + sqr (X.data[k]));
	if (norm_x < EPS) {
		// вектор уже повернут
		return ALREADY_ROTATED;
	}
	// X = X / norm_x
	for (row = k; row < n; row++) {
		X.data[row] /= norm_x;
	}
	return SUCCESS;
}

void matrix_triagonalize (Matrix A, Vector V) {
	if (i_am_the_master && A.n_rows != V.size) {
		errno = EINVAL;
		perror ("matrix_triagonalize: размеры не совпадают!");
		exit (ERROR);
	}
	int n = A.n_rows;
	int k;
	int row;
	int col;
	int code;
	Vector X = vector_new (n);

	for (k = 0; k < n-1; k++) {
		// считаем вектор x матрицы отражения
		// x = (a - ||a||*e) / ||a - ||a||*e||

		int rank_has_column = indxg2pc (k);

		if (g_grid_proc.mycol == rank_has_column) {

			// его считает только процесс, у которого есть этот столбец
			code = calculate_x (n, k, X, A);
		}

		// этот процесс рассылает код возврата остальным
		MPI_Bcast (&code, 1, MPI_INT, rank_has_column, MPI_COMM_WORLD);
		if (code == ALREADY_ROTATED)
			continue;
		else if (code == ALREADY_TRIAGONAL) {
			vector_free (X);
			return;
		}

		// далее этот процесс рассылает вектор X
		MPI_Bcast (X.data, X.size, MPI_DOUBLE, rank_has_column, MPI_COMM_WORLD);

		// считаем преобразование матрицы A
		// a = (E - 2xx*)a = a - (2x*a)x
		double dot_product;
		for (col = k; col < n; col++) {
			rank_has_column = indxg2pc (col);
			if (g_grid_proc.mycol == rank_has_column) {
				int column_start = INDEX(0,col,n);

				// dot_product = x*a
				dot_product = 0;
				for (row = k; row < n; row++) {
					dot_product += A.data[column_start+row]*X.data[row];
				}

				// 2x*a
				dot_product *= 2;

				// a = a - (2x*a)x
				for (row = k; row < n; row++) {
					int index = column_start+row;
					A.data[index] -= dot_product*X.data[row];
				}
			}
		}
		if (i_am_the_master) {
			// считаем преобразование вектора b
			// a = (E - 2xx*)a = a - (2x*a)x
			// dot_product = x*a
			dot_product = 0;
			for (row = k; row < n; row++) {
				dot_product += V.data[row]*X.data[row];
			}
			// 2x*a
			dot_product *= 2;
			// a = a - (2x*a)x
			for (row = k; row < n; row++) {
				V.data[row] -= dot_product*X.data[row];
			}	
		}
	}
	vector_free (X);
}

Vector gaussian_elimination (Matrix A, Vector V) {
	if (i_am_the_master && (
			A.n_rows != V.size ||
			A.n_rows != A.n_cols
		)
	) {
		errno = EINVAL;
		perror ("gaussian_elimination: размеры не совпадают!");
		exit (ERROR);
	}

	int n = A.n_rows;
	int k;
	int row;
	int code;
	Vector X, A_k;

	X = vector_new (n); // вектор решения
	if (i_am_the_master) {
		A_k = vector_new (n); // дополнительный для хранения столбцов матрицы
		memset (A_k.data, 0, n*sizeof(double));
	}

	for (k = n-1; k >= 0; k--) {
		MPI_Barrier (MPI_COMM_WORLD);

		// находим процесс, у которого есть k-ый столбец матрицы
		int process_col = indxg2pc (k);
		if (g_grid_proc.mycol == process_col) {

			// посылаем на MASTER k-ый столбец матрицы с 0 по k-ый элементы включительно
			MPI_Request request;
			MPI_Isend (
				A.data+INDEX(0,k,n), // первый элемент в k-ом столбце
				k+1, 
				MPI_DOUBLE, 
				MASTER, 
				0, MPI_COMM_WORLD, 
				&request
			);
		}

		if (i_am_the_master) {
			// принимаем k-ый столбец матрицы
			MPI_Status status;
			MPI_Recv (A_k.data, k+1, MPI_DOUBLE, process_col, 0, MPI_COMM_WORLD, &status);

			// MASTER вычисляет X[k]
			if (fabs (A_k.data[k]) < EPS) {
				printf ("\n\n--- Матрица вырожденная\n\n");
				vector_free (A_k);
				code = ZERO_COEFFICIENT;
			} else {
				// находим очередное решение
				X.data[k] = V.data[k] / A_k.data[k];

				// и вычитаем из правой части k-ый столбец, умноженный на это решение
				for (row = 0; row <= k; row++)
					V.data[row] -= X.data[k] * A_k.data[row];

				code = SUCCESS;
			}
		}

		MPI_Bcast (&code, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
		if (code == ZERO_COEFFICIENT) {
			vector_free (X);
			X.data = NULL;
			X.size = 0;
			return X;
		}
		MPI_Barrier (MPI_COMM_WORLD);
	}
	MPI_Bcast (X.data, n, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	return X;
}

Vector matrix_vector_mult (Matrix A, Vector X) {
	if (i_am_the_master && (
			A.n_cols != X.size ||
			A.n_rows != A.n_cols
		)
	) {
		errno = EINVAL;
		perror ("gaussian_elimination: размеры не совпадают!");
		exit (ERROR);
	}

	int n = A.n_rows;
	int row;
	int col;
	int k;
	Vector V;

	if (i_am_the_master) {
		V = vector_new (n); // результат
	} else {
		V.data = NULL;
		V.size = 0;
	}

	// некоторые процессы могли не получить ни одного столбца
	Vector Zero;
	int is_zero_size = A.local_rows*A.local_cols == 0;
	if (is_zero_size) {
		Zero = vector_new (n);
		// они отправляют нулевые векторы
		memset (Zero.data, 0, n*sizeof(double));
	}


	for (k = 0; k < n; k++) {		
		if (indxg2pc (k) == g_grid_proc.mycol) {
			// каждый процесс умножает свою колонку в матрице на число X[k]
			int row;
			for (row = 0; row < n; row++) {
				A.data[INDEX(row,k,n)] *= X.data[k];
			}
		}
	}

	// каждый процесс суммирует все локальные столбцы и кладет результат в первый
	for (row = 0; row < A.local_rows; row++) {
		for (col = 1; col < A.local_cols; col++) {
			A.data[row] += A.data[col*A.local_rows+row];
		}
	}

	// все процессы посылают MASTER-у свой первый столбец, они все суммируются
	if (is_zero_size) 
		MPI_Reduce (
			Zero.data,
			V.data, n,
			MPI_DOUBLE,
			MPI_SUM,
			MASTER,
			MPI_COMM_WORLD
		);
	else
		MPI_Reduce (
			A.data,
			V.data, n,
			MPI_DOUBLE,
			MPI_SUM,
			MASTER,
			MPI_COMM_WORLD
		);
	if (is_zero_size)
		vector_free (Zero);
	return V;
}

double vector_norm (Vector X) {
	int row;
	double sum = 0;
	for (row = 0; row < X.size; row++)
		sum += sqr (X.data[row]);
	return sum;
}

double matrix_slau_difference (Matrix A, Vector X, Vector B) {
	// невязка = ||Ax-b||

	int n = A.n_rows;
	int row;
	double norm;

	// Diff = Ax
	Vector Diff = matrix_vector_mult (A, X);
	if (i_am_the_master) {
		// Diff = Ax-b
		for (row = 0; row < n; row++)
			Diff.data[row] -= B.data[row];
		norm = vector_norm (Diff);
		vector_free (Diff);
	}
	MPI_Bcast (&norm, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	// возвращаем норму
	return norm;
}