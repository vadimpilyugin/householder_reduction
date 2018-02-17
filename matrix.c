// cpp header
#include "matrix.h"

// Standard headers
#include <stdio.h>
// memcpy
#include <string.h>
// errno
#include <errno.h>
// lseek, read
#include <sys/types.h>
#include <unistd.h>
// open
#include <sys/stat.h>
#include <fcntl.h>
// sqrt
#include <math.h>
// exit
#include <stdlib.h>
// int usleep(useconds_t usec);
#include <unistd.h>


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

static void guess_process_grid_dimensions (int *n_proc_rows_, int *n_proc_cols_) {
	int nproc;
	int myrank; // не используется
	int n_proc_rows;

	// Cblacs_pinfo(&myrank, &nproc);
	MPI_Comm_size (MPI_COMM_WORLD, &nproc);
	//for (n_proc_rows = (int) sqrt(nproc); (nproc % n_proc_rows != 0); n_proc_rows--);
	*n_proc_rows_ = 1;//n_proc_rows;
	*n_proc_cols_ = nproc;//nproc / n_proc_rows;
}

static GridProcess new_grid (int n_proc_rows, int n_proc_cols, int myrank, int nproc) {
	int context;

	// Cblacs_get(NO_PARAM, GET_DEFAULT_CONTEXT, &context);
	// Cblacs_gridinit(&context, (char *) "Row-major", n_proc_rows, n_proc_cols);

	GridProcess grid_proc;
	// grid_proc.context = context;
	// Cblacs_gridinfo(
	// 	context, 
	// 	&(grid_proc.n_proc_rows),
	// 	&(grid_proc.n_proc_cols),
	// 	&(grid_proc.myrow),
	// 	&(grid_proc.mycol)
	// );
	grid_proc.n_proc_rows = n_proc_rows;
	grid_proc.n_proc_cols = n_proc_cols;
	grid_proc.myrow = 0;
	grid_proc.mycol = myrank;

	grid_proc.nproc = nproc; //grid_proc.n_proc_rows*grid_proc.n_proc_cols;
	grid_proc.myrank = myrank; // grid_proc.myrow*grid_proc.n_proc_cols+grid_proc.mycol;
	return grid_proc;
}

void matrix_init (const int block_size) {

	int myrank, nproc;
	MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
	MPI_Comm_size (MPI_COMM_WORLD, &nproc);
	i_am_the_master = myrank == MASTER;

	g_block_size = block_size;
	int n_proc_rows = 1;
	int n_proc_cols = nproc;
	// guess_process_grid_dimensions (&n_proc_rows, &n_proc_cols);
	g_grid_proc = new_grid (n_proc_rows, n_proc_cols, myrank, nproc);
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
}

void matrix_exit () {
	// Cblacs_gridexit (g_grid_proc.context);
	// Cblacs_exit (1);
}

Matrix matrix_new(const int dim_) {
	Matrix A;
	int master = MASTER;
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
					snprintf (g_str, STR_LEN, "%7.3f", local_elem);
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

// void matrix_print_diag (Matrix A) {
// 	int global_row;
// 	int global_col;

// 	for (global_row = 0; global_row < A.n_rows; global_row++) {
// 		MPI_Barrier (MPI_COMM_WORLD);


// 		global_col = global_row;
// 		int process_row = (global_row/g_block_size)%g_grid_proc.n_proc_rows;
// 		int process_col = (global_col/g_block_size)%g_grid_proc.n_proc_cols;
		
// 		if (process_row == g_grid_proc.myrow && process_col == g_grid_proc.mycol) {
// 			{
// 				int local_row = indxg2lr (global_row);
// 				int local_col = indxg2lc (global_col);
// 				double local_elem = A.data[local_col*A.local_rows+local_row];
// 				snprintf (g_str, STR_LEN, "%8.4f", cabs (local_elem));
// 			}
// 			MPI_Request request;
// 			MPI_Isend (g_str,STR_LEN,MPI_CHAR,MASTER,0,MPI_COMM_WORLD, &request);
// 		}
// 		if (i_am_the_master) {
// 			int source = process_row*g_grid_proc.n_proc_cols+process_col;
// 			MPI_Status status;

// 			MPI_Recv (g_str,STR_LEN,MPI_CHAR,source,0,MPI_COMM_WORLD,&status);
// 			printf ("%s",g_str);
// 			fflush (stdout);
// 			if (global_col == A.n_cols-1)
// 				printf ("\n\n");
// 		}

// 		MPI_Barrier (MPI_COMM_WORLD);
// 	}
// }

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
