// Standard headers
#include <stdio.h>
// exit, atoi
#include <stdlib.h>
// time
#include <time.h>
// rand
#include <stdlib.h>

// Libraries
#include "mpi.h"
#include "matrix.h"

// Definitions
#define N_ARGS_FORMULA 2
#define ARG_PROG 0
#define ARG_SIZE 1

#define MASTER 0
#define EPS 1e-7
#define ERROR 2
#define DAY 3600*24

typedef unsigned int uint;
int i_am_the_master;
uint randr_seed;
int g_n;
// Измерение времени
double sum_time = 0, tmp_time;

void go()
{
	tmp_time = MPI_Wtime();
}
double stop()
{
	tmp_time = MPI_Wtime() - tmp_time;
	sum_time += tmp_time;
	double t_sec = sum_time;
	sum_time = 0;
	return t_sec;
}

double matrix_filler (int row, int col) {
  // row и col от 0
  return 1/(row+col+1+0.0);
}
double vector_filler (int row) {
  double sum = 0;
  for (int col = 0; col < g_n; col++)
    if (col % 2 == 1)
      sum += matrix_filler (row, col);
  return sum;
}
double rand_d () {
	return (rand_r (&randr_seed) - RAND_MAX/2)/((double)RAND_MAX)*100;
}
double vector_filler_random (int row) {
	return rand_d();
}
double matrix_filler_random (int row, int col) {
	return rand_d();
}

int main (int argc, char * argv[]) {

  MPI_Init (&argc, &argv);
  matrix_init ();

  int rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  randr_seed = time(NULL)+rank*DAY;
  int nproc;
  MPI_Comm_size (MPI_COMM_WORLD, &nproc);
  i_am_the_master = rank == MASTER;
	Matrix A;
	Vector V;

  // ввод матрицы
  if (argc != N_ARGS_FORMULA) {
  	if (i_am_the_master)
  		fprintf (stderr, "\n\nUsage: %s <size>\n", argv[0]);
		exit (ERROR);
  } else {
  	int size = atoi (argv[ARG_SIZE]);
    g_n = size;
		A = matrix_new_and_fill (size, matrix_filler_random);
		if (i_am_the_master)
			V = vector_new_and_fill (size, vector_filler_random);
  }
  MPI_Barrier (MPI_COMM_WORLD);
  if (i_am_the_master)
  	fprintf (stderr,"--- Число процессов: %d\n", nproc);
  if (i_am_the_master)
  	fprintf (stderr,"--- Размер системы: %d\n", A.n_rows);

	// приведение к треугольному виду
	go ();
	matrix_triagonalize (A, V);
	double triag_time = stop ();
	if (i_am_the_master)
		fprintf (stderr,"--- Приведение к треугольному виду: %lf сек\n", triag_time);

	// вектор правой части портится при обратном ходе
	Vector V_copy;
	if (i_am_the_master)
		V_copy = vector_new_copy (V); 

	// решение обратным ходом
	go ();
	Vector X = gaussian_elimination (A, V);
	double gauss_time = stop ();
	if (i_am_the_master && X.data != NULL)
		fprintf (stderr,"--- Обратный ход метода Гаусса: %lf сек\n", gauss_time);


	if (X.data != NULL) {
		// матрица невырожденная

		// считаем невязку
		double diff = matrix_slau_difference (A, X, V_copy);
		if (i_am_the_master)
			fprintf (stderr,"--- Невязка = %f\r\n\r\n", diff);
		if (diff > EPS) {
			if (i_am_the_master)
				fprintf (stderr, "Невязка не нулевая\n");
			exit (ERROR);
		}

		if (i_am_the_master) {
			vector_print (X);
		}
	}
	matrix_free (A);
	vector_free (X);
	if (i_am_the_master) {
		vector_free (V);
		vector_free (V_copy);
	}
  matrix_exit ();
  MPI_Finalize ();
  return 0;
}

