
// Libraries
#include "mpi.h"
#include "matrix.h"

int main (int argc, char * argv[]) {
  MPI_Init (&argc, &argv);
  matrix_init (1);
  matrix_exit ();
  MPI_Finalize ();
  return 0;
}
