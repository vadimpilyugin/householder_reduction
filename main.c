
// Libraries
#include "mpi.h"
#include "matrix.h"

int main (int argc, char * argv[]) {
  MPI_Init (&argc, &argv);
  placeholder ();
  MPI_Finalize ();
  return 0;
}
