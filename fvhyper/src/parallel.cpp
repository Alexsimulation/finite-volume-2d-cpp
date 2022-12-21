#include <fvhyper/parallel.h>


namespace fvhyper {


mpi_wrapper::mpi_wrapper() {
    MPI_Init(NULL, NULL);
    // Find out rank, size
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
}

int mpi_wrapper::exit() {
    return MPI_Finalize();
}




}
