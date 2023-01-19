/*

       ___     __                    
      / _/  __/ /  __ _____  ___ ____
     / _/ |/ / _ \/ // / _ \/ -_) __/
    /_/ |___/_//_/\_, / .__/\__/_/   
                 /___/_/             

    Finite Volumes for High Performance

    - Description : Testing module sources
    - Author : Alexis Angers
    - Contact : alexis.angers@polymtl.ca

*/
#include <fvhyper/test.h>
#include <fvhyper/parallel.h>



namespace fvhyper {



tester::tester(std::string name_in, status (*to_test_in)(mpi_wrapper&), mpi_wrapper& pool_in) {
    name = name_in;
    to_test = to_test_in;
    pool = &pool_in;
}

void tester::operator()() {
    auto s = to_test(*pool);

    if (pool->rank == 0) {
        std::cout << "Test " << name << " (" << std::flush;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    
    for (uint i=0; i<pool->size; ++i) {
        int a;
        if (i == pool->rank) {
            a = s.success;
        }
        MPI_Bcast(&a, 1, MPI_INT, i, MPI_COMM_WORLD);
        if (pool->rank == 0) {
            if (i!=0) {std::cout << ", ";}
            std::cout << std::to_string(i) + ": ";
            if (a == STATUS_SUCCESS) {
                std::cout << "\033[1;32msuccess\033[0m" << std::flush;
            } else {
                std::cout << "\033[1;31mfailure(" + std::to_string(a) + ")\033[0m" << std::flush;
            }
        }
    }
    if (pool->rank == 0) {std::cout << ")" << std::endl;}

    MPI_Barrier(MPI_COMM_WORLD);

}



}


