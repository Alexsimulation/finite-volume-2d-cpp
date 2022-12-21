#pragma once

#include <mpi.h>


namespace fvhyper {


class mpi_wrapper {

public:
    int rank;
    int size;

    mpi_wrapper();

    int exit();
};



}
