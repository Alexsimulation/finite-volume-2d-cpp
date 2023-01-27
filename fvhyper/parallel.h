/*

       ___     __                    
      / _/  __/ /  __ _____  ___ ____
     / _/ |/ / _ \/ // / _ \/ -_) __/
    /_/ |___/_//_/\_, / .__/\__/_/   
                 /___/_/             

    Finite Volumes for High Performance

    - Description : Parallel mpi wrapper header
    - Author : Alexis Angers
    - Contact : alexis.angers@polymtl.ca

*/
#pragma once

#include <mpi.h>
#include <vector>


namespace fvhyper {


class mpi_wrapper {

public:
    int rank;
    int size;

    mpi_wrapper();

    int exit();
};


class mpi_comm_cells {
public:
    std::vector<uint> snd_indices;
    std::vector<uint> rec_indices;
    std::vector<double> snd_q;
    std::vector<double> rec_q;

    uint own_rank;
    uint out_rank;

    uint start_index;
    uint length;
};



}
