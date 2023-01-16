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


namespace fvhyper {


class mpi_wrapper {

public:
    int rank;
    int size;

    mpi_wrapper();

    int exit();
};



}
