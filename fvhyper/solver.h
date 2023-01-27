/*

       ___     __                    
      / _/  __/ /  __ _____  ___ ____
     / _/ |/ / _ \/ // / _ \/ -_) __/
    /_/ |___/_//_/\_, / .__/\__/_/   
                 /___/_/             

    Finite Volumes for High Performance

    - Description : Solver header
    - Author : Alexis Angers
    - Contact : alexis.angers@polymtl.ca

*/
#pragma once

#include <fvhyper/mesh.h>
#include <fvhyper/physics.h>
#include <fvhyper/explicit.h>
#include <iostream>



namespace fvhyper {




class solver {
public:
    std::vector<physics> phys;

    // MPI
    mpi_wrapper* pool;

    // Stage coefficients
    std::vector<double> alpha = {1.0};

    // Options
    solverOptions opt;

    // Variables
    uint step_counter = 0;
    

protected:


    virtual void print();

    virtual void dt();

    virtual void step();

public:
    virtual void run();

}



}




