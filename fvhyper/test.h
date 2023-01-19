/*

       ___     __                    
      / _/  __/ /  __ _____  ___ ____
     / _/ |/ / _ \/ // / _ \/ -_) __/
    /_/ |___/_//_/\_, / .__/\__/_/   
                 /___/_/             

    Finite Volumes for High Performance

    - Description : Testing module headers
    - Author : Alexis Angers
    - Contact : alexis.angers@polymtl.ca

*/
#pragma once

#include <iostream>
#include <fvhyper/parallel.h>


namespace fvhyper {

const int STATUS_SUCCESS = 1982615;


class status {
public:
    int success;
};


class tester {

    status (*to_test)(mpi_wrapper&);
    std::string name;
    mpi_wrapper *pool;

public:
    tester(std::string name_in, status (*to_test_in)(mpi_wrapper&), mpi_wrapper& pool_in);

    void operator()();

};



}

