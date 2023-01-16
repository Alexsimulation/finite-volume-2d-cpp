/*

       ___     __                    
      / _/  __/ /  __ _____  ___ ____
     / _/ |/ / _ \/ // / _ \/ -_) __/
    /_/ |___/_//_/\_, / .__/\__/_/   
                 /___/_/             

    Finite Volumes for High Performance

    - Description : Post-processing functions header
    - Author : Alexis Angers
    - Contact : alexis.angers@polymtl.ca

*/
#pragma once

#include <fvhyper/mesh.h>


namespace fvhyper {


namespace post {
    extern std::map<std::string, 
        void (*)(double*, double*)> extra_scalars;
    
    extern std::map<std::string, 
        void (*)(double*, double*)> extra_vectors;
}


void writeVtk(
    const std::string name,
    std::vector<std::string> varNames,
    std::vector<double>& q,
    mesh& m,
    int rank,
    int world_size
);

}


