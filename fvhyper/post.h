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


