#pragma once

#include <fvhyper/mesh.h>


namespace fvhyper {

void writeVtk(
    const std::string name,
    std::vector<std::string> varNames,
    std::vector<double>& q,
    mesh& m,
    int rank,
    int world_size
);

}


