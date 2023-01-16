/*

       ___     __                    
      / _/  __/ /  __ _____  ___ ____
     / _/ |/ / _ \/ // / _ \/ -_) __/
    /_/ |___/_//_/\_, / .__/\__/_/   
                 /___/_/             

    Finite Volumes for High Performance

    - Description : Mesh class and functions header
    - Author : Alexis Angers
    - Contact : alexis.angers@polymtl.ca

*/
#pragma once

#include <mpi.h>
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <fvhyper/parallel.h>


namespace fvhyper {



extern const int vars;


namespace boundaries {
    extern std::map<std::string, 
        void (*)(double*, double*, double*)> bounds;
}




class mpi_comm_cells {
public:
    std::vector<uint> snd_indices;
    std::vector<uint> rec_indices;
    std::vector<double> snd_q;
    std::vector<double> rec_q;

    uint own_rank;
    uint out_rank;
};


template<uint N>
class meshArray {
private:
    std::vector<uint> nodes;
    uint n;
public:
    meshArray();
    uint& operator()(const uint& i, const uint& j);
    void push_back(const std::vector<uint> v);
    uint cols() const;
    uint rows() const;
    void dump();
    void move_to_end(const uint& i);
};

template<uint N>
meshArray<N>::meshArray() {n=0;}

template<uint N>
uint& meshArray<N>::operator()(const uint& i, const uint& j) {
    return nodes[i*N + j];
}

template<uint N>
void meshArray<N>::push_back(const std::vector<uint> v) {
    if (v.size() != N) {
        throw std::invalid_argument( "size of input array " + std::to_string(v.size()) + " not equal to " + std::to_string(N) );
    }
    for (uint i=0; i<v.size(); ++i) {
        nodes.push_back(v[i]);
    }
    n += 1;
}

template<uint N>
uint meshArray<N>::cols() const {return n;};

template<uint N>
uint meshArray<N>::rows() const {return N;};

template<uint N>
void meshArray<N>::dump() {
    for (uint i=0; i<n; ++i) {
        for (uint j=0; j<N; ++j) {
            std::cout << this->operator()(i, j) << " ";
        }
        std::cout << "\n";
    }
}

template<uint N>
void meshArray<N>::move_to_end(const uint& i) {
    // Move row i to the end
    std::rotate(nodes.begin() + N*i, nodes.begin() + N*i + N, nodes.end());
}



// Class for a partitionned mesh domain
class mesh {
public:

    std::map<std::string, std::string> meshFormat;

    std::vector<double> nodesX;
    std::vector<double> nodesY;
    std::unordered_map<uint, uint> originalNodesRef;

    meshArray<2> edgesNodes;
    meshArray<2> edgesCells;
    std::vector<double> edgesLengths;
    std::vector<double> edgesNormalsX;
    std::vector<double> edgesNormalsY;
    std::vector<double> edgesCentersX;
    std::vector<double> edgesCentersY;

    std::vector<uint> boundaryEdges;
    std::vector<void (*)(double*, double*, double*)> boundaryFuncs;

    std::map<std::tuple<uint, uint>, uint> edgesRef;

    meshArray<4> cellsNodes;
    std::vector<bool> cellsIsTriangle;
    std::vector<double> cellsAreas;
    std::vector<double> cellsCentersX;
    std::vector<double> cellsCentersY;
    std::vector<bool> cellsIsGhost;

    std::unordered_map<uint, uint> originalCellsRef;
    std::unordered_map<uint, uint> originalCellsRefInv;

    std::vector<uint> ghostCellsOriginalIndices;
    std::vector<uint> ghostCellsCurrentIndices;
    std::vector<uint> ghostCellsOwners;

    std::vector<mpi_comm_cells> comms;

    uint nRealCells;    // Number of real cells (start of ghost cells)
    uint nNonBoundCells;    // Number of non-boundary cells (start of boundary cells)

    void read_file(std::string filename, mpi_wrapper& pool);

    void make_comms(uint rank);

    int find_edge_with_nodes(const uint n0, const uint n1);
    uint find_bound_with_nodes(const uint n0, const uint n1);
    bool find_if_edge_in_mesh(const uint n0, const uint n1);

    void convert_node_face_info();
    void compute_mesh();
    void add_cell_edges(uint cell_id);

    void send_mesh_info();
};

std::vector<std::string> cut_str(const std::string s, const char c);
std::vector<uint> str_to_ints(std::string s);





}
