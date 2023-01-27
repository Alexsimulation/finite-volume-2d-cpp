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


template<class T>
void move_to_end(std::vector<T>& v, const uint& i) {
    // Move element i in vector to the end
    std::rotate(v.begin() + i,  v.begin() + i + 1, v.end());
}


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
    void swap(const uint& i, const uint& j);
    void move_to_end(const uint& i);
    const std::vector<uint>& get_vector();
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
void meshArray<N>::swap(const uint& i, const uint& j) {
    // Swap row i and row j
    for (uint k=0; k<N; ++k) std::iter_swap(nodes.begin()+N*i+k, nodes.begin()+N*j+k);
}

template<uint N>
void meshArray<N>::move_to_end(const uint& i) {
    // Move row i to the end
    std::rotate(nodes.begin() + N*i, nodes.begin() + N*i + N, nodes.end());
}

template<uint N>
const std::vector<uint>& meshArray<N>::get_vector() {
    return nodes;
}



// Class for a partitionned mesh domain
class mesh {
public:

    std::string filename;

    bool do_compute_wall_dist = false;

    std::map<std::string, std::string> meshFormat;

    std::map<uint, std::string> physicalNames;
    std::map<uint, uint> entityTagToPhysicalTag;
    std::map<std::string, uint> entitiesNumber;

    std::vector<double> nodesX;
    std::vector<double> nodesY;
    std::map<uint, uint> originalNodesRef;

    meshArray<2> edgesNodes;
    meshArray<2> edgesCells;
    std::vector<double> edgesLengths;
    std::vector<double> edgesNormalsX;
    std::vector<double> edgesNormalsY;
    std::vector<double> edgesCentersX;
    std::vector<double> edgesCentersY;

    std::vector<uint> boundaryEdges;
    std::vector<uint> boundaryEdges0;
    std::vector<uint> boundaryEdges1;
    std::vector<uint> boundaryEdgesIntTag;
    std::vector<void (*)(double*, double*, double*)> boundaryFuncs;

    std::map<std::tuple<uint, uint>, uint> edgesRef;

    meshArray<4> cellsNodes;
    std::vector<bool> cellsIsTriangle;
    std::vector<double> cellsAreas;
    std::vector<double> cellsCentersX;
    std::vector<double> cellsCentersY;
    std::vector<bool> cellsIsGhost;

    std::map<uint, uint> originalToCurrentCells;    // maps original index -> current index
    std::map<uint, uint> currentToOriginalCells;    // maps original index -> current index

    std::vector<uint> ghostCellsOriginalIndices;
    std::vector<uint> ghostCellsCurrentIndices;
    std::vector<uint> ghostCellsOwners;

    std::vector<double> wall_dist;

    uint nRealCells;

    std::vector<mpi_comm_cells> comms;

    void read_entities();
    void read_nodes();
    void read_boundaries();
    void read_elements();
    void read_ghost_elements();

    void add_boundary_cells();

    void read_file(std::string filename, mpi_wrapper& pool);

    void make_comms(uint rank);

    int find_edge_with_nodes(const uint n0, const uint n1);
    uint find_bound_with_nodes(const uint n0, const uint n1);
    bool find_if_edge_in_mesh(const uint n0, const uint n1);

    void convert_node_face_info();
    void compute_mesh();
    void add_cell_edges(uint cell_id);

    void compute_wall_dist(mpi_wrapper& pool);

    void send_mesh_info();
};

std::vector<std::string> cut_str(const std::string s, const char c);
std::vector<uint> str_to_ints(std::string s);





}
