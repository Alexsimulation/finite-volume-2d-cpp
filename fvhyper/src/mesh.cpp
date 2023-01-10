#include <fvhyper/mesh.h>
#include <fstream>
#include <map>
#include <algorithm>
#include <stdexcept>



namespace fvhyper {





// Cut string by character c
std::vector<std::string> cut_str(const std::string s, const char c) {
	std::vector<std::string> out;
	std::vector<int> I;
	// find zeros
	for (int i=0; i < s.size(); ++i) {
		if (s[i] == c) {
			I.push_back(i);
		}
	}
	// get each integer value
	if (I.size() == 0) {
		out.push_back( s );
	} else {
		out.push_back( s.substr(0, I[0]) );
		for (int j=0; j < (I.size()-1); ++j) {
			out.push_back( s.substr(I[j]+1, I[j+1]-I[j]-1) );
		}
		if (I[I.size()-1] != s.size()-1) {
			out.push_back( s.substr(I[I.size()-1]+1, s.size()-I[I.size()-1]-1) );
		}
	}
	return out;
}

std::vector<uint> str_to_ints(std::string s) {
	std::vector<uint> out;
	std::vector<uint> I;
	// find zeros
	for (uint i=0; i < s.size(); ++i) {
		if (s[i] == ' ') {
			I.push_back(i);
		}
	}
	// get each integer value
	if (I.size() == 0) {
		out.push_back( std::stoi(s) );
	} else {
		out.push_back( std::stoi(s.substr(0, I[0])) );
		for (uint j=0; j < (I.size()-1); ++j) {
			out.push_back( std::stoi(s.substr(I[j]+1, I[j+1]-I[j])) );
		}
		if (I[I.size()-1] != s.size()-1) {
			out.push_back( std::stoi(s.substr(I[I.size()-1]+1, s.size()-I[I.size()-1]-1)) );
		}
	}
	return out;
}

std::vector<double> str_to_floats(const std::string s) {
	std::vector<double> out;
	std::vector<int> I;
	// find zeros
	for (int i=0; i < s.size(); ++i) {
		if (s[i] == ' ') {
			I.push_back(i);
		}
	}
	// get each integer value
	if (I.size() == 0) {
		out.push_back( std::stod(s) );
	} else {
		out.push_back( std::stod(s.substr(0, I[0])) );
		for (int j=0; j < (I.size()-1); ++j) {
			out.push_back( std::stod(s.substr(I[j]+1, I[j+1]-I[j])) );
		}
		if (I[I.size()-1] != s.size()-1) {
			out.push_back( std::stod(s.substr(I[I.size()-1]+1, s.size()-I[I.size()-1]-1)) );
		}
	}
	return out;
}

int mesh::find_edge_with_nodes(const uint n0, const uint n1) {
    const uint nmin = std::min(n0, n1);
    const uint nmax = std::max(n0, n1);
    // Also check in bounds
    const auto tp = std::make_tuple(nmin, nmax);
    if (edgesRef.find(tp) != edgesRef.end()) {
        // Found
        return edgesRef.at(tp);
    } else {
        return -1;
    }
}


bool mesh::find_if_edge_in_mesh(const uint n0, const uint n1) {
    const uint nmin = std::min(n0, n1);
    const uint nmax = std::max(n0, n1);
    const auto tp = std::make_tuple(nmin, nmax);
    if (edgesRef.find(tp) != edgesRef.end()) {
        // Found
        return true;
    }
    return false;
}



void mesh::add_last_cell_edges(uint size) {

    uint nfaces = cellsAreas.size()-1;
    for (uint i=0; i<size; ++i) {
        uint j =  (i<(size-1)) ? (i+1) : 0;

        if (!find_if_edge_in_mesh(cellsNodes(nfaces, i), cellsNodes(nfaces, j))) {
            edgesNormalsX.push_back(0.);
            edgesNormalsY.push_back(0.);
            edgesCentersX.push_back(0.);
            edgesCentersY.push_back(0.);

            std::vector<uint> efi = {nfaces, 0};
            edgesCells.push_back(efi);

            std::vector<uint> eni = {cellsNodes(nfaces, i), cellsNodes(nfaces, j)};
            edgesNodes.push_back(eni);
            edgesLengths.push_back(0.);

            const uint nmin = std::min(cellsNodes(nfaces, i), cellsNodes(nfaces, j));
            const uint nmax = std::max(cellsNodes(nfaces, i), cellsNodes(nfaces, j));
            edgesRef.insert({ std::make_tuple(nmin, nmax), edgesLengths.size()-1 });
        }
    }
}


void mesh::convert_node_face_info() {
    // Take a mesh input, and convert its info to correct fvm data
    // Currently,
    //      edges contain no face connectivity data
    //      bounds contain no face connectivity data

    // Loop over all faces
    std::vector<int> cellEdges(4);
    for (int i=0; i<cellsAreas.size(); ++i) {

		// Find edges
        cellEdges[0] = find_edge_with_nodes(cellsNodes(i, 0), cellsNodes(i, 1));
        cellEdges[1] = find_edge_with_nodes(cellsNodes(i, 1), cellsNodes(i, 2));
        if (cellsIsTriangle[i]) {
            cellEdges[2] = find_edge_with_nodes(cellsNodes(i, 2), cellsNodes(i, 0));
            cellEdges[3] = -1;
        } else {
            cellEdges[2] = find_edge_with_nodes(cellsNodes(i, 2), cellsNodes(i, 3));
            cellEdges[3] = find_edge_with_nodes(cellsNodes(i, 3), cellsNodes(i, 0));
        }

        // Now update these edges with face connectivity info
        for (int ei : cellEdges) {
            if (ei != -1) {
                if (i != edgesCells(ei, 0)) {
                    edgesCells(ei, 1) = i;
                }
            }
        }
    }
}


void mesh::compute_mesh() {

	// Compute mesh normals, areas, lengths
    // Loop over each face
	for (int i=0; i<cellsAreas.size(); ++i) {
		// Evaluate face centroid position
        double cellC[2];
        cellC[0] = 0.; cellC[1] = 0.;

        uint n_cells = cellsIsTriangle[i] ? 3 : 4;
        for (uint j=0; j<n_cells; ++j) {
            cellC[0] += nodesX[cellsNodes(i, j)]/((double) n_cells);
            cellC[1] += nodesY[cellsNodes(i, j)]/((double) n_cells);
        }
        cellsCentersX[i] = cellC[0];
        cellsCentersY[i] = cellC[1];
    }

    // Evaluate edge normals and edge centers
	// Loop over each edge
	for (int i=0; i<edgesLengths.size(); ++i) {
        edgesCentersX[i] = (nodesX[edgesNodes(i, 1)] + nodesX[edgesNodes(i, 0)])*0.5;
        edgesCentersY[i] = (nodesY[edgesNodes(i, 1)] + nodesY[edgesNodes(i, 0)])*0.5;

        double dex = nodesX[edgesNodes(i, 1)] - nodesX[edgesNodes(i, 0)];
        double dey = nodesY[edgesNodes(i, 1)] - nodesY[edgesNodes(i, 0)];
        double l = sqrt(dex*dex + dey*dey);
		edgesLengths[i] = l;
		edgesNormalsX[i] = -dey/l;
        edgesNormalsY[i] =  dex/l;

        // Check normal direction
        double dot_f_f = 
              edgesNormalsX[i]*(edgesCentersX[i] -  cellsCentersX[edgesCells(i, 0)])
            + edgesNormalsY[i]*(edgesCentersY[i] -  cellsCentersY[edgesCells(i, 0)]);
        if (dot_f_f < 0) {
            edgesNormalsX[i] *= -1.;
            edgesNormalsY[i] *= -1.;
        }
	}

    // Compute face area
    for (int i=0; i<cellsAreas.size(); ++i) {
        double& x1 = nodesX[cellsNodes(i, 0)];
        double& x2 = nodesX[cellsNodes(i, 1)];
        double& x3 = nodesX[cellsNodes(i, 2)];

        double& y1 = nodesY[cellsNodes(i, 0)];
        double& y2 = nodesY[cellsNodes(i, 1)];
        double& y3 = nodesY[cellsNodes(i, 2)];
        
        if (cellsIsTriangle[i]) {
            cellsAreas[i] = 0.5*abs(
                x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
            );
        } else {
            double& x4 = nodesX[cellsNodes(i, 3)];
            double& y4 = nodesY[cellsNodes(i, 3)];

            cellsAreas[i] = 0.5*abs(
                x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
            ) + 0.5*abs(
                x1*(y3-y4) + x3*(y4-y1) + x4*(y1-y3)
            );
        }

        if (cellsAreas[i] < 1e-10) {
            std::cout << "Null cell area for node " << i << " ";
        }
    }
}


void mesh::read_file(std::string name, mpi_wrapper& pool) {
    uint rank = pool.rank;

    std::string filename = "";
    filename += (pool.size > 1) ? (name + "_" + std::to_string(pool.rank + 1) + ".msh") : name + ".msh";

    std::ifstream infile(filename);

    std::string line;

    std::string entitiesSection;

    uint n = 0;
    uint ns = 0;
    uint nss = 0;
    uint start_offset = 0;

    std::map<uint, std::string> physicalNames;
    std::map<uint, uint> entityTagToPhysicalTag;
    std::map<std::string, uint> entitiesNumber;

    std::vector<uint> boundaryEdges0;
    std::vector<uint> boundaryEdges1;
    std::vector<uint> boundaryEdgesIntTag;

    uint n_in_block = 0;
    std::vector<uint> block_tags;

    uint blockPhysicalTag;
    uint nGhostEntities;

    std::string currentSection = "start";

    uint counter = 0;
    while (std::getline(infile, line)) {
        //std::cout << pool.rank << " got to line " << counter << std::endl;
        counter += 1;

        if (line[0] == '$') {
            // This line tells us a section number
            currentSection = line.substr(1, line.size());
            // Restart current section line
            ns = 0;
        } else if (currentSection == "MeshFormat") {
            // Read mesh format
            if (ns == 1) {
                meshFormat["format"] = line;
            }
        } else if (currentSection == "PhysicalNames") {
            if (ns > 1) {
                auto lineS = cut_str(line, ' ');
                std::string pnamei = lineS[2];
                pnamei.erase(remove(pnamei.begin(), pnamei.end(), '"'), pnamei.end());
                physicalNames[std::stoi(lineS[1])] = pnamei;
            }
        } else if (currentSection == "Entities") {
            // Entities are nodes on border elements
            if (ns == 1) {
                // First line, contains number of entities
                auto l = str_to_ints(line);
                entitiesNumber["points"] = l[0];
                entitiesNumber["curves"] = l[1] + l[0];
                entitiesNumber["surfaces"] = l[2] + l[1] + l[0];
                entitiesNumber["volumes"] = l[3] + l[2] + l[1] + l[0];

                nss = 0;
            } else if (nss <= entitiesNumber["points"]) {
                // Reading points
            } else if (nss <= entitiesNumber["curves"]) {
                // Reading curves
                auto l = str_to_ints(line);
                entityTagToPhysicalTag[l[0]] = l[8];
            } else if (nss <= entitiesNumber["surfaces"]) {
                // Reading surfaces
                auto l = str_to_ints(line);
            }
        } else if (currentSection == "PartitionedEntities") {
            if (ns < 2) {
            } else if (ns == 2) {
                auto l = str_to_ints(line);
                nGhostEntities = l[0];
                nss = 524288;
            } else if (ns == (3 + nGhostEntities)) {
                auto l = str_to_ints(line);
                entitiesNumber["points"] = l[0];
                entitiesNumber["curves"] = l[1] + l[0];
                entitiesNumber["surfaces"] = l[2] + l[1] + l[0];
                entitiesNumber["volumes"] = l[3] + l[2] + l[1] + l[0];

                nss = 0;
            } else if (nss <= entitiesNumber["points"]) {
                // Reading points
            } else if (nss <= entitiesNumber["curves"]) {
                // Reading curves
                auto l = str_to_ints(line);
                if (l[1] == 1) 
                    entityTagToPhysicalTag[l[0]] =
                        entityTagToPhysicalTag.at(l[2]);
            } else if (nss <= entitiesNumber["surfaces"]) {
                // Reading surfaces
            }
        } else if (currentSection == "Nodes") {
            // Read nodes
            if (ns == 1) {
                // Read number of node blocks
                start_offset = 1;
                nss = 0;
            } else if (ns > (start_offset + 2*n_in_block)) {
                // Read block
                auto l = str_to_ints(line);
                n_in_block = l[3];
                start_offset = ns;
                block_tags.resize(n_in_block);
                nss = 0;
            } else if (nss <= n_in_block) {
                // Read current node tag
                block_tags[nss-1] = std::stoi(line) - 1;
            } else {
                // Read current node
                uint tag_i = block_tags[nss - 1 - n_in_block];
                auto l = str_to_floats(line);
                originalNodesRef[tag_i] = nodesX.size();
                nodesX.push_back(l[0]);
                nodesY.push_back(l[1]);
            }
        } else if (currentSection == "Elements") {
            // Read elements
            if (ns == 1) {
                // Read number of node blocks
                start_offset = 1;
                n_in_block = 0;
                nss = 0;
            } else if (ns > (start_offset + n_in_block)) {
                // Read block
                auto l = str_to_ints(line);
                if (l[0] == 1) {
                    blockPhysicalTag = entityTagToPhysicalTag.at(l[1]);
                }
                n_in_block = l[3];
                start_offset = ns;
                nss = 0;
            } else {
                // Read current tag
                auto li = str_to_ints(line);
                uint tag_i = li[0] - 1;
                std::vector<uint> l(li.size()-1);
                for (uint i=1; i<li.size(); ++i) {
                    l[i-1] = li[i] - 1;
                }
                if (l.size() == 2) {
                    // Boundary edge
                    boundaryEdges0.push_back(originalNodesRef.at(l[0]));
                    boundaryEdges1.push_back(originalNodesRef.at(l[1]));
                    boundaryEdgesIntTag.push_back(blockPhysicalTag);
                } else {
                    // Triangle or quad cell
                    originalCellsRef[cellsIsTriangle.size()] = tag_i;
                    originalCellsRefInv[tag_i] = cellsIsTriangle.size();

                    uint l_size = l.size();

                    for (uint i=0; i<l.size(); ++i) {
                        l[i] = originalNodesRef.at(l[i]);
                    }

                    if (l_size == 3) {
                        l.push_back(0);
                        cellsIsTriangle.push_back(true);
                    } else {
                        cellsIsTriangle.push_back(false);
                    }

                    cellsNodes.push_back(l);
                    cellsAreas.push_back(0.);
                    cellsCentersX.push_back(0.);
                    cellsCentersY.push_back(0.);
                    cellsIsGhost.push_back(false);
                    
                    add_last_cell_edges(l_size);
                }
            }
        } else if (currentSection == "GhostElements") {
            // Read ghost elements
            if (ns == 1) {
                // Read number of ghost elements
                uint n_ghost_cells = std::stoi(line);
                ghostCellsOriginalIndices.reserve(n_ghost_cells);
                ghostCellsCurrentIndices.reserve(n_ghost_cells);
                ghostCellsOwners.reserve(n_ghost_cells);
                nss = 0;
            } else {
                // Read block
                auto l = str_to_ints(line);
                ghostCellsOriginalIndices.push_back(l[0] - 1);
                ghostCellsCurrentIndices.push_back(originalCellsRefInv.at(l[0] - 1));
                ghostCellsOwners.push_back(l[1] - 1);
                cellsIsGhost[originalCellsRefInv.at(l[0] - 1)] = true;
            }
        }

        n += 1;
        ns += 1;
        nss += 1;
    }


    // Fix edges part of ghost cells
    for (auto& i : ghostCellsCurrentIndices) {
        const uint cell_size = cellsIsTriangle[i] ? 3 : 4;
        for (uint j=0; j<cell_size; ++j) {
            uint k = (j==(cell_size-1)) ? 0 : j+1;
            uint n0 = cellsNodes(i, j);
            uint n1 = cellsNodes(i, k);

            const uint nmin = std::min(n0, n1);
            const uint nmax = std::max(n0, n1);
            const auto tp = std::make_tuple(nmin, nmax);
            const uint e = edgesRef.at(tp);

            edgesCells(e, 1) = edgesCells(e, 0);
        }
    }


    // Convert nodes and faces info
    convert_node_face_info();

    // Compute mesh metrics
    compute_mesh();

    // Generate communicators
    make_comms(rank);


    // Add boundary cells
    nRealCells = cellsAreas.size();
    for (uint i=0; i<boundaryEdges0.size(); ++i) {

        // Find edge index on boundary
        const uint nmin = std::min(boundaryEdges0[i], boundaryEdges1[i]);
        const uint nmax = std::max(boundaryEdges0[i], boundaryEdges1[i]);
        std::tuple<uint, uint> tp = std::make_tuple(nmin, nmax);
        uint e = 0;
        if (edgesRef.find(tp) != edgesRef.end()) {
            // Found
            e = edgesRef.at(tp);
        } else {
            throw std::invalid_argument( "invalid edge ref (" + std::to_string(nmin) + ", " + std::to_string(nmax) + ")" );
        }

        uint cell = edgesCells(e, 0);

        // Compute virtual cell properties
        double area = cellsAreas[cell];
        double dist  = area / edgesLengths[e];
        double cx = edgesCentersX[e] + dist*edgesNormalsX[e];
        double cy = edgesCentersY[e] + dist*edgesNormalsY[e];
        std::vector<uint> thisCellNodes = {nmin, nmax, 0, 0};

        // Add virtual cell
        cellsAreas.push_back(area);
        cellsCentersX.push_back(cx);
        cellsCentersY.push_back(cy);
        cellsIsTriangle.push_back(true);
        cellsNodes.push_back(thisCellNodes);

        // Add cell index to edge
        edgesCells(e, 1) = cellsAreas.size() - 1;

        // Add edge index to boundary edges
        boundaryEdges.push_back(e);
        if (boundaries::bounds.find(physicalNames[boundaryEdgesIntTag[i]]) == boundaries::bounds.end()) {
            throw std::invalid_argument("Physical name " + physicalNames[boundaryEdgesIntTag[i]] + " of tag " + std::to_string(boundaryEdgesIntTag[i]) + " not in mesh");
        }
        boundaryFuncs.push_back(
            boundaries::bounds.at(physicalNames.at(boundaryEdgesIntTag[i]))
        );
    }

}



void mesh::make_comms(uint rank) {
    // Make communicators

    // Compute all ranks of connected nodes
    std::vector<uint> connected_nodes;

    for (auto i : ghostCellsOwners) {
        bool is_in = false;
        for (auto& k : connected_nodes) {
            if (k == i) {
                is_in = true;
            }
        }
        if (!is_in) {
            connected_nodes.push_back(i);
        }
    }

    // Set number of communicators and ranks
    comms.resize(connected_nodes.size());
    for (uint i=0; i<connected_nodes.size(); ++i) {
        comms[i].own_rank = rank;
        comms[i].out_rank = connected_nodes[i];
    }
    // Set snd and rec sizes
    std::vector<uint> sizes(connected_nodes.size());
    for (auto i=0; i<sizes.size(); ++i) {
        sizes[i] = 0;
    }
    for (auto i : ghostCellsOwners) {
        for (uint j=0; j<comms.size(); ++j) {
            if (i == comms[j].out_rank) {
                sizes[j] += 1;
            }
        }
    }
    for (uint i=0; i<comms.size(); ++i) {
        comms[i].snd_indices.reserve(sizes[i]);
        comms[i].rec_indices.reserve(sizes[i]);
        comms[i].snd_q.resize(vars*sizes[i]);
        comms[i].rec_q.resize(vars*sizes[i]);
    }

    // Set rec indices
    for (uint i=0; i<ghostCellsCurrentIndices.size(); ++i) {
        for (uint j=0; j<comms.size(); ++j) {
            if (ghostCellsOwners[i] == comms[j].out_rank) {
                comms[j].rec_indices.push_back(
                    ghostCellsOriginalIndices[i]
                );
                comms[j].snd_indices.push_back(0);
            }
        }
    }

    // Send to other partition the number of cells we are sending
    for (auto& comm : comms) {
        // Send values
        uint this_size = comm.rec_indices.size();
        MPI_Send(
        /* data         = */ &this_size, 
        /* count        = */ 1, 
        /* datatype     = */ MPI_UNSIGNED,
        /* destination  = */ comm.out_rank, 
        /* tag          = */ 0,
        /* communicator = */ MPI_COMM_WORLD
        );
    }
    // Recieve from other partitions the number of cells they're sending
    for (auto& comm : comms) {
        uint this_size = 0;
        // Recieve values
        MPI_Recv(
        /* data         = */ &this_size, 
        /* count        = */ 1, 
        /* datatype     = */ MPI_UNSIGNED, 
        /* source       = */ comm.out_rank, 
        /* tag          = */ 0,
        /* communicator = */ MPI_COMM_WORLD,
        /* status       = */ MPI_STATUS_IGNORE
        );
        comm.snd_indices.resize(this_size);
        comm.snd_q.resize(vars*this_size);
    }



    // Send to other partitions the cells we want to recieve
    for (auto& comm : comms) {
        // Send values
        MPI_Send(
        /* data         = */ &comm.rec_indices[0], 
        /* count        = */ comm.rec_indices.size(), 
        /* datatype     = */ MPI_UNSIGNED,
        /* destination  = */ comm.out_rank, 
        /* tag          = */ 0,
        /* communicator = */ MPI_COMM_WORLD
        );
    }
    // Recieve other partitions cells
    for (auto& comm : comms) {
        
        // Recieve values
        MPI_Recv(
        /* data         = */ &comm.snd_indices[0], 
        /* count        = */ comm.snd_indices.size(), 
        /* datatype     = */ MPI_UNSIGNED, 
        /* source       = */ comm.out_rank, 
        /* tag          = */ 0,
        /* communicator = */ MPI_COMM_WORLD,
        /* status       = */ MPI_STATUS_IGNORE
        );
    }

    // These are the original indices, we must map them to current indices
    for (auto& comm : comms) {
        for (uint i=0; i<comm.rec_indices.size(); ++i) {
            comm.rec_indices[i] = originalCellsRefInv.at(comm.rec_indices[i]);
        }
        for (uint i=0; i<comm.snd_indices.size(); ++i) {
            comm.snd_indices[i] = originalCellsRefInv.at(comm.snd_indices[i]);
        }
    }
}


}

