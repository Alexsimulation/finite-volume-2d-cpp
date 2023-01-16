#include <fstream>
#include <fvhyper/post.h>

namespace fvhyper {

void writeVtk(
    const std::string name,
    std::vector<std::string> varNames,
    std::vector<double>& q,
    mesh& m,
    int rank,
    int world_size
) {

    std::string filename =
        (world_size>1)
        ? name + "_" + std::to_string(rank) + ".vtu"
        : name + ".vtu";

    if ((rank == 0)&(world_size > 1)) {
        std::string coreFileName = name + "_parallel.pvtu";
        // Write vtk header file
        std::string core = "";
        core += "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n";
        core += "<PUnstructuredGrid GhostLevel=\"1\">\n";
        core += "  <PPoints>\n";
        core += "    <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n";
        core += "  </PPoints>\n";
        core += "  <PCells>\n";
        core += "    <PDataArray type=\"Int32\" Name=\"connectivity\"/>\n";
        core += "    <PDataArray type=\"Int32\" Name=\"offsets\"/>\n";
        core += "    <PDataArray type=\"Int32\" Name=\"types\"/>\n";
        core += "  </PCells>\n";
        core += "  <PCellData Scalars=\"scalars\">\n";
        core += "    <PDataArray type=\"Float32\" Name=\"rank\"/>\n";
        for (auto varname : varNames) {
            core += "    <PDataArray type=\"Float32\" Name=\"" + varname + "\"/>\n";
        }
        for (auto& keyval : post::extra_scalars) {
            core += "    <PDataArray type=\"Float32\" Name=\"" + keyval.first + "\"/>\n";
        }
        if (post::extra_vectors.size() > 0) {
            for (auto& keyval : post::extra_vectors) {
                core += "    <PDataArray type=\"Float32\" Name=\"" + keyval.first + "\" NumberOfComponents=\"3\"/>\n";
            }
        }
        core += "  </PCellData>\n";
        for (uint i=0; i<world_size; ++i) {
            core += "  <Piece Source=\"" + name + "_" + std::to_string(i) + ".vtu\"/>\n";
        }
        core += "</PUnstructuredGrid>\n";
        core += "</VTKFile>";

        std::ofstream out(coreFileName);
        out << core;
        out.close();
    }
    // Write solution q to .vtk file filename
    std::string s = "";

    s += "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n";
    s += "  <UnstructuredGrid>\n";
    s += "    <Piece NumberOfPoints=\"" + std::to_string(m.nodesX.size()) + "\" NumberOfCells=\"" + std::to_string(m.nRealCells) + "\">\n";
    s += "      <Points>\n";
    s += "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
    for (uint i=0; i<m.nodesX.size(); ++i) {
        s += "          " + std::to_string(m.nodesX[i]) + " " + std::to_string(m.nodesY[i]) + " 0.0\n";
    }
    s += "        </DataArray>\n";
    s += "      </Points>\n";

    // Figure out the total number of integer tags in the cells
    s += "      <Cells>\n";
    s += "        <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n";

    // Add each cell
    s += "          ";
    for (int fi=0; fi<m.nRealCells; ++fi) {

        std::vector<uint> ns = {
            m.cellsNodes(fi, 0),
            m.cellsNodes(fi, 1),
            m.cellsNodes(fi, 2)
        };
        if (!m.cellsIsTriangle[fi]) {
            ns.push_back(m.cellsNodes(fi, 3));
        }
        for (auto ni : ns) {
            s += std::to_string(ni) + "  ";
        }
        s += "\n          ";
    }
    s += "        </DataArray>\n";

    // Save cell offsets
    s += "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n";
    s += "          ";
    uint current_offset = 0;
    for (uint i=0; i<m.nRealCells; ++i) {
        current_offset += m.cellsIsTriangle[i] ? 3 : 4;
        s += std::to_string(current_offset) + "  ";
    }
    s += "\n";
    s += "        </DataArray>\n";

    // Write cell types
    s += "        <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n";
    s += "          ";
    for (uint i=0; i<m.nRealCells; ++i) {
        bool isTriangle = m.cellsIsTriangle[i];
        if (isTriangle) {
            // Tri cell
            s += "5  ";
        } else {
            // Quad cell
            s += "9  ";
        }
    }
    s += "\n";
    s += "        </DataArray>\n";
    s += "      </Cells>\n";

    // Save scalar data
    s += "      <CellData Scalars=\"scalars\">\n";

    // Save cell rank
    s += "        <DataArray type=\"Float32\" Name=\"rank\" Format=\"ascii\">\n";
    for (int j=0; j<m.nRealCells; ++j) {
        s += "          " + std::to_string(rank) + "\n";
    }
    s += "        </DataArray>\n";

    // Save variables
    for (uint i=0; i<varNames.size(); ++i) {
        s += "        <DataArray type=\"Float32\" Name=\"" + varNames[i] + "\" Format=\"ascii\">\n";
        for (int j=0; j<m.nRealCells; ++j) {
            s += "          " + std::to_string(q[varNames.size()*j + i]) + "\n";
        }
        s += "        </DataArray>\n";
    }
    for (auto& keyval : post::extra_scalars) {
        s += "        <DataArray type=\"Float32\" Name=\"" + keyval.first + "\" Format=\"ascii\">\n";
        for (int j=0; j<m.nRealCells; ++j) {
            double scal_i[1];
            keyval.second(scal_i, &q[varNames.size()*j]);
            s += "          " + std::to_string(scal_i[0]) + "\n";
        }
        s += "        </DataArray>\n";
    }

    // Save vectors
    if (post::extra_vectors.size() > 0) {
        for (auto& keyval : post::extra_vectors) {
            s += "        <DataArray type=\"Float32\" Name=\"" + keyval.first + "\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
            for (int j=0; j<m.nRealCells; ++j) {
                double vec_i[2];
                keyval.second(vec_i, &q[varNames.size()*j]);
                s += "          " + std::to_string(vec_i[0]) + " ";
                s += std::to_string(vec_i[1]) + " 0.0\n";
            }
            s += "        </DataArray>\n";
        }
    }

    s += "      </CellData>\n";

    s += "    </Piece>\n";
    s += "  </UnstructuredGrid>\n";
    s += "</VTKFile>";

    // Save data to file
    std::ofstream out(filename);
    out << s;
    out.close();

}



}
