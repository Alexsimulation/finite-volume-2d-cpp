/*

       ___     __                    
      / _/  __/ /  __ _____  ___ ____
     / _/ |/ / _ \/ // / _ \/ -_) __/
    /_/ |___/_//_/\_, / .__/\__/_/   
                 /___/_/             

    Finite Volumes for High Performance

    - Description : Post-processing sources
    - Author : Alexis Angers
    - Contact : alexis.angers@polymtl.ca

*/
#include <fstream>
#include <fvhyper/post.h>

#include <iomanip>
#include <sstream>

namespace fvhyper {


std::string double2string(const double& x, const int precision) {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(precision) << x;
    return stream.str();
}


void writeVtk(
    const std::string name,
    std::vector<double>& q,
    mesh& m,
    int rank,
    int world_size,
    const std::string time
) {
    std::string dash_time = (time == "") ? "" : "_" + time;

    std::string filename =
        (world_size>1)
        ? name + "_" + std::to_string(rank) + dash_time + ".vtu"
        : name + dash_time + ".vtu";
    
    if (time != "") {
        filename = "times/" + filename;
    }

    if ((rank == 0)&(world_size > 1)) {
        std::string coreFileName = name + "_parallel" + dash_time + ".pvtu";
        if (time != "") {
            coreFileName = "times/" + coreFileName;
        }
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
        for (auto varname : var_names) {
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
        if (m.do_compute_wall_dist) {
            core += "    <PDataArray type=\"Float32\" Name=\"Wall Distance\"/>\n";
        }
        core += "  </PCellData>\n";
        for (uint i=0; i<world_size; ++i) {
            std::string filename_i = name + "_" + std::to_string(i) + dash_time + ".vtu";
            core += "  <Piece Source=\"" + filename_i + "\"/>\n";
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

    // Save wall distance
    if (m.do_compute_wall_dist) {
        s += "        <DataArray type=\"Float32\" Name=\"Wall Distance\" Format=\"ascii\">\n";
        for (int j=0; j<m.nRealCells; ++j) {
            s += "          " + double2string(m.wall_dist[j], 16) + "\n";
        }
        s += "        </DataArray>\n";
    }

    // Save variables
    for (uint i=0; i<var_names.size(); ++i) {
        s += "        <DataArray type=\"Float32\" Name=\"" + var_names[i] + "\" Format=\"ascii\">\n";
        for (int j=0; j<m.nRealCells; ++j) {
            s += "          " + double2string(q[var_names.size()*j + i], 16) + "\n";
        }
        s += "        </DataArray>\n";
    }
    for (auto& keyval : post::extra_scalars) {
        s += "        <DataArray type=\"Float32\" Name=\"" + keyval.first + "\" Format=\"ascii\">\n";
        for (int j=0; j<m.nRealCells; ++j) {
            double scal_i[1];
            keyval.second(scal_i, &q[var_names.size()*j]);
            s += "          " + double2string(scal_i[0], 16) + "\n";
        }
        s += "        </DataArray>\n";
    }

    // Save vectors
    if (post::extra_vectors.size() > 0) {
        for (auto& keyval : post::extra_vectors) {
            s += "        <DataArray type=\"Float32\" Name=\"" + keyval.first + "\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
            for (int j=0; j<m.nRealCells; ++j) {
                double vec_i[2];
                keyval.second(vec_i, &q[var_names.size()*j]);
                s += "          " + double2string(vec_i[0], 16) + " ";
                s += double2string(vec_i[1], 16) + " 0.0\n";
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
