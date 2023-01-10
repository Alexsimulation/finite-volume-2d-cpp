
/* Gen mesh with:
    First convert to new version:
        gmsh naca0012-129.msh convert.geo -parse_and_exit
    Then partition:
        gmsh naca.msh convert.geo -setnumber n 4 -parse_and_exit
*/



// Should we create the boundary representation of the partition entities?
Mesh.PartitionCreateTopology = 0;

// Should we create ghost cells?
Mesh.PartitionCreateGhostCells = 1;

// Should we automatically create new physical groups on the partition entities?
Mesh.PartitionCreatePhysicals = 0;

// Should we save one mesh file per partition?
Mesh.PartitionSplitMeshFiles = 1;

PartitionMesh n;

Save "naca.msh";
