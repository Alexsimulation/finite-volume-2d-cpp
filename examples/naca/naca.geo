
/*
    first convert with:
        gmsh naca0012_129.msh -o naca.msh -save
    Then partition:
        gmsh naca.msh naca.geo -setnumber n 6 -parse_and_exit
*/


// Should we create the boundary representation of the partition entities?
Mesh.PartitionCreateTopology = 0;

// Should we create ghost cells?
Mesh.PartitionCreateGhostCells = 1;

// Should we automatically create new physical groups on the partition entities?
Mesh.PartitionCreatePhysicals = 0;

// Should we save one mesh file per partition?
Mesh.PartitionSplitMeshFiles = 1;

//PartitionMesh n;

Save "naca.msh";

