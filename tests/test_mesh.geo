// Gen mesh with:
// gmsh test_mesh.geo -setnumber n 4 -setstring outfile test_mesh.msh -parse_and_exit
SetFactory("OpenCASCADE");

lc = 0.5;
Point (1) = {0, 0, 0, lc};
Point (2) = {1, 0, 0, lc};
Point (3) = {1, 1, 0, lc};

Circle (1) = {1, 1, 0, 1, 0, 2*Pi};

Curve Loop (1) = {1};
Plane Surface (1) = {1};

Physical Curve("wall") = {1};
Physical Surface("internal") = {1};

// Mesh the surface
Mesh 2;

// Should we create the boundary representation of the partition entities?
Mesh.PartitionCreateTopology = 1;

// Should we create ghost cells?
Mesh.PartitionCreateGhostCells = 1;

// Should we automatically create new physical groups on the partition entities?
Mesh.PartitionCreatePhysicals = 1;

// Should we save one mesh file per partition?
Mesh.PartitionSplitMeshFiles = 1;

PartitionMesh n;

Save Str(outfile);
