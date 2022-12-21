
// Gen mesh with:
// gmsh test.geo -setnumber n 6 -setstring outfile square.msh -parse_and_exit



Point (1) = {0, 0, 0, 0.025};
Point (2) = {1, 0, 0, 0.025};
Point (3) = {1, 1, 0, 0.025};
Point (4) = {0, 1, 0, 0.025};

Line (1) = {1, 2};
Line (2) = {2, 3};
Line (3) = {3, 4};
Line (4) = {4, 1};

Curve Loop (1) = {1, 2, 3, 4};

Plane Surface (1) = {1};

Physical Curve("wall") = {1, 2, 3, 4};
Physical Surface("internal") = {1};

Transfinite Surface {1};

Recombine Surface (1);

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
