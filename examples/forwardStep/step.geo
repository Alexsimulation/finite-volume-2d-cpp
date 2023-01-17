
// Gen mesh with:
// gmsh step.geo -setnumber n 6 -setstring outfile step.msh -parse_and_exit

//lc = 0.02;
lc = 0.02;

Point (1) = {-0.6, 0, 0, lc};
Point (2) = {0, 0, 0, lc*0.5};
Point (3) = {0, 0.2, 0, lc*0.5};
Point (4) = {2.4, 0.2, 0, lc};
Point (5) = {2.4, 1, 0, lc};
Point (6) = {0, 1, 0, lc};
Point (7) = {-0.6, 1, 0, lc};
Point (8) = {-0.6, 0.2, 0, lc};

Line (1) = {1, 2};
Line (2) = {2, 3};
Line (3) = {3, 4};
Line (4) = {4, 5};
Line (5) = {5, 6};
Line (6) = {6, 7};
Line (7) = {7, 8};
Line (8) = {8, 1};
//Line (9) = {3, 8};
//Line (10) = {3, 6};

Curve Loop (1) = {1, 2, 3, 4, 5, 6, 7, 8};
//Curve Loop (1) = {1, 2, 9, 8};
//Curve Loop (2) = {-9, 10, 6, 7};
//Curve Loop (3) = {3, 4, 5, -10};

Plane Surface (1) = {1};
//Plane Surface (2) = {2};
//Plane Surface (3) = {3};

Physical Curve("wall") = {1, 2, 3, 5, 6};
Physical Curve("inlet") = {7, 8};
Physical Curve("outlet") = {4};
Physical Surface("internal") = {1};

/*
Transfinite Curve {1, -9, -6} = 30 Using Progression 0.975;
Transfinite Curve {7, -10, -4} = 60 Using Progression 0.98;
Transfinite Curve {-3, 5} = 70 Using Progression 0.98;
Transfinite Curve {2, -8} = 20 Using Progression 0.98;

Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};

Recombine Surface (1);
Recombine Surface (2);
Recombine Surface (3);
*/

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
