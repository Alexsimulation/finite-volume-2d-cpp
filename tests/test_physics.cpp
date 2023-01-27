/*

       ___     __                    
      / _/  __/ /  __ _____  ___ ____
     / _/ |/ / _ \/ // / _ \/ -_) __/
    /_/ |___/_//_/\_, / .__/\__/_/   
                 /___/_/             

    Finite Volumes for High Performance

    - Description : Tests for solver edits validation
    - Author : Alexis Angers
    - Contact : alexis.angers@polymtl.ca

*/
#include <fvhyper/test.h>
#include <fvhyper/mesh.h>
#include <fvhyper/physics.h>
#include <fvhyper/solver.h>
#include <fvhyper/explicit.h>
#include <fvhyper/parallel.h>
#include <fvhyper/post.h>


#include <map>
#include <string>



class heat : public fvhyper::physics {

public:
    heat() : physics(1) {

        // Set a fixed dt
        set_fixed_dt(1e-4);

        // Set the variable names
        var_names = {"T"};

        // Set extra scalars
        extra_scalars = {"2T"};

        // Enable the source term calculation
        do_source_term = true;
    }

    void flux(
        double* flux,
        const double* qi,
        const double* qj,
        const double* gx,
        const double* gy,
        const double* n
    ) {
        flux[0] = -1.0*(gx[0]*n[0] + gy[0]*n[1]);
    }

    void source(
        double* source,
        const double* q
    ) {
        source[0] = 1.0;
    }

    void boundary(
        std::string& name,
        double* b,
        double* q,
        double* n
    ) {
        // Check which boundary we're at
        if (name == "wall") {
            b[0] = 0.;
        }
    }

    void init(
        double* q,
        double& x,
        double& y
    ) {
        double r2 = x*x + y*y;
        q[0] = 0;
    }


    void calculate_extra_scalars(double* s, double* q, std::string& name) {
        if (name == "2T") {
            s[0] = 2.0*q[0];
        }
    }

};





class burgers : public fvhyper::physics {

public:
    burgers() : physics(2) {

        // Set a fixed dt
        set_fixed_dt(1e-4);

        // Set the variable names
        var_names = {"u", "v"};

        // Set extra vectors
        extra_vectors = {"U"};

        // Enable the source term calculation
        do_source_term = false;
    }

    void flux(
        double* flux,
        const double* qi,
        const double* qj,
        const double* gx,
        const double* gy,
        const double* n
    ) {
        double Vi = qi[0]*n[0] + qi[1]*n[1];
        double Vj = qj[0]*n[0] + qj[1]*n[1];

        flux[0] = 0.5*(Vi*qi[0] + Vj*qj[0]);
        flux[1] = 0.5*(Vi*qi[1] + Vj*qj[1]);

        double V = std::max(abs(Vi), abs(Vj));

        flux[0] -= 0.5*V*(qj[0] - qi[0]);
        flux[1] -= 0.5*V*(qj[1] - qi[1]);
    }

    void boundary(
        uint& id,
        double* b,
        double* q,
        double* n
    ) {
        // Check which boundary we're at
        if (id == 0) {
            b[0] = 1.;
            b[1] = 0.;
        }
    }

    void init(
        double* q,
        double& x,
        double& y
    ) {
        q[0] = 0;
        q[1] = 0;
    }


    void calculate_extra_vectors(double* v, double* q, std::string& name) {
        if (name == "U") {
            v[0] = q[0];
            v[1] = q[1];
        }
    }

};






int main() {

    heat pde_0;
    burgers pde_1;

    // Problem files name
    std::string name = "square";

    // Initialize the MPI environment
    fvhyper::mpi_wrapper pool;
    
    // Read the file
    pde_0.m.read_file(name, pool, pde_0);

    // Solver options
    fvhyper::solverOptions options;
    options.max_step = 5000;
    options.max_time = 10.;
    options.print_interval = 500;

    fvhyper::solver solver;
    solver.opt = options;

    // Run solver
    solver.run();

    // Save file
    fvhyper::writeVtk(name, pde_0, pde_0.s.q, pde_0.m, pool.rank, pool.size);

    return pool.exit();
}



