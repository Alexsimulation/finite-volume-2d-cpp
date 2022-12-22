#include <fvhyper/mesh.h>
#include <fvhyper/explicit.h>
#include <fvhyper/parallel.h>
#include <fvhyper/post.h>


// Implement a heat tranfer problem using fvhyper
namespace fvhyper {

    // Define global constants
    int vars = 1;
    namespace solver {
        bool do_calc_gradients = false;
        bool do_calc_limiters = false;
    }

    /* 
        Define initial solution
        Here, we have a non-uniform initial condition
            u = 1   in the circle with center (0.5, 0.5) and radius 0.4
            u = 0   elsewhere
    */
    void generate_initial_solution(
        std::vector<double>& v,
        const mesh& m
    ) {
        for (uint i=0; i<m.cellsAreas.size(); ++i) {
            double x = m.cellsCentersX[i] - 0.5;
            double y = m.cellsCentersY[i] - 0.5;
            if (sqrt(x*x + y*y) < 0.4*0.4) {
                v[1*i] = 1.;
            } else {
                v[1*i] = 0.;
            }
        }
    }

    // Here the limiter function is not used, default implementation
    void limiter_func(double* l, const double* r) {
        l[0] = 1.;
    }

    /*
        Define flux function
        The pde to solve is:
            u_t = div(grad(u))
        In finite volume form:
            u_t - div(grad(u)) = 0
            int(u_t)dV - int(grad(u)*n)dS
        So the flux function is:
            F(u) = -grad(u)*n
        With n the normal vector.
    */
    void calc_flux(
        double* f,
        const double* qi,
        const double* qj,
        const double* gxi,
        const double* gyi,
        const double* gxj,
        const double* gyj,
        const double* l,
        const double* n,
        const double* di,
        const double* dj,
        const double area,
        const double len
    ) {
        // Diffusion equation
        double grad[2];
        grad[0] = (qj[0] - qi[0])*n[0]*len/area;
        grad[1] = (qj[0] - qi[0])*n[1]*len/area;

        f[0] = -(grad[0]*n[0] + grad[1]*n[1]);   
    }

    /*
        Define the time step. Here, time step is constant
    */
    void calc_dt(
        double& dt,
        const std::vector<double>& q,
        mesh& m
    ) {
        // Constant time step
        dt = 1e-5;
    }


    namespace boundaries {

        /*
            Define boundary conditions
            We apply a neumann zero flux condition to all boundaries
        */
        void zero_flux(double* b, double* q, double* n) {
            b[0] = q[0];
        }
        std::map<std::string, void (*)(double*, double*, double*)> 
        bounds = {
            {"wall", zero_flux}
        };
    }


}



int main() {
    // Initialize the MPI environment
    fvhyper::mpi_wrapper pool;

    // Create mesh object m
    fvhyper::mesh m;

    /*
        Define problem files name
        mesh files must be name according to:
            name_{rank+1}.msh
        where rank is the mpi rank associated with this mesh process
        per example, for 3 mpi ranks, we would have the files:
            name_1.msh  name_2.msh  name_3.msh
    */
    std::string name = "square";

    // Read the file
    m.read_file(name, pool);

    fvhyper::solverOptions options;
    options.max_step = 1000;
    options.max_time = 0.01;
    options.print_interval = 100;

    // Run solver
    std::vector<double> q;
    fvhyper::run(q, pool, m, options);

    // Save file
    fvhyper::writeVtk(name, {"u"}, q, m, pool.rank, pool.size);

    return pool.exit();
}

