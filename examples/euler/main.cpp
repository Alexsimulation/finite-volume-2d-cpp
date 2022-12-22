#include <fvhyper/mesh.h>
#include <fvhyper/explicit.h>
#include <fvhyper/parallel.h>
#include <fvhyper/post.h>


// Implement the sod shock tube problem using fvhyper
namespace fvhyper {


    // Define global constants
    int vars = 4;
    namespace solver {
        bool do_calc_gradients = false;
        bool do_calc_limiters = false;
    }

    namespace consts {
        double gamma = 1.4;
    }

    /*
        Define initial solution
        Here, we have a non-uniform initial condition corresponding to
        the sod shock tube problem
        https://en.wikipedia.org/wiki/Sod_shock_tube
    */
    void generate_initial_solution(
        std::vector<double>& v,
        const mesh& m
    ) {
        for (uint i=0; i<m.cellsAreas.size(); ++i) {
            if (m.cellsCentersX[i] < 0.5) {
                v[4*i] = 1.;
                v[4*i+1] = 0.;
                v[4*i+2] = 0.;
                v[4*i+3] = 1./(consts::gamma-1);
            } else {
                v[4*i] = 0.125;
                v[4*i+1] = 0.;
                v[4*i+2] = 0.;
                v[4*i+3] = 0.1/(consts::gamma-1);
            }
        }
    }

    // Van albada 1 limiter function
    void limiter_func(double* l, const double* r) {
        for (uint i=0; i<4; ++i) {
            l[i] = (r[i]*r[i] + r[i])/(r[i]*r[i] + 1.);
        }
    }

    // Helper function for entropy correction
    inline double entropy_correction(double l, double d) {
        if (l > d) { return l; }
        else { return (l*l + d*d)/(2.*d); }
    }

    /*
        Define flux function for the euler equations
        Roe flux vector differencing
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
        // Central flux
        double pi, pj, Vi, Vj;
        pi = (consts::gamma - 1)*(qi[3] - 0.5/qi[0]*(qi[1]*qi[1] + qi[2]*qi[2]));
        pj = (consts::gamma - 1)*(qj[3] - 0.5/qj[0]*(qj[1]*qj[1] + qj[2]*qj[2]));
        Vi = (qi[1]*n[0] + qi[2]*n[1])/qi[0];
        Vj = (qj[1]*n[0] + qj[2]*n[1])/qj[0];

        f[0] = qi[0]*Vi;
        f[1] = qi[1]*Vi + pi*n[0];
        f[2] = qi[2]*Vi + pi*n[1];
        f[3] = (qi[3] + pi)*Vi;

        f[0] += qj[0]*Vj;
        f[1] += qj[1]*Vj + pj*n[0];
        f[2] += qj[2]*Vj + pj*n[1];
        f[3] += (qj[3] + pj)*Vj;

        for (uint i=0; i<4; ++i) f[i] *= 0.5;

        // Upwind flux
        const double pL = pi;
        const double pR = pj;

        // Roe variables
        const double uL = qi[1]/qi[0];
        const double uR = qj[1]/qj[0];
        const double vL = qi[2]/qi[0];
        const double vR = qj[2]/qj[0];

        const double srhoL = sqrt(qi[0]);
        const double srhoR = sqrt(qj[0]);
        const double rho = srhoR*srhoL;
        const double u = (uL*srhoL + uR*srhoR)/(srhoL + srhoR);
        const double v = (vL*srhoL + vR*srhoR)/(srhoL + srhoR);
        const double h = ((qi[3] + pL)/qi[0]*srhoL + (qj[3] + pR)/qj[0]*srhoR)/(srhoL + srhoR);
        const double q2 = u*u + v*v;
        const double c = sqrt( (consts::gamma - 1.) * (h - 0.5*q2) );
        const double V = u*n[0] + v*n[1];
        const double VR = uR*n[0] + vR*n[1];
        const double VL = uL*n[0] + vL*n[1];

        const double delta = 0.05*c;
        double lambda_cm = entropy_correction(std::abs(V-c), delta);
        double lambda_c  = entropy_correction(std::abs(V), delta);
        double lambda_cp = entropy_correction(std::abs(V+c), delta);

        const double kF1 = lambda_cm*((pR-pL) - rho*c*(VR-VL))/(2.*c*c);
        const double kF234_0 = lambda_c*((qj[0] - qi[0]) - (pR-pL)/(c*c));
        const double kF234_1 = lambda_c*rho;
        const double kF5 = lambda_cp*((pR-pL) + rho*c*(VR-VL))/(2*c*c);

        // Roe flux

        f[0] -= 0.5*(kF1            + kF234_0                                                       + kF5);
        f[1] -= 0.5*(kF1*(u-c*n[0]) + kF234_0*u      + kF234_1*(uR - uL - (VR-VL)*n[0])             + kF5*(u+c*n[0]));
        f[2] -= 0.5*(kF1*(v-c*n[1]) + kF234_0*v      + kF234_1*(vR - vL - (VR-VL)*n[1])             + kF5*(v+c*n[1]));
        f[3] -= 0.5*(kF1*(h-c*V)    + kF234_0*q2*0.5 + kF234_1*(u*(uR-uL) + v*(vR-vL) - V*(VR-VL))  + kF5*(h+c*V)); 
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
        dt = 2e-5;
    }


    namespace boundaries {

        /*
            Define boundary conditions
            We apply a neumann zero flux condition to all boundaries
        */
        void zero_flux(double* b, double* q, double* n) {
            b[0] = q[0];
            b[1] = q[1];
            b[2] = q[2];
            b[3] = q[3];
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
    options.max_step = 10000;
    options.max_time = 0.2;
    options.print_interval = 100;

    // Run solver
    std::vector<double> q;
    fvhyper::run(q, pool, m, options);

    // Save file
    fvhyper::writeVtk(name, {"rho", "rhou", "rhov", "rhoe"}, q, m, pool.rank, pool.size);

    return pool.exit();
}

