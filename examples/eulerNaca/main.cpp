#include <fvhyper/mesh.h>
#include <fvhyper/explicit.h>
#include <fvhyper/parallel.h>
#include <fvhyper/post.h>


// Implement the sod shock tube problem using fvhyper
namespace fvhyper {


    // Define global constants
    const int vars = 4;
    namespace solver {
        const bool do_calc_gradients = true;
        const bool do_calc_limiters = true;
        const bool linear_interpolate = true;
        const bool diffusive_gradients = false;
        const bool global_dt = false;
    }

    namespace consts {
        double gamma = 1.4;
        double cfl = 1.25;
    }

    // Helper function for pressure calc
    inline double calc_p(const double* q) {
        return (consts::gamma - 1)*(q[3] - 0.5/q[0]*(q[1]*q[1] + q[2]*q[2]));
    }

    /*
        Define initial solution
    */
    void generate_initial_solution(
        std::vector<double>& v,
        const mesh& m
    ) {
        for (uint i=0; i<m.cellsAreas.size(); ++i) {
            v[4*i] = 1.4;
            v[4*i+1] = 1.4 * 0.8;
            v[4*i+2] = 1.4 * 0.0175;
            v[4*i+3] = 1.0/(consts::gamma-1) + 0.5*1.4*(0.8*0.8 + 0.0175*0.0175);
        }
    }

    // Michalak limiter function
    double limiter_func(const double& y) {
        const double yt = 2.0;
        if (y >= yt) {
            return 1.0;
        } else {
            const double a = 1.0/(yt*yt) - 2.0/(yt*yt*yt);
            const double b = -3.0/2.0*a*yt - 0.5/yt;
            return a*y*y*y + b*y*y + y;
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
        const double* gx,
        const double* gy,
        const double* n
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
        std::vector<double>& dt,
        const std::vector<double>& q,
        mesh& m
    ) {

        double cfl = consts::cfl;

        // Set max dt
        for (uint i=0; i<dt.size(); ++i) {
            dt[i] = 1.0;
        }
        // Calculate time step
        for (uint e=0; e<m.edgesNodes.cols(); ++e) {
            double n[2];

            const auto& i = m.edgesCells(e, 0);
            const auto& j = m.edgesCells(e, 1);
            const auto& le = m.edgesLengths[e];
            
            n[0] = m.edgesNormalsX[e];
            n[1] = m.edgesNormalsY[e];

            // From cell i to cell j
            // Compute max eigenvalue
            const double* qi = &q[4*i];
            const double* qj = &q[4*j];
            double ci = sqrt(calc_p(qi)*consts::gamma / qi[0]);
            double cj = sqrt(calc_p(qj)*consts::gamma / qj[0]);

            double max_eig_i = ci + abs(n[0]*qi[1]/qi[0] + n[1]*qi[2]/qi[0]);
            double max_eig_j = cj + abs(n[0]*qj[1]/qj[0] + n[1]*qj[2]/qj[0]);

            double dt_i = cfl * (m.cellsAreas[i] / le) / max_eig_i;
            double dt_j = cfl * (m.cellsAreas[j] / le) / max_eig_j;

            for (uint k=0; k<vars; ++k) {
                dt[vars*i + k] = std::min(dt[4*i + k], dt_i);
                dt[vars*j + k] = std::min(dt[4*j + k], dt_j);
            }
        }
    }


    namespace boundaries {

        /*
            Define boundary conditions
            We apply a neumann zero flux condition to all boundaries
        */
        void farfield(double* b, double* q, double* n) {
            
            double b_pressure = 1.0;
            double bv[4];
            bv[0] = 1.4;
            bv[1] = bv[0] * 0.8;
            bv[2] = bv[0] * 0.0175;
            bv[3] = b_pressure / (consts::gamma - 1) + 0.5/bv[0]*(bv[1]*bv[1] + bv[2]*bv[2]);

            double U[2];
            U[0] = q[1]/q[0];
            U[1] = q[2]/q[0];

            double u_norm = sqrt(U[0]*U[0] + U[1]*U[1]);
            double u_dot_n = U[0]*n[0] + U[1]*n[1];
            double p = calc_p(q);
            double c = sqrt(consts::gamma*p/q[0]);
            double mach = u_norm / c;

            if (mach > 1) {
                // Supersonic
                if (u_dot_n < 0) {
                    // Inlet
                    b[0] = bv[0];
                    b[1] = bv[1];
                    b[2] = bv[2];
                    b[3] = bv[3];
                } else {
                    // Outlet
                    b[0] = q[0];
                    b[1] = q[1];
                    b[2] = q[2];
                    b[3] = q[3];
                }
            } else {
                // Subsonic
                // a variables -> outside
                double pa = calc_p(bv);
                double rhoa = bv[0];
                double ua = bv[1] / bv[0];
                double va = bv[2] / bv[0];
                
                // d variables -> inside
                double pd = p;
                double rhod = q[0];
                double ud = q[1] / q[0];
                double vd = q[2] / q[0];

                double rho0 = q[0];
                double c0 = c;

                if (u_dot_n < 0) {
                    // Inlet
                    double pb = 0.5*(pa + pd - rho0*c0*(n[0]*(ua - ud) + n[1]*(va - vd)));
                    b[0] = rhoa + (pb - pa)/(c0*c0);
                    b[1] = b[0] * (ua - n[0]*(pa - pb)/(rho0*c0));
                    b[2] = b[1] * (va - n[1]*(pa - pb)/(rho0*c0));
                    b[3] = pb / (consts::gamma - 1) + 0.5/b[0]*(b[1]*b[1] + b[2]*b[2]);
                } else {
                    // Outlet
                    double pb = pa;
                    b[0] = rhod + (pb - pd)/(c0*c0);
                    b[1] = b[0] * (ud + n[0]*(pd - pb)/(rho0*c0));
                    b[2] = b[1] * (va + n[1]*(pd - pb)/(rho0*c0));
                    b[3] = pb / (consts::gamma - 1) + 0.5/b[0]*(b[1]*b[1] + b[2]*b[2]);
                }
            }
        }
        void wall(double* b, double* q, double* n) {

            double U[2];
            U[0] = q[1]/q[0];
            U[1] = q[2]/q[0];

            // Flip velocity, doesn't change norm of velocity
            U[0] = U[0] - 2.0 * n[0] * (n[0]*U[0] + n[1]*U[1]);
            U[1] = U[1] - 2.0 * n[1] * (n[0]*U[0] + n[1]*U[1]);

            b[0] = q[0];
            b[1] = U[0] * q[0];
            b[2] = U[1] * q[0];
            b[3] = q[3];
        }
        std::map<std::string, void (*)(double*, double*, double*)> 
        bounds = {
            {"airfoil", wall},
            {"farfield", farfield}
        };
    }


    /*
        Define extra output variables
    */
    namespace post {
        void calc_output_u(double* u, double* q) {
            // Compute vector u
            u[0] = q[1] / q[0];
            u[1] = q[2] / q[0];
        }

        void calc_output_p(double* p, double* q) {
            // Compute pressure p
            p[0] = calc_p(q);
        }

        std::map<std::string, void (*)(double*, double*)> 
            extra_scalars = {
                {"p", calc_output_p}
            };
        
        std::map<std::string, void (*)(double*, double*)> 
            extra_vectors = {
                {"U", calc_output_u}
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
    std::string name = "naca";

    // Read the file
    m.read_file(name, pool);

    fvhyper::solverOptions options;
    options.max_step = 15000;
    options.print_interval = 10;
    options.tolerance = 1e-12;

    // Run solver
    std::vector<double> q;
    fvhyper::run(q, pool, m, options);

    // Save file
    fvhyper::writeVtk(name, {"rho", "rhou", "rhov", "rhoe"}, q, m, pool.rank, pool.size);

    return pool.exit();
}

