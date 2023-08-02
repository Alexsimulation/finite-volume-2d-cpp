#include <fvhyper/mesh.h>
#include <fvhyper/explicit.h>
#include <fvhyper/parallel.h>
#include <fvhyper/post.h>

/*
    Implementation of the flat plate laminar boundary layer using fvhyper
*/
namespace fvhyper {


    // Define global constants
    const int vars = 5;
    const std::vector<std::string> var_names = {
        "rho",
        "rhou",
        "rhov",
        "rhoe",
        "nu_scaled"
    };
    const std::vector<double> vars_limiters = {
        1., 1., 1., 1., 0.
    };
    namespace solver {
        const bool do_calc_gradients = true;
        const bool do_calc_limiters = true;
        const bool linear_interpolate = true;
        const bool diffusive_gradients = true;
        const bool global_dt = false;
        const bool smooth_residuals = false;
        const bool source_term = true;
    }

    namespace consts {
        double gamma = 1.4;
        double cfl = 2.0;
        double mu = 2.8e-7;
        double pr = 0.72;
        double cp = 1;
        double r = 8.3145;

        double pr_T = 0.9;

        double nu_scale = 1e-6;
    }

    // Helper function for pressure calc
    inline double calc_p(const double* q) {
        return (consts::gamma - 1)*(q[3] - 0.5/q[0]*(q[1]*q[1] + q[2]*q[2]));
    }
    // Helper function for pressure gradient calc
    void calc_grad_p(double* gp, const double* q, const double* gx, const double* gy) {
        // p = (g-1)*rhoe - (g-1)*0.5/rho*(rhou*rhou + rhov*rhov))
        // p = (g-1)*rhoe - (g-1)*0.5*rhou*rhou/rho - (g-1)*0.5/rho*rhov*rhov
        // dp = 0.5*(g-1)*(2*drhoe - rhou*rhou/rho - rhov*rhov/rho)
        gp[0] = 2.*gx[3];
        gp[0] -= gx[1]*q[1]/q[0] + q[1]*(gx[1]*q[0] - gx[0]*q[1])/(q[0]*q[0]);
        gp[0] -= gx[2]*q[2]/q[0] + q[2]*(gx[2]*q[0] - gx[0]*q[2])/(q[0]*q[0]);
        gp[0] *= 0.5*(consts::gamma - 1);

        gp[1] = 2.*gy[3];
        gp[1] -= gy[1]*q[1]/q[0] + q[1]*(gy[1]*q[0] - gy[0]*q[1])/(q[0]*q[0]);
        gp[1] -= gy[2]*q[2]/q[0] + q[2]*(gy[2]*q[0] - gy[0]*q[2])/(q[0]*q[0]);
        gp[1] *= 0.5*(consts::gamma - 1);
    }
    // Helper function for temperature calc
    inline double calc_t(const double* q) {
        // p = rho * r * t
        // t = (1/r) * p/rho
        // t = (consts::gamma - 1)/r*(q[3]/q[0] - 0.5/(q[0]*q[0])*(q[1]*q[1] + q[2]*q[2]))
        return calc_p(q)/(q[0]*consts::r);
    }
    // Helper function for temperature gradient calc
    void calc_grad_t(double* gt, const double* q, const double* gx, const double* gy) {
        // t = (1/r) * p/rho
        // dt = (1/r)*( (dp*rho - drho*p)/(rho*rho) )
        double gp[2];
        const double p = calc_p(q);
        calc_grad_p(gp, q, gx, gy);

        gt[0] = (1./consts::r)*(
            (gp[0]*q[0] - gx[0]*p)/(q[0]*q[0])
        );
        gt[1] = (1./consts::r)*(
            (gp[1]*q[0] - gy[0]*p)/(q[0]*q[0])
        );
    }
    // Compute gradient of velocity u
    void calc_grad_u(double* gu, const double* q, const double* gx, const double* gy) {
        // u = rhou/rho
        // du = (rho*drhou - rhou*drho)/(rho*rho)
        gu[0] = (q[0]*gx[1] - q[1]*gx[0])/(q[0]*q[0]);
        gu[1] = (q[0]*gy[1] - q[1]*gy[0])/(q[0]*q[0]);
    }
    // Compute gradient of velocity v
    void calc_grad_v(double* gv, const double* q, const double* gx, const double* gy) {
        // v = rhov/rho
        // dv = (rho*drhov - rhov*drho)/(rho*rho)
        gv[0] = (q[0]*gx[2] - q[2]*gx[0])/(q[0]*q[0]);
        gv[1] = (q[0]*gy[2] - q[2]*gy[0])/(q[0]*q[0]);
    }


    /*
        Define initial solution
    */
    void generate_initial_solution(
        std::vector<double>& q,
        const mesh& m
    ) {
        for (uint i=0; i<m.cellsAreas.size(); ++i) {
            q[vars*i] = 1.4;
            q[vars*i+1] = 1.4 * 0.2;
            q[vars*i+2] = 1.4 * 0.;
            q[vars*i+3] = 1.0/(consts::gamma-1) + 0.5*1.4*(0.2*0.2 + 0.0*0.0);
            q[vars*i+4] = (consts::mu / 1.4) * 0.1 / consts::nu_scale;
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

    /*
        Define source function for the ns-sa model
    */
    void calc_source(
        double* s,
        const double* q,
        const double* gx,
        const double* gy,
        const double& wall_dist
    ) {

        // Turbulence source term
        double gradu[2];
        double gradv[2];
        calc_grad_u(gradu, q, gx, gy);
        calc_grad_v(gradv, q, gx, gy);

        // Model constants
        const double Cb1 = 0.1355;
        const double Cb2 = 0.622;
        const double Cv1 = 7.1;
        const double Cv2 = 5;
        const double sigma = 2./3.;
        const double kappa = 0.41;
        const double Cw1 = Cb1/(kappa*kappa) + (1. + Cb2)/sigma;
        const double Cw2 = 0.3;
        const double Cw3 = 2;
        const double Ct1 = 1;
        const double Ct2 = 2;
        const double Ct3 = 1.3;
        const double Ct4 = 0.5;

        // Velocity gradient tensor
        const double S_xx = 0.5*(gradu[0] + gradu[0]);
        const double S_xy = 0.5*(gradu[1] + gradv[0]);
        const double S_yy = 0.5*(gradv[1] + gradv[1]);

        const double O_xx = 0.5*(gradu[0] - gradu[0]);
        const double O_xy = 0.5*(gradu[1] - gradv[0]);
        const double O_yy = 0.5*(gradv[1] - gradv[1]);

        // Other terms
        const double nu_L = consts::mu / q[0];
        const double nu_tilda = q[4] * consts::nu_scale;

        const double X = nu_tilda / nu_L;
        const double X3 = X*X*X;

        const double fv1 = (X3)/(X3 + Cv1*Cv1*Cv1);
        
        double fv2 = 1 + X/Cv2;
        fv2 = 1/(fv2*fv2*fv2);

        const double fv3 = (1 + X*fv1)*(1 - fv2)/std::max(X, 0.001);

        double S_tilda = fv3*sqrt(2*O_xy*O_xy) + nu_tilda/(kappa*kappa * wall_dist*wall_dist) * fv2;

        const double Cw3_6 = Cw3*Cw3*Cw3 * Cw3*Cw3*Cw3;
        const double r = nu_tilda/(S_tilda*kappa*kappa*wall_dist*wall_dist);
        const double g = r + Cw2*(r*r*r*r*r*r - r);
        double fw = (1 + Cw3_6)/(g*g*g * g*g*g + Cw3_6);
        fw = g * std::pow(fw, 1./6.);

        const double ft2 = Ct3*std::exp(-Ct4*X*X);
        
        // Source terms
        const double term_0 = Cb1*(1 - ft2)*S_tilda*nu_tilda;
        const double term_1 = Cb2/sigma*(gx[4]*gx[4] + gy[4]*gy[4])*consts::nu_scale*consts::nu_scale;
        const double term_2 = -(Cw1*fw - Cb1/(kappa*kappa)*ft2)*(nu_tilda*nu_tilda/(wall_dist*wall_dist));
        // Ignore production from the trip point

        s[0] = 0;
        s[1] = 0;
        s[2] = 0;
        s[3] = 0;
        s[4] = (term_0 + term_1 + term_2)/consts::nu_scale;
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
        f[4] = qi[4]*Vi;

        f[0] += qj[0]*Vj;
        f[1] += qj[1]*Vj + pj*n[0];
        f[2] += qj[2]*Vj + pj*n[1];
        f[3] += (qj[3] + pj)*Vj;
        f[4] += qj[4]*Vj;

        for (uint i=0; i<vars; ++i) f[i] *= 0.5;

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

        // Roe correction
        // From https://www.researchgate.net/publication/305638346_Cures_for_the_Expansion_Shock_and_the_Shock_Instability_of_the_Roe_Scheme
        const double lambda_cm = abs(std::min(V-c, VL-c));
        const double lambda_c  = abs(V);
        const double lambda_cp = abs(std::max(V+c, VR+c));

        const double kF1 = lambda_cm*((pR-pL) - rho*c*(VR-VL))/(2.*c*c);
        const double kF234_0 = lambda_c*((qj[0] - qi[0]) - (pR-pL)/(c*c));
        const double kF234_1 = lambda_c*rho;
        const double kF5 = lambda_cp*((pR-pL) + rho*c*(VR-VL))/(2*c*c);

        // Roe flux
        f[0] -= 0.5*(kF1            + kF234_0                                                       + kF5);
        f[1] -= 0.5*(kF1*(u-c*n[0]) + kF234_0*u      + kF234_1*(uR - uL - (VR-VL)*n[0])             + kF5*(u+c*n[0]));
        f[2] -= 0.5*(kF1*(v-c*n[1]) + kF234_0*v      + kF234_1*(vR - vL - (VR-VL)*n[1])             + kF5*(v+c*n[1]));
        f[3] -= 0.5*(kF1*(h-c*V)    + kF234_0*q2*0.5 + kF234_1*(u*(uR-uL) + v*(vR-vL) - V*(VR-VL))  + kF5*(h+c*V)); 
        f[4] -= 0.5*(abs(V)+c)*(qj[4] - qi[4]);

        // Gradients for viscous fluxes
        // Diffusion equation

        // Central values
        double qc[vars];
        for (uint i=0; i<vars; ++i) {
            qc[i] = 0.5*(qi[i] + qj[i]);
        }

        // SA things
        const double nu_L = consts::mu / qc[0];
        const double Cv1 = 7.1;

        const double nu_tilda = qc[4] * consts::nu_scale;
        const double X = nu_tilda / nu_L;
        const double X3 = X*X*X;

        const double fv1 = (X3)/(X3 + Cv1*Cv1*Cv1);
        const double nu_T = fv1 * nu_tilda;

        // Temperature gradient, null for now
        double gradT[2];
        calc_grad_t(gradT, qc, gx, gy);
        double gradu[2];
        double gradv[2];
        calc_grad_u(gradu, qc, gx, gy);
        calc_grad_v(gradv, qc, gx, gy);

        // Total viscosity
        const double mu_T = qc[0] * nu_T;
        const double mu_tot = consts::mu + mu_T;

        // Viscous stresses
        const double div_v = gradu[0] + gradv[1];
        const double tau_xx = 2. * mu_tot * (gradu[0] - div_v/3.);
        const double tau_yy = 2. * mu_tot * (gradv[1] - div_v/3.);
        const double tau_xy = mu_tot * (gradu[1] + gradv[0]);

        double phi[2];
        const double k = consts::cp * (consts::mu / consts::pr + mu_T / consts::pr_T);
        phi[0] = qc[1]/qc[0]*tau_xx + qc[2]/qc[0]*tau_xy + k*gradT[0];
        phi[1] = qc[1]/qc[0]*tau_xy + qc[2]/qc[0]*tau_yy + k*gradT[1];

        // Turbulence
        const double sigma = 2./3.;

        const double tau_xx_T = 1/sigma*(nu_L + nu_tilda)*(gx[4]*consts::nu_scale);
        const double tau_yy_T = 1/sigma*(nu_L + nu_tilda)*(gy[4]*consts::nu_scale);

        // Viscous fluxes
        f[1] -= n[0]*tau_xx   + n[1]*tau_xy;
        f[2] -= n[0]*tau_xy   + n[1]*tau_yy;
        f[3] -= n[0]*phi[0]   + n[1]*phi[1];
        f[4] -= (n[0]*tau_xx_T + n[1]*tau_xx_T)/consts::nu_scale;
    }

    /*
        Define the time step. Here, time step is constant
    */
    void calc_dt(
        std::vector<double>& dt,
        const std::vector<double>& q,
        const std::vector<double>& gx,
        const std::vector<double>& gy,
        mesh& m
    ) {

        double cfl = consts::cfl;
        const double C = 2;

        // Set max dt
        for (uint i=0; i<dt.size(); ++i) {
            dt[i] = 0.0;
        }
        // Calculate time step scale
        for (uint e=0; e<m.edgesNodes.cols(); ++e) {
            double n[2];

            const auto& i = m.edgesCells(e, 0);
            const auto& j = m.edgesCells(e, 1);
            const auto& le = m.edgesLengths[e];
            
            n[0] = m.edgesNormalsX[e];
            n[1] = m.edgesNormalsY[e];

            // From cell i to cell j
            // Compute max eigenvalue
            const double* qi = &q[vars*i];
            const double* qj = &q[vars*j];
            double ci = sqrt(calc_p(qi)*consts::gamma / qi[0]);
            double cj = sqrt(calc_p(qj)*consts::gamma / qj[0]);

            double eig_ci = (ci + abs(n[0]*qi[1]/qi[0]) + abs(n[1]*qi[2]/qi[0]))*le;
            double eig_cj = (cj + abs(n[0]*qj[1]/qj[0]) + abs(n[1]*qj[2]/qj[0]))*le;

            double eig_vi = std::max(4./(3.*qi[0]), consts::gamma/qi[0])
                            * (consts::mu/consts::pr) * le*le / m.cellsAreas[i];
            double eig_vj = std::max(4./(3.*qj[0]), consts::gamma/qj[0])
                            * (consts::mu/consts::pr) * le*le / m.cellsAreas[j];

            for (uint k=0; k<vars; ++k) {
                dt[vars*i + k] += eig_ci + C*eig_vi;
                dt[vars*j + k] += eig_cj + C*eig_vj;
            }
        }
        // Compute resulting dt from cfl
        for (uint i=0; i<m.cellsAreas.size(); ++i) {
            // Calc source term jacobian
            double source[vars];
            double source_p[vars];
            double q_p[vars];
            for (int k=0; k<vars; ++k)
                q_p[k] = q[vars*i+k];
            q_p[4] += 1e-4;
            calc_source(source, &q[vars*i], &gx[vars*i], &gy[vars*i], m.wall_dist[i]);
            calc_source(source_p, q_p, &gx[vars*i], &gy[vars*i], m.wall_dist[i]);

            double eig_source = abs(source_p[4] - source[4])/1e-4;

            // Compute dt
            for (uint k=0; k<vars-1; ++k) {
                dt[vars*i + k] = cfl * m.cellsAreas[i]/dt[vars*i + k];
            }
            dt[vars*i + vars-1] = cfl * std::min(m.cellsAreas[i]/dt[vars*i + vars-1], 1./eig_source);
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
            bv[1] = bv[0] * 0.2;
            bv[2] = bv[0] * 0.0;
            bv[3] = b_pressure / (consts::gamma - 1) + 0.5/bv[0]*(bv[1]*bv[1] + bv[2]*bv[2]);
            bv[4] = (consts::mu / bv[0]) * 0.1 / consts::nu_scale;

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
                    b[4] = bv[4];
                } else {
                    // Outlet
                    b[0] = q[0];
                    b[1] = q[1];
                    b[2] = q[2];
                    b[3] = q[3];
                    b[4] = q[4];
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
                    b[2] = b[0] * (va - n[1]*(pa - pb)/(rho0*c0));
                    b[3] = pb / (consts::gamma - 1) + 0.5/b[0]*(b[1]*b[1] + b[2]*b[2]);
                    b[4] = bv[4];
                } else {
                    // Outlet
                    double pb = pa;
                    b[0] = rhod + (pb - pd)/(c0*c0);
                    b[1] = b[0] * (ud + n[0]*(pd - pb)/(rho0*c0));
                    b[2] = b[0] * (va + n[1]*(pd - pb)/(rho0*c0));
                    b[3] = pb / (consts::gamma - 1) + 0.5/b[0]*(b[1]*b[1] + b[2]*b[2]);
                    b[4] = q[4];
                }
            }
        }
        void slip_wall(double* b, double* q, double* n) {
            // Flip velocity, doesn't change norm of velocity
            b[0] = q[0];
            b[1] = q[1] - 2.0 * n[0] * (n[0]*q[1] + n[1]*q[2]);
            b[2] = q[2] - 2.0 * n[1] * (n[0]*q[1] + n[1]*q[2]);
            b[3] = q[3];
            b[4] = q[4];
        }
        void wall(double* b, double* q, double* n) {
            b[0] = q[0];
            b[1] = -q[1];
            b[2] = -q[2];
            b[3] = q[3];
            b[4] = 0;
        }
        std::map<std::string, void (*)(double*, double*, double*)> 
        bounds = {
            {"top", farfield},
            {"bot0", slip_wall},
            {"wall", wall},
            {"left", farfield},
            {"right", farfield}
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

        void calc_output_muT(double* muT, double* q) {
            // Compute muT, turbulent viscosity
            const double nu_L = consts::mu / q[0];
            const double Cv1 = 7.1;

            const double nu_tilda = q[4] * consts::nu_scale;
            const double X = nu_tilda / nu_L;
            const double X3 = X*X*X;

            const double fv1 = (X3)/(X3 + Cv1*Cv1*Cv1);
            const double nu_T = fv1 * nu_tilda;
            
            muT[0] = nu_T*q[0];
        }

        std::map<std::string, void (*)(double*, double*)> 
        extra_scalars = {
            {"p", calc_output_p},
            {"muT", calc_output_muT}
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
    std::string name = "square";

    // Read the file
    m.read_file(name, pool);

    fvhyper::solverOptions options;
    options.max_step = 10000;
    options.print_interval = 20;
    options.tolerance = 1e-20;

    // Run solver
    std::vector<double> q;
    fvhyper::run(name, q, pool, m, options);

    // Save file
    fvhyper::writeVtk(name, q, m, pool.rank, pool.size);

    return pool.exit();
}

