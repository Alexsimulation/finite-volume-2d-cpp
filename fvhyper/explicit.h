/*

       ___     __                    
      / _/  __/ /  __ _____  ___ ____
     / _/ |/ / _ \/ // / _ \/ -_) __/
    /_/ |___/_//_/\_, / .__/\__/_/   
                 /___/_/             

    Finite Volumes for High Performance

    - Description : Explicit solver header
    - Author : Alexis Angers
    - Contact : alexis.angers@polymtl.ca

*/
#pragma once

#include <fvhyper/mesh.h>
#include <fvhyper/parallel.h>
#include <mpi.h>
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <math.h>


namespace fvhyper {


extern const int vars;
extern const std::vector<std::string> var_names;

namespace solver {
    extern const bool do_calc_gradients;
    extern const bool do_calc_limiters;
    extern const bool linear_interpolate;
    extern const bool diffusive_gradients;
    extern const bool global_dt;
    extern const bool smooth_residuals;
    extern const double limiter_k_value;
}


void smooth_residuals(
    std::vector<double>& qt_,
    std::vector<double>& smoother_qt,
    std::vector<double>& smoother,
    mesh& m
);


double limiter_func(const double& r);


void gradient_for_diffusion(
    double* gradx, double* grady, 
    const double* gxi, const double* gyi,
    const double* gxj, const double* gyj,
    const double* qi, const double* qj,
    const double* di, const double* dj
);


void generate_initial_solution(
    std::vector<double>& v,
    const mesh& m
);


void calc_flux(
    double* f,
    const double* qi,
    const double* qj,
    const double* gx,
    const double* gy,
    const double* n
);


void calc_dt(
    std::vector<double>& dt,
    const std::vector<double>& q,
    mesh& m
);


void calc_gradients(
    std::vector<double>& gx,
    std::vector<double>& gy,
    const std::vector<double>& q,
    const mesh& m
);


void calc_limiters(
    std::vector<double>& limiters,
    const std::vector<double>& q,
    const std::vector<double>& gx,
    const std::vector<double>& gy,
    const mesh& m
);


void calc_time_derivatives(
    std::vector<double>& qt,
    const std::vector<double>& q,
    const std::vector<double>& gx,
    const std::vector<double>& gy,
    const std::vector<double>& limiters,
    mesh& m
);


void update_cells(
    std::vector<double>& q,
    std::vector<double>& ql,
    const std::vector<double>& qt,
    const std::vector<double>& dt,
    const double v
);

void update_bounds(
    std::vector<double>& q,
    std::vector<double>& gx,
    std::vector<double>& gy,
    std::vector<double>& limiters,
    mesh& m
);


void update_comms(
    std::vector<double>& q,
    mesh& m
);


void calc_residuals(
    double* R,
    std::vector<double>& qt,
    mesh& m,
    mpi_wrapper& pool
);


void min_dt(std::vector<double>& dt, mesh& m);


void validate_dt(std::vector<double>& dt, mpi_wrapper& pool);


struct solverOptions {
    double max_time = 1e10;
    uint max_step = 1e8;
    uint print_interval = 1;
    double tolerance = 1e-16;
    bool save_time_series = false;
    double time_series_interval = 0.2;
};

void complete_calc_qt(
    std::vector<double>& qt,
    std::vector<double>& q,
    std::vector<double>& gx,
    std::vector<double>& gy,
    std::vector<double>& qmin,
    std::vector<double>& qmax,
    std::vector<double>& limiters,
    mesh& m,
    mpi_wrapper& pool
);


void run(
    const std::string name,
    std::vector<double>& q,
    mpi_wrapper& pool,
    mesh& m,
    solverOptions& opt
);



}


