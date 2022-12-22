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


extern int vars;

namespace solver {
    extern bool do_calc_gradients;
    extern bool do_calc_limiters;
}


void limiter_func(double* l, const double* r);


void generate_initial_solution(
    std::vector<double>& v,
    const mesh& m
);


void calc_flux(
    double* f,
    const double* vi,
    const double* vj,
    const double* gxi,
    const double* gyi,
    const double* gxj,
    const double* gyj,
    const double* n,
    const double* l,
    const double* di,
    const double* dj,
    const double area,
    const double len
);


void calc_dt(
    double& dt,
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
    const double dt
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


void validate_dt(double& dt, mpi_wrapper& pool);


struct solverOptions {
    double max_time = 1e10;
    uint max_step = 1e8;
    uint print_interval = 1;
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
    std::vector<double>& q,
    mpi_wrapper& pool,
    mesh& m,
    solverOptions& opt
);



}


