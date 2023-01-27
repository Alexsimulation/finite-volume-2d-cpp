/*

       ___     __                    
      / _/  __/ /  __ _____  ___ ____
     / _/ |/ / _ \/ // / _ \/ -_) __/
    /_/ |___/_//_/\_, / .__/\__/_/   
                 /___/_/             

    Finite Volumes for High Performance

    - Description : Solver sources
    - Author : Alexis Angers
    - Contact : alexis.angers@polymtl.ca

*/
#pragma once

#include <fvhyper/solver.h>




namespace fvhyper {



void solver::print() {
    if (pool.rank ==0) {
        std::cout << "step = " << step_counter << std::endl;
    }
}


void solver::dt() {
    // Compute dt
    for (uint i=0; i<phys.size(); ++i) {
        auto& p = phys[i];
        auto& sol = p.s;

        calculate_dt(sol.dt, sol.q, p.m, p, opt.cfl);
        if (pool.size > 1) update_comms(sol.dt, p.m, p);
    }
    // Minimize dt over all solutions
    double dt_minimal = 1e10;
    for (uint i=0; i<phys.size(); ++i) {
        auto& p = phys[i];
        auto& sol = p.s;

        for (uint j=0; j<sol.dt.size(); ++j) {
            dt_minimal = std::min(dt_minimal, sol.dt[j]);
        }
    }
    for (uint i=0; i<phys.size(); ++i) {
        auto& p = phys[i];
        auto& sol = p.s;

        sol.dt[0] = dt_minimal;
        min_dt(sol.dt, p.m, p);

        if (pool.size > 1) validate_dt(sol.dt, *pool, p);
    }
}

void solver::step() {
    // Handle the dt 
    dt();

    // Runge kutta iterations
    for (uint k=0; k<phys.size(); ++k) {
        auto& p = phys[i];
        auto& sol = p.s;

        for (uint i=0; i<sol.q.size(); ++i) sol.qk[i] = sol.q[i];
    }
    for (const double& a : alpha) {
        for (uint k=0; k<phys.size(); ++k) {
            auto& p = phys[i];
            auto& sol = p.s;

            complete_calc_qt(sol.qt, sol.qk, sol.gx, sol.gy, sol.qmin, sol.qmax, sol.l, p.m, *pool, p);
            // Ignore residual smoothing
            update_cells(sol.qk, sol.q, sol.qt, sol.dt, a, p);
            update_bounds(sol.qk, sol.gx, sol.gy, sol.l, p.m, p);
            if (pool.size > 1) update_comms(sol.qk, p.m, p);
        }
    }
    // Get back qk values into q
    for (uint k=0; k<phys.size(); ++k) {
        auto& p = phys[i];
        auto& sol = p.s;

        for (uint i=0; i<sol.q.size(); ++i) sol.q[i] = sol.qk[i];
    }
}


void solver::run() {
    step_counter = 0;
    while (step_counter < opt.max_step) {
        step();
        step_counter += 1;
    }
}




}