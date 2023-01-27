/*

       ___     __                    
      / _/  __/ /  __ _____  ___ ____
     / _/ |/ / _ \/ // / _ \/ -_) __/
    /_/ |___/_//_/\_, / .__/\__/_/   
                 /___/_/             

    Finite Volumes for High Performance

    - Description : Physics header
    - Author : Alexis Angers
    - Contact : alexis.angers@polymtl.ca

*/
#pragma once

#include <vector>
#include <string>
#include <fvhyper/parallel.h>
#include <fvhyper/mesh.h>


namespace fvhyper {


class solution {
public:
    std::vector<double> q;
    std::vector<double> gx;
    std::vector<double> gy;

    std::vector<double> qmin;
    std::vector<double> qmax;
    std::vector<double> l;

    std::vector<double> dt;

    std::vector<double> qk;

    std::vector<mpi_comm_cells> comms;
};


/*
    Physics class

    Allows user to defines a given boundary value problem on a mesh domain
*/
class physics {

public:
    const int vars;
    std::vector<std::string> var_names;

    std::vector<std::string> extra_scalars;
    std::vector<std::string> extra_vectors;

    // Default 
    bool do_calc_gradients = true;
    bool do_calc_limiters = false;
    bool do_linear_interpolate = false;
    bool do_diffusive_gradients = true;
    bool do_smooth_residuals;
    double limiter_k_value = 10.;

    bool do_source_term = true;

    bool fixed_dt = false;
    double dt = -1;   // fixed dt value. for a variable dt, keep fixed_dt < 0

    // The container for the solution variables
    solution s;

    mesh m;


    physics(const int& n) : vars(n) {

    }

    void set_names(std::vector<std::string> names);

    void set_fixed_dt(double dt_in);


    // Physics definition members

    // Required members

    virtual void flux(
        double* flux,
        const double* qi,
        const double* qj,
        const double* gx,
        const double* gy,
        const double* n
    );

    virtual void boundary(
        uint& id,
        double* b,
        double* q,
        double* n
    );

    virtual void init(
        double* q,
        double& x,
        double& y
    );

    // Optional members

    virtual void source(
        double* source,
        const double* q
    );

    virtual double eigenvalue_for_cfl(
        double* q,
        double* n,
        const double& dx
    );

    virtual double limiter(const double& r);

    void calculate_extra_scalars(double* s, double* q, std::string& name);
    void calculate_extra_vectors(double* v, double* q, std::string& name);

};







}


