/*

       ___     __                    
      / _/  __/ /  __ _____  ___ ____
     / _/ |/ / _ \/ // / _ \/ -_) __/
    /_/ |___/_//_/\_, / .__/\__/_/   
                 /___/_/             

    Finite Volumes for High Performance

    - Description : Physics sources
    - Author : Alexis Angers
    - Contact : alexis.angers@polymtl.ca

*/
#include <fvhyper/physics.h>
#include <stdexcept>



namespace fvhyper {


// Default members
void physics::set_names(std::vector<std::string> names) {
    var_names = names;
}

void physics::set_fixed_dt(double dt_in) {
    fixed_dt = true;
    dt = dt_in;
}

// Required members
void physics::flux(
    double* flux,
    const double* qi,
    const double* qj,
    const double* gx,
    const double* gy,
    const double* n
) {}

void physics::boundary(
    std::string& name,
    double* b,
    double* q,
    double* n
) {}

void physics::init(
    double* q,
    double& x,
    double& y
) {}

// Optional members
void physics::source(
    double* source,
    const double* q
) {
    for (uint i=0; i<vars; ++i) {source[i] = 0;}
}


double physics::eigenvalue_for_cfl(
    double* q,
    double* n,
    const double& dx
) {
    if (!fixed_dt) {
        throw std::logic_error( "Tried to call fvhyper::physics::eigenvalue_for_cfl() without implementation." );
    }
    return 1;
}

double physics::limiter(const double& r) {
    return 1.;
}


void physics::calculate_extra_scalars(double* s, double* q, std::string& name) {}
void physics::calculate_extra_vectors(double* v, double* q, std::string& name) {}


}


