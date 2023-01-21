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

void physics::set_names(std::vector<std::string> names) {
    var_names = names;
}

void physics::set_fixed_dt(double dt_in) {
    fixed_dt = true;
    dt = dt_in;
}


virtual void physics::flux(
    double* flux,
    const double* qi,
    const double* qj,
    const double* gx,
    const double* gy,
    const double* n
) {}


virtual void physics::source(
    double* source,
    const double* q,
) {}


virtual void boundary(
    std::string& name,
    double* b,
    double* q,
    double* n
) {}


virtual void init(
    double* q,
    double& x,
    double& y
) {}


void physics::eigenvalue_for_cfl(
    double* q,
    double* n,
    double& length
    double& area
) {
    if (!fixed_dt) {
        throw std::logic_error( "Tried to call fvhyper::physics::eigenvalue_for_cfl() without implementation." );
    }
}

double physics::limiter(const double& r) {
    return 1.;
}



}


