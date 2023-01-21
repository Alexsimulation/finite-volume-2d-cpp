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
#include <fvhyper/mesh.h>
#include <vector>


namespace fvhyper {


/*
    Physics class

    Allows user to defines a given boundary value problem on a mesh domain
*/
class physics {

protected:
    const int vars;
    std::vector<std::string> var_names;

    // Default 
    bool do_calc_gradients = true;
    bool do_calc_limiters = false;
    bool linear_interpolate = false;
    bool diffusive_gradients = true;
    double limiter_k_value = 10.;

public:

    physics(const int& n) : vars(n) {}

    void set_names(std::vector<std::string> names);


    void operator()(
        double* flux,
        double* source,
        const double* qi,
        const double* qj,
        const double* gx,
        const double* gy,
        const double* n
    );

    void dt(
        std::vector<double>& dt,
        const std::vector<double>& q,
        mesh& m
    );

    void initial_q(
        double* q,
        double& x,
        double& y
    );

    double limiter(const double& r);

};







// Implementations

void physics::set_names(std::vector<std::string> names) {
    var_names = names;
}





}


