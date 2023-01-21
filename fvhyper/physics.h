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
    bool fixed_dt = false;
    double limiter_k_value = 10.;

    bool flux_term = true;
    bool source_term = true;

    bool fixed_dt = false;
    double dt = -1;   // fixed dt value. for a variable dt, keep fixed_dt < 0

public:

    physics(const int& n) : vars(n) {}

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

    virtual void source(
        double* source,
        const double* q,
    );

    virtual void boundary(
        std::string& name,
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

    void eigenvalue_for_cfl(
        double* q,
        double* n,
        double& length
        double& area
    );

    double limiter(const double& r);

};







}


