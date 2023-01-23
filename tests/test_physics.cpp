/*

       ___     __                    
      / _/  __/ /  __ _____  ___ ____
     / _/ |/ / _ \/ // / _ \/ -_) __/
    /_/ |___/_//_/\_, / .__/\__/_/   
                 /___/_/             

    Finite Volumes for High Performance

    - Description : Tests for solver edits validation
    - Author : Alexis Angers
    - Contact : alexis.angers@polymtl.ca

*/
#include <fvhyper/test.h>
#include <fvhyper/mesh.h>
#include <fvhyper/physics.h>
#include <fvhyper/explicit.h>
#include <fvhyper/parallel.h>




class heat : public fvhyper::physics {

public:
    heat() : physics(1) {

        // Set a fixed dt
        set_fixed_dt(1e-5);

        // Set the variable names
        var_names = {"T"};

        // Set extra scalars
        extra_scalars = {"2T"};

        // Disable the source term calculation
        do_source_term = false;
    }

    void flux(
        double* flux,
        const double* qi,
        const double* qj,
        const double* gx,
        const double* gy,
        const double* n
    ) {
        flux[0] = -1.0*(gx[0]*n[0] + gy[0]*n[1]);
    }

    void boundary(
        std::string& name,
        double* b,
        double* q,
        double* n
    ) {
        // Check which boundary we're at
        if (name == "wall") {
            b[0] = 1.;
        }
    }

    void init(
        double* q,
        double& x,
        double& y
    ) {
        double r2 = x*x + y*y;
        if (r2 < 0.5*0.5) {
            q[0] = 1.;
        } else {
            q[0] = 0.;
        }
    }


    void calculate_extra_scalars(double* s, double* q, std::string& name) {
        if (name == "2T") {
            s[0] = 2.0*q[0];
        }
    }

};




int main() {

    heat physics;

    return 0;
}



