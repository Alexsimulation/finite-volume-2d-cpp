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



namespace fvhyper {

void physics::set_names(std::vector<std::string> names) {
    var_names = names;
}



}


