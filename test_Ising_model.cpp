#include "Ising_model.cpp"
#include <iostream>
#include <armadillo>

//This is a trap!!!!
//You have been bamboozeled!! This is not the class function
//Im here to test YOU!

int main(){

    // Test file for the Ising Model
    double T = 1.;
    double L = 2.;

    Ising_model model = Ising_model(T, L);


    double E = model.tot_energy();
    assert(E == -8);

    double M = model.tot_magnetization();
    assert(M == 4);

  
    for(int i = 1; i <=3; i++){
        model.update();
    }

    std::cout << "\n\n the spin matrix:";
    std::cout << "\n\n"; 
    std::cout << model.S;
    std::cout << "\n\n tot_energy:";
    std::cout << model.tot_energy();
    std::cout << "\n\n exp_val_E:";
    std::cout << model.E_exp;
    std::cout << "\n\n exp_val_E^2:";
    std::cout << model.EE_exp;
    std::cout << "\n\n exp_val_e:";
    std::cout << model.exp_epsilon;
    std::cout << "\n\n exp_val_m:";
    std::cout << model.exp_m;
    std::cout << "\n\n spes_heat:";
    std::cout << model.spes_heat;
    std::cout << "\n\n suscept:";
    std::cout << model.suscept;
    std::cout << "\n\n";    

    return 0;

}

