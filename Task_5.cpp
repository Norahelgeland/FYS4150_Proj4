#include "Ising_model.cpp"
#include <iostream>
#include <armadillo>

int main(){

    double T = 1.;
    double L = 20.;

    Ising_model model = Ising_model(T, L);

    for(int i = 1; i <=100; i++){
        model.MCMC();
    }
    std::cout << model.S;
    std::cout << "###########################";
    std::cout << model.exp_val_E;
    std::cout << "###########################";
    std::cout << model.exp_val_e;
    std::cout << "###########################";
    std::cout << model.exp_val_m;
    std::cout << "###########################";
    std::cout << model.tot_energy(model.S)/model.N;
    std::cout << "###########################";
    std::cout << model.tot_magnetization(model.S)/model.N;
    std::cout << "###########################";
    std::cout << model.spes_heat;
    std::cout << "###########################";
    std::cout << model.suscept;
    std::cout << "###########################";
}