#include "Ising_model.cpp"
#include <iostream>
#include <armadillo>
#include "random_model.hpp"
#include <math.h>

int main(){

    // Define temperature and lattice size
    double T = 2.4;
    double L = 20.; 

    // Creating an instance of the class
    Ising_model model = Ising_model(T, L);

    // Decide if initial s is random or not
    // model.S = random_model(model.S, L);


    double cycles = 100000.; // Number of cycles

    // Create vectors to store the expectation values
    arma::vec avg_val_e = arma::vec(cycles).fill(0);
    arma::vec avg_val_m = arma::vec(cycles).fill(0);

    // Create doubles to store the sum
    double Esum = 0;
    double Msum = 0;

    for (int n = 0; n < cycles; n++){
        
    
        model.update();

        avg_val_e(n) = model.epsilon; // add the energy per spin afther one MCMC cycle
        avg_val_m(n) = abs(model.m); // add the magnitization per spin afther one MCMC cycle

  }
 
    // Write the vectors to files
    std::string filename = "Test_E.txt";
    std::ofstream ofile;
    ofile.open(filename);
    int width = 12;
    int prec  = 4;

    // Loop over steps
    for (int i = 0; i < avg_val_e.size(); i++){
    ofile << std::setw(width) << std::setprecision(prec) << std::scientific << avg_val_e[i]
            << std::setw(width) << std::setprecision(prec) << std::scientific << avg_val_m[i]
            << std::endl; 
    }   


}