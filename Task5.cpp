#include "Ising_model.cpp"
#include <iostream>
#include <armadillo>
#include "random_model.hpp"
#include <math.h>

int main(){

    // Define temperature and lattice size
    double T = 1.;
    double L = 20.; 

    // Creating an instance of the class
    Ising_model model = Ising_model(T, L);

    // Decide if initial s is random or not
    model.S = random_model(model.S, L);
    model.E_tot = model.tot_energy();
    model.M_tot = model.tot_magnetization();


    double cycles = 10000.; // number of cycles

    // Create vectors to store the expectation values
    arma::vec avg_exp_val_e = arma::vec(cycles).fill(0);
    arma::vec avg_exp_val_m = arma::vec(cycles).fill(0);

    // Create doubles to store the sum
    double Esum = 0;
    double Msum = 0;

    // Run the model 
    for (int n = 0; n < cycles; n++){
    
        model.update();
        Esum += model.epsilon;
        Msum += abs(model.m);

        double avg_e = Esum / (n); // Calculate the average energy per spin
        double avg_m= Msum / (n); // Calculate the average magnitization per spin
        avg_exp_val_e(n) = avg_e; // Add to vector
        avg_exp_val_m(n) = avg_m; // Add to vector

  }
 
    // Write the vectors to files
    std::string filename = "Exp_e_m_random_1_TEST2.txt";
    std::ofstream ofile;
    ofile.open(filename);
    int width = 12;
    int prec  = 4;

    // Loop over steps
    for (int i = 0; i < avg_exp_val_e.size(); i++){
    ofile << std::setw(width) << std::setprecision(prec) << std::scientific << avg_exp_val_e[i]
            << std::setw(width) << std::setprecision(prec) << std::scientific << avg_exp_val_m[i]
            << std::endl; 
    }  
    ofile.close();  


}