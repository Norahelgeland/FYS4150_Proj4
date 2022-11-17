#include "Ising_model.cpp"
#include <iostream>
#include <armadillo>
#include "random_model.hpp"
#include <math.h>

int main(){

    double T = 1.;
    // Burde ikke dette være større siden det er et latice og ikke input?
    double L = 2.; 

    Ising_model model = Ising_model(T, L);

    // Decide if initial s is random or not
   // model.S = random_model(model.S, L);
    model.E_tot = model.tot_energy();
    model.M_tot = model.tot_magnetization();


    double cycles = 5.;
    arma::vec avg_exp_val_e = arma::vec(cycles).fill(0);
    arma::vec avg_exp_val_m = arma::vec(cycles).fill(0);
    double Esum = 0;
    double Msum = 0;

    for (int n = 0; n < cycles; n++){
    
        model.update();
        Esum += model.epsilon;
        Msum += abs(model.m);

        double avg_e = Esum / (n);
        double avg_m= Msum / (n);
        avg_exp_val_e(n) = avg_e;
        avg_exp_val_m(n) = avg_m;

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