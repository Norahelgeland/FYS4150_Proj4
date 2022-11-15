#include "Ising_model.cpp"
#include <iostream>
#include <armadillo>
#include "random_model.hpp"
#include <math.h>

int main(){

    double T = 2.4;
    // Burde ikke dette være større siden det er et latice og ikke input?
    double L = 20.; 

    Ising_model model = Ising_model(T, L);

    // Decide if initial s is random or not
    // model.S = random_model(model.S, L);


    double cycles = 100000.;

    //arma::vec e_TEST =  arma::vec(1).fill(0);
    arma::vec avg_val_e = arma::vec(cycles).fill(0);
    arma::vec avg_val_m = arma::vec(cycles).fill(0);
    double Esum = 0;
    double Msum = 0;

    for (int n = 0; n < cycles; n++){
        
    
        model.update();
        Esum += model.epsilon;
        Msum += abs(model.m);

        //e_TEST = join_cols(E_func, e_TEST);

        double avg_e = Esum / (n);
        double avg_m= Msum / (n);

        avg_val_e(n) = avg_e;
        avg_val_m(n) = avg_m;

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