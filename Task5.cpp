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

    //model.S = random_model(model.S, L);


    double cycles = 1000000.;
    arma::vec avg_exp_val_e = arma::vec(cycles).fill(0);
    arma::vec avg_exp_val_m = arma::vec(cycles).fill(0);
    double Esum = 0;
    double Msum = 0;

    for (int n = 0; n < cycles; n++){
    
        model.update();
        Esum += model.E_col;
        Msum += abs(model.B_col);

        double avg_e = Esum / (n * model.N);
        double avg_m= Esum / (n * model.N);
        avg_exp_val_e(n) = avg_e;
        avg_exp_val_m(n) = avg_m;

    //avgEps_sqrd = sumEE / (nCycles * N * N);
    //avgM = sumM / (nCycles * N);
    //avgM_sqrd = sumMM / (nCycles * N * N);
  }
 
    // Write the vectors to files
    std::string filename = "Exp_e_m_2.4.txt";
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
    /* 
    std::cout << "\n\n exp_val_E:";
    std::cout << model.exp_val_E;
    std::cout << "\n\n exp_val_e:";
    std::cout << model.exp_val_e;
    std::cout << "\n\n exp_val_m:";
    std::cout << model.exp_val_m;
    std::cout << "\n\n tot_energy:";
    std::cout << model.tot_energy(model.S)/model.N;
    std::cout << "\n\n tot_magnetization:";
    std::cout << model.tot_magnetization(model.S)/model.N;
    std::cout << "\n\n spes_heat:";
    std::cout << model.spes_heat;
    std::cout << "\n\n suscept:";
    std::cout << model.suscept;
    std::cout << "####";    */


}