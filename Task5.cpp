#include "Ising_model.cpp"
#include <iostream>
#include <armadillo>
#include "random_model.hpp"

int main(){

    double T = 1;
    // Burde ikke dette være større siden det er et latice og ikke input?
    double L = 20.; // correct (???)

    Ising_model model = Ising_model(T, L);
    //std::cout << model.possible_E(); 

    // model.S = random_model(model.S, L);
    

   /*  arma::mat S_new = model.S;
    //S_new(1,1) = -1;
    S_new.col(0) = S_new.col(L);
    S_new.col(L+1) = S_new.col(1);
    S_new.row(0) = S_new.row(L);
    S_new.row(L+1) = S_new.row(1); */

    //std::cout << S_new;

    //std::cout << model.tot_energy(S_new); 
 

    //arma::vec p = model.Possible_p();
    double cycles = 1000.;
    arma::vec exp_val_e = arma::vec(cycles).fill(0);
    arma::vec exp_val_m = arma::vec(cycles).fill(0);

    for(int i = 0; i < cycles; i++){
        model.update();
        exp_val_e(i) = model.exp_val_e;
        exp_val_m(i) = model.exp_val_m;
    }

    // Write the vectors to files
    std::string filename = "Exp_e_m_1.txt";
    std::ofstream ofile;
    ofile.open(filename);
    int width = 12;
    int prec  = 4;

    // Loop over steps
    for (int i = 0; i < exp_val_m.size(); i++){
    ofile << std::setw(width) << std::setprecision(prec) << std::scientific << exp_val_e[i]
            << std::setw(width) << std::setprecision(prec) << std::scientific << exp_val_m[i]
            << std::endl; 
    }  
    ofile.close();  
    
    
    std::cout << model.S;
    std::cout << "\n\n tot_energy:";
    std::cout << model.possible_E();
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
    std::cout << "####";   


}