#include "Ising_model.cpp"
#include <iostream>
#include <armadillo>

int main(){

    // Test file for the Ising Model
    double T = 1.;
    double L = 2.;

    Ising_model model = Ising_model(T, L);

    //std::cout << model.S;

    //arma::mat S = arma::mat(L+2, L+2).fill(1);
    //S.col(0).fill(0);
    //S.col(L+1).fill(0);
    //S.row(0).fill(0);
    //S.row(L+1).fill(0);

    //std::cout << "###########################";

    //double E_tot = model.tot_energy(S_new);
    //double E_tot = model.tot_energy(model.S);
    //std::cout << E_tot;


    arma::vec possible_E = model.possible_E();

    assert(possible_E(0) == -8.);
    assert(possible_E(1) == 0.);
    assert(possible_E(2) == 0.);
    assert(possible_E(3) == 0.);
    assert(possible_E(4) == -8.);

    

    arma::vec possible_M = model.possible_M();

    assert(possible_M(0) == 4.);
    assert(possible_M(1) == 2.);
    assert(possible_M(2) == 0.);
    assert(possible_M(3) == -2.);
    assert(possible_M(4) == -4.); 



    double Z = model.Z_fun();
    assert(floor(Z) == 5964.);
   /*

    double Boltz = model.boltzmann_dist(model.S);
    //assert(Boltz == 0.146745);

    arma::vec P_poss = model.Possible_p();
    std::cout << P_poss;

    //Testing the expectation value function for the possible energies
    double exp_value = model.Exp_value(possible_E);

   // assert(floor(exp_value)==-20);
   */
 
    for(int i = 1; i <=1; i++){
        model.update();
    }
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

    return 0;

}

