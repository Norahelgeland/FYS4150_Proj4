#include "Ising_model.cpp"
#include <iostream>
#include <armadillo>

int main(){

    // Test file for the Ising Model
    double T = 1.;
    double L = 3.;

    Ising_model model = Ising_model(T, L);

    //std::cout << model.S;

    //arma::mat S = arma::mat(L+2, L+2).fill(1);
    //S.col(0).fill(0);
    //S.col(L+1).fill(0);
    //S.row(0).fill(0);
    //S.row(L+1).fill(0);

    //std::cout << "###########################";

    double E_tot = model.tot_energy(model.S);
    //std::cout << E_tot;

    arma::vec possible_E = model.possible_E();

    assert(possible_E(0) == -36.);
    assert(possible_E(9) == -36.);
    assert(possible_E(1) == -20.);
    assert(possible_E(8) == -20.);
    assert(possible_E(2) == -12.);
    assert(possible_E(3) == -12.);
    assert(possible_E(7) == -12.);
    assert(possible_E(6) == -12.);
    assert(possible_E(4) == -4.);
    assert(possible_E(5) == -4.);

    arma::vec possible_M = model.possible_M();

    assert(possible_M(0) == 9.);
    assert(possible_M(1) == 7.);
    assert(possible_M(2) == 5.);
    assert(possible_M(3) == 3.);
    assert(possible_M(4) == 1.);
    assert(possible_M(5) == -1.);
    assert(possible_M(6) == -3.);
    assert(possible_M(7) == -5.);
    assert(possible_M(8) == -7.);
    assert(possible_M(9) == -9.);


    double Z = model.Z_fun();
   // assert(floor(Z) == 14.);

    double Boltz = model.boltzmann_dist(model.S);
    //assert(Boltz == 0.146745);

    arma::vec P_poss = model.Possible_p();
    std::cout << P_poss;

    //Testing the expectation value function for the possible energies
    double exp_value = model.Exp_value(possible_E);

   // assert(floor(exp_value)==-20);

    for(int i = 1; i <=9; i++){
        model.MCMC();
    }
    std::cout<< model.S;

    return 0;

}

