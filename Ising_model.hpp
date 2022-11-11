#ifndef __filename_hpp__
#define __filename_hpp__
#include <armadillo>
#include <iostream>

class Ising_model{

    private:

    public:
        double T;
        int L;
        double N;
        arma::mat S;
        double Z;
        double exp_val_E;
        double exp_val_e;
        double exp_val_m;
        double spes_heat;
        double suscept;
        std::map <double, double> boltzmann_factor;


        Ising_model(){};
        // Constructor
        Ising_model(double T, int L);

        
        double tot_energy(arma::mat S);
        
        std::map <double, double>  make_boltzmann_factors();

        //arma::vec possible_E();

        //double Z_fun();

        double tot_magnetization(arma::mat S);

        //arma::vec possible_M();

        //double boltzmann_dist(arma::mat S);

        //arma::vec Possible_p();
        
        // double Exp_value(arma::vec input);    

        void MCMC();

        void update();

};

#endif
