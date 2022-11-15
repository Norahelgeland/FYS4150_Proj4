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
        double E;
        double E_exp;
        double EE_exp;
        double epsilon;
        double exp_epsilon;
        double M;
        double m;
        double exp_m;
        double spes_heat;
        double suscept;
        std::map <double, double> boltzmann_factor;
        double E_tot;


        Ising_model(){};
        // Constructor
        Ising_model(double T, int L);

        
        double tot_energy();
        
        std::map <double, double>  make_boltzmann_factors();

        //arma::vec possible_E();

        //double Z_fun();

        double tot_magnetization();

        //arma::vec possible_M();

        //double boltzmann_dist(arma::mat S);

        //arma::vec Possible_p();
        
        // double Exp_value(arma::vec input);    

        void MCMC();

        void update();

};

#endif
