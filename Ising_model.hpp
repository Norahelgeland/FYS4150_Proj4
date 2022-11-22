#ifndef __filename_hpp__
#define __filename_hpp__
#include <armadillo>
#include <iostream>

class Ising_model{

    private:

    public:
        double T; // Temperature
        int L; // Lattice length
        double N; // Lattice size
        arma::mat S; // Stores the spins
        double E_tot; // total energy in the system
        double epsilon; // Energy per spin
        double m; // Magnitization per spin
        double M_tot; // Total magnitization in the system
        std::map <double, double> boltzmann_factor; // Stores boltzmann factors


        Ising_model(){};
        
        // Constructor
        Ising_model(double T, int L);

        // Function that calculates the total energy of the system
        double tot_energy();
        
        // Function that calculates the possible Boltzmann factors for the system
        std::map <double, double>  make_boltzmann_factors();

        // Function that calculates the total magnitization of the system
        double tot_magnetization();

        // Function that does one iteration trough an MCMC cycle
        void MCMC();

        // Function that does one MCMC cycle
        void update();

};

#endif
