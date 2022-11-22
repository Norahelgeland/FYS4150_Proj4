

#include "Ising_model.hpp"
#include <armadillo>
#include <iostream>
#include <math.h>

// Creating random generator
std::mt19937 generator (123);
std::uniform_int_distribution<int> disINT(0.0, 99.);
std::uniform_real_distribution<double> disDOUBLE(0.0, 1.0);


Ising_model::Ising_model(double T_in, int L_in){
    
    // Setting values
    T = T_in; // Temprature
    L = L_in; // Lattice length
    N = L*L; // Lattice area
    S = arma::mat(L, L).fill(1); // The spins
    boltzmann_factor = make_boltzmann_factors(); // Creating and storing the boltzmann factors
    E_tot = tot_energy(); // Total energy of the system
    M_tot = tot_magnetization(); // Total magnitization of the system
    
}
    // Function that calculates the total energy of the system
    double Ising_model::tot_energy(){
        // Create initial variable
        double E_sys = 0;
        // Go trough every spin in the lattice and add the result
        for(int i = 0; i < L; i++){
            for(int j = 0; j < L; j++){
                E_sys += S(i,j)*(S((i+1)%L,j));
                E_sys += S(i,j)*(S(i,(j+1)%L));
            }
         }
        
        return -E_sys;
    }

    // Function that calculates the possible Boltzmann factors for the system
    std::map <double, double> Ising_model::make_boltzmann_factors(){
        // Creating a map where the values wil be stored
        std::map <double, double> boltzmann_factor;
        // Calculating the values
        boltzmann_factor[-8.] = exp((-1./T)*-8.);
        boltzmann_factor[-4.] = exp((-1./T)*-4.);
        boltzmann_factor[0.] = exp((-1./T)*0.);
        boltzmann_factor[4.] = exp((-1./T)*4.);
        boltzmann_factor[8.] = exp((-1./T)*8.);
        return boltzmann_factor;
    }

    // Function that calculates the total magnitization of the system
    double Ising_model::tot_magnetization(){
        // Create initial variable
        int Msum = 0;
        // Go trough every spin in the lattice and add to the sum
        for(int i = 0; i < L; i++){
            for(int j = 0; j < L; j++){
                Msum += S(i,j);
            } 
        } 
    return Msum;
    }

    // Function that does one iteration trough an MCMC cycle
    void Ising_model::MCMC(){

        // Pick random spinn position
        int num_i = disINT(generator);
        int num_j = disINT(generator);

        // Calculate delta E if the spinn is turned
        double delta_E = 2.*S(num_i, num_j)*(S((num_i-1+L)%L, num_j) + S(num_i, (num_j-1+L)%L) + S((num_i+1)%L, num_j) + S(num_i, (num_j+1)%L));
        
        // Pick random number between 0 and 1
        double r = disDOUBLE(generator);

        // Check if it accept the change, and flips the spinn
        if(delta_E < 0 || r <= boltzmann_factor[delta_E]){
            S(num_i, num_j) = S(num_i, num_j)*(-1); // Flipp the spinn
            E_tot += delta_E; // Updates E
            M_tot += 2.*S(num_i, num_j); // Updates M
        }

    }

    // Function that does one MCMC cycle
    void Ising_model::update(){

            // Run one MCMC cycle
            for(int j = 0; j < N; j++){
                MCMC();
                }  
            // Store the epsilon and m value afther the cycle
            epsilon = E_tot/N;
            m = M_tot/N;

    }
    






    



    


