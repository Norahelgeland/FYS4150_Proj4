

#include "Ising_model.hpp"
#include <armadillo>
#include <iostream>
#include <math.h>

std::mt19937 generator (123);
std::uniform_int_distribution<int> disINT(0.0, 1.);
std::uniform_real_distribution<double> disDOUBLE(0.0, 1.0);


Ising_model::Ising_model(double T_in, int L_in){
    
    // Setting values
    T = T_in;
    L = L_in;
    N = L*L;
    S = arma::mat(L, L).fill(1);
    boltzmann_factor = make_boltzmann_factors();
    E_tot = tot_energy();
    M_tot = tot_magnetization();
    
}
    double Ising_model::tot_energy(){
        double E_sys = 0;

        for(int i = 0; i < L; i++){
            for(int j = 0; j < L; j++){
                E_sys += S(i,j)*(S((i+1)%L,j));
                E_sys += S(i,j)*(S(i,(j+1)%L));
            }
         }
        
        return -E_sys;
    }

    std::map <double, double> Ising_model::make_boltzmann_factors(){
        std::map <double, double> boltzmann_factor;
        boltzmann_factor[-8.] = exp((-1/T)*-8.);
        boltzmann_factor[-4.] = exp((-1/T)*-4.);
        boltzmann_factor[0.] = exp((-1/T)*0.);
        boltzmann_factor[4.] = exp((-1/T)*4.);
        boltzmann_factor[8.] = exp((-1/T)*8.);
        return boltzmann_factor;
    }


    double Ising_model::tot_magnetization(){
        int Msum = 0;
        for(int i = 0; i < L; i++){
            for(int j = 0; j < L; j++){
                Msum += S(i,j);
            } 
        } 
    return Msum;
    }


    void Ising_model::MCMC(){

        int num_i = disINT(generator);
        int num_j = disINT(generator);

        // Pick random spinn position and turn the spinn
       // int old_sum =  S(num_i, num_j)*(S((num_i-1+L)%L, num_j) + S(num_i, (num_j-1+L)%L) + S((num_i+1)%L, num_j) + S(num_i, (num_j+1)%L));
        
       // S(num_i, num_j) = S(num_i, num_j)*(-1);
        //int new_sum =  S(num_i, num_j)*(S((num_i-1+L)%L, num_j) + S(num_i, (num_j-1+L)%L) + S((num_i+1)%L, num_j) + S(num_i, (num_j+1)%L));
        
        //S(num_i, num_j) = S(num_i, num_j)*(-1);

        //double delta_E = new_sum-old_sum;
        double delta_E = -2.*S(num_i, num_j)*(S((num_i-1+L)%L, num_j) + S(num_i, (num_j-1+L)%L) + S((num_i+1)%L, num_j) + S(num_i, (num_j+1)%L));
        std::cout << "\n\n delta_E:";
        std::cout << delta_E;

        double r = disDOUBLE(generator);

        // Check if it accept the change, and flips the spinn
        if(delta_E < 0 || r <= boltzmann_factor[delta_E]){
            S(num_i, num_j) = S(num_i, num_j)*(-1);
            E_tot += delta_E;
            M_tot += 2.*S(num_i, num_j);
        }

    }

    void Ising_model::update(){

/*             std::cout << "\n\n Energy:";
            std::cout << tot_energy(); 
            std::cout << "\n\n Energy:";
            std::cout << E_tot;
            std::cout << "\n\n magnetization:";
            std::cout << tot_magnetization(); 
            std::cout << M_tot;
            //M = tot_magnetization();
 */
            for(int j = 0; j < N; j++){
                MCMC();
                //E += E_tot; 
                //M += tot_magnetization();
                std::cout << "\n\n Energy_funk:";
                std::cout << tot_energy(); 
                std::cout << "\n\n Energy:";
                std::cout << E_tot;
                std::cout << "\n\n magnetization:";
                std::cout << tot_magnetization(); 
                std::cout << M_tot;
                }  

            //E_exp = E/N;
           // M_exp = M/N;
            //EE_exp = (E*E)/(N*N);
            epsilon = E_tot/N;
            m = M_tot/N;
            
            //exp_epsilon = E/(N*N);
            //exp_m = M/(N*N);
            //spes_heat = 1./N*1./pow(T,2)*(((E*E)/N)-(pow((E/N),2)));
            //suscept = 1./N*1./T*((abs(M*M)/(N))-(pow(abs(M)/N, 2)));

    }
    






    



    


