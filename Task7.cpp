#include "omp.h"  // OpenMP header
#include "Ising_model.cpp"
#include <iostream>
#include <armadillo>
#include "random_model.hpp"
#include <math.h>

int main(int argc, const char* argv[]){

    // Check number of command line arguments
    assert(argc == 6);

    // Read command line arguments
    const double T_min = atof(argv[1]);
    const double T_max = atof(argv[2]);
    const int n_T = atoi(argv[3]);
    const int n_cycles_per_thread = atoi(argv[4]);
    const string output_file_name = argv[5];
    
    const double delta_T = (T_max - T_min) / (n_T - 1);  // n_A points correspond to (n_A - 1) intervals

    #pragma omp parallel // Start parallel region
    {

    //double T = 1.;
    // Burde ikke dette være større siden det er et latice og ikke input? 

    // Prepare for file output
    static int print_prec = 10;
    // Each thread will get its own output file name
    const int my_thread = omp_get_thread_num();
    ofstream ofile;
    string my_output_file_name = output_file_name + ".thread_" + to_string(my_thread);
    ofile.open(my_output_file_name.c_str(), ofstream::trunc);  // ofstream::trunc makes sure old content is deleted



    // Here we start the parallelized loop over A
    #pragma omp for
    for (int i = 0; i < n_T; ++i)
    {

        double T = T_min + i * delta_T;
        double L = 20.;

        Ising_model model = Ising_model(T, L);

        // Decide if initial s is random or not
        // model.S = random_model(model.S, L);


        double cycles = 100000.;
        arma::vec avg_val_e = arma::vec(cycles).fill(0);
        arma::vec avg_val_m = arma::vec(cycles).fill(0);
        double Esum = 0;
        double Msum = 0;

        for (int n = 0; n < cycles; n++){
        
            model.update();
            Esum += model.epsilon;
            Msum += abs(model.m);

            double avg_e = Esum / (n);
            double avg_m= Msum / (n);

            avg_val_e(n) = avg_e;
            avg_val_m(n) = avg_m;
            }
    
        // Write the vectors to files
        //std::string filename = "e_m_1.txt";
        std::ofstream ofile;
        ofile.open(output_file_name);
        int width = 12;
        int prec  = 4;

        // Loop over steps

        // Write results for this A value

        for (int i = 0; i < avg_val_e.size(); i++){
        ofile << std::setw(width) << std::setprecision(prec) << std::scientific << avg_val_e[i]
                << std::setw(width) << std::setprecision(prec) << std::scientific << avg_val_m[i]
                << std::endl; 
        }  
        ofile.close();  

     } // End parallelized loop over A

  } // End entire parallel region


}