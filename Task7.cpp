#include "omp.h"  // OpenMP header
//#include "/usr/local/opt/libomp/include/omp.h"
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
    const int cycles = atoi(argv[4]);
    const std::string output_file_name = argv[5];

    arma::mat results = arma::mat(n_T, 5, arma::fill::zeros); // Creates matrix that wil store the results
    
    const double delta_T = (T_max - T_min) / (n_T - 1);  // n_T points correspond to (n_T - 1) intervals

    #pragma omp parallel // Start parallel region
    {

    // Here we start the parallelized loop over A
    #pragma omp for
    for (int i = 0; i < n_T; ++i)
    {

        double T = T_min + i * delta_T;
        double L = 100.;

        Ising_model model = Ising_model(T, L);

        // Decide if initial s is random or not
        // model.S = random_model(model.S, L);


    // Create doubles where the values wil be stored
    double e_sum = 0; // Energy per spin sum
    double m_sum = 0; // Magnitization per spin sum
    double EE_sum = 0; // Energy^2 sum
    double E_sum = 0; // Energy sum
    double M_sum = 0; // Magnitization sum
    double MM_sum = 0; // Magnitization^2 sum

    for (int n = 0; n < cycles; n++){
    
        model.update();
        e_sum += model.epsilon;
        m_sum += abs(model.m);
        EE_sum += pow(model.E_tot, 2);
        E_sum += model.E_tot;
        M_sum += model.M_tot;
        MM_sum += pow(model.M_tot, 2);

        // Printing the status
        std::cout << "\n\n I'm thread";
        std::cout << i;
        std::cout << ", and I'm working on cycle:";
        std::cout << n;
    }

    // Calculate the values
    double avg_e = e_sum / cycles; // Average energy per spin
    double avg_m= m_sum / cycles; // Average magnitization per spin
    double C_v = 1./cycles*1./pow(T,2)*(((EE_sum)/cycles)-(pow((E_sum/cycles),2))); // Heat capacity
    double X = 1./cycles*1./T*(((abs(MM_sum)/(cycles)))-(pow(abs(M_sum)/cycles, 2))); // Susceptebility
    
    // Add the values to the matrix
    results(i, 0) = T;
    results(i, 1) = avg_e;
    results(i, 2) = avg_m;
    results(i, 3) = C_v;
    results(i, 4) = X;
    


     } // End parallelized loop over A

  } // End entire parallel region


    std::string my_output_file_name = output_file_name + ".csv";
    results.save(my_output_file_name, arma::csv_ascii);

    return 0;
}