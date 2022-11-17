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

    arma::mat results = arma::mat(n_T, 5, arma::fill::zeros);
    
    const double delta_T = (T_max - T_min) / (n_T - 1);  // n_A points correspond to (n_A - 1) intervals

    #pragma omp parallel // Start parallel region
    {

    //double T = 1.;
    // Burde ikke dette være større siden det er et latice og ikke input? 

    // Prepare for file output
    //tatic int print_prec = 10;
    // Each thread will get its own output file name
    //const int my_thread = omp_get_thread_num();
    //std::ofstream ofile;
    //std::string my_output_file_name = output_file_name + "_thread_" + std::to_string(my_thread) + ".txt";
    //ofile.open(my_output_file_name.c_str(), std::ofstream::trunc);  // ofstream::trunc makes sure old content is deleted
    //ofile.open(my_output_file_name);


    // Here we start the parallelized loop over A
    #pragma omp for
    for (int i = 0; i < n_T; ++i)
    {

        double T = T_min + i * delta_T;
        double L = 100.;

        Ising_model model = Ising_model(T, L);

        // Decide if initial s is random or not
        // model.S = random_model(model.S, L);



     //double cycles = 10000.;
    //arma::vec avg_exp_val_e = arma::vec(cycles).fill(0);
    //arma::vec avg_exp_val_m = arma::vec(cycles).fill(0);
    double e_sum = 0;
    double m_sum = 0;
    double EE_sum = 0;
    double E_sum = 0;
    double M_sum = 0;
    double MM_sum = 0;

    for (int n = 0; n < cycles; n++){
    
        model.update();
        e_sum += model.exp_epsilon;
        m_sum += abs(model.exp_m);
        EE_sum += pow(model.E_exp, 2);
        E_sum += model.E_exp;
        M_sum += model.M_exp;
        MM_sum += pow(model.M_exp, 2);

        //double avg_e = Esum / (n);
        //double avg_m= Msum / (n);
        //avg_exp_val_e(n) = avg_e;
        //avg_exp_val_m(n) = avg_m;

        std::cout << "\n\n I'm thread";
        std::cout << i;
        std::cout << ", and I'm working on cycle:";
        std::cout << n;

    }
    double avg_e = e_sum / cycles;
    double avg_m= m_sum / cycles;
    //double C_v = model.spes_heat;
    //double X = model.suscept;
    double C_v = 1./cycles*1./pow(T,2)*(((EE_sum)/cycles)-(pow((E_sum/cycles),2)));
    double X = 1./cycles*1./T*(((abs(MM_sum)/(cycles)))-(pow(abs(M_sum)/cycles, 2)));

    results(i, 0) = T;
    results(i, 1) = avg_e;
    results(i, 2) = avg_m;
    results(i, 3) = C_v;
    results(i, 4) = X;
    

        // Write the vectors to files
        //std::string filename = "e_m_1.txt";
        //std::ofstream ofile;
        //ofile.open(output_file_name);
        int width = 12;
        int prec  = 4;

        // Loop over steps

        // Write results for this A value

    

     } // End parallelized loop over A

  } // End entire parallel region

/*         ofile << std::setw(width) << std::setprecision(prec) << std::scientific << T
                << std::setw(width) << std::setprecision(prec) << std::scientific << avg_e
                << std::setw(width) << std::setprecision(prec) << std::scientific << avg_m
                << std::setw(width) << std::setprecision(prec) << std::scientific << C_v
                << std::setw(width) << std::setprecision(prec) << std::scientific << avg_m
                << std::endl; 
        
        ofile.close();   */

    std::string my_output_file_name = output_file_name + ".csv";
    results.save(my_output_file_name, arma::csv_ascii);

    return 0;
}