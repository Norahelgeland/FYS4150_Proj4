# FYS4150_Proj4

## Our files

- In the source file Ising_model.cpp we have our Ising Model class. In this file we have also included the header file Ising_model.hpp containg the class declaration.
  
- In the sourcefile Task5.cpp we use the Ising Model class from Ising_model.cpp. It also contains the header file random_model.hpp. Here we run the class function that run the monte Carlo method N times for 1000000 cycles. This we do for the 20 x 20 Lattice. Returns expectation values for the energy.
  
- In the source file Task6.cpp we use the Ising Model class from Ising_model.cpp. Here we do the same as in Task5.cpp but it returns energy per spin for the two temperatures instead of expectation values.
  
- In the sourcefile Task7.cpp we we use the Ising Model class from Ising_model.cpp.  Here have the monte carlo iterations that returns different values for different temperatures. Parallelization is used here.

- The header file random_model.hpp contains a function that generates random numbers and places them in our lattice. 

- Our tests
    - test_Ising_model.cpp

## How to Compile and run the files

- Task5.cpp should be compiled and run for the two temperatures, 1 and 2.4. Then for each temperature there has to be run one with, and one without the random function being commented out. The name of the files that the values are written to should be changed each time.
- Task6.cpp should be compiled and run to produce the energies per spins after each cycle and it gets saved to a file. It should be run for the same temperatures as in Task5.cpp   and the filenames should be changed accordingly.
- Task7.cpp should be compiled and rund for the Lattice sizes 0, 60, 80 and 100.  When changing the lattice size in Task7.cpp, the random number generator (line 9) has to be changed to correspond to the lattice size with max being 39, 59, 79 and 99. In the terminal each lattice size should be run three times. One time for the temperature intervall of 2.1 to 2.4, 2.23 to 2.35 and 2.25 to 2.3. For the first two use 12 threads and for the last set threads to 15. All of them should have 1000 000 Markov Chain Monte Carlo Cycles. The  name of the file that the values gets written to has to be changed in the terminal each time its run.
- To test the Ising Model class run and compile test_Ising_model.cpp
