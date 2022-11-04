


class ising_model{

    private:

    public:

    double T;
    int L;
    int N;

    Ising_model(){};
        // Constructor
    PenningTrap(double T, int L);


    double tot_energy();

    double tot_magnetization();

    double energy_spin();

    double magnetization_spin();

    double boltzmann_dist();

}
