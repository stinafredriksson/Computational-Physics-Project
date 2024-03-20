//#################################################################
// HEAD

#include <random> // for RNG
#include <cmath> // for abs, exp
#include <iostream> // for writing in the terminal
#include <fstream> // for writing in files
#include <string> // imports the string class

using namespace std;

//#################################################################
// GLOBALS

const double k = 1.380649e-23; // Boltzmann constant

//#################################################################
// SYSTEM CLASS

class System {
public:
    int N; // lattice width
    double J; // spin coupling strength
    double B; // magnetic field interaction strength
    double T; // temperature

    bool bath_flag; // to use heat bath or not

    int **lattice;

    // constructors
    System();
    System(int N_in, double J_in, double B_in, double T_in, bool bflag);

    // pseudoconstructor
    void setup(int N_in, double J_in, double B_in, double T_in, bool bflag);

    double magnetization();
    double hamiltonian();
    void sweep(int Nsweeps);
    void step(int Nsteps);
    
private:
    double calc_dE(int i, int j); // energy change for a spin-flip at i,j
    double calc_p_i(int i, int j); // heat bath probability
    int calc_Snn(int i, int j);
};

//#################################################################
// SYSTEM CLASS PUBLIC METHODS

System::System(){}


System::System(int N_in, double J_in, double B_in, double T_in, bool bflag = false) 
{

    N = N_in;
    J = J_in;
    B = B_in;
    T = T_in;

    bath_flag = bflag;

    double frac;

    if (J_in > 0)
    {
        frac = 0.75; // faster approach to equilibrium
    }
    else
    {
        frac = 0.5;
    }

    // initializes the 2D lattice
    lattice = new int*[N];
    for (int i = 0; i < N; ++i)
    {
        lattice[i] = new int[N];
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0., 1.);

    // populates 2D lattice
    for (int i = 0; i < N; ++i) 
    {
        for (int j = 0; j < N; ++j) 
        {
            double rand = dist(gen);
            int spin = (rand < frac) ? 1 : -1;

            lattice[i][j] = spin;
        }
    }

}


void System::setup(int N_in, double J_in, double B_in, double T_in, bool bflag = false)
{
    N = N_in;
    J = J_in;
    B = B_in;
    T = T_in;

    bath_flag = bflag;

    double frac;

    if (J_in > 0)
    {
        frac = 0.75; // faster approach to equilibrium
    }
    else
    {
        frac = 0.5;
    }

    // initializes 2D lattice
    lattice = new int*[N];
    for (int i = 0; i < N; ++i)
    {
        lattice[i] = new int[N];
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0., 1.);

    // populates 2D lattice
    for (int i = 0; i < N; ++i) 
    {
        for (int j = 0; j < N; ++j) 
        {
            double rand = dist(gen);
            int spin = (rand < frac) ? 1 : -1;

            lattice[i][j] = spin;
        }
    }
}


double System::magnetization()
{
    int total_M = 0;

    for (int i = 0; i < N; ++i) 
    {
        for (int j = 0; j < N; ++j) 
        {
            total_M += lattice[i][j];
        }
    }

    return static_cast<double>(total_M);
}


double System::hamiltonian()
{
    double H = 0.;
    double tot_spin = 0;


    for (int i = 0; i < N; ++i) 
    {
        for (int j = 0; j < N; ++j) 
        {
            int Snn = calc_Snn(i,j);

            H += - J*lattice[i][j]*Snn;
            tot_spin += 1.*lattice[i][j];
        }
    }
    
    return H - B*tot_spin;
}


void System::step(int Nsteps)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, N-1);
    std::uniform_real_distribution<> dist2(0., 1.);

    for (int n = 0; n < Nsteps; n++)
    {
        int i = dist(gen);
        int j = dist(gen);

        if (bath_flag)
        { 
            double p_i = calc_p_i(i,j);

            if (dist2(gen) < p_i)
            {
                lattice[i][j] = 1;
            }
            else
            {
                lattice[i][j] = -1;
            }
        }
        else
        {
            double dE = calc_dE(i,j);

            if (dE < 0 || dist2(gen)<std::exp(-dE/T/k))
            {
                lattice[i][j] = -lattice[i][j];
            }
        }
    }
}


void System::sweep(int Nsweeps)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist2(0., 1.);

    for (int i = 1; i <= Nsweeps; i++)
    {
        for (int i = 0; i < N; ++i) 
        {
            for (int j = 0; j < N; ++j) 
            {
                if (bath_flag)
                {
                    double p_i = calc_p_i(i,j);

                    if (dist2(gen) < p_i)
                    {
                        lattice[i][j] = 1;
                    }
                    else
                    {
                        lattice[i][j] = -1;
                    }
                }

                else
                {
                    double dE = calc_dE(i,j);

                    if (dE < 0 || dist2(gen)<std::exp(-dE/k/T))
                    {
                        lattice[i][j] = -lattice[i][j];
                    }
                }
            }
        }
    }
}

//#################################################################
// SYSTEM CLASS PRIVATE METHODS

double System::calc_dE(int i, int j)
{
    int Snn = calc_Snn(i,j);

    return 2.*J*static_cast<double>(lattice[i][j])*static_cast<double>(Snn) + 2.*B*static_cast<double>(lattice[i][j]);
}


double System::calc_p_i(int i, int j)
{
    int Snn = calc_Snn(i,j);

    double exponential = std::exp(2*J/(k*T)*Snn);

    return exponential/(1+exponential);
}


int System::calc_Snn(int i, int j)
{
    int i_down = i-1;
    int j_down = j-1;

    if (i == 0)
    {
        i_down = N-1; // goes to end of lattice
    }
    if (j == 0)
    {
        j_down = N-1; // goes to end of lattice
    }

    return lattice[(i+1)%N][j]+lattice[i_down][j]+lattice[i][(j+1)%N]+lattice[i][j_down];
}

//#################################################################
// FUNCTIONS

double* linspace(double low, double high, int steps)
// this generates a linear array between low-high inclusive with
// steps+1 points
{
    double* arr = new double[steps+1];

    int i = 0;
    for (double T = low; T<high; T = T+(high-low)/static_cast<double>(steps))
    {
        arr[i] = T;
        i++;
    }

    arr[steps] = high;

    return arr;
}


void loading(int j, int total)
{
    double frac = static_cast<double>(j)/static_cast<double>(total);

    int percentage = frac*100 + 0.5;

    std::cout << "[";

    for (int i = 0; i < frac*20; i++)
    {
        std::cout << "#";
    }
    for (int i = 0; i <= (1-frac)*20; i++)
    {
        std::cout << "-";
    }

    std::cout << "] [" << percentage << "%]\r";
    std::cout.flush();
}


void LLsizes(double J = k, double B = 0., int n_real = 10, int n_samples = 1000, int T_steps = 40, bool bflag = false)
{
    int sizes[] = {4,8,16,32};
    int num_size = 4;

    double low = 1.;
    double high = 4.;

    if (B!=0)
    {
        high = 10.;
    }

    double* Ts = linspace(low,high,T_steps);
    T_steps = T_steps + 1;

    // these are the arrays that will store the final values
    double mag_tot[num_size][T_steps];
    double sus_tot[num_size][T_steps];
    double energy_tot[num_size][T_steps];
    double specific_tot[num_size][T_steps];
    double cumul_tot[num_size][T_steps];

    for (int i = 0; i < T_steps; i++)
    {
        // these arrays will average over the realizations
        double mags[num_size];
        double sus[num_size];
        double energy[num_size];
        double specific[num_size];
        double cumul[num_size];

        for (int s = 0; s<num_size;s++)
        {
            mags[s] = 0.;
            sus[s] = 0.;
            energy[s] = 0.;
            specific[s] = 0.;
            cumul[s] = 0.;
        }

        for (int j = 0; j < n_real; j++)
        {
            // these arrays will average over samples
            double temp_mag[num_size];
            double temp_mag2[num_size];
            double temp_E[num_size];
            double temp_spec[num_size];
            double temp_mag4[num_size];

            loading(i*n_real+j,T_steps*n_real);

            System systems[num_size];

            for (int s = 0; s<num_size;s++)
            {
                systems[s].setup(sizes[s],J,B,Ts[i], bflag);
                systems[s].sweep(20);

                temp_mag[s] = 0.;
                temp_mag2[s] = 0.;
                temp_mag4[s] = 0.;
                temp_E[s] = 0.;
                temp_spec[s] = 0.;
            }

            for (int l = 0; l < n_samples; l++)
            {
                for (int s = 0; s<num_size;s++)
                {
                    double M = systems[s].magnetization();
                    double E = systems[s].hamiltonian();
                    double E2 = std::pow(E,2);
                    double M2 = std::pow(M,2);
                    double M4 = std::pow(M,4);

                    if (J > 0)
                    {
                        temp_mag[s] += std::abs(M); // ferromagnet, look at |M|
                    }
                    else
                    {
                        temp_mag[s] += M; // antiferromagnet, M is about 0
                    }
                    
                    temp_mag2[s] += M2;
                    temp_E[s] += E;
                    temp_spec[s] += E2;
                    temp_mag4[s] += M4;

                    systems[s].sweep(1);
                    systems[s].step(systems[s].N*systems[s].N*5);
                }
            }

            // averages over samples
            for (int s = 0; s<num_size;s++)
            {
                double M_mean = temp_mag[s]/static_cast<double>(n_samples);
                double M2_mean = temp_mag2[s]/static_cast<double>(n_samples);
                double M4_mean = temp_mag4[s]/static_cast<double>(n_samples);
                double E_mean = temp_E[s]/static_cast<double>(n_samples);
                double E2_mean = temp_spec[s]/static_cast<double>(n_samples);

                mags[s] += M_mean/(static_cast<double>(systems[s].N)*static_cast<double>(systems[s].N));
                sus[s] += (M2_mean-std::pow(M_mean,2))/(static_cast<double>(systems[s].N)*static_cast<double>(systems[s].N)*Ts[i]);
                energy[s] += E_mean/k/static_cast<double>(systems[s].N)/static_cast<double>(systems[s].N);
                specific[s] += (E2_mean-std::pow(E_mean,2))/(k*k*static_cast<double>(systems[s].N)*static_cast<double>(systems[s].N)*Ts[i]*Ts[i]);
                cumul[s] += 1.-M4_mean/(3.*std::pow(M2_mean,2));
            }

        }

        // averages over realizations
        for (int s = 0; s<num_size;s++)
        {
            mag_tot[s][i] = mags[s]/static_cast<double>(n_real);
            sus_tot[s][i] = sus[s]/static_cast<double>(n_real);
            energy_tot[s][i] = energy[s]/static_cast<double>(n_real);
            specific_tot[s][i] = specific[s]/static_cast<double>(n_real);
            cumul_tot[s][i] = cumul[s]/static_cast<double>(n_real);
        }
    }

    // writes the results to csv files for plotting in Python
    for (int s = 0; s < num_size; s++)
    {
        std::string filename = "output_N" + std::to_string(sizes[s]) + ".csv";

        std::ofstream file(filename);

        file << "kT/J,M,X,E,Cv,U4\n";

        for (int i=0; i<T_steps;i++)
        {   
            file << k*Ts[i]/std::abs(J);
            file << ","; 
            file << mag_tot[s][i];
            file << ","; 
            file << sus_tot[s][i];
            file << ","; 
            file << energy_tot[s][i];
            file << ","; 
            file << specific_tot[s][i];
            file << ","; 
            file << cumul_tot[s][i];
            file << "\n";
        }
        file.close();
    }    
}

//#################################################################
// MAIN

int main()
{   
    int n_real = 10; // standard 10
    int n_samples = 1000; // standard 1000
    int T_steps = 40; // standard 40

    double B = 0; // standard 0
    double J = k; // standard +k (ferro) -k (anti-ferro)

    bool bflag = true; // standard false

    LLsizes(J,B,n_real,n_samples,T_steps,bflag);
    
    return 0;
}
