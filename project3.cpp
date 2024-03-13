#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

const double k = 1.3806503e-23;

class System {
public:
    // Member variables
    int N;
    double J;
    double B;
    double T;

    int **lattice;

    // Constructor
    System(); // Default constructor
    System(int N_in, double J_in, double B_in, double T_in);

    double magnetization();
    double hamiltonian();
    void sweep(int Nsweeps);
    void step();
    void setup(int N_in, double J_in, double B_in, double T_in);
    
private:
    double calc_dE(int i, int j);
};

System::System(){}

System::System(int N_in, double J_in, double B_in, double T_in) {

    N = N_in;
    J = J_in;
    B = B_in;
    T = T_in;

    lattice = new int*[N];
    for (int i = 0; i < N; ++i) {
        lattice[i] = new int[N];
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0, 1);

    // Initialize 2D lattice
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double rand = dist(gen);
            int spin = (rand < 0.75) ? 1 : -1;

            lattice[i][j] = spin;
        }
    }

}

void System::setup(int N_in, double J_in, double B_in, double T_in)
{
    N = N_in;
    J = J_in;
    B = B_in;
    T = T_in;

    double frac;

    if (J_in > 0)
    {
        frac = 0.75;
    }
    else
    {
        frac = 0.5;
    }

    lattice = new int*[N];
    for (int i = 0; i < N; ++i) {
        lattice[i] = new int[N];
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0, 1);

    // Initialize 2D lattice
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            // int spin = 1;
            double rand = dist(gen);
            int spin = (rand < frac) ? 1 : -1;

            // int spin = dist(gen);
            // if (spin == 0)
            // {
            //     spin = -1;
            // }

            lattice[i][j] = spin;
        }
    }
}

double System::magnetization()
{
    int total_M = 0;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            total_M = total_M + lattice[i][j];
        }
    }

    return total_M/(static_cast<double>(N)*static_cast<double>(N));
}

double System::hamiltonian()
{
    double H = 0.;
    

    for (int i = 0; i < N; ++i) 
    {
        for (int j = 0; j < N; ++j) 
        {
            int i_down = i;
            int j_down = j;
            
            if (i == 0)
            {
                i_down = N-1;
            }
            if (j == 0)
            {
                j_down = N-1;
            }

            H = H - J*lattice[i][j]*(lattice[(i+1)%N][j]+lattice[i_down][j]+lattice[i][(j+1)%N]+lattice[i][j_down]);
        }
    }

    if (B != 0)
    {
        double tot_spin = 0;

        for (int i = 0; i < N; ++i) 
        {
            for (int j = 0; j < N; ++j) 
            {
                tot_spin += 1.*lattice[i][j];
            }
        }
        H = H-B*tot_spin;
    }
    
    return H;
}

double System::calc_dE(int i, int j)
{

    int i_down = i;
    int j_down = j;

    if (i == 0)
    {
        i_down = N-1;
    }
    if (j == 0)
    {
        j_down = N-1;
    }

    return 2.*J*lattice[i][j]*(lattice[(i+1)%N][j]+lattice[i_down][j]+lattice[i][(j+1)%N]+lattice[i][j_down]) + 2.*B*lattice[i][j];
}

void System::step()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dist(0, N-1);
    std::uniform_real_distribution<double> dist2(0, 1);

    int i = dist(gen);
    int j = dist(gen);

    double dE = calc_dE(i,j);

    if (dE < 0 || dist2(gen)<std::exp(-dE/k/T))
    {
        lattice[i][j] = -lattice[i][j];
    }
}

void System::sweep(int Nsweeps)
{

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist2(0, 1);

    for (int i = 1; i <= Nsweeps; i++)
    {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                
                double dE = calc_dE(i,j);

                if (dE < 0 || dist2(gen)<std::exp(-dE/k/T))
                    {
                        lattice[i][j] = -lattice[i][j];
                    }
            }
        }
    }
}

double* linspace(double low, double high, int steps)
{
    double* arr = new double[steps];

    int i = 0;
    for (double T = low; T<=high; T = T+(high-low)/static_cast<double>(steps))
    {
        arr[i] = T;
        i++;
    }

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

void LLsizes(double J = k, double B = 0., int n_real = 10, int n_samples = 10, int T_steps = 100)
{
    int sizes[] = {4,8,16,32};
    int num_size = 4;

    double low = 1.;
    double high = 4.;

    double* Ts = linspace(low,high,T_steps);

    double mag_tot[num_size][T_steps];
    double sus_tot[num_size][T_steps];
    double energy_tot[num_size][T_steps];
    double specific_tot[num_size][T_steps];
    double cumul_tot[num_size][T_steps];

    for (int i = 0; i < T_steps; i++)
    {
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

            loading(i*n_real+j,T_steps*n_real);

            double temp_mag[num_size];
            double temp_mag2[num_size];
            double temp_E[num_size];
            double temp_spec[num_size];
            double temp_mag4[num_size];

            // System systems[] = {System(4,k,0,Ts[i]),System(8,k,0,Ts[i]),System(16,k,0,Ts[i]),System(32,k,0,Ts[i])};

            System systems[num_size];

            for (int s = 0; s<num_size;s++)
            {
                systems[s].setup(sizes[s],J,B,Ts[i]);
                systems[s].sweep(200);
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
                    double M = std::abs(systems[s].magnetization());
                    double E = std::abs(systems[s].hamiltonian());
                    // systems[s].sweep(2);
                    double E2 = std::pow(systems[s].hamiltonian(),2);
                    double M2 = std::abs(std::pow(systems[s].magnetization(),2));
                    // systems[s].sweep(2);
                    double M4 = std::abs(std::pow(systems[s].magnetization(),4));

                    temp_mag[s] += M;
                    temp_mag2[s] += M2;
                    temp_E[s] += E;
                    temp_spec[s] += E2;
                    temp_mag4[s] += M4;

                    systems[s].sweep(2);
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

                mags[s] += M_mean;
                sus[s] += (M2_mean-std::pow(M_mean,2))/(Ts[i]);
                energy[s] += E_mean/k;
                specific[s] += (E2_mean-std::pow(E_mean,2))/(k*k);
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

    for (int s = 0; s < num_size; s++)
    {
        std::string filename = "output_N" + std::to_string(sizes[s]) + ".csv";

        std::ofstream file(filename);

        file << "T,M,X,E,Cv,U4\n";

        for (int i=0; i<T_steps;i++)
        {   
            file << Ts[i];
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

int main()
{
    
    int n_real = 20;
    int n_samples = 1000;
    int T_steps = 40;

    double B = 0;
    double J = -k;

    LLsizes(J,B,n_real,n_samples,T_steps);
    
    return 0;
}