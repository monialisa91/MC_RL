#include <armadillo>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "observables.h" // calculation of observables
#include "some_functions.h" // helping functions

using namespace arma;
using namespace std;

int main() {

    arma_rng::set_seed_random();

    // PARAMETERS
    int L = 10; // linear size of the lattice
    int MC_steps = 10;
    int t = 1;
    int n = 1; // the multiple of space between the independent configurations
    int space = L * L * n; // space between the independent configurations
    int mes = 10000; // number of measurements to average over
    int n_conf;
    double U = 4.0; // potential
    double cp = U / 2; //chemical potential
    double T = 0.150;
    double beta = 1.0 / T;
    double r, delta, E_new, E0, corr_new, corr, corr2, ipr, n_electr;
    double delta_sub, delta_sub2, energies, energies2, cv; // heat capacity
    double c_sub, cc, E_norm;
    mat new_hamiltonian;
    const string filename = "L=16_U=4_T=0.145.txt";
    cout << "cos" << endl;
    //mat initial_hamiltonian = LastHamiltonian(L, t, U, filename);
    mat initial_hamiltonian = random_hamiltonian(U, L, t);
    E0 = energy_conf(initial_hamiltonian, beta, cp);
    ostringstream ss1;
    ostringstream ss2;

    ss1 << L;
    ss2 << U;

    string Size = ss1.str();
    string Pot = ss2.str();

    const string name_cdw = "L=" + Size + "_U=" + Pot + "_" + "cdw.txt";
    const string name_heat = "L=" + Size + "_U=" + Pot + "_" + "heat.txt";


    ofstream heat_capacity;
    ofstream cdw;

    ofstream histogram_energies;
    ofstream histogram_cdw;
    ofstream conf;

    cout << "L= " << Size << endl;
    cout << "U = " << Pot << endl;

    // thermalisation

    for (int k = 16; k > 14; k -=2) {


        // opening file to put the configurations into
        T = 0.001 * k;
        beta = 1.0/T;

        cout << "temperature= " << T << endl;
        ostringstream ss;
        ss << T;

        string Temp = ss.str();

        const string name = "L=" + Size + "_U=" + Pot + "_T=" + Temp + ".txt";
        const string name1 = "L=" + Size + "_U=" + Pot + "_T=" + Temp + "cdw_hist.txt";
        const string name2 = "L=" + Size + "_U=" + Pot + "_T=" + Temp + "energy_hist.txt";




        int acc = 0;
        //printf("conf przed\n");
        new_hamiltonian = MC(initial_hamiltonian, MC_steps, beta, cp);
        //printf("conf po\n");
       // print_conf(new_hamiltonian);

//        for (int i = 0; i < MC_steps; i++) {
//            new_hamiltonian = Swap_sites(initial_hamiltonian);
//            eig_sym(eigval, eigvec, new_hamiltonian);
//            E_new = sum_eigvals(eigval, beta, cp);
//            delta = E_new - E0;
//            //cout<<delta<<endl;
//            r = ((double) rand() / (RAND_MAX));
//            if (delta < 0 || exp(-delta * beta) > r) {
//                initial_hamiltonian = new_hamiltonian;
//                E0 = E_new;
//                acc++;
//            }
//
//        }
        cout<<"koniec termalizacji"<<endl;





// measurement
        corr = 0;
        corr2 = 0;
        //ipr = 0;
        //n_electr = 0;
        //delta_sub = 0;
        //delta_sub2 = 0;
        energies = 0;
        energies2 = 0;


        for (int k=0; k<mes; k++) {
            if(k%100 == 0) {
                cout << k << " konfiguracji" << endl;
                //print_conf(new_hamiltonian);
            }
            new_hamiltonian = MC(new_hamiltonian, space, beta, cp);
            E_new = energy_conf(new_hamiltonian, beta, cp);
            //cout<<"E_new: " <<E_new<<endl;
            corr_new = correlation(new_hamiltonian, L, U);
            //cout<<"corr_new:" <<corr_new<<endl;
            // correlation
            corr += corr_new;
            corr2 += corr_new * corr_new;
            E_norm = E_new/ (L * L);
            energies += E_norm;
            energies2 += E_norm * E_norm;
            conf.open(name.c_str(), ios_base::app);
            save_conf(new_hamiltonian, conf);
            conf.close();
            histogram_cdw.open(name1.c_str() , ios_base::app);
            histogram_cdw << left << setw(1) << corr_new << '\n';
            histogram_cdw.close();
            histogram_energies.open(name2.c_str(), ios_base::app);
            histogram_energies << left << setw(3) << E_norm << '\n';
            histogram_energies.close();


        }
//        for (int k = 0; k < mes; k++) {
//            n_conf = 0; // iterator over independent configurations
//            while (n_conf < L*L) {
//                new_hamiltonian = Swap_sites(initial_hamiltonian);
//                eig_sym(eigval, eigvec, new_hamiltonian);
//                E_new = sum_eigvals(eigval, beta, cp);
//                delta = E_new - E0;
//                r = ((double) rand() / (RAND_MAX));
//                if (delta < 0 || exp(-delta * beta) > r) {
//                    initial_hamiltonian = new_hamiltonian;
//                    E0 = E_new;
//                    n_conf++;
//                }
//
//            }
//
//
//            // IPR
//            //ipr += ipr(eigvec);
//            // NFERM
//            //n_electr += nferm(eigval, beta, cp);
//            //sub_diff
//            //delta_sub += sub_lattice(new_hamiltonian);
//            //delta_sub2 += sub_lattice(new_hamiltonian) * sub_lattice(new_hamiltonian);
//            // energies
//
//
//        }
        // AVERAGES

        energies = energies / mes;
        energies2 = energies2 / mes;
        cv = energies2 - energies * energies;// heat capacity
        cv = cv * beta * beta;
        //cout<<"cv="<<cv<<endl;
        heat_capacity.open(name_heat.c_str(), ios_base::app);
        heat_capacity << left << setw(3) << T << right << setw(12) << cv << '\n';
        heat_capacity.close();

        //n_electr = n_electr/mes;

//	delta_sub = delta_sub/mes;
//	delta_sub2 = delta_sub2/mes;
//	c_sub = delta_sub2 - delta_sub * delta_sub;
//	c_sub = beta * c_sub;

        corr = corr / mes;
        corr2 = corr2 / mes;
        cc = corr2 - corr * corr;
        cc = beta * cc;
        //cout<<cdw<<endl;
        cdw.open(name_cdw.c_str(), ios_base::app);
        cdw << left << setw(3) << T << right << setw(12) << cc << '\n';
        cdw.close();


    }

    return 0;

}
