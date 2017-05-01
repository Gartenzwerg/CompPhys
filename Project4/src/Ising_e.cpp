#include<iostream>
#include<stdlib.h>
#include<cmath>
#include<time.h>
#include<fstream>
#include<iomanip>

void LatticeInitialization(int** &Lattice, int dim, bool R);
void PrintLattice(int** &Lattice, int dim);
double TotalEnergy(int** &Lattice, int dim, double J);
int MagneticMoment(int** &Lattice, int dim);
void RandomSpinFlips(int** &Lattice, int dim, int cycleNum, double &energy, double J, int &MagMoment, double temp);

using namespace std;
int main(int argc, char* argv[]) {
	
	int dim;
	if(argc < 2) {
		cout<<"ERROR!"<<endl;
		cout<<"Please input an integer as the dimension of the lattice.\n";
		return 1;
	} else {
		dim = atoi(argv[1]);
	}
	
	int **Lattice;
	bool R = false; // Set R to 'true' if a random initialization is needed.

	LatticeInitialization(Lattice, dim, R);
	//PrintLattice(Lattice, dim);

	double temperature = 2.2, J = 1.; // unit of temperature: kT/J (?); J -- coupling constant;
	double delta_T = 0.01, T_max = 2.35; 

	//cout<<TotalEnergy(Lattice, dim, J)<<" "<<MagneticMoment(Lattice, dim)<<" "<<endl;

	long cycleNum = 1000000*dim*dim;
	double energy = TotalEnergy(Lattice, dim, J);
	int MagMoment = MagneticMoment(Lattice, dim);

	//cout<<" "<<"T"<<"    "<<"<E>"<<"     "<<"C_v"<<"       "<<"|M|"<<"   "<<"\\Chi"<<endl;
	// for(double T = 1.; T < temperature; T += 0.1) {
	// 	RandomSpinFlips(Lattice, dim, cycleNum, energy, J, MagMoment, T);
	// }
	while (temperature <= T_max+delta_T/10.) {
		RandomSpinFlips(Lattice, dim, cycleNum, energy, J, MagMoment, temperature);
		temperature += delta_T;
	}
	// for (double T = T_max+0.1; T < 3.01; T += 0.1) {
	// 	RandomSpinFlips(Lattice, dim, cycleNum, energy, J, MagMoment, T);
	// }
	// PrintLattice(Lattice, dim);
	// cout<<TotalEnergy(Lattice, dim, J)<<" "<<MagneticMoment(Lattice, dim)<<" "<<endl;

	return 0;
}

void LatticeInitialization(int** &Lattice, int dim, bool R) {

	Lattice = new int *[dim];
	for(int i = 0; i < dim; i++) {
		Lattice[i] = new int [dim];
	}

	if(R) {
		srand(time(NULL));
		for(int i = 0; i < dim; i++) {
			for(int j = 0; j < dim; j++){
				if(((double) rand()/RAND_MAX) >= 0.5) {
					Lattice[i][j] = 1;
				} else {
					Lattice[i][j] = -1;
				}
			}
		}

	} else {
		// Setting spins at all lattice points to be spin-up(1); incidentally, spin-down(-1)
		for(int i = 0; i < dim; i++) {
			for(int j = i; j < dim; j++) {
				Lattice[i][j] = 1;
				Lattice[j][i] = 1;
			}
		}
	}

}

void PrintLattice(int** &Lattice, int dim) {
	if (Lattice) {
		for(int i = 0; i < dim; i++) {
			for(int j = 0; j < dim; j++) {
				cout<<Lattice[i][j]<<"  ";
			}
			cout<<endl;
		}
	} else {
		cout<<"ERROR!"<<endl;
		cout<<"The Lattice has not yet been initialized.\n";		
	}
}

double TotalEnergy(int** &Lattice, int dim, double J) {
	double E = 0.;
	int d = dim-1;

	//First ignore the periodic boundary conditions and the "right" and "bottom" boundaries
	for (int i = 0; i < d; i++) {
		for (int j = 0; j < d; j++) {
			E -= Lattice[i][j]*(Lattice[i+1][j] + Lattice[i][j+1]);
		}
	}

	//Take care of the boundaries
	for (int i = 0; i < d; i++) {
		//"bottom" boundary
		E -= Lattice[d][i]*(Lattice[0][i] + Lattice[d][i+1]);
		//"right" boundary
		E -= Lattice[i][d]*(Lattice[i][0] + Lattice[i+1][d]);
	}

	//Special care for the bottom right corner
	E -= Lattice[d][d]*(Lattice[d][0] + Lattice[0][d]);

	E *= J;
	return E;
}

int MagneticMoment(int** &Lattice, int dim) {
	int M = 0;

	for(int i = 0; i < dim; i++) {
		for(int j = 0; j < dim; j++) {
			M += Lattice[i][j];
		}
	}
	return M;
}

void RandomSpinFlips(int** &Lattice, int dim, int cycleNum, double &energy, double J, int &MagMoment, double temp) {
	//pre-calculated Boltzmann factors for energy change ~ 4J and 8J
	double fourJ, eightJ;
	fourJ = exp(-4*J/temp);
	eightJ = fourJ*fourJ;

	int x, y, neighbor_sum, change;
	double acceptance_prop, random_prop;
	bool accept = false;

	double mean_energy, mean_energy2, mean_m, mean_m_, mean_m2, specific_heat, susceptibility;
	int criterion, count;

	mean_energy = 0.;
	mean_energy2 = 0.;
	mean_m = 0.;
	mean_m_ = 0.;
	mean_m2 = 0.;

	criterion = 4*dim*dim-1;  //start collecting data after this number of cycles
	count = cycleNum - criterion - 1; //total count of "data" collected

	for(long i = 0; i < cycleNum; i++) {
		if(((i >> 31) << 31) == i) srand(time(NULL)+(i>>abs(change))%31);  //regenerate the seed if i is divisible by 2^31
		accept = false;
		x = rand()%dim;
		y = rand()%dim;
		neighbor_sum = Lattice[(x-1+dim)%dim][y] + Lattice[(x+1)%dim][y] + Lattice[x][(y-1+dim)%dim] + Lattice[x][(y+1)%dim];
		change = 2*Lattice[x][y]*neighbor_sum;
		if(change <= 0) {
			accept = true;
		} else if(change > 0) {
			if(change == 4) {
				acceptance_prop = fourJ;
			} else {
				acceptance_prop = eightJ;
			}
			random_prop = ((double) rand()/RAND_MAX);  //RAND_MAX = 2147483647 = 2^31-1
			if (acceptance_prop > random_prop) accept = true;
		}
		if(accept) {
			Lattice[x][y] = -Lattice[x][y];
			energy += J*change;
			MagMoment += 2*Lattice[x][y];
		}
		if(i > criterion) {
			mean_energy += energy;
			mean_energy2 += energy*energy;
			mean_m += abs(MagMoment);
			mean_m_ += MagMoment;
			mean_m2 += MagMoment*MagMoment;
		}
	}

	mean_energy /= count;
	mean_energy2 /= count;
	mean_m /= (double) count;
	mean_m2 /= (double) count;
	mean_m_ /= (double) count;
	specific_heat = (mean_energy2 - mean_energy*mean_energy)/temp/temp;
	susceptibility = (mean_m2 - mean_m*mean_m)/temp;

	//per spin
	mean_energy /= (dim*dim);
	specific_heat /= (dim*dim);
	mean_m /= (dim*dim);
	susceptibility /= (dim*dim);

	cout<<setfill(' ')<<setw(8)<<temp<<setfill(' ')<<setw(12)<<mean_energy<<setfill(' ')<<setw(12)<<specific_heat<<setfill(' ')<<setw(12)<<mean_m<<setfill(' ')<<setw(12)<<susceptibility<<endl;
}