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
void RandomSpinFlips(int** &Lattice, int dim, int cycleNum, double &energy, double J, int &MagMoment, double temp, double P_E[], int size);

using namespace std;
int main(int argc, char* argv[]) {
	
	int dim;
	if(argc < 2) {
		cout<<"ERROR!"<<endl;
		cout<<"Please input an integer as the dimension of the lattice.\n";
	} else {
		dim = atoi(argv[1]);
	}
	
	int **Lattice;
	bool R = false; // Set R to 'true' if a random initialization is needed.

	LatticeInitialization(Lattice, dim, R);
	//PrintLattice(Lattice, dim);

	double temperature = 1., J = 1.; // unit of temperature: kT/J (?); J -- coupling constant; 

	//cout<<TotalEnergy(Lattice, dim, J)<<" "<<MagneticMoment(Lattice, dim)<<" "<<endl;

	long cycleNum = 100000*dim*dim;
	double energy = TotalEnergy(Lattice, dim, J);
	int MagMoment = MagneticMoment(Lattice, dim);
	double *P_E;
	int size;
	size = ((int) 2*temperature/0.01)*2 + 1;
	P_E = new double[size];

	for(int i = 0; i < size; i++) P_E[i] = 0;

	//cout<<" "<<"T"<<"    "<<"<E>"<<"     "<<"C_v"<<"       "<<"|M|"<<"   "<<"\\Chi"<<endl;
	RandomSpinFlips(Lattice, dim, cycleNum, energy, J, MagMoment, temperature, P_E, size);
	// PrintLattice(Lattice, dim);
	// cout<<TotalEnergy(Lattice, dim, J)<<" "<<MagneticMoment(Lattice, dim)<<" "<<endl;
	
	double average = 0., var = 0.;
	for(int i = 0; i < size; i++) {
		cout<<setw(8)<<(i-(size-1)/2)*0.01<<setfill(' ')<<setw(10)<<setprecision(4)<<P_E[i]<<endl;
		average += ((i-(size-1)/2)*0.01)*P_E[i];
		var += ((i-(size-1)/2)*0.01)*((i-(size-1)/2)*0.01)*P_E[i];
	}
	cout <<"energy_variance (from P_E) = "<<(var - average*average)<<endl;

	delete []P_E;

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

void RandomSpinFlips(int** &Lattice, int dim, int cycleNum, double &energy, double J, int &MagMoment, double temp, double P_E[], int size) {
	//pre-calculated Boltzmann factors for energy change ~ 4J and 8J
	double fourJ, eightJ;
	fourJ = exp(-4*J/temp);
	eightJ = fourJ*fourJ;

	int x, y, neighbor_sum, change;
	double acceptance_prop, random_prop;
	bool accept = false;

	double mean_energy, mean_energy2, mean_m, mean_m_, mean_m2, specific_heat, susceptibility;
	int criterion, count, index;
	double mE;

	mean_energy = 0.;
	mean_energy2 = 0.;
	mean_m = 0.;
	mean_m_ = 0.;
	mean_m2 = 0.;

	criterion = 4*(dim*dim-1);  //start collecting data after this number of cycles
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
			mE = energy/(dim*dim);
			mean_energy += mE;
			mean_energy2 += mE*mE;
			if(mE < 0) {
				index = ((int) (mE/0.01-0.5))+(size-1)/2;
			} else {
				index = ((int) (mE/0.01+0.5))+(size-1)/2;
			}
			P_E[index] += 1;
		}
	}

	for(int i = 0; i<size; i++) {
		P_E[i] /= (double) count;
	}

	mean_energy /= count;
	mean_energy2 /= count;

	cout<<"energy_variance = "<<(mean_energy2 - mean_energy*mean_energy)<<endl;

	// mean_m /= (double) count;
	// mean_m2 /= (double) count;
	// specific_heat = (mean_energy2 - mean_energy*mean_energy)/temp/temp;
	// susceptibility = mean_m2/temp;

	//per spin
	// mean_energy /= (dim*dim);
	// specific_heat /= (dim*dim);
	// mean_m /= (dim*dim);
	// susceptibility /= (dim*dim);

}