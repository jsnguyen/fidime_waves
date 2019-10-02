// waves.cpp : Defines the entry point for the application.
//

#include "waves.h"

int main()
{
	const double length = 100;
	const double spacestep = 1;
	double max_time = 100;
	double timestep = 0.5;
	
	double c = 1;

	double eta = c * timestep / spacestep;

	int n_iter = (int) (max_time / timestep);
	const int n_nodes = (const int) (length/spacestep);

	double init_rod[n_nodes];
	double init_rod_v[n_nodes];

	for (int i = 0; i < n_nodes; i++) {
		init_rod[i] = GaussianWavePacket(i,50,c,0,5);
		init_rod_v[i] = GaussianWavePacketDerivative(i, 50, c, 0,5);
	}

	double prev_rod[n_nodes];
	double rod[n_nodes];
	double new_rod[n_nodes];

	for (int i = 0; i < n_nodes; i++) {
		rod[i] = init_rod[i];
		new_rod[i] = 0;
	}

	std::ofstream file;
	file.open("out.txt");

	for (int j = 0; j < n_nodes; j++) {
		file << init_rod[j] << ' ';
	}
	file << std::endl;


	
	for (int j = 2; j < n_nodes-2; j++) {
		double f[3] = { rod[j - 1],rod[j],rod[j + 1] };
		new_rod[j] = WaveEquationExplicitInitial(f , eta, timestep, init_rod_v[j]);
	}

	for (int j = 0; j < n_nodes; j++) {
		file << new_rod[j] << ' ';
		prev_rod[j] = rod[j];
		rod[j] = new_rod[j];
	}
	file << std::endl;

	std::cout << "computing..." << std::endl;

	for (int i = 0; i < n_iter; i++) {
		for (int j = 2; j < n_nodes-2; j++) {
			double f[3] = { rod[j - 1],rod[j],rod[j + 1] };
			new_rod[j] = WaveEquationExplicit(f , eta, prev_rod[j]);
		}

		for (int j = 0; j < n_nodes; j++) {
			file << new_rod[j] << ' ';
			prev_rod[j] = rod[j];
			rod[j] = new_rod[j];
		}
		file << std::endl;
	}

	file.close();
	
	std::cout << "done." << std::endl;
	
	return 0;
}

double WaveEquationExplicitInitial(double f[3], double eta, double timestep, double init_v)
{
	return eta * eta * f[0] + (2 - 2 * eta * eta - 1) * f[1] + eta * eta * f[2] + 2*timestep*init_v;
}

double WaveEquationExplicit(double f[3], double eta, double prev)
{
	return eta * eta * f[0] + (2 - 2 * eta * eta) * f[1] + eta * eta * f[2] - prev;
}

double GaussianWavePacket(double x, double start, double c, double t, double sigma)
{
	return exp(-(x-start - c*t)*(x-start -c*t)/(2*sigma*sigma));
}

double GaussianWavePacketDerivative(double x, double start, double c, double t, double sigma)
{
	//return c*exp(-(x - start - c * t) * (x - start - c * t) / (2*sigma*sigma)) * (x-start -c*t) / (sigma*sigma);
	return 0;
}
