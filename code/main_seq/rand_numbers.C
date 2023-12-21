#include "aux_func.h"
#include "sobol.h"
#include <iostream>
#include <random>
#include <iterator>
#include <fstream>
#include <vector>
using namespace std;

int main (int argc, char *argv[]) {

	int M = 1000;
	double integral;
	double phi = 1.6180339887498948482;
	double alpha = 1.0/phi;
	double seed = 0;
	vector<int> samples(50);

	for (int it = 0; it < 50; it++) {
		samples[it] = (int)(modf(seed+alpha*it, &integral)*M)+1;
	}

	string filename_error = "outputs/seq/halton_rows.txt";
	ofstream output_file_1(filename_error);
	ostream_iterator<int> output_iterator_1(output_file_1, "\n");
	copy(samples.begin(), samples.end(), output_iterator_1);

	for (int it = 0; it < 50; it++) {
		samples[it] = (int)(sobol::sample(it, 0)*M)+1;
	}

	filename_error = "outputs/seq/sobol_rows.txt";
	ofstream output_file_2(filename_error);
	ostream_iterator<int> output_iterator_2(output_file_2, "\n");
	copy(samples.begin(), samples.end(), output_iterator_2);

	mt19937_64 gen(0);
	uniform_int_distribution<> dis(0, M-1);

	for (int it = 0; it < 50; it++) {
		samples[it] = dis(gen)+1;
	}

	filename_error = "outputs/seq/rand_rows.txt";
	ofstream output_file_3(filename_error);
	ostream_iterator<int> output_iterator_3(output_file_3, "\n");
	copy(samples.begin(), samples.end(), output_iterator_3);

	return 0;
}