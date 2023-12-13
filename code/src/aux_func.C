#include "aux_func.h"
#include "cgls.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <algorithm>
#include <iomanip>
using namespace std;

double sqrNorm(double * vec, int size) {
	double norm = 0;
	for (int i = 0; i < size; i++)
		norm += vec[i]*vec[i];
	return norm;
}

double sqrNormDiff(double * vec1, double * vec2, int size) {
	double norm = 0;
	double dif = 0;
	for (int i = 0; i < size; i++) {
		dif = vec1[i] - vec2[i];
		norm += dif*dif;
	}
	return norm;
}

double sqrNormMatrixCol(double ** matrix, int col, int size) {
	double norm = 0;
	for (int i = 0; i < size; i++)
		norm += matrix[i][col]*matrix[i][col];
	return norm;
}

double dotProduct(double * vec1, double * vec2, int size) {
	double dot_p = 0;
	for (int i = 0; i < size; i++)
		dot_p += vec1[i]*vec2[i];
	return dot_p;
}

double dotProductCol(double ** matrix, int col, double * vec2, int size) {
	double dot_p = 0;
	for (int i = 0; i < size; i++)
		dot_p += matrix[i][col]*vec2[i];
	return dot_p;
}

void importDenseSystemTXT(int M, int N, string A_filename, string b_filename, string x_filename, double**& A, double*& b, double*& x) {
	
	ifstream file_A;
	file_A.open(A_filename);
	if (!file_A.is_open()) {
		cout << "Error in opening the matrix file" << endl;
		exit(1);
	}

	double * aux;
	aux = new double[(long)M*(long)N];
	A = new double*[(long)M];
	for (long i = 0; i < M; i++)
		A[i] = &aux[i * N];

	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			file_A >> A[i][j];

	file_A.close();

	ifstream file_b;
	file_b.open(b_filename);
	if (!file_b.is_open()) {
		cout << "Error in opening the b vector file" << endl;
		delete[] A[0];
		delete[] A;
		exit(1);
	}

	b = new double[M];

	for (int i = 0; i < M; i++)
		file_b >> b[i];

	file_b.close();

	ifstream file_x;
	file_x.open(x_filename);
	if (!file_x.is_open()) {
		cout << "Error in opening the x vector file" << endl;
		delete[] A[0];
		delete[] A;
		delete[] b;
		exit(1);
	}

	x = new double[N];

	for (int i = 0; i < N; i++)
		file_x >> x[i];

	file_x.close();
}

void importDenseSystemBIN(int M, int N, string A_filename, string b_filename, string x_filename, double**& A, double*& b, double*& x) {
	
	ifstream file_A;
	file_A.open(A_filename, ios::binary);
	if (!file_A.is_open()) {
		cout << "Error in opening the matrix file" << endl;
		exit(1);
	}

	double * aux;
	aux = new double[(long)M*(long)N];
	A = new double*[(long)M];
	for (long i = 0; i < M; i++)
		A[i] = &aux[i * N];

	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			file_A.read(reinterpret_cast<char *>(&A[i][j]), sizeof(A[i][j]));

	file_A.close();

	ifstream file_b;
	file_b.open(b_filename, ios::binary);
	if (!file_b.is_open()) {
		cout << "Error in opening the b vector file" << endl;
		delete[] A[0];
		delete[] A;
		exit(1);
	}

	b = new double[M];

	for (int i = 0; i < M; i++)
		file_b.read(reinterpret_cast<char *>(&b[i]), sizeof(b[i]));

	file_b.close();

	ifstream file_x;
	file_x.open(x_filename, ios::binary);
	if (!file_x.is_open()) {
		cout << "Error in opening the x vector file" << endl;
		delete[] A[0];
		delete[] A;
		delete[] b;
		exit(1);
	}

	x = new double[N];

	for (int i = 0; i < N; i++)
		file_x.read(reinterpret_cast<char *>(&x[i]), sizeof(x[i]));

	file_x.close();
}

void importDenseMatrixTXT(int M, int N, string A_filename, double**& A) {

	ifstream file_A;
	file_A.open(A_filename);
	if (!file_A.is_open()) {
		cout << "Error in opening the matrix file" << endl;
		exit(1);
	}

	double * aux = new double[M*N];
	A = new double*[M];
	for (int i = 0; i < M; i++)
		A[i] = &aux[i * N];

	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			file_A >> A[i][j];

	file_A.close();
}

void importDenseMatrixBIN(int M, int N, string A_filename, double**& A) {

	ifstream file_A;
	file_A.open(A_filename, ios::binary);
	if (!file_A.is_open()) {
		cout << "Error in opening the matrix file" << endl;
		exit(1);
	}

	double * aux;
	aux = new double[(long)M*(long)N];
	A = new double*[(long)M];
	for (long i = 0; i < M; i++)
		A[i] = &aux[i * N];

	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			file_A.read(reinterpret_cast<char *>(&A[i][j]), sizeof(A[i][j]));

	file_A.close();
}

void importbVectorTXT(int M, string b_filename, double*& b) {
	ifstream file;
	file.open(b_filename);
	if (!file.is_open())
		exit(1);

	b = new double[M];

	for (int i = 0; i < M; i++)
		file >> b[i];
	file.close();
}

void importxVectorTXT(int N, string x_filename, double*& x) {
	ifstream file;
	file.open(x_filename);
	if (!file.is_open())
		exit(1);

	x = new double[N];

	for (int i = 0; i < N; i++)
		file >> x[i];
	file.close();
}

void importbVectorBIN(int M, string b_filename, double*& b) {

	b = new double[M];

	ifstream file;
	file.open(b_filename, ios::binary);
	if (file.is_open())	{
		for (int i = 0; i < M; i++) {
			file.read(reinterpret_cast<char *>(&b[i]), sizeof(b[i]));
		}
		file.close();
	}
	else {
		cout << "ERROR: Invalid input file for b vector." << endl;
		delete[] b;
		exit(1);
	}
}

void importxVectorBIN(int N, string x_filename, double*& x) {

	x = new double[N];

	ifstream file;
	file.open(x_filename, ios::binary);
	if (file.is_open())	{
		for (int i = 0; i < N; i++) {
			file.read(reinterpret_cast<char *>(&x[i]), sizeof(x[i]));
		}
		file.close();
	}
	else {
		cout << "ERROR: Invalid input file for x vector." << endl;
		delete[] x;
		exit(1);
	}
}

void convertMatrixTXT(int M, int N, string A_filename) {

	string input_etx = ".bin";
	string output_etx = ".txt";

	string input_name = A_filename;
	string output_name = input_name.substr(0, input_name.length() - input_etx.length());

	output_name = output_name + output_etx;

	ifstream file_in(input_name, ios::binary);
	ofstream file_out(output_name);
	file_out << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);
	if (file_in.is_open() && file_out.is_open()) {
		double num;
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < N-1; j++) {
				file_in.read(reinterpret_cast<char *>(&num), sizeof(num));
				file_out << num << " ";
			}
			for (int j = N-1; j < N; j++) {
				file_in.read(reinterpret_cast<char *>(&num), sizeof(num));
				file_out << num << endl;
			}
		}
		file_in.close();
		file_out.close();
	}	
	else {
		cout << "ERROR: Invalid input file for matrix." << endl;
		exit(1);
	}

	return;	
}

void convertMatrixBIN(int M, int N, string A_filename) {

	string input_etx = ".txt";
	string output_etx = ".bin";

	string input_name = A_filename;
	string output_name = input_name.substr(0, input_name.length() - input_etx.length());

	output_name = output_name + output_etx;

	ifstream file_in(input_name);
	ofstream file_out(output_name, ios::binary);
	if (file_in.is_open() && file_out.is_open()) {
		double num;
		for (int i = 0; i < M; i++) {
			for (int j = 0; j < N; j++) {
				file_in >> num;
				file_out.write(reinterpret_cast<char*>(&num), sizeof(num));
			}
		}
		file_in.close();
		file_out.close();
	}	
	else {
		cout << "ERROR: Invalid input file for matrix." << endl;
		exit(1);
	}

	return;	
}

void convertxVectorBIN(int N, string x_filename) {

	string input_etx = ".txt";
	string output_etx = ".bin";

	string input_name = x_filename;
	string output_name = input_name.substr(0, input_name.length() - input_etx.length());

	output_name = output_name + output_etx;

	ifstream file_in(input_name);
	ofstream file_out(output_name, ios::binary);
	if (file_in.is_open() && file_out.is_open()) {
		double num;
		for (int i = 0; i < N; i++) {
			file_in >> num;
			file_out.write(reinterpret_cast<char*>(&num), sizeof(num));
		}
		file_in.close();
		file_out.close();
	}	
	else {
		cout << "ERROR: Invalid input file for x vector." << endl;
		exit(1);
	}

	return;
}

void convertbVectorBIN(int M, string b_filename) {

	string input_etx = ".txt";
	string output_etx = ".bin";

	string input_name = b_filename;
	string output_name = input_name.substr(0, input_name.length() - input_etx.length());

	output_name = output_name + output_etx;

	ifstream file_in(input_name);
	ofstream file_out(output_name, ios::binary);
	if (file_in.is_open() && file_out.is_open()) {
		double num;
		for (int i = 0; i < M; i++) {
			file_in >> num;
			file_out.write(reinterpret_cast<char*>(&num), sizeof(num));
		}
		file_in.close();
		file_out.close();
	}	
	else {
		cout << "ERROR: Invalid input file for b vector." << endl;
		exit(1);
	}

	return;
}

void convertxVectorTXT(int N, string x_filename) {

	string input_etx = ".bin";
	string output_etx = ".txt";

	string input_name = x_filename;
	string output_name = input_name.substr(0, input_name.length() - input_etx.length());

	output_name = output_name + output_etx;

	ifstream file_in(input_name, ios::binary);
	ofstream file_out(output_name);
	file_out << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);
	if (file_in.is_open() && file_out.is_open()) {
		double num;
		for (int i = 0; i < N; i++) {
			file_in.read(reinterpret_cast<char *>(&num), sizeof(num));
			file_out << num << endl;
		}
		file_in.close();
		file_out.close();
	}	
	else {
		cout << "ERROR: Invalid input file for x vector." << endl;
		exit(1);
	}

	return;
}

void convertbVectorTXT(int M, string b_filename) {

	string input_etx = ".bin";
	string output_etx = ".txt";

	string input_name = b_filename;
	string output_name = input_name.substr(0, input_name.length() - input_etx.length());

	output_name = output_name + output_etx;

	ifstream file_in(input_name, ios::binary);
	ofstream file_out(output_name);
	file_out << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);
	if (file_in.is_open() && file_out.is_open()) {
		double num;
		for (int i = 0; i < M; i++) {
			file_in.read(reinterpret_cast<char *>(&num), sizeof(num));
			file_out << num << endl;
		}
		file_in.close();
		file_out.close();
	}	
	else {
		cout << "ERROR: Invalid input file for b vector." << endl;
		exit(1);
	}

	return;
}

void genConsistDenseSystems() {

	int M_max = 160000;
	int N_max = 20000;

	vector<int> N_values{50, 100, 200, 500, 750, 1000, 2000, 4000, 10000, 20000};
	vector<int> M_values{2000, 4000, 20000, 40000, 80000, 160000};

	default_random_engine gen_mu;
	gen_mu.seed(1);
	uniform_int_distribution<int> dist_mu(-5, 5);

	default_random_engine gen_sigma;
	gen_sigma.seed(2);
	uniform_int_distribution<int> dist_sigma(1, 20);

	int mu;
	int sigma;

	double* aux = new double[(long)M_max*(long)N_max];
	double** A = new double*[(long)M_max];
	for (long i = 0; i < M_max; i++)
		A[i] = &aux[i * N_max];

	default_random_engine gen;
	gen.seed(3);

	for (int i = 0; i < M_max; i++) {
		mu = dist_mu(gen_mu);
		sigma = dist_sigma(gen_sigma);
		normal_distribution<double> dist(mu,sigma);
		for (int j = 0; j < N_max; j++) {
			A[i][j] = dist(gen);
		}
	}

	mu = dist_mu(gen_mu);
	sigma = dist_sigma(gen_sigma);
	normal_distribution<double> dist(mu,sigma);
	double* x = new double[N_max];
	for (int i = 0; i < N_max; i++)
		x[i] = dist(gen);

	cout << "Full matrix A and solution generated!" << endl;

	string filename_A;
	string filename_x;
	string filename_b;

	int M;
	int N;
	double* b = new double[M_max];
	for (int i = 0; i < M_values.size(); i++)
		for (int j = 0; j < N_values.size(); j++) {
			M = M_values[i];
			N = N_values[j];
			if ((double)M/N > 2.0) {

				for (int k = 0; k < M; k++) {
					b[k] = 0;
					for (int l = 0; l < N; l++)
						b[k] += A[k][l]*x[l];
				}

				filename_A = "../data/dense/A_" + to_string(M) + "_" + to_string(N) + ".bin";
				filename_x = "../data/dense/x_" + to_string(M) + "_" + to_string(N) + ".bin";
				filename_b = "../data/dense/b_" + to_string(M) + "_" + to_string(N) + ".bin";

				ofstream file_A(filename_A, ios::binary);
				if(file_A.is_open()) {
					for(int k = 0; k < M; k++)
						for (int l = 0; l < N; l++)
							file_A.write(reinterpret_cast<char*>(&A[k][l]), sizeof(A[k][l]));
					file_A.close();
				}
				else {
					cout << "ERROR: Invalid output file for matrix A." << endl;
					delete[] A[0];
					delete[] A;
					delete[] x;
					delete[] b;
					exit(1);
				}

				ofstream file_x(filename_x, ios::binary);
				if(file_x.is_open()) {
					for (int l = 0; l < N; l++)
						file_x.write(reinterpret_cast<char*>(&x[l]), sizeof(x[l]));
					file_x.close();
				}
				else {
					cout << "ERROR: Invalid output file for vector x." << endl;
					delete[] A[0];
					delete[] A;
					delete[] x;
					delete[] b;
					exit(1);
				}

				ofstream file_b(filename_b, ios::binary);
				if(file_b.is_open()) {
					for(int k = 0; k < M; k++)
						file_b.write(reinterpret_cast<char*>(&b[k]), sizeof(b[k]));
					file_b.close();
				}
				else {
					cout << "ERROR: Invalid output file for vector b." << endl;
					delete[] A[0];
					delete[] A;
					delete[] x;
					delete[] b;
					exit(1);
				}

				cout << "M = " << M << " and N = " << N << " completed!" << endl;

			}
		}

	delete[] A[0];
	delete[] A;
	delete[] x;
	delete[] b;

	return;
}

void genConsistDenseSystemsNorm() {

	int M_max = 160000;
	int N_max = 20000;

	vector<int> N_values{50, 100, 200, 500, 750, 1000, 2000, 4000, 10000, 20000};
	vector<int> M_values{2000, 4000, 20000, 40000, 80000, 160000};

	int mu = 0;
	int sigma = 1;

	default_random_engine gen;
	gen.seed(3);
	normal_distribution<double> dist(mu,sigma);

	double* aux = new double[(long)M_max*(long)N_max];
	double** A = new double*[(long)M_max];
	for (long i = 0; i < M_max; i++)
		A[i] = &aux[i * N_max];

	for (int i = 0; i < M_max; i++)
		for (int j = 0; j < N_max; j++) {
			A[i][j] = dist(gen);
		}

	double* x = new double[N_max];
	for (int i = 0; i < N_max; i++)
		x[i] = dist(gen);

	cout << "Full matrix A and solution generated!" << endl;

	string filename_A;
	string filename_x;
	string filename_b;

	int M;
	int N;
	double* b = new double[M_max];
	for (int i = 0; i < M_values.size(); i++)
		for (int j = 0; j < N_values.size(); j++) {
			M = M_values[i];
			N = N_values[j];
			if ((double)M/N > 2.0) {

				for (int k = 0; k < M; k++) {
					b[k] = 0;
					for (int l = 0; l < N; l++)
						b[k] += A[k][l]*x[l];
				}

				filename_A = "../data/dense_norm/A_" + to_string(M) + "_" + to_string(N) + ".bin";
				filename_x = "../data/dense_norm/x_" + to_string(M) + "_" + to_string(N) + ".bin";
				filename_b = "../data/dense_norm/b_" + to_string(M) + "_" + to_string(N) + ".bin";

				ofstream file_A(filename_A, ios::binary);
				if(file_A.is_open()) {
					for(int k = 0; k < M; k++)
						for (int l = 0; l < N; l++)
							file_A.write(reinterpret_cast<char*>(&A[k][l]), sizeof(A[k][l]));
					file_A.close();
				}
				else {
					cout << "ERROR: Invalid output file for matrix A." << endl;
					delete[] A[0];
					delete[] A;
					delete[] x;
					delete[] b;
					exit(1);
				}

				ofstream file_x(filename_x, ios::binary);
				if(file_x.is_open()) {
					for (int l = 0; l < N; l++)
						file_x.write(reinterpret_cast<char*>(&x[l]), sizeof(x[l]));
					file_x.close();
				}
				else {
					cout << "ERROR: Invalid output file for vector x." << endl;
					delete[] A[0];
					delete[] A;
					delete[] x;
					delete[] b;
					exit(1);
				}

				ofstream file_b(filename_b, ios::binary);
				if(file_b.is_open()) {
					for(int k = 0; k < M; k++)
						file_b.write(reinterpret_cast<char*>(&b[k]), sizeof(b[k]));
					file_b.close();
				}
				else {
					cout << "ERROR: Invalid output file for vector b." << endl;
					delete[] A[0];
					delete[] A;
					delete[] x;
					delete[] b;
					exit(1);
				}

				cout << "M = " << M << " and N = " << N << " completed!" << endl;

			}
		}

	delete[] A[0];
	delete[] A;
	delete[] x;
	delete[] b;

	return;
}

void genCoherentDenseSystems() {

	int M_max = 160000;
	int N_max = 1000;

	vector<int> N_values{1000};
	vector<int> M_values{4000, 20000, 40000, 80000, 160000};

	default_random_engine gen_col_idx;
	gen_col_idx.seed(1);
	uniform_int_distribution<int> dist_col_idx(0, N_max-1);

	default_random_engine gen_col_val;
	gen_col_val.seed(2);
	uniform_int_distribution<int> dist_col_val(2, 20);

	double* aux = new double[(long)M_max*(long)N_max];
	double** A = new double*[(long)M_max];
	for (long i = 0; i < M_max; i++)
		A[i] = &aux[i * N_max];

	int mu = 0;
	int sigma = 1;
	default_random_engine gen;
	gen.seed(3);
	normal_distribution<double> dist(mu,sigma);

	for (int j = 0; j < N_max; j++) {
		A[0][j] = dist(gen);
	}

	for (int i = 1; i < M_max; i++) {
		for (int j = 0; j < N_max; j++) {
			A[i][j] = A[i-1][j];
		}
		for (int k = 0; k < 5; k++)
			A[i][dist_col_idx(gen_col_idx)] = dist_col_val(gen_col_val);
	}

	double* angles = new double[M_max-1];
	for (int i = 1; i < M_max; i++) {
		angles[i] = dotProduct(A[i], A[0], N_max)/sqrt(sqrNorm(A[i], N_max))/sqrt(sqrNorm(A[0], N_max));
		angles[i] = acos(angles[i]);
	}

	string filename_angles = "../data/dense_coherent/angles_" + to_string(M_max) + "_" + to_string(N_max) + ".bin";

	ofstream file_angles(filename_angles, ios::binary);
	if(file_angles.is_open()) {
		for (int l = 0; l < N_max; l++)
			file_angles.write(reinterpret_cast<char*>(&angles[l]), sizeof(angles[l]));
		file_angles.close();
	}
	else {
		cout << "ERROR: Invalid output file for angles." << endl;
		delete[] A[0];
		delete[] A;
		delete[] angles;
		exit(1);
	}

	delete[] angles;

	double* x = new double[N_max];
	for (int i = 0; i < N_max; i++)
		x[i] = dist(gen);

	cout << "Full matrix A and vector b generated!" << endl;

	string filename_A;
	string filename_x;
	string filename_b;

	int M;
	int N;
	double* b = new double[M_max];
	for (int i = 0; i < M_values.size(); i++)
		for (int j = 0; j < N_values.size(); j++) {
			M = M_values[i];
			N = N_values[j];
			if ((double)M/N > 2.0) {

				for (int k = 0; k < M; k++) {
					b[k] = 0;
					for (int l = 0; l < N; l++)
						b[k] += A[k][l]*x[l];
				}

				filename_A = "../data/dense_coherent/A_" + to_string(M) + "_" + to_string(N) + ".bin";
				filename_x = "../data/dense_coherent/x_" + to_string(M) + "_" + to_string(N) + ".bin";
				filename_b = "../data/dense_coherent/b_" + to_string(M) + "_" + to_string(N) + ".bin";

				ofstream file_A(filename_A, ios::binary);
				if(file_A.is_open()) {
					for(int k = 0; k < M; k++)
						for (int l = 0; l < N; l++)
							file_A.write(reinterpret_cast<char*>(&A[k][l]), sizeof(A[k][l]));
					file_A.close();
				}
				else {
					cout << "ERROR: Invalid output file for matrix A." << endl;
					delete[] A[0];
					delete[] A;
					delete[] x;
					delete[] b;
					exit(1);
				}

				ofstream file_x(filename_x, ios::binary);
				if(file_x.is_open()) {
					for (int l = 0; l < N; l++)
						file_x.write(reinterpret_cast<char*>(&x[l]), sizeof(x[l]));
					file_x.close();
				}
				else {
					cout << "ERROR: Invalid output file for vector x." << endl;
					delete[] A[0];
					delete[] A;
					delete[] x;
					delete[] b;
					exit(1);
				}

				ofstream file_b(filename_b, ios::binary);
				if(file_b.is_open()) {
					for(int k = 0; k < M; k++)
						file_b.write(reinterpret_cast<char*>(&b[k]), sizeof(b[k]));
					file_b.close();
				}
				else {
					cout << "ERROR: Invalid output file for vector b." << endl;
					delete[] A[0];
					delete[] A;
					delete[] x;
					delete[] b;
					exit(1);
				}

				cout << "M = " << M << " and N = " << N << " completed!" << endl;

			}
		}

	delete[] A[0];
	delete[] A;
	delete[] x;
	delete[] b;

	return;
}

void genLSDenseSystems() {

	int M_max = 160000;
	int N_max = 20000;

	vector<int> N_values{50, 100, 200, 500, 750, 1000, 2000, 4000, 10000, 20000};
	vector<int> M_values{2000, 4000, 20000, 40000, 80000, 160000};

	double* aux = new double[(long)M_max*(long)N_max];
	double** A = new double*[(long)M_max];
	for (long i = 0; i < M_max; i++)
		A[i] = &aux[i * N_max];

	default_random_engine gen_mu;
	gen_mu.seed(1);
	uniform_int_distribution<int> dist_mu(-5, 5);

	default_random_engine gen_sigma;
	gen_sigma.seed(2);
	uniform_int_distribution<int> dist_sigma(1, 20);

	int mu;
	int sigma;

	default_random_engine gen;
	gen.seed(3);

	for (int i = 0; i < M_max; i++) {
		mu = dist_mu(gen_mu);
		sigma = dist_sigma(gen_sigma);
		normal_distribution<double> dist(mu,sigma);
		for (int j = 0; j < N_max; j++) {
			A[i][j] = dist(gen);
		}
	}

	mu = dist_mu(gen_mu);
	sigma = dist_sigma(gen_sigma);
	normal_distribution<double> dist(mu,sigma);
	double* x_temp = new double[N_max];
	for (int i = 0; i < N_max; i++)
		x_temp[i] = dist(gen);

	double* b = new double[M_max];

	int M;
	int N;
	double* x = new double[N_max];
	for (int i = 0; i < M_values.size(); i++)
		for (int j = 0; j < N_values.size(); j++) {

			M = M_values[i];
			N = N_values[j];

			if ((double)M/N > 2.0) {

				mu = 0;
				sigma = 1;
				normal_distribution<double> dist(mu,sigma);

				for (int k = 0; k < M; k++) {
					b[k] = 0;
					for (int l = 0; l < N; l++)
						b[k] += A[k][l]*x_temp[l];
					b[k] += dist(gen);
				}

				int it;
				double* x_sol = cglsSolve(M, N, it, A, b, 1E-15);

				for (int i = 0; i < N; i++)
					x[i] = x_sol[i];

				string filename_A = "../data/dense_ls/A_" + to_string(M) + "_" + to_string(N) + ".bin";
				string filename_x = "../data/dense_ls/x_" + to_string(M) + "_" + to_string(N) + ".bin";
				string filename_b = "../data/dense_ls/b_" + to_string(M) + "_" + to_string(N) + ".bin";

				ofstream file_A(filename_A, ios::binary);
				if(file_A.is_open()) {
					for(int k = 0; k < M; k++)
						for (int l = 0; l < N; l++)
							file_A.write(reinterpret_cast<char*>(&A[k][l]), sizeof(A[k][l]));
					file_A.close();
				}
				else {
					cout << "ERROR: Invalid output file for matrix A." << endl;
					delete[] A[0];
					delete[] A;
					delete[] x;
					delete[] b;
					exit(1);
				}

				ofstream file_x(filename_x, ios::binary);
				if(file_x.is_open()) {
					for (int l = 0; l < N; l++)
						file_x.write(reinterpret_cast<char*>(&x[l]), sizeof(x[l]));
					file_x.close();
				}
				else {
					cout << "ERROR: Invalid output file for vector x." << endl;
					delete[] A[0];
					delete[] A;
					delete[] x;
					delete[] b;
					exit(1);
				}

				ofstream file_b(filename_b, ios::binary);
				if(file_b.is_open()) {
					for(int k = 0; k < M; k++)
						file_b.write(reinterpret_cast<char*>(&b[k]), sizeof(b[k]));
					file_b.close();
				}
				else {
					cout << "ERROR: Invalid output file for vector b." << endl;
					delete[] A[0];
					delete[] A;
					delete[] x;
					delete[] b;
					exit(1);
				}

				cout << "M = " << M << " and N = " << N << " compleated!" << endl;

				delete[] x_sol;

			}
		}

	delete[] x_temp;

	delete[] A[0];
	delete[] A;
	delete[] x;
	delete[] b;

	return;
}