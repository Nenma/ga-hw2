// Tema2.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <iterator>
#include <algorithm>
#include <random>
#include <chrono>
#include <thread>
#include <future>
#include <iomanip>

using namespace std;

ofstream f1("file1.txt");
ofstream f2("file2.txt");
ofstream f3("file3.txt");
ofstream f4("file4.txt");

int COUNTER = 1;

constexpr auto PI = 3.14;
constexpr auto PRECISION = 3;
constexpr auto POP_SIZE = 100;
constexpr auto GEN_MAX = 250;
constexpr auto MUT_PROB = 0.01;
constexpr auto CROSS_PROB = 0.3;
constexpr auto ELITE = 5;

typedef chrono::high_resolution_clock Clock;

#define FIXED_FLOAT(x) std::fixed << std::setprecision(4) << (x)

double lowerBound;
double upperBound;
int dim;

double minVal = DBL_MAX;
double maxVal = DBL_MIN;
double avgVal = 0;
double stdDev = 0;

double minTime = DBL_MAX;
double maxTime = DBL_MIN;
double avgTime = 0;


double randomDouble(double lowerBound, double upperBound) {
	uniform_real_distribution<double> unif(lowerBound, upperBound);
	random_device rd;
	return unif(rd);
}

int randomInt(int lowerBound, int upperBound) {
	uniform_int_distribution<int> unif(lowerBound, upperBound);
	random_device rd;
	return unif(rd);
}

template<typename T>
void printVector(vector<T> vector, ofstream &stream) {
	for (unsigned int i = 0; i < vector.size(); ++i) {
		stream << vector[i] << " ";
	}
	//cout << endl;
}


vector<bool> randomBitstring(int size) {
	vector<bool> bitstring;
	for (int i = 0; i < size; ++i) {
		bitstring.push_back(randomInt(0, 1));
	}
	return bitstring;
}


// 1. Rastrigin's function (f : [-5.12, 5.12]) with min at 0
double Rastrigin(vector<double> x) {
	double res = 0;
	for (unsigned int i = 0; i < x.size(); ++i) {
		res += x[i] * x[i] - 10 * cos(2 * PI * x[i]);
	}
	res += (double)10 * x.size();
	return res;
}

// 2. Rosenbrock's function (f : [-5, 10]) with min at 0
double Rosenbrock(vector<double> x) {
	double res = 0;
	for (unsigned int i = 0; i < x.size() - 1; ++i) {
		res += 100 * pow(x[i] * x[i] - x[i + 1], 2) + (x[i] - 1) * (x[i] - 1);
	}
	return res;
}

// 3. De Jong's (Sphere) function (f : [-5.12, 5.12]) with min at 0
double Sphere(vector<double> x) {
	double res = 0;
	for (unsigned int i = 0; i < x.size(); ++i) {
		res += x[i] * x[i];
	}
	return res;
}

// 4. Michalewicz's function (f : [0, pi]) with min at -4.687(dim = 5), -9.66(dim = 10)
double Michalewicz(vector<double> x) {
	double res = 0;
	for (unsigned int i = 0; i < x.size(); ++i) {
		res += sin(x[i]) * pow(sin(((i + 1) * x[i] * x[i]) / PI), 20);
	}
	return -res;
}



vector<bool> generatePopulation(int cromozomeSize) {
	vector<bool> population;
	for (int i = 0; i < cromozomeSize * POP_SIZE; ++i) {
		population.push_back(randomInt(0, 1));
	}
	return population;
}

// Each bit has a chance of being inverted, except the first 5 representing the elite
void mutation(vector<bool> &population) {
	for (unsigned int i = ELITE * population.size() / POP_SIZE; i < population.size(); ++i) {
		if (randomDouble(0, 1) < MUT_PROB) {
			population[i] = !population[i];
		}
	}
}

// ===StackOverflow answer===
pair<double, int> flip_pair(const pair<int, double>& p)
{
	return pair<double, int>(p.second, p.first);
}

multimap<double, int> flip_map(const map<int, double>& src)
{
	multimap<double, int> dst;
	transform(src.begin(), src.end(), inserter(dst, dst.begin()), flip_pair);
	return dst;
}
// ==========================

void crossover(vector<bool> &population) {
	// Give each individual a chance to participate in crossover
	map<int, double> crossoverChance;
	for (int i = 0; i < POP_SIZE; ++i) {
		crossoverChance.insert({ i, randomDouble(0, 1) });
	}
	// Sort them by chance
	multimap<double, int> crossoverCandidates = flip_map(crossoverChance);
	
	// Have crossovers between pairs of individuals who got their chance
	int length = population.size() / POP_SIZE;
	for (auto it = crossoverCandidates.begin(); it != crossoverCandidates.end(), it->first < CROSS_PROB; ++it) {
		vector<bool> firstParent;
		for (int i = 0; i < length; ++i) {
			firstParent.push_back(population[length * it->second + i]);
		}

		it++;

		vector<bool> secondParent;
		for (int i = 0; i < length; ++i) {
			secondParent.push_back(population[length * it->second + i]);
		}

		int cutPoint = randomInt(0, length - 1);
		vector<bool> firstChild = firstParent;
		vector<bool> secondChild = secondParent;
		for (int i = cutPoint; i < length; ++i) {
			firstChild[i] = secondParent[i];
			secondChild[i] = firstParent[i];
		}


		for (int i = 0; i < length; ++i) {
			population.push_back(firstChild[i]);
		}
		for (int i = 0; i < length; ++i) {
			population.push_back(secondChild[i]);
		}

		// Deallocating memory
		vector<bool>().swap(firstParent);
		vector<bool>().swap(secondParent);
		vector<bool>().swap(firstChild);
		vector<bool>().swap(secondChild);
	}
}

int decimal(vector<bool> cromozome) {
	int x = 0;
	for (unsigned int i = 0; i < cromozome.size(); ++i) {
		x *= 2;
		x += cromozome[i];
	}
	return x;
}

vector<double> decode(vector<bool> cromozome) {
	vector<double> candidate;
	int size = ceil(log2((upperBound - lowerBound) * pow(10, PRECISION)));
	int length = size * dim;
	for (int i = 0; i < length; i += size) {
		vector<bool> temporary;
		for (int j = 0; j < size; ++j) {
			temporary.push_back(cromozome[j + i]);
		}
		double component = lowerBound + decimal(temporary) * (upperBound - lowerBound) / (pow(2, size) - 1);
		candidate.push_back(component);
	}
	return candidate;
}

vector<double> evaluate(vector<bool> population, double f(vector<double>)) {
	double max = DBL_MIN;
	vector<double> value;
	int size = ceil(log2((upperBound - lowerBound) * pow(10, PRECISION)));
	int length = size * dim;
	
	for (unsigned int i = 0; i < population.size(); i += length) {
		vector<bool> cromozome;
		for (int j = 0; j < length; ++j) {
			cromozome.push_back(population[i + j]);
		}
		double eval = f(decode(cromozome));
		value.push_back(eval);
		if (eval > max) {
			max = eval;
		}
	}

	vector<double> fitness;
	for (unsigned int i = 0; i < value.size(); ++i) {
		fitness.push_back(1.1 * max - value[i]);
	}

	return fitness;
}

int chooseOne(vector<double> fitSum) {
	double pos = randomDouble(0, fitSum[fitSum.size() - 1]);
	for (unsigned int i = 0; i < fitSum.size(); ++i) {
		if (pos <= fitSum[i]) {
			return i;
		}
	}
	return 0;
}

vector<bool> selection(vector<bool> population, double f(vector<double>)) {
	int size = ceil(log2((upperBound - lowerBound) * pow(10, PRECISION)));
	int length = size * dim;
	
	vector<double> fitness = evaluate(population, f);


	vector<bool> newPopulation;

	// Making sure 5 of the best cromozomes make it into the next generation
	map<int, double> elites;
	for (unsigned int i = 0; i < fitness.size(); ++i) {
		elites.insert({ i, fitness[i] });
	}

	multimap<double, int> ranking = flip_map(elites);
	auto it = ranking.rbegin(); // reverse iterator
	for (int i = 0; it != ranking.rend(), i < ELITE; ++it, ++i) {
		for (int j = 0; j < length; ++j) {
			newPopulation.push_back(population[it->second * length + j]);
		}
	}

	// Choosing the rest according to the Roulette Wheel
	vector<double> fitSum;
	fitSum.push_back(fitness[0]);

	for (unsigned int i = 1; i < fitness.size(); ++i) {
		fitSum.push_back(fitSum[i - 1] + fitness[i]);
	}
	
	for (int i = ELITE; i < POP_SIZE; ++i) {
		for (int j = 0; j < length; ++j) {
			newPopulation.push_back(population[chooseOne(fitSum) * length + j]);
		}
	}

	return newPopulation;
}

void GA(double f(vector<double>), ofstream &stream) {
	int size = ceil(log2((upperBound - lowerBound) * pow(10, PRECISION)));
	int length = size * dim;
	vector<bool> population = generatePopulation(length);
	
	auto start = Clock::now();

	int generation = 0;

	double min = DBL_MAX;
	double max = DBL_MIN;
	double avg = 0;

	double localMin = DBL_MAX;
	double localMax = DBL_MIN;

	while (generation != GEN_MAX) { //floor(localMin) != floor(localMax) (?)
		mutation(population);
		crossover(population);
		population = selection(population, f);

		generation++;
		
		localMin = DBL_MAX;
		localMax = DBL_MIN;

		for (unsigned int i = 0; i < population.size(); i += length) {
			stream << i / length + 1 << ". ";
			vector<bool> cromozome;
			for (int j = 0; j < length; ++j) {
				cromozome.push_back(population[i + j]);
			}
			printVector(decode(cromozome), stream);
			double value = f(decode(cromozome));
			stream << " -----> " << FIXED_FLOAT(value) << endl;

			if (value > localMax) localMax = value;
			if (value < localMin) localMin = value;

			if (value > max) max = value;
			if (value < min) min = value;
			avg += value;
		}
		stream << "===============" << endl;
		stream << "Min: " << FIXED_FLOAT(min) << "    Max: " << FIXED_FLOAT(max) << endl;
		stream << "=======" << generation << "=======" << endl;

	}

	//stream.close();

	cout << "DONE " << COUNTER <<"   MIN: " << FIXED_FLOAT(min) << endl;
	COUNTER++;

	if (min < minVal) minVal = min;
	if (min > maxVal) maxVal = min;
	avgVal += min;

	auto finish = Clock::now();
	int time = chrono::duration_cast<chrono::seconds>(finish - start).count();

	cout << "Found in: " << time << "s" << endl << endl;

	if (time < minTime) minTime = time;
	if (time > maxTime) maxTime = time;
	avgTime += time;
}



int main()
{
	lowerBound = -5.12; // -5 | -5.12 | 0
	upperBound = 5.12; // 10 | 5.12 | pi
	dim = 30;
	GA(Sphere, f1);  // Rosenbrock | Griewangk | Michalewicz

	for (int i = 0; i < 7; ++i) {
		auto future1 = async(GA, Sphere, ref(f1));
		auto future2 = async(GA, Sphere, ref(f2));
	}

	cout << "MIN: " << minVal << endl
		<< "MAX: " << maxVal << endl
		<< "AVG: " << avgVal << endl << endl
		<< "MIN: " << minTime << "s" << endl
		<< "MAX: " << maxTime << "s" << endl
		<< "AVG:" << avgTime << "s" << endl;

	return 0;
}
