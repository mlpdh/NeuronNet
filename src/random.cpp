#include "random.h"


 RandomNumbers::RandomNumbers(unsigned long int s) {
	 if (s == 0) {
		 std::random_device rd;
		 seed = rd();
	}	
	 
	rng = std::mt19937(seed);
} 

void RandomNumbers::uniform_double(std::vector<double>& vec, double lower, double upper){
	
	std::uniform_real_distribution<double> d(lower, upper); 
	
	for(auto& x : vec)
	  x = d(rng); 
	
}
	
double RandomNumbers::uniform_double(double lower, double upper){
	
	std::uniform_real_distribution<double> d(lower, upper); 
	return d(rng);
	
}
	
void RandomNumbers::normal(std::vector<double>& vec, double mean, double sd){
	
	  std::normal_distribution<double> d(mean, sd);
	  
	  for(auto& x : vec)
		x = d(rng); 
	  
}
	
double RandomNumbers::normal(double mean, double sd){
	
	std::normal_distribution<double> d(mean, sd);
	return d(rng);
	
}
	
void RandomNumbers::poisson(std::vector<int>& vec, double mean){
	
	std::poisson_distribution<int> d(mean);
	  
	for(auto& x : vec)
	  x = d(rng); 
	   
}

int RandomNumbers::poisson(double mean){
	
	std::poisson_distribution<int> d(mean);
	return d(rng);
	
}
