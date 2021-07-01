#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>

std::mt19937 g1(1); //random number generator
std::uniform_real_distribution<double> RNG(0.0,1.0); 

int width;
double R;
double fractionFree;
double dx;
double D0;
double koff;
double khop;
double kon;

int ID; //for seeding
int Seed;//


#include "Binding.h"
#include "AddParticle.h"
#include "MoveParticle.h"
#include "Output.h"
#include "Input.h"

int main(int argc, char *argv[]) {

    //input parameters and seed RNG
    Input(argc, argv);
    g1.seed(Seed);

	//output current parameters
    std::cout << "R: " << R << " D0: " << D0 << " dx: " << dx
    	      << " kon: " << kon << " koff: " << koff 
	      << " khop: " << khop << " width: " << width << std::endl;
    
    int step = 0;
    double time = 0;
    int numberFree = 0; //Number of free particles
    int numberBound = 0; // number of bound particles

    std::vector<int> freeParticles; //positions of free particles
    std::vector<int> boundParticles; //positions of bound particles
    freeParticles.reserve(width+2);
    boundParticles.reserve(width+2);
	
    std::vector<int> positions(width, 0); //0 if empty, 1 if particle
	
    //Initial particles
    for(int pos = 0; pos < width; pos++) {
		
		if(RNG(g1) < 1.0 - pos/double(width))
			CreateParticle(pos, numberFree, numberBound, positions, freeParticles, boundParticles);
        
    }
	
    
    double nextTime = 1.0; //time to record at
    double tFactor = 2.0; //factors of time at which to record at
    int scale = 0; 
	
    ///START EVOLUTION
    while(scale < 9) { //run for 2^scale steps
		
		std::vector<double> averageDensity(width+2, 0); //moving average of density over intervals
		
		double thisTime = nextTime;
		double leftFlux = 0;
		double rightFlux = 0;
		
		
		while(thisTime > 0) {
			
			//add to left at infinite rate
			if(positions[0] == 0)
				AddParticle(numberFree, positions, freeParticles, leftFlux);
			
			//kmc 
			double totalRate = numberFree * (kon + khop) + numberBound * koff;
			double deltaT = -1.0/totalRate * log(1.0 - RNG(g1));
			double event = RNG(g1) * totalRate;
			
			//save density
			for(int i = 0; i < width; i++) 
				averageDensity[i] += positions[i] * deltaT;
		    
			averageDensity[width] += deltaT;
			
			//bind
			if(event < numberFree * kon) { 
				
				int particle = floor(RNG(g1) * numberFree);
				Bind(particle, numberFree, numberBound, positions, freeParticles, boundParticles);
				
			}
			
			//move
			else if(event < numberFree * (kon + khop)) {
				
				int particle = floor(RNG(g1) * numberFree);
				MoveParticle(particle, numberFree, positions, freeParticles, rightFlux);
				
			}
			
			//unbind
			else {
				
				int particle = floor(RNG(g1) * numberBound);
				Unbind(particle, numberFree, numberBound, positions, freeParticles, boundParticles);
				
			}
			
			time += deltaT;
			
			step++;
			thisTime -= deltaT;
			
		}
		
		leftFlux /= 1.0*nextTime;
		rightFlux /= 1.0*nextTime;
		averageDensity[width + 1] = (leftFlux + rightFlux)/2.0;
		
		OutputDensityProfile(averageDensity, scale);
		OutputBound(numberFree, numberBound, scale);
		
		scale += 1;
		
		// output to monitor simulation
		std::cout << nextTime << "\t" << averageDensity[width+1] << "\t" << leftFlux << "\t" << rightFlux << "\t" << numberFree << "\t" << numberBound << "\t" << numberFree + numberBound << std::endl;
		
		nextTime *= tFactor;
		
    }
	
}
