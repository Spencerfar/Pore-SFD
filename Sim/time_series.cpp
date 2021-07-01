#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <iomanip>

std::mt19937 g1(1); //Random number generator
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
int scale_series;
int scale_numtime;
std::string Folder = "DataMar9/";

#include "Binding.h"
#include "AddParticle.h"
#include "MoveParticle.h"
#include "Output.h"
#include "Input.h"


int main(int argc, char *argv[]) {

    //input parameters and seed RNG
    Input(argc, argv);
    g1.seed(Seed);
    
    std::cout << "R: " << R << " D0: " << D0 << " dx: " << dx
    	      << " kon: " << kon << " koff: " << koff 
			  << " khop: " << khop << " width: " << width << " ID: "
			  << ID << " seed: " << Seed << std::endl;
    ///
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
		
		if(RNG(g1) < 1.0 - pos/float(width))
			CreateParticle(pos, numberFree, numberBound, positions, freeParticles, boundParticles);
        
    }
	
    double nextTime = 1.0;
    double tFactor = 2.0;
    int scale = 0;
	
    std::ofstream timeSeries_output;
    std::string timeSeries_name = Folder + "TimeSeries_Width" + std::to_string(width) + "R" + to_string(R) + "koff" + to_string(koff) + "ID" + std::to_string(ID);
    timeSeries_output.open(timeSeries_name.c_str());
	
    ///START EVOLUTION
    while(scale < 100) {
		
		std::vector<double> averageDensity(width+2, 0); //moving average of density over intervals
		
		double thisTime = nextTime;
		double leftFlux = 0;
		double rightFlux = 0;
		
		bool time_series = false;
		if (scale == scale_series)
			time_series = true;
		
		int num_series = 0;
		
		while(thisTime > 0) {
			
			int fluxCheck = rightFlux;
			int boundCheck = numberBound;
			
			
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
			
			if(time_series) {
				
				if(num_series < scale_numtime)  {
					
					thisTime = 1000000;
					
					if(boundCheck != numberBound || fluxCheck < rightFlux) {
						
						
						if(boundCheck != numberBound)
							timeSeries_output << time << "\t" << step << "\t" << 0 << "\t" << numberBound << "\t" << numberBound + numberFree << "\n";
						
						if(fluxCheck < rightFlux)
							timeSeries_output << std::setprecision(15) << time << "\t" << step << "\t" << 1 << "\t" << numberBound << "\t" << numberBound + numberFree << "\n";
						
						num_series += 1;
					}
					
				}
				
				else
					return 0;
				
			}
			
		}
		
		double Q = (leftFlux + rightFlux)/2.0;
		leftFlux /= 1.0*nextTime;
		rightFlux /= 1.0*nextTime;
		averageDensity[width + 1] = (leftFlux + rightFlux)/2.0;
		
		timeSeries_output << std::flush;
		std::cout << scale << "\t" << time << "\t" << averageDensity[width + 1] << "\t" << numberBound << std::endl;
		
		scale += 1;
		nextTime *= tFactor;
		
    }
	
    timeSeries_output.close();
	
}
