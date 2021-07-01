//make a particle at a given position with random bound/unbound (for initial conditions)
void CreateParticle(int position, int &numberFree, int &numberBound, std::vector<int> &positions, std::vector<int> &freeParticles, std::vector<int> &boundParticles) { //add particle at this position
    
	positions[position] = 1;
	
	if(RNG(g1) < fractionFree) {
		
		freeParticles.push_back(position);
		numberFree++;
		
	}
	
	else {
		
		boundParticles.push_back(position);
		numberBound++;
		
	}
	
}

//add particle at the entrance
void AddParticle(int &numberFree, std::vector<int> &positions, std::vector<int> &freeParticles, double &leftFlux) { //add particle at this position
	
    positions[0] = 1;
    freeParticles.push_back(0);
    numberFree++;
    leftFlux++;
	
}
