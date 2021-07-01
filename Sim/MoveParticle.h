void MoveParticle(int particle, int &numberFree, std::vector<int> &positions, std::vector<int> &freeParticles, double &rightFlux) {
	
    int position = freeParticles[particle];
    
    //move to left
    if (RNG(g1) < 0.5) { 
		
		//check if possible
		if(position > 0 && positions[position - 1] == 0) {
			
			positions[position] = 0;
			positions[position - 1] = 1;
			freeParticles[particle] = position - 1;
			
		}
		
    }
    
    //move to right
    else {
		
		//delete particle at edge
		if (position == width-1) {
			
			rightFlux++;
			positions[position] = 0;
			
			//remove particle
			freeParticles.erase(freeParticles.begin() + particle);
			numberFree--;
			
		}
		
		//check if possible to move right
		else if (positions[position + 1] == 0) {
			
			positions[position] = 0;
			positions[position + 1] = 1;
			freeParticles[particle] = position + 1;
			
		}
		
    }
    
}
