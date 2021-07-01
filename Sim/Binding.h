//bind particle
void Bind(int particle, int &numberFree, int &numberBound, std::vector<int> &positions, std::vector<int> &freeParticles, std::vector<int> &boundParticles) {
	
    //get position of particle
    int position = freeParticles[particle];
	
    //add this particle to list of bound particles
    boundParticles.push_back(position);
    //remove it from free particles
    freeParticles.erase(freeParticles.begin() + particle);
	
    //update particle #s
    numberFree -=1;
    numberBound++;
    

}


//unbind particle
void Unbind(int particle, int &numberFree, int &numberBound, std::vector<int> &positions, std::vector<int> &freeParticles, std::vector<int> &boundParticles) {

    //get position of particle
    int position = boundParticles[particle];

    //add this particle to list of free particles
    freeParticles.push_back(position);
    //remove it from free particles
    boundParticles.erase(boundParticles.begin() + particle);

    //update particle #s
    numberFree++;
    numberBound -= 1;

}
