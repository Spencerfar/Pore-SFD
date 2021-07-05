//turn doubles into strings without the trailling zeros
std::string to_string(double x) {
	
	std::string temp = std::to_string(x);
	temp.erase ( temp.find_last_not_of('0') + 1, std::string::npos );
	if(*temp.rbegin() == '.') temp = temp + '0';
	
	return temp;
	
}

//output densities
void OutputDensityProfile(const std::vector<double> &density, int scale) {
	
    std::ofstream output;
	
    std::string name = "Data/DensityProfile_Width" + std::to_string(width) + "Scale" + std::to_string(scale) + "R" + to_string(R) + "koff" + to_string(koff) + "ID" + std::to_string(ID);
    
    output.open(name.c_str());
	
    for(int i = 0; i < width; i++) {
		
		output << i << "\t" << density[i]/density[width] << "\n";
		
    }
	
    //second last is total time for average
    output << width << "\t" << density[width] << "\n";
    //last is flux
    output << width+1 << "\t" << density[width+1] << "\n";
	
    output.close();
	
}

//output number bound
void OutputBound(int numberFree, int numberBound, int scale) {
	
    std::ofstream output;
	
    std::string name = "Data/NumberBound_Width" + std::to_string(width) + "Scale" + std::to_string(scale) + "R" + to_string(R) + "koff" + to_string(koff) + "ID" + std::to_string(ID);
    
    output.open(name.c_str());

	
    output << numberFree << "\t" << numberBound << "\n";

    output.close();


}
