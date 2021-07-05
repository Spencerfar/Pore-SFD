//read input parameters
void Input(int argc, char *argv[]) {
  
	// command line parameter args
	width = atoi(argv[1]);
	R = atof(argv[2]);
	dx = atof(argv[3]);
	D0 = atof(argv[4]);
	koff = atof(argv[5]);
	
	// create parameters
	fractionFree = R/(1.0 + R);
	khop = 2.0 * D0/(dx*dx)/fractionFree;
	kon = koff/R;
	
	
	//get seed
	ID = atoi(argv[6]);
	
	std::ifstream SeedFile;
	SeedFile.open("SeedFile");
	int lineCounter = 0;
	SeedFile >> Seed;
	while(lineCounter < ID) {
		
		SeedFile >> Seed;
		lineCounter++;
		
    }
    std::cout << Seed << std::endl; //output seed
    
}
