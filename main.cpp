#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include "cblas.h"
#include "mpi.h"
#include "SPH.h"
using namespace std;

int main(int argc, char* argv[]) {
    
    int     numOfProcessors, currRank, localN, nInput;
    double  dtInput, TInput, hInput;
    
    try {
        for (int i = 1; i < argc; ++i) {
            if (string(argv[i]) == "--ic-dam-break") {
                nInput = 400;
            }
            if (string(argv[i]) == "--ic-block-drop") {
                nInput = 650;
            }
            if (string(argv[i]) == "--ic-droplet") {
                nInput = 360;
            }
            if (string(argv[i]) == "--ic-one-particle") {
                nInput = 1;
            }
            if (string(argv[i]) == "--ic-two-particles") {
                nInput = 2;
            }
            if (string(argv[i]) == "--ic-three-particles") {
                nInput = 3;
            }
            if (string(argv[i]) == "--ic-four-particles") {
                nInput = 4;
            }
            if (string(argv[i]) == "--dt") {
                dtInput = atof(argv[i+1]);
            }
            if (string(argv[i]) == "--T") {
                TInput = atof(argv[i+1]);
            }
            if (string(argv[i]) == "--h") {
                hInput = atof(argv[i+1]);
            }
        }
    } catch (const logic_error& e) {
        cout << "Error found: " << e.what() << endl;
    } catch (const overflow_error& e) {
        cout << "Range exceeded: " << e.what() << endl;
    } catch (const bad_alloc& e) {
        cout << "Storage allocation error: " << e.what() << endl;
    }
    
    // Initialise problem as SPH object
    // SPH(const unsigned int& numOfParticles, const double& timeStep, const double& finalT, const double& radOfInfl)
    try {
        SPH obj1(nInput, dtInput, TInput, hInput);
        
            
        MPI_Init(&argc, &argv);
            
        MPI_Comm_size(MPI_COMM_WORLD, &numOfProcessors);
        MPI_Comm_rank(MPI_COMM_WORLD, &currRank);
        
        ofstream xOut("output.txt", ios::out | ios::trunc);
        xOut.precision(10);
        ofstream energyOut("energy.txt", ios::out | ios::trunc);
        energyOut.precision(10);
                
        localN = obj1.getN() / numOfProcessors;
                
        obj1.iterate(xOut, energyOut, localN, numOfProcessors, currRank);
            
        MPI_Finalize();
            
        obj1.closePPOutputFile(xOut);
        obj1.closePPOutputFile(energyOut);
    } catch (const logic_error& e) {
        cout << "Error found: " << e.what() << endl;
    }

    return 0;
    
}