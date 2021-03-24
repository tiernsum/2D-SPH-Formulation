#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "cblas.h"
#include "mpi.h"
#include "SPH.h"
using namespace std;

int main(int argc, char* argv[]) {
    
    struct timespec beginReal, endReal, beginCPU, endCPU;
    clock_gettime(CLOCK_REALTIME, &beginReal);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &beginCPU);
        
    int     numOfProcessors, currRank, localN, nInput;
    double  dtInput, TInput, hInput;
    
    try {
        for (int i = 1; i < argc; ++i) {
            if (string(argv[i]) == "--ic-dam-break") {
                nInput = 400;
            }
            if (string(argv[i]) == "--ic-block-drop") {
                nInput = 652;
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
        ofstream dataOut("data.txt", ios::out | ios::trunc);
        energyOut.precision(10);
                
        localN = obj1.getN() / numOfProcessors;
                
        obj1.iterate(xOut, energyOut, dataOut, localN, numOfProcessors, currRank);
            
        MPI_Finalize();
            
        obj1.closePPOutputFile(xOut);
        obj1.closePPOutputFile(energyOut);
    } catch (const logic_error& e) {
        cout << "Error found: " << e.what() << endl;
    }
    
    clock_gettime(CLOCK_REALTIME, &endReal);
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &endCPU);
            
    long secondsReal = endReal.tv_sec - beginReal.tv_sec;
    long nanosecondsReal = endReal.tv_nsec - beginReal.tv_nsec;
    double elapsedReal = secondsReal + nanosecondsReal*1e-9;
            
    long secondsCPU = endCPU.tv_sec - beginCPU.tv_sec;
    long nanosecondsCPU = endCPU.tv_nsec - beginCPU.tv_nsec;
    double elapsedCPU = secondsCPU + nanosecondsCPU*1e-9;
            
    cout << "Elapsed Time (Real): " << elapsedReal << endl;
    cout << "Elapsed Time (CPU): " << elapsedCPU << endl;

    return 0;
    
}