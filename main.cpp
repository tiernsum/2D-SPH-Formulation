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
    
    int numOfProcessors, currRank, localN, nInput;
    double dtInput, TInput, hInput;
    
    for (unsigned int i = 1; i < argc; ++i) {
        
        if (string(argv[i]) == "--ic-dam-break") {
            
            nInput = 900;
            
        }
        
        if (string(argv[i]) == "--ic-block-drop") {
            
            nInput = 1200;
            
        }
        
        if (string(argv[i]) == "--ic-droplet") {
            
            nInput = 360;
            
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
    
    MPI_Init(&argc, &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProcessors);
    MPI_Comm_rank(MPI_COMM_WORLD, &currRank);
    
    // cout << "" << endl;
    // cout << "Number of Processors: " << numOfProcessors << endl;
    // cout << "Current Rank: " << currRank << endl;
    
    // cout << "MPI Process Initiated..." << endl;
    
    // Initialise problem as SPH object
    // SPH(const unsigned int& numOfParticles, const double& timeStep, const double& finalT, const double& radOfInfl)
    SPH obj1(nInput, dtInput, TInput, hInput);
    
    // ofstream outputPP("data.txt", ios::out | ios::trunc);
    ofstream xOut("output.txt", ios::out | ios::trunc);
    xOut.precision(10);
    ofstream energyOut("energy.txt", ios::out | ios::trunc);
    energyOut.precision(10);
    ofstream dataOut("data.txt", ios::out | ios::trunc);
    energyOut.precision(10);
    
    localN = obj1.getN() / numOfProcessors;
    
    // cout << "Local N: " << localN << endl;
    
    /*obj1.calcDensityWithMPI(localN, numOfProcessors, currRank);
    obj1.calcPressureWithMPI(localN, numOfProcessors, currRank);
    obj1.calcPressureForceWithMPI(localN, numOfProcessors, currRank);
    obj1.calcViscousForceWithMPI(localN, numOfProcessors, currRank);
    obj1.calcGravityForceWithMPI(localN, numOfProcessors, currRank);    
    obj1.calcAcceleration();
    obj1.generateVInit();
    obj1.getNextParticlePos();*/
    
    obj1.iterate(xOut, energyOut, dataOut, localN, numOfProcessors, currRank);
    
    MPI_Finalize();
    
    // cout << "MPI Process Terminated..." << endl;
    
    // obj1.writeToPPOutputFile(outputPP);
    
    obj1.closePPOutputFile(xOut);
    obj1.closePPOutputFile(energyOut);
    
    /*
    SPH obj1(12, 0.001, 10, 0.01);
    // SPH obj1(2, 0.0001, 0.001, 0.01);
    
    ofstream outputPP("data.txt", ios::out | ios::trunc);
    // ofstream outputPP("allEnergy.txt", ios::out | ios::trunc);
    outputPP.precision(10);
    
    obj1.writeToPPOutputFile(outputPP);
    
    // Iterate from t = dt to t = T  
    obj1.iterate(outputPP);
    
    while (obj1.getCurrT() <= obj1.getTotalIntTime()) {
        
        obj1.calcDensity();
        obj1.calcPressure();
        obj1.calcPressureForce();
        obj1.calcViscousForce();
        obj1.calcGravityForce();
        obj1.calcAcceleration();
        obj1.getNextParticleVel();
        
        // cout << "Particle Location (@t = " << obj1.getCurrT() << "): " << endl;
        
        obj1.getNextParticlePos();
        
        obj1.writeToPPOutputFile(outputPP);
        
        // obj1.applyBC();
        
        // cout << obj1.getCurrT() << " " << obj1.calcKineticEnergy() << " " << obj1.calcPotentialEnergy() << " " << obj1.calcTotalEnergy() << endl;
        
        obj1.setCurrT();
        
    }
    
    
    obj1.writeToPPOutputFile(outputPP);
    
    obj1.closePPOutputFile(outputPP);
    */

    return 0;
    
}