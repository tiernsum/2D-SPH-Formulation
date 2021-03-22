#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "cblas.h"
#include "mpi.h"
#include "SPH.h"
// #include "SPHwithMPI.h"
using namespace std;

int main(int argc, char* argv[]) {
    
    /*
    int numOfProcessors, currRank, localN, startN, endN;
    
    MPI_Init(&argc, &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &numOfProcessors);
    MPI_Comm_rank(MPI_COMM_WORLD, &currRank);
    
    cout << "" << endl;
    cout << "Number of Processors: " << numOfProcessors << endl;
    cout << "Current Rank: " << currRank << endl;
    
    // cout << argc << endl;
    
    // for (int i = 0; i < argc; ++i) {
        
    //     cout << argv[i] << endl;
    // }
    
    cout << "MPI Process Initiated..." << endl;
    
    // Initialise problem as SPH object
    // SPH(const unsigned int& numOfParticles, const double& timeStep, const double& finalT, const double& radOfInfl)
    SPHwithMPI obj1(4, 0.0001, 1, 0.01);
    
    localN = obj1.getN() / numOfProcessors;
    
    cout << "Local N: " << localN << endl;
    
    obj1.iterWithMPI(localN, numOfProcessors, currRank);
    
    MPI_Finalize();
    
    cout << "MPI Process Terminated..." << endl;
    */
    
    SPH obj1(12, 0.001, 10, 0.01);
    // SPH obj1(2, 0.0001, 0.001, 0.01);
    
    ofstream outputPP("data.txt", ios::out | ios::trunc);
    // ofstream outputPP("allEnergy.txt", ios::out | ios::trunc);
    outputPP.precision(10);
    
    obj1.writeToPPOutputFile(outputPP);
    
    // Iterate from t = dt to t = T  
    obj1.iterate(outputPP);
    /*
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
    */
    
    obj1.writeToPPOutputFile(outputPP);
    
    obj1.closePPOutputFile(outputPP);
    
    
    return 0;
    
}