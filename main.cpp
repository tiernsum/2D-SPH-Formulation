#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "cblas.h"
#include "mpi.h"
#include "SPH.h"
using namespace std;

int main(int argc, char* argv[]) {
    
    /*
    cout << argc << endl;
    
    for (int i = 0; i < argc; ++i) {
        
        cout << argv[i] << endl;
    }
    */
    
    // double x[16] = {0.00, 0.00, 0.10, 0.00, 0.20, 0.00, 0.00, 0.10, 0.10, 0.10, 0.20, 0.10, 0.00, 0.20, 0.20, 0.20};
    
    // Initialise problem as SPH object
    // SPH(const unsigned int& numOfParticles, const double& timeStep, const double& finalT, const double& radOfInfl)
    SPH obj1(4, 0.0001, 10, 0.01);
    
    // ofstream outputPP("output.txt", ios::out | ios::trunc);
    ofstream outputPP("data.txt", ios::out | ios::trunc);
    outputPP.precision(10);
    
    obj1.writeToPPOutputFile(outputPP);
    
    // Iterate from t = dt to t = T    
    while (obj1.getCurrT() <= obj1.getTotalIntTime()) {
        
        obj1.calcDensity();
        obj1.calcPressure();
        obj1.calcPressureForce();
        obj1.calcViscousForce();
        obj1.calcGravityForce();
        obj1.calcAcceleration();
        obj1.getNextParticleVel();
        
        cout << "Particle Location (@t = " << obj1.getCurrT() << "): " << endl;
        
        obj1.getNextParticlePos();
        
        obj1.writeToPPOutputFile(outputPP);
        
        obj1.applyBC();
        
        obj1.setCurrT();
        
    }
    
    // obj1.writeToPPOutputFile(outputPP);
    
    obj1.closePPOutputFile(outputPP);
    
    return 0;
    
}