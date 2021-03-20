#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "cblas.h"
#include "mpi.h"
#include "SPH.h"
using namespace std;

SPH::SPH(const unsigned int& numOfParticles, const double& timeStep, const double& finalT, const double& radOfInfl) {
    
    cout << "Starting computation..." << endl;
    
    h = radOfInfl;             // Coefficient of restitution
    dt = timeStep;             // Time-step delta_t = 10^-4
    T = finalT;
    N = numOfParticles;
    
    t = dt;
    
    x = new double[2 * N]();
    rho = new double[N]();     // Particle density array [rhox0, rhoy0, rhox1, rhoy1]
    p = new double[N]();       // Particle pressure array
    v = new double[2 * N]();   // Particle velocity array
    a = new double[2 * N]();   // Particle velocity array
    F_p = new double[2 * N](); // Pressure force array
    F_v = new double[2 * N](); // Viscous force array
    F_g = new double[2 * N](); // Gravity force array
    
    getExecCase(N);
    
    cout << "Initialised particle location: " << endl;
    for (unsigned int i = 0; i < N; ++i) {
        cout << x[2*i] << " " << x[2*i+1] << endl;
    }
    
    calcDensityInit();
    m = scaleMass();
    calcDensity();
    calcPressure();
    calcPressureForce();
    calcViscousForce();
    calcGravityForce();
    calcAcceleration();
    generateVInit();
    getNextParticlePos();
    applyBC();
    
    cout << "Particle Location (@t = dt): " << endl;
    for (unsigned int i = 0; i < N; ++i) {
        cout << x[2*i] << " " << x[2*i+1] << endl;
    }
    
}

SPH::~SPH() {
    
    cout << "Concluding computation..." << endl;
    
    delete[] rho;
    delete[] x;
    delete[] p;
    delete[] v;
    delete[] a;
    delete[] F_p;
    delete[] F_v;
    delete[] F_g;
    
    cout << "Program Terminated..." << endl;
    
}
/*
void SPH::createPPOutputFile() {
    
    ofstream vOut("data.txt", ios::out | ios::trunc);
    vOut.precision(10);
    
}

void SPH::writeToPPOutputFile() {
    
    vOut << 
    for (unsigned int i = 0; i < N; ++i) {
        
        x[i] << " " <<
        
    }
    endl;
    
}
    
void SPH::closePPOutputFile() {
    
    vOut.close();
    
}
*/
int SPH::getN() {
    
    return N;
    
}

double SPH::getCurrT() {
    
    return t;
    
}

void SPH::setCurrT() {
    
    t += dt;
    
}

double SPH::getDeltaT() {
    
    return dt;
    
}

double SPH::getTotalIntTime() {
    
    return T;
    
}

void SPH::getExecCase(unsigned int caseID) {
    
    switch(caseID) {
        
        case 1:
            x[0] = 0.5;
            x[1] = 0.5;
            break;
            
        case 2:
            x[0] = 0.5;
            x[1] = 0.5;
            x[2] = 0.5;
            x[3] = h;
            break;
            
        case 3:
            x[0] = 0.5;
            x[1] = 0.5;
            x[2] = 0.495;
            x[3] = h;
            x[4] = 0.505;
            x[5] = h;
            break;
            
        case 4:
            x[0] = 0.505;
            x[1] = 0.5;
            x[3] = 0.515;
            x[4] = 0.5;
            x[5] = 0.51;
            x[6] = 0.45;
            x[7] = 0.5;
            x[8] = 0.45;
            break;
        
    }
    
}

void SPH::calcDensity() {
    
    double* r_ij = new double[2*N];
    double* q = new double[N];
    
    for (unsigned int i = 0; i < N; ++i) {
        
        // Compute r_ij
        for (unsigned int j = 0; j < N; ++j) {
            
            r_ij[2*j] = x[2*i] - x[2*j];
        
            r_ij[2*j+1] = x[2*i+1] - x[2*j+1];
            
        }
        
        // Compute q
        for (unsigned int k = 0; k < N; ++k) {
            
            q[k] = sqrt((r_ij[2*k]*r_ij[2*k] + r_ij[2*k+1]*r_ij[2*k+1]) / h);
            
        }
        
        // Compute rho_i
        for (unsigned int j = 0; j < N; ++j) {
            
            if (q[j] < 1.0) {
                
                rho[i] += ((4.0 * m) / (M_PI * h * h)) * ((1 - (q[j] * q[j])) * (1 - (q[j] * q[j])) * (1 - (q[j] * q[j])));
                
            } else {
                
                rho[i] += 0.0;
            }
        }
        
        // cout << rho[i] << endl;
        
    }
    
    delete[] r_ij;
    delete[] q;
    
}

void SPH::calcDensityInit() {
    
    double* r_ij = new double[2*N];
    double* q = new double[N];
    
    for (unsigned int i = 0; i < N; ++i) {
        
        // Compute r_ij
        for (unsigned int j = 0; j < N; ++j) {
            
            r_ij[2*j] = x[2*i] - x[2*j];
        
            r_ij[2*j+1] = x[2*i+1] - x[2*j+1];
            
        }
        
        // Compute q
        for (unsigned int k = 0; k < N; ++k) {
            
            q[k] = sqrt((r_ij[2*k]*r_ij[2*k] + r_ij[2*k+1]*r_ij[2*k+1]) / h);
            
        }
        
        // Compute rho_i
        for (unsigned int j = 0; j < N; ++j) {
            
            if (q[j] < 1.0) {
                
                rho[i] += ((4.0 * m_init) / (M_PI * h * h)) * ((1 - (q[j] * q[j])) * (1 - (q[j] * q[j])) * (1 - (q[j] * q[j])));
                
            } else {
                
                rho[i] += 0.0;
            }
        }
        
        // cout << rho[i] << endl;
        
    }
    
    delete[] r_ij;
    delete[] q;
    
}

void SPH::calcPressure() {
    
    for (unsigned int i = 0; i < N; ++i) {
      
        p[i] = k * (rho[i] - rho_0);
        
    } 
    
    // cout << "Pressure Calculation Complete..." << endl;
   
}

void SPH::calcPressureForce() {
    
    double* r_ij = new double[2*N];
    double* q = new double[N];
    
    for (unsigned int i = 0; i < N; ++i) {
        
        // Compute r_ij
        for (unsigned int j = 0; j < N; ++j) {
            
            r_ij[2*j] = x[2*i] - x[2*j];
        
            r_ij[2*j+1] = x[2*i+1] - x[2*j+1];
            
        }
        
        // Compute q
        for (unsigned int k = 0; k < N; ++k) {
            
            q[k] = sqrt((r_ij[2*k]*r_ij[2*k] + r_ij[2*k+1]*r_ij[2*k+1]) / h);
            
        }
        
        // Compute F^p_i
        for (unsigned int j = 0; j < N; ++j) {
            
            if (q[j] < 1.0 && i != j) {
                
                F_p[2*i] += -((m * (p[i] + p[j])) / (rho[j] * 2.0)) * (((-30.0 * r_ij[2*j]) / (M_PI * h * h * h)) * (((1.0 - q[j]) * (1.0 - q[j])) / q[j]));
                F_p[2*i+1] += -((m * (p[i] + p[j])) / (rho[j] * 2.0)) * (((-30.0 * r_ij[2*j+1]) / (M_PI * h * h * h)) * (((1.0 - q[j]) * (1.0 - q[j])) / q[j]));
            
            } else {
            
                F_p[2*i] += 0.0;
                F_p[2*i+1] += 0.0;
            }
        }
        
    }
    
    delete[] r_ij;
    delete[] q;
    
}

void SPH::calcViscousForce() {
    
    double* r_ij = new double[2*N];
    double* v_ij = new double[2*N];
    double* q = new double[N];
    
    for (unsigned int i = 0; i < N; ++i) {
        
        // Compute r_ij
        for (unsigned int j = 0; j < N; ++j) {
            
            r_ij[2*j] = x[2*i] - x[2*j];
        
            r_ij[2*j+1] = x[2*i+1] - x[2*j+1];
            
        }
        
        // Compute v_ij
        for (unsigned int l = 0; l < N; ++l) {
            
            v_ij[2*l] = v[2*i] - v[2*l];
        
            v_ij[2*l+1] = v[2*i+1] - v[2*l+1];
            
        }
        
        // Compute q
        for (unsigned int k = 0; k < N; ++k) {
            
            q[k] = sqrt((r_ij[2*k]*r_ij[2*k] + r_ij[2*k+1]*r_ij[2*k+1]) / h);
            
        }
        
        // Compute F^v_i
        for (unsigned int j = 0; j < N; ++j) {
            
            if (q[j] < 1.0 && i != j) {
                
                F_v[2*i] += ((m * v_ij[2*j]) / rho[j]) * ((40 / (M_PI * h * h * h * h)) * (1 - q[j]));
                F_v[2*i+1] += ((m * v_ij[2*j+1]) / rho[j]) * ((40 / (M_PI * h * h * h * h)) * (1 - q[j]));;
            
            } else {
            
                F_v[2*i] += 0.0;
                F_v[2*i+1] += 0.0;
            }
        }
        
        F_v[2*i] *= -mu;
        F_v[2*i+1] *= -mu;
        
    }
    
    delete[] r_ij;
    delete[] v_ij;
    delete[] q;
    
}

void SPH::calcGravityForce() {
    
    for (unsigned int i = 0; i < N; ++i) {
        
        F_g[2*i] = 0.0;
        
        F_g[2*i+1] = -rho[i] * g;
        
    }
    
}

void SPH::calcAcceleration() {
    
    for (unsigned int i = 0; i < N; ++i) {
        
        a[2*i] = (F_p[2*i] + F_v[2*i] + F_g[2*i]) / rho[i];
        
        a[2*i+1] = (F_p[2*i+1] + F_v[2*i+1] + F_g[2*i+1]) / rho[i];
        
    }
    
}

void SPH::generateVInit() {
    
    for (unsigned int i = 0; i < N; ++i) {
        
        v[2*i] = 0 + (a[2*i] * (dt / 2));
        
        v[2*i+1] = 0 + (a[2*i+1] * (dt / 2));
        
    }
    
}

void SPH::getNextParticleVel() {
    
    for (unsigned int i = 0; i < N; ++i) {
        
        v[2*i] += a[2*i] * dt;
        
        v[2*i+1] += a[2*i+1] * dt;
        
    }
    
}

void SPH::getNextParticlePos() {
    
    for (unsigned int i = 0; i < N; ++i) {
        
        x[2*i] += v[2*i] * dt;
        
        x[2*i+1] += v[2*i+1] * dt;
        
    }
    
    for (unsigned int i = 0; i < N; ++i) {
        cout << x[2*i] << " " << x[2*i+1] << endl;
    }
    
}

void SPH::applyBC() {
    
    for (unsigned int i = 0; i < N; ++i) {
        
        if (x[2*i] < h) {
            
            x[2*i] = h;
            
            v[2*i] *= -e; 
            
        } else if (x[2*i] > (1 - h)) {
            
            x[2*i] = 1 - h;
            
            v[2*i] *= -e;
            
        }
        
        if (x[2*i+1] < h) {
            
            x[2*i+1] = h;
            
            v[2*i+1] *= -e; 
            
        } else if (x[2*i+1] > (1 - h)) {
            
            x[2*i+1] = 1 - h;
            
            v[2*i+1] *= -e;
            
        }
        
    }
    
}

double SPH::scaleMass() {
    
    double sumRho = 0.0;
    double scaledMass = 0.0;
    
    for (unsigned int i = 0; i < N; ++i) {
        
        sumRho += rho[i];
        
    }
    
    scaledMass = (N * rho_0) / sumRho;
    
    return scaledMass;
    
}