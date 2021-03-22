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
    rhoInit = new double[N]();
    p = new double[N]();       // Particle pressure array
    v = new double[2 * N]();   // Particle velocity array
    a = new double[2 * N]();   // Particle velocity array
    F_p = new double[2 * N](); // Pressure force array
    F_v = new double[2 * N](); // Viscous force array
    F_g = new double[2 * N](); // Gravity force array
    r_ij = new double[2 * N]();
    v_ij = new double[2 * N]();
    
    getExecCase(N);
    /*
    cout << "Initialised particle location: " << endl;
    for (unsigned int i = 0; i < N; ++i) {
        cout << x[2*i] << " " << x[2*i+1] << endl;
    }
    */
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
    // applyBC();
    /*
    cout << "Particle Location (@t = dt): " << endl;
    for (unsigned int i = 0; i < N; ++i) {
        cout << x[2*i] << " " << x[2*i+1] << endl;
    }
    */
    
    // cout << t << " " << calcKineticEnergy() << " " << calcPotentialEnergy() << " " << calcTotalEnergy() << endl;
    
}

SPH::~SPH() {
    
    delete[] rho;
    delete[] rhoInit;
    delete[] x;
    delete[] p;
    delete[] v;
    delete[] a;
    delete[] F_p;
    delete[] F_v;
    delete[] F_g;
    delete[] r_ij;
    delete[] v_ij;
    
    cout << "Program Terminated..." << endl;
    
}
/*
void SPH::createPPOutputFile(ofstream& fileName) {
    
    ofstream fileName("particlePosition.txt", ios::out | ios::trunc);
    fileName.precision(10);
    
}
*/
void SPH::writeToPPOutputFile(ofstream& fileName) {
    
    if (fileName.is_open()) {
        /*
        for (unsigned int i = 0; i < N; ++i) {
            
            fileName << x[2*i] << " " << x[2*i+1] << endl;
            
        }
        */
        
        // fileName << x[0] << " " << x[1] << endl;
        // fileName << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << endl;
        // fileName << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4] << " " << x[5] << endl;
        // fileName << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4] << " " << x[5] << " " << x[6] << " " << x[7] << endl;
        // fileName << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4] << " " << x[5] << " " << x[6] << " " << x[7] << " " << x[8] << " " << x[9] << " " << x[10] << " " << x[11] << " " << 
        // x[12] << " " << x[13] << " " << x[14] << " " << x[15] << endl;
        fileName << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4] << " " << x[5] << " " << x[6] << " " << x[7] << " " << x[8] << " " << x[9] << " " << x[10] << " " << x[11] << " " << 
        x[12] << " " << x[13] << " " << x[14] << " " << x[15] << " " << x[16] << " " << x[17] << " " << x[18] << " " << x[19] << " " << x[20] << " " << x[21] << " " << x[22] << " " << x[23] << endl;
        
        // fileName << getCurrT() << " " << calcKineticEnergy() << " " << calcPotentialEnergy() << " " << calcTotalEnergy() << endl;
        
    }
    
}
  
void SPH::closePPOutputFile(ofstream& fileName) {
    
    fileName.close();
    
}

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
            x[2] = 0.515;
            x[3] = 0.5;
            x[4] = 0.51;
            x[5] = 0.45;
            x[6] = 0.5;
            x[7] = 0.45;
            break;
        
        /*
        case 8:
            x[0] = 0.075;
            x[1] = 0.05;
            x[2] = 0.1;
            x[3] = 0.05;
            x[4] = 0.125;
            x[5] = 0.05;
            x[6] = 0.15;
            x[7] = 0.05;
            x[8] = 0.075;
            x[9] = 0.15;
            x[10] = 0.1;
            x[11] = 0.15;
            x[12] = 0.125;
            x[13] = 0.15;
            x[14] = 0.15;
            x[15] = 0.15;
            break;
        */
        
        case 8:
            x[0] = 0.075;
            x[1] = 0.05;
            x[2] = 0.1;
            x[3] = 0.05;
            x[4] = 0.12;
            x[5] = 0.05;
            x[6] = 0.15;
            x[7] = 0.05;
            x[8] = 0.075;
            x[9] = 0.15;
            x[10] = 0.1;
            x[11] = 0.15;
            x[12] = 0.125;
            x[13] = 0.15;
            x[14] = 0.15;
            x[15] = 0.15;
            break;
        /*    
        case 12:
            x[0] = 0.10;
            x[1] = 0.30;
            x[2] = 0.20;
            x[3] = 0.30;
            x[4] = 0.30;
            x[5] = 0.30;
            x[6] = 0.10;
            x[7] = 0.40;
            x[8] = 0.20;
            x[9] = 0.40;
            x[10] = 0.30;
            x[11] = 0.40;
            x[12] = 0.10;
            x[13] = 0.49;
            x[14] = 0.20;
            x[15] = 0.50;
            x[16] = 0.30;
            x[17] = 0.50;
            x[18] = 0.10;
            x[19] = 0.60;
            x[20] = 0.20;
            x[21] = 0.60;
            x[22] = 0.33;
            x[23] = 0.60;
            break;
        */
        
        case 12:
            srand(time(0));
            for (unsigned int i = 0; i < 12; ++i) {
                
                x[2*i] = (double)rand()/RAND_MAX * 0.2 + 0.1;
        
                x[2*i+1] = (double)rand()/RAND_MAX * 0.3 + 0.3;
                
            }
        
    }
    
}
/*
void SPH::calcRij(unsigned int i, unsigned int j) {
    
    for (unsigned int j = 0; j < N; ++j) {
        
        r_ij[2*j] = x[2*i] - x[2*j];
        
        r_ij[2*j+1] = x[2*i+1] - x[2*j+1];
        
    }
    
}
*/

void SPH::calcRij(unsigned int i, unsigned int j) {
    
    r_ij[2*j] = x[2*i] - x[2*j];
    
    r_ij[2*j+1] = x[2*i+1] - x[2*j+1];
    
}

void SPH::calcVij(unsigned int i, unsigned int j) {
    
    v_ij[2*j] = v[2*i] - v[2*j];
    
    v_ij[2*j+1] = v[2*i+1] - v[2*j+1];
    
}
/*
void SPH::calcDensity() {
    
    // double* r_ij = new double[2*N];
    double* q = new double[N];
    
    for (unsigned int i = 0; i < N; ++i) {
        
        // Compute r_ij
        for (unsigned int j = 0; j < N; ++j) {
            
            r_ij[2*j] = x[2*i] - x[2*j];
        
            r_ij[2*j+1] = x[2*i+1] - x[2*j+1];
            
        }
        
        // calcRij(i);
    
        // Compute q
        for (unsigned int k = 0; k < N; ++k) {
            
            q[k] = sqrt(r_ij[2*k]*r_ij[2*k] + r_ij[2*k+1]*r_ij[2*k+1]) / h;
            
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
    
    
    cout << "Density at t = " << t << endl;
    for (unsigned int i = 0; i < N; ++i) {
        cout << rho[i] << endl;
    }
    
    
    // delete[] r_ij;
    delete[] q;
    
}

void SPH::calcDensityInit() {
    
    // double* r_ij = new double[2*N];
    double* q = new double[N];
    
    for (unsigned int i = 0; i < N; ++i) {
        
        // Compute r_ij
        for (unsigned int j = 0; j < N; ++j) {
            
            r_ij[2*j] = x[2*i] - x[2*j];
        
            r_ij[2*j+1] = x[2*i+1] - x[2*j+1];
            
        }
      
        // calcRij(i);
 
        // Compute q
        for (unsigned int k = 0; k < N; ++k) {
            
            q[k] = sqrt(r_ij[2*k]*r_ij[2*k] + r_ij[2*k+1]*r_ij[2*k+1]) / h;
            
        }
        
        // Compute rho_i
        for (unsigned int j = 0; j < N; ++j) {
            
            if (q[j] < 1.0) {
                
                rhoInit[i] += ((4.0 * m_init) / (M_PI * h * h)) * ((1 - (q[j] * q[j])) * (1 - (q[j] * q[j])) * (1 - (q[j] * q[j])));
                
            } else {
                
                rhoInit[i] += 0.0;
            }
        }
        
        // cout << rho[i] << endl;
        
    }
    
    // delete[] r_ij;
    delete[] q;
    
}
*/

void SPH::calcDensity() {
    
    double q;
    
    for (unsigned int i = 0; i < N; ++i) {
        
        for (unsigned int j = 0; j < N; ++j) {
            
            calcRij(i, j);
            
            q = sqrt(r_ij[2*j]*r_ij[2*j] + r_ij[2*j+1]*r_ij[2*j+1]) / h;
            
            if (q < 1) {
                
                rho[i] += ((4 * m) / (M_PI * h * h)) * ((1 - (q * q)) * (1 - (q * q)) * (1 - (q * q)));
                
            } else {
                
                rho[i] += 0.0;
                
            }
        
        }
        
    }
    /*
    cout << "Density at t = " << t << endl;
    for (unsigned int i = 0; i < N; ++i) {
        cout << rho[i] << endl;
    }
    */
}

void SPH::calcDensityInit() {
    
    double q;
    
    for (unsigned int i = 0; i < N; ++i) {
        
        for (unsigned int j = 0; j < N; ++j) {
            
            calcRij(i, j);
            
            q = sqrt(r_ij[2*j]*r_ij[2*j] + r_ij[2*j+1]*r_ij[2*j+1]) / h;
            
            if (q < 1) {
                
                rhoInit[i] += ((4 * m_init) / (M_PI * h * h)) * ((1 - (q * q)) * (1 - (q * q)) * (1 - (q * q)));
                
            } else {
                
                rhoInit[i] += 0.0;
                
            }
        
        }
        
    }
    
}

void SPH::calcPressure() {
    
    for (unsigned int i = 0; i < N; ++i) {
      
        p[i] = k * (rho[i] - rho_0);
        
    } 
    /*
    cout << "Pressure at t = " << t << endl;
    for (unsigned int i = 0; i < N; ++i) {
        cout << p[i] << endl;
    }
    */
    // cout << "Pressure Calculation Complete..." << endl;
   
}
/*
void SPH::calcPressureForce() {
    
    // double* r_ij = new double[2*N];
    // double* q = new double[N];
    double q;
    
    for (unsigned int i = 0; i < N; ++i) {
        
        // Compute r_ij
        for (unsigned int j = 0; j < N; ++j) {
            
            r_ij[2*j] = x[2*i] - x[2*j];
        
            r_ij[2*j+1] = x[2*i+1] - x[2*j+1];
            
            cout << "R_ij: " << endl;
            cout << r_ij[2*j] << " " << r_ij[2*j+1] << endl;
            
        }
        
        // calcRij(i);
        
        // Compute q
        for (unsigned int k = 0; k < N; ++k) {
            
            q[k] = sqrt(r_ij[2*k]*r_ij[2*k] + r_ij[2*k+1]*r_ij[2*k+1]) / h;
            
            cout << "q: " << q[k] << endl;
            
        }
        
        // Compute F^p_i
        for (unsigned int j = 0; j < N; ++j) {
            
            r_ij[2*j] = x[2*i] - x[2*j];
        
            r_ij[2*j+1] = x[2*i+1] - x[2*j+1];
            
            q = sqrt(r_ij[2*j]*r_ij[2*j] + r_ij[2*j+1]*r_ij[2*j+1]) / h;
            
            if (q < 1.0 && i != j) {
                
                F_p[2*i] += ((m * (p[i] + p[j])) / (rho[j] * 2.0)) * (((-30.0 * r_ij[2*j]) / (M_PI * h * h * h)) * (((1.0 - q) * (1.0 - q)) / q));
                F_p[2*i+1] += ((m * (p[i] + p[j])) / (rho[j] * 2.0)) * (((-30.0 * r_ij[2*j+1]) / (M_PI * h * h * h)) * (((1.0 - q) * (1.0 - q)) / q));
            
            
            if (q[j] < 1.0 && i != j) {
                
                F_p[2*i] += ((m * (p[i] + p[j])) / (rho[j] * 2.0)) * (((-30.0 * r_ij[2*j]) / (M_PI * h * h * h)) * (((1.0 - q[j]) * (1.0 - q[j])) / q[j]));
                F_p[2*i+1] += ((m * (p[i] + p[j])) / (rho[j] * 2.0)) * (((-30.0 * r_ij[2*j+1]) / (M_PI * h * h * h)) * (((1.0 - q[j]) * (1.0 - q[j])) / q[j]));
            
            } else {
            
                F_p[2*i] += 0.0;
                F_p[2*i+1] += 0.0;
            }
        }
        
        F_p[2*i] *= -1.0;
        F_p[2*i+1] *= -1.0;
        
    }
    
    
    cout << "Pressure Force at t = " << t << endl;
    for (unsigned int i = 0; i < N; ++i) {
        cout << F_p[2*i] << " " << F_p[2*i+1] << endl;
    }
    
    
    // delete[] r_ij;
    // delete[] q;
    
}
*/

void SPH::calcPressureForce() {
    
    double q;
    
    for (unsigned int i = 0; i < N; ++i) {
        
        for (unsigned int j = 0; j < N; ++j) {
        
            calcRij(i, j);
            
            // cout << "R_ij at t = " << t << endl;
            // cout << r_ij[2*j] << " " << r_ij[2*j+1] << endl;
            
            q = sqrt(r_ij[2*j]*r_ij[2*j] + r_ij[2*j+1]*r_ij[2*j+1]) / h;
            
            // cout << "q = " << q << endl;
            
            if (q < 1 && i != j) {
                
                F_p[2*i] += ((m / rho[j]) * ((p[i] + p[j]) / 2.0)) * ((-30 / (M_PI * h * h * h)) * (r_ij[2*j]) * (((1 - q) * (1 - q)) / q));
            
                F_p[2*i+1] += ((m / rho[j]) * ((p[i] + p[j]) / 2.0)) * ((-30 / (M_PI * h * h * h)) * (r_ij[2*j+1]) * (((1 - q) * (1 - q)) / q));
                
            } else {
                
                F_p[2*i] += 0.0;
            
                F_p[2*i+1] += 0.0;
                
            }
        
        }
        
        F_p[2*i] *= -1.0;
        F_p[2*i+1] *= -1.0;
        
    }
    /*
    cout << "Pressure Force at t = " << t << endl;
    for (unsigned int i = 0; i < N; ++i) {
        cout << F_p[2*i] << " " << F_p[2*i+1] << endl;
    }
    */
}
/*
void SPH::calcViscousForce() {
    
    // double* r_ij = new double[2*N];
    double* v_ij = new double[2*N];
    double* q = new double[N];
    
    for (unsigned int i = 0; i < N; ++i) {
        
        // Compute r_ij
        for (unsigned int j = 0; j < N; ++j) {
            
            r_ij[2*j] = x[2*i] - x[2*j];
        
            r_ij[2*j+1] = x[2*i+1] - x[2*j+1];
            
        }
       
        // calcRij(i);
        
        // Compute v_ij
        for (unsigned int l = 0; l < N; ++l) {
            
            v_ij[2*l] = v[2*i] - v[2*l];
        
            v_ij[2*l+1] = v[2*i+1] - v[2*l+1];
            
        }
        
        // Compute q
        for (unsigned int k = 0; k < N; ++k) {
            
            q[k] = sqrt(r_ij[2*k]*r_ij[2*k] + r_ij[2*k+1]*r_ij[2*k+1]) / h;
            
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
    
    
    cout << "Viscous Force at t = " << t << endl;
    for (unsigned int i = 0; i < N; ++i) {
        cout << F_v[2*i] << " " << F_v[2*i+1] << endl;
    }
    
    
    // delete[] r_ij;
    delete[] v_ij;
    delete[] q;
    
}
*/

void SPH::calcViscousForce() {
    
    double q;
    
    for (unsigned int i = 0; i < N; ++i) {
        
        for (unsigned int j = 0; j < N; ++j) {
        
            calcVij(i, j);
            
            calcRij(i, j);
            
            q = sqrt(r_ij[2*j]*r_ij[2*j] + r_ij[2*j+1]*r_ij[2*j+1]) / h;
            
            if (q < 1 && i != j) {
                
                F_v[2*i] += (m / rho[j]) * (v_ij[2*j]) * ((40 / (M_PI * h * h * h * h)) * (1 - q));
                
                F_v[2*i+1] += (m / rho[j]) * (v_ij[2*j+1]) * ((40 / (M_PI * h * h * h * h)) * (1 - q));
                
            } else {
                
                F_v[2*i] += 0.0;
                
                F_v[2*i+1] += 0.0;
                
            }
        
        }
        
        F_v[2*i] *= -mu;
                
        F_v[2*i+1] *= -mu;
        
    }
    /*
    cout << "Viscous Force at t = " << t << endl;
    for (unsigned int i = 0; i < N; ++i) {
        cout << F_v[2*i] << " " << F_v[2*i+1] << endl;
    }
    */
}

void SPH::calcGravityForce() {
    
    for (unsigned int i = 0; i < N; ++i) {
        
        F_g[2*i] = 0.0;
        
        F_g[2*i+1] = -rho[i] * g;
        
    }
    
    /*
    cout << "Gravity Force at t = " << t << endl;
    for (unsigned int i = 0; i < N; ++i) {
        cout << F_g[2*i] << " " << F_g[2*i+1] << endl;
    }
    */
    
}

void SPH::calcAcceleration() {
    
    for (unsigned int i = 0; i < N; ++i) {
        
        a[2*i] = (F_p[2*i] + F_v[2*i] + F_g[2*i]) / rho[i];
        
        a[2*i+1] = (F_p[2*i+1] + F_v[2*i+1] + F_g[2*i+1]) / rho[i];
        
    }
    
    /*
    cout << "Acceleration at t = " << t << endl;
    for (unsigned int i = 0; i < N; ++i) {
        cout << a[2*i] << " " << a[2*i+1] << endl;
    }
    */
    
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
    
    /*
    cout << "Velocity at t = " << t << endl;
    for (unsigned int i = 0; i < N; ++i) {
        cout << v[2*i] << " " << v[2*i+1] << endl;
    }
    */
    
}

void SPH::getNextParticlePos() {
    
    for (unsigned int i = 0; i < N; ++i) {
        
        x[2*i] += v[2*i] * dt;
        
        x[2*i+1] += v[2*i+1] * dt;
        
        applyBC(i);
        
    }
    /*
    cout << "Position at t = " << t << endl;
    for (unsigned int i = 0; i < N; ++i) {
        cout << x[2*i] << " " << x[2*i+1] << endl;
    }
    */
}

void SPH::applyBC(unsigned int i) {
    
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

double SPH::scaleMass() {
    
    double sumRho = 0.0;
    double scaledMass = 0.0;
    
    for (unsigned int i = 0; i < N; ++i) {
        
        sumRho += rhoInit[i];
        
    }
    
    scaledMass = (N * rho_0) / sumRho;
    
    return scaledMass;
    
}

double SPH::calcKineticEnergy() {
    
    double E_k = 0.0;
    
    for (unsigned int i = 0; i < N; ++i) {
        
        E_k += (v[2*i] * v[2*i]) + v[2*i+1] * v[2*i+1];
        
    }
    
    return 0.5 * m * E_k;
    
}

double SPH::calcPotentialEnergy() {
    
    double E_p = 0.0;
    
    for (unsigned int i = 0; i < N; ++i) {
        
        E_p += x[2*i+1];
        
    }
    
    return m * g * E_p;
    
}

double SPH::calcTotalEnergy() {
    
    return calcKineticEnergy() + calcPotentialEnergy();
    
}

void SPH::iterate(ofstream& fileName) {
    
    while (getCurrT() <= getTotalIntTime()) {
        
        for (unsigned int i = 0; i < N; ++i) {
            
            rho[i] = 0.0;
            F_p[2*i] = 0.0;
            F_p[2*i+1] = 0.0;
            F_v[2*i] = 0.0;
            F_v[2*i+1] = 0.0;
            
        }
        
        calcDensity();
        calcPressure();
        calcPressureForce();
        calcViscousForce();
        calcGravityForce();
        calcAcceleration();
        getNextParticleVel();
        
        // cout << "Particle Location (@t = " << obj1.getCurrT() << "): " << endl;
        
        getNextParticlePos();
        
        writeToPPOutputFile(fileName);
        
        // obj1.applyBC();
        
        // cout << obj1.getCurrT() << " " << obj1.calcKineticEnergy() << " " << obj1.calcPotentialEnergy() << " " << obj1.calcTotalEnergy() << endl;
        
        setCurrT();
        
    }
    
    cout << "Position at t = " << t - dt << endl;
    for (unsigned int i = 0; i < N; ++i) {
        cout << x[2*i] << " " << x[2*i+1] << endl;
    }
    
}