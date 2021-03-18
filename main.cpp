#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;

double k = 2000.0; // Gas constant
double rho_0 = 1000.0; // Resting density
double mu = 1.0; // Viscosity
double g = 9.81; // Acceleration due to gravity
double h = 0.01; // Radius of influence
double e = 0.8; // Coefficient of restitution
double dt = 0.0001; // Time-step delta_t = 10^-4

void calcDensity(double* rho, double* x, double h, double m, int N);
void calcPressure(double* p, double k, double* rho, double rho_0, int N);
void calcPressureForce(double* F_p, double* x, double* p, double* rho, double h, double m, int N);
void calcViscousForce(double* F_v, double* x, double* v, double* rho, double mu, double h, double m, int N);
void calcGravityForce(double* F_g, double* rho, double g, int N);
void calcAcceleration(double* a, double* F_p, double* F_v, double* F_g, double* rho, int N);
void generateVInit(double* v, double* a, double dt, int N);
void getNextParticleVel(double* v, double* a, double dt, int N);
void getNextParticlePos(double* x, double* v, double dt, int N);
void applyBC(double* x, double* v, double h, double e, int N);
double scaleMass(double* rho, double rho_0, int N);

int main() {
    
    int N = 8; // Number of particles
    double m_init = 1.0; // Initialised particle mass
    double* x = new double[2 * N](); // Particle velocity array
    double* rho = new double[N](); // Particle density array [rhox0, rhoy0, rhox1, rhoy1]
    double* p = new double[N](); // Particle pressure array
    double* v = new double[2 * N](); // Particle velocity array
    double* a = new double[2 * N](); // Particle velocity array
    double* F_p = new double[2 * N](); // Pressure force array
    double* F_v = new double[2 * N](); // Viscous force array
    double* F_g = new double[2 * N](); // Gravity force array
    
    
    srand(time(0));
    
    for (int i = 0; i < N; ++i) {
        
        x[2*i] = (double)rand()/RAND_MAX * 0.1 + 0.5;
        
        x[2*i+1] = (double)rand()/RAND_MAX * 0.1 + 0.5;
        
    }
    
    
    // double x[8] = {0.505, 0.5, 0.515, 0.5, 0.51, 0.45, 0.5, 0.45};

    // double x[8] = {0.50, 0.50, 0.50, 0.501, 0.50, 0.502, 0.50, 0.53};
    
    // double x[8] = {0.501, 0.50, 0.503, 0.50, 0.505, 0.50, 0.507, 0.50};
    
    // double x[16] = {0.10, 0.50, 0.20, 0.50, 0.30, 0.50, 0.40, 0.50, 0.50, 0.50, 0.60, 0.50, 0.70, 0.50, 0.80, 0.50};
    
    // double x[2] = {0.50, 0.50};
    
    // double x[4] = {0.5, 0.5, 0.5, h};
    
    cout << "Particle Position Array: " << endl;
    for (int i = 0; i < N; ++i) {
        cout << x[2*i] << " " << x[2*i+1] << endl;
    }
    
    ofstream vOut("data.txt", ios::out | ios::trunc);
    vOut.precision(10);
    
    
    vOut << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4] << " " << x[5] << " " << x[6] << " " << x[7] << " " <<
    x[8] << " " << x[9] << " " << x[10] << " " << x[11] << " " << x[12] << " " << x[13] << " " << x[14] << " " << x[15] << " " << endl;
    
    
    // vOut << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << endl;
    
    // vOut << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4] << " " << x[5] << " " << x[6] << " " << x[7] << " " << endl;
    
    // First time-step
    calcDensity(rho, x, h, m_init, N);
    
    double m = scaleMass(rho, rho_0, N);
        
    calcPressure(p, k, rho, rho_0, N);
    
    calcPressureForce(F_p, x, p, rho, h, m, N);
    
    calcViscousForce(F_v, x, v, rho, mu, h, m, N);
    
    calcGravityForce(F_g, rho, g, N);
    
    calcAcceleration(a, F_p, F_v, F_g, rho, N);
    
    generateVInit(v, a, dt, N);
    
    getNextParticleVel(v, a, dt, N);
    
    getNextParticlePos(x, v, dt, N);
    
    
    vOut << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4] << " " << x[5] << " " << x[6] << " " << x[7] << " " <<
    x[8] << " " << x[9] << " " << x[10] << " " << x[11] << " " << x[12] << " " << x[13] << " " << x[14] << " " << x[15] << " " << endl;
    
    
    // vOut << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4] << " " << x[5] << " " << x[6] << " " << x[7] << " " << endl;
    
    // Remaining time-step
    for (int t = 1; t < 1899; ++t) {
        
        calcDensity(rho, x, h, m, N);
    
        m = scaleMass(rho, rho_0, N);
            
        calcPressure(p, k, rho, rho_0, N);
        
        calcPressureForce(F_p, x, p, rho, h, m, N);
        /*
        cout << "Pressure Force Array: " << endl;
        for (int i = 0; i < N; ++i) {
            cout << F_p[2*i] << " " << F_p[2*i+1] << endl;
        }
        */
        calcViscousForce(F_v, x, v, rho, mu, h, m, N);
        /*
        cout << "Viscous Force Array: " << endl;
        for (int i = 0; i < N; ++i) {
            cout << F_v[2*i] << " " << F_v[2*i+1] << endl;
        }
        */
        calcGravityForce(F_g, rho, g, N);
        
        calcAcceleration(a, F_p, F_v, F_g, rho, N);
        /*
        cout << "Acceleration (@ next time-step) Array: " << endl;
        for (int i = 0; i < N; ++i) {
            cout << a[2*i] << " " << a[2*i+1] << endl;
        }
        */
        getNextParticleVel(v, a, dt, N);
        
        cout << "Velocity (@ next time-step) Array: " << endl;
        for (int i = 0; i < N; ++i) {
            cout << v[2*i] << " " << v[2*i+1] << endl;
        }
        
        getNextParticlePos(x, v, dt, N);
        
        cout << "Position (@ next time-step) Array: " << endl;
        for (int i = 0; i < N; ++i) {
            cout << x[2*i] << " " << x[2*i+1] << endl;
        }
        
        applyBC(x, v, h, e, N);
        
        
        vOut << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4] << " " << x[5] << " " << x[6] << " " << x[7] << " " <<
        x[8] << " " << x[9] << " " << x[10] << " " << x[11] << " " << x[12] << " " << x[13] << " " << x[14] << " " << x[15] << " " << endl;
        
        
        // vOut << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << endl;
        
        // vOut << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4] << " " << x[5] << " " << x[6] << " " << x[7] << " " << endl;
        
    }
    
    vOut.close();
    
    /*
    cout << "Particle Density Array: " << endl;
    for (int i = 0; i < N; ++i) {
        cout << rho[i] << endl;
    }
    cout << "Particle Pressure Array: " << endl;
    for (int i = 0; i < N; ++i) {
        cout << p[i] << endl;
    }
    cout << "Pressure Force Array: " << endl;
    for (int i = 0; i < N; ++i) {
        cout << F_p[2*i] << " " << F_p[2*i+1] << endl;
    }
    cout << "Viscous Force Array: " << endl;
    for (int i = 0; i < N; ++i) {
        cout << F_v[2*i] << " " << F_v[2*i+1] << endl;
    }
    cout << "Gravity Force Array: " << endl;
    for (int i = 0; i < N; ++i) {
        cout << F_g[2*i] << " " << F_g[2*i+1] << endl;
    }
    cout << "Acceleration Array: " << endl;
    for (int i = 0; i < N; ++i) {
        cout << a[2*i] << " " << a[2*i+1] << endl;
    }
    cout << "Initial Velocity Array: " << endl;
    for (int i = 0; i < N; ++i) {
        cout << v[2*i] << " " << v[2*i+1] << endl;
    }
    cout << "Velocity (@ next time-step) Array: " << endl;
    for (int i = 0; i < N; ++i) {
        cout << v[2*i] << " " << v[2*i+1] << endl;
    }
    
    cout << "Position (@ next time-step) Array: " << endl;
    for (int i = 0; i < N; ++i) {
        cout << x[2*i] << " " << x[2*i+1] << endl;
    }
    */
    
    delete[] rho;
    delete[] p;
    delete[] v;
    delete[] a;
    delete[] F_p;
    delete[] F_v;
    delete[] F_g;
    
    return 0;
    
}

void calcRij(double* r_ij, double* x, int i, int N) {
    
    for (int j = 0; j < N; ++j) {
        
        r_ij[2*j] = x[2*i] - x[2*j];
        
        r_ij[2*j+1] = x[2*i+1] - x[2*j+1];
        
    }
    
}

void calcQ(double* r_ij, double* q, double h, int N) {
    
    for (int i = 0; i < N; ++i) {
        
        q[i] = sqrt((r_ij[2*i]*r_ij[2*i] + r_ij[2*i+1]*r_ij[2*i+1]) / h);
        
    }
    
}

void calcDensity(double* rho, double* x, double h, double m, int N) {
    
    double* r_ij = new double[2*N];
    double* q = new double[N];
    
    for (int i = 0; i < N; ++i) {
        
        // Compute r_ij
        for (int j = 0; j < N; ++j) {
            
            r_ij[2*j] = x[2*i] - x[2*j];
        
            r_ij[2*j+1] = x[2*i+1] - x[2*j+1];
            
        }
        
        // Compute q
        for (int k = 0; k < N; ++k) {
            
            q[k] = sqrt((r_ij[2*k]*r_ij[2*k] + r_ij[2*k+1]*r_ij[2*k+1]) / h);
            
        }
        
        // Compute rho_i
        for (int j = 0; j < N; ++j) {
            
            if (q[j] < 1.0) {
                
                rho[i] += ((4.0 * m) / (M_PI * h * h)) * ((1 - (q[j] * q[j])) * (1 - (q[j] * q[j])) * (1 - (q[j] * q[j])));
                
            } else {
                
                rho[i] += 0.0;
            }
        }
        
    }
    
    delete[] r_ij;
    delete[] q;
    
}

void calcPressure(double* p, double k, double* rho, double rho_0, int N) {
    
    for (int i = 0; i < N; ++i) {
      
        p[i] = k * (rho[i] - rho_0);
        
    } 
   
}

void calcPressureForce(double* F_p, double* x, double* p, double* rho, double h, double m, int N) {
    
    double* r_ij = new double[2*N];
    double* q = new double[N];
    
    for (int i = 0; i < N; ++i) {
        
        // Compute r_ij
        for (int j = 0; j < N; ++j) {
            
            r_ij[2*j] = x[2*i] - x[2*j];
        
            r_ij[2*j+1] = x[2*i+1] - x[2*j+1];
            
        }
        
        // Compute q
        for (int k = 0; k < N; ++k) {
            
            q[k] = sqrt((r_ij[2*k]*r_ij[2*k] + r_ij[2*k+1]*r_ij[2*k+1]) / h);
            
        }
        
        // Compute F^p_i
        for (int j = 0; j < N; ++j) {
            
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

void calcViscousForce(double* F_v, double* x, double* v, double* rho, double mu, double h, double m, int N) {
    
    double* r_ij = new double[2*N];
    double* v_ij = new double[2*N];
    double* q = new double[N];
    
    for (int i = 0; i < N; ++i) {
        
        // Compute r_ij
        for (int j = 0; j < N; ++j) {
            
            r_ij[2*j] = x[2*i] - x[2*j];
        
            r_ij[2*j+1] = x[2*i+1] - x[2*j+1];
            
        }
        
        // Compute v_ij
        for (int l = 0; l < N; ++l) {
            
            v_ij[2*l] = v[2*i] - v[2*l];
        
            v_ij[2*l+1] = v[2*i+1] - v[2*l+1];
            
        }
        
        // Compute q
        for (int k = 0; k < N; ++k) {
            
            q[k] = sqrt((r_ij[2*k]*r_ij[2*k] + r_ij[2*k+1]*r_ij[2*k+1]) / h);
            
        }
        
        // Compute F^v_i
        for (int j = 0; j < N; ++j) {
            
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

void calcGravityForce(double* F_g, double* rho, double g, int N) {
    
    for (int i = 0; i < N; ++i) {
        
        F_g[2*i] = 0.0;
        
        F_g[2*i+1] = -rho[i] * g;
        
    }
    
}

void calcAcceleration(double* a, double* F_p, double* F_v, double* F_g, double* rho, int N) {
    
    for (int i = 0; i < N; ++i) {
        
        a[2*i] = (F_p[2*i] + F_v[2*i] + F_g[2*i]) / rho[i];
        
        a[2*i+1] = (F_p[2*i+1] + F_v[2*i+1] + F_g[2*i+1]) / rho[i];
        
    }
    
}

void generateVInit(double* v, double* a, double dt, int N) {
    
    for (int i = 0; i < N; ++i) {
        
        v[2*i] = 0 + (a[2*i] * (dt / 2));
        
        v[2*i+1] = 0 + (a[2*i+1] * (dt / 2));
        
    }
    
}

void getNextParticleVel(double* v, double* a, double dt, int N) {
    
    for (int i = 0; i < N; ++i) {
        
        v[2*i] += a[2*i] * dt;
        
        v[2*i+1] += a[2*i+1] * dt;
        
    }
    
}

void getNextParticlePos(double* x, double* v, double dt, int N) {
    
    for (int i = 0; i < N; ++i) {
        
        x[2*i] += v[2*i] * dt;
        
        x[2*i+1] += v[2*i+1] * dt;
        
    }
    
}

void applyBC(double* x, double* v, double h, double e, int N) {
    
    for (int i = 0; i < N; ++i) {
        
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

double scaleMass(double* rho, double rho_0, int N) {
    
    double sumRho = 0.0;
    
    for (int i = 0; i < N; ++i) {
        
        sumRho += rho[i];
        
    }
    
    return (N * rho_0) / sumRho;
    
}
