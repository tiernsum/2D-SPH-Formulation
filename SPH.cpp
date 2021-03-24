#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "cblas.h"
#include "mpi.h"
#include "SPH.h"
using namespace std;

/**
 * @brief Initialises an SPH object, the relevant matrices used to in the algorithm to solve the problem, getting the correct test/validation 
 * case, initialising particle density and scaling particle mass
 * @param numOfParticles
 * @param timeStep
 * @param finalT
 * @param radOfInfl
 */
SPH::SPH(const unsigned int& numOfParticles, const double& timeStep, const double& finalT, const double& radOfInfl) {
    
    h       = radOfInfl;            // Radius of influence
    dt      = timeStep;             // Time-step delta_t
    T       = finalT;               // Final time
    N       = numOfParticles;       // Number of particles based on the test/validation case selected
    
    x       = new double[2 * N]();  // Initialising particle position vector [x_0, y_0, x_1, y_1, ...]
    v       = new double[2 * N]();  // Initialising particle velocity vector [vx_0, vy_0, vx_1, vy_1, ...]
    a       = new double[2 * N]();  // Initialising particle acceleration vector [ax_0, ay_0, ax_1, ay_1, ...]
    r_ij    = new double[2 * N]();  // Initialising particle distance vector
    v_ij    = new double[2 * N]();  // Initialising particle velocity difference vector
    rho     = new double[N]();      // Initialising particle density vector [rho_0, rho_1, rho_2, ..., rho_N]
    rhoInit = new double[N]();      // Initialising particle density initalisation vector [rho_0, rho_1, rho_2, ..., rho_N]
    p       = new double[N]();      // Initialising particle pressure vector [p_0, p_1, p_2, ..., p_N]
    F_p     = new double[2 * N]();  // Initialising pressure force vector
    F_v     = new double[2 * N]();  // Initialising viscous force vector
    F_g     = new double[2 * N]();  // Initialising gravity force vector
    
    t       = 0.0;                  // Tracker for current time-step within the iteration
    
    getExecCase(N);                 // Obtaining the particle position vector x from seclected test/validation case input
    
    calcDensityInit();              // Initialising particle density based on intial particle position (t = 0)
    m = scaleMass();                // Scaling the mass based on initial particle position such that resting density equals rho_0
    
}

/**
 * @brief Removes arrays declared on heap upon completion of usage to prevent memory leakage
 */
SPH::~SPH() {
    
    delete[] x;
    delete[] v;
    delete[] a;
    delete[] r_ij;
    delete[] v_ij;
    delete[] rho;
    delete[] rhoInit;
    delete[] p;
    delete[] F_p;
    delete[] F_v;
    delete[] F_g;
    
}

/**
 * @brief Selecting the particle position at (t = 0) based on test/validation case selected
 * @param caseID
 */
void SPH::getExecCase(unsigned int caseID) {
    
    switch(caseID) {
        
        // One-particle validation case
        case 1:
            x[0] = 0.5;
            x[1] = 0.5;
            break;
        
        // Two-particles validation case    
        case 2:
            x[0] = 0.5;
            x[1] = 0.5;
            x[2] = 0.5;
            x[3] = h;
            break;
        
        // Three-particles validation case
        case 3:
            x[0] = 0.5;
            x[1] = 0.5;
            x[2] = 0.495;
            x[3] = h;
            x[4] = 0.505;
            x[5] = h;
            break;
        
        // Four-particles validation case
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
        
        // Dam break test case
        case 400:
            {
                srand(time(0));
                
                xCoor = 0.0;
                yCoor = 0.01;
                
                for (unsigned int i = 0; i < 400; ++i) {
                        
                    xCoor += 0.01;
                    
                    if (xCoor > 0.21) {
                        
                        xCoor = 0.01;
                        x[2*i] = xCoor;
                        
                        yCoor += 0.01;
                        x[2*i+1] = yCoor;
                        
                    } else {
                        
                        x[2*i] = xCoor;
                        x[2*i+1] = yCoor;
                        
                    }
                    
                    // Adding noise to the bottom right of the "dam"
                    if (i % 10 == 0 && yCoor == 0.01) {
                    
                        x[2*i] = (double)rand()/RAND_MAX * 0.02 + 0.18;
                        x[2*i+1] = (double)rand()/RAND_MAX * 0.02;
                        
                    }
                        
                }
                break;
            }
        
        // Block drop test case
        case 650:
            {
                xCoor = 0.09;
                yCoor = 0.3;
                
                for (unsigned int i = 0; i < 650; ++i) {
                        
                    xCoor += 0.01;
                    
                    if (xCoor > 0.31) {
                        
                        xCoor = 0.1;
                        x[2*i] = xCoor;
                        
                        yCoor += 0.01;
                        x[2*i+1] = yCoor;
                        
                    } else {
                        
                        x[2*i] = xCoor;
                        x[2*i+1] = yCoor;
                        
                    }
                        
                }
                break;
            }
        
        // Droplet Test Case
        case 360:
            {
                counter = 1;
                for (unsigned int i = 0; i < 360; ++i) {
                    
                    if (i % 36 == 0 && i != 0) {
                        counter++;
                    }
                
                    x[2*i] = counter * 0.01 * cos(i % 36 * 10 * M_PI / 180) + 0.5;
                    x[2*i+1] = counter * 0.01 * sin(i % 36 * 10 * M_PI / 180) + 0.7;
                
                }
                break;
            }
    }
    
}

/**
 * @brief Calculating the distance between the current particle [i] and the rest of the particles [j]
 * @param i
 * @param j
 */
void SPH::calcRij(unsigned int i, unsigned int j) {
    
    r_ij[2*j] = xLocal[2*i] - x[2*j];
    r_ij[2*j+1] = xLocal[2*i+1] - x[2*j+1];
    
}

/**
 * @brief Calculating the velocity difference between the current particle [i] and the rest of the particles [j]
 * @param i
 * @param j
 */
void SPH::calcVij(unsigned int i, unsigned int j) {
    
    v_ij[2*j] = vLocal[2*i] - v[2*j];
    v_ij[2*j+1] = vLocal[2*i+1] - v[2*j+1];
    
}

/**
 * @brief Calculates the initial density of the particles at (t = 0) for subsequent mass scaling
 */
void SPH::calcDensityInit() {
    
    // q is the normalised distance between particles (||r_ij||/h)
    double q;
    
    for (unsigned int i = 0; i < N; ++i) {
        
        for (unsigned int j = 0; j < N; ++j) {
            
            r_ij[2*j] = x[2*i] - x[2*j];
            r_ij[2*j+1] = x[2*i+1] - x[2*j+1];
            
            q = sqrt(r_ij[2*j]*r_ij[2*j] + r_ij[2*j+1]*r_ij[2*j+1]) / h;
            
            if (q < 1) {
                rhoInit[i] += ((4 * m_init) / (M_PI * h * h)) * ((1 - (q * q)) * (1 - (q * q)) * (1 - (q * q)));
            } else {
                rhoInit[i] += 0.0;
            }
        
        }
        
    }
    
}

/**
 * @brief Scales the mass based on the initial density calculated to ensure the particles' resting density equals rho_0
 * @return Scaled mass of particles
 */
double SPH::scaleMass() {
    
    double sumRho = 0.0;
    
    for (unsigned int i = 0; i < N; ++i) {
        
        sumRho += rhoInit[i];
        
    }
    
    return (N * rho_0) / sumRho;
    
}

/**
 * @brief Iterates the algorithm to solve the SPH Formulation to final time, T, with MPI parallelisation
 * @param xCoor     - Writes the particle positions at final time, T, to "output.txt"
 * @param energyTxt - Writes the system's kinetic, potential, and overall energy at every time-step to "energy.txt"
 * @param dataTxt
 * @param localN    - Allocates the number of particles iterated per processor
 * @param nProc     - The number of processors available to iterate
 * @param currRank  - The rank of the current processor
 */
void SPH::iterate(ofstream& xCoor, ofstream& energyTxt, unsigned int localN, unsigned int nProc, unsigned int currRank) {
    
    // Iterates when the current time-step is less than the final time-step specified
    while (getCurrT() <= getTotalIntTime()) {
        
        // These operations are carried out with MPI parallelisation across various processors
        calcDensityWithMPI(localN, nProc, currRank);        // Calculating the density of particles
        calcPressureWithMPI(localN, nProc, currRank);       // Calculating the pressure of particles
        calcPressureForceWithMPI(localN, nProc, currRank);  // Calculating the pressure force acting on each particle
        calcViscousForceWithMPI(localN, nProc, currRank);   // Calculating the viscous force acting on each particle
        calcGravityForceWithMPI(localN, nProc, currRank);   // Calculating the gravitational force acting on each particle
        
        // These operations are carried out within the root processor
        if (currRank == 0) {
            calcAcceleration();                             // Calculating the particles' acceleration as a result of the force interactions
            getNextParticleVel();                           // Calculating the particles' velocity as a result of its acceleration
            getNextParticlePos();                           // Calculating the particles' resulting position
            
            writeEnergy(energyTxt);
        }
        
        setCurrT();                                         // Updating the current time-step
        
    }
    
    if (currRank == 0) {
        writeParticlePosition(xCoor);                       // Writing the particles' coordinates at final time-step to "output.txt"
    }
    
}

/**
 * @brief Calculating the particle's density based on its relative distance to other particles with parallelisation
 * @param localN    - Allocates the number of particles iterated per processor
 * @param nProc     - The number of processors available to iterate
 * @param currRank  - The rank of the current processor
 */
void SPH::calcDensityWithMPI(unsigned int localN, unsigned int nProc, unsigned int currRank) {
    
    // Ensuring all processors have the global set of particles' location
    MPI_Bcast(x, 2*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // Operations within the root
    if (currRank == 0){
            
        for (unsigned int i = 1; i < nProc; ++i) {
                
            startN = i * localN;
            endN = (i+1) * localN;
            
            // Sending the corresponding allocation of particles to different processors
            MPI_Send(&startN, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(&endN, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
             
        }
        
        xLocal = new double[2 * localN]();  // Initialising the local copy of the particles' location
        rhoLocal = new double[localN]();    // Initialising the local copy of the particles' density
            
        startN = currRank * localN;         // First particle this processor should iterate
        endN = (currRank + 1) * localN;     // Last particle this processor should iterate
        
        // Copying the relevant particle coordinates for this processor
        for (unsigned int i = 0; i < localN; ++i) {
                
            xLocal[2*i] = x[2*(i + startN)];
            xLocal[2*i+1] = x[2*(i + startN) + 1];
                
        }
            
        double q;
    
        for (unsigned int i = 0; i < localN; ++i) {
            
            for (unsigned int j = 0; j < N; ++j) {
                
                calcRij(i, j);
                q = sqrt(r_ij[2*j]*r_ij[2*j] + r_ij[2*j+1]*r_ij[2*j+1]) / h;
                
                if (q < 1) {
                    rhoLocal[i] += ((4 * m) / (M_PI * h * h)) * ((1 - (q * q)) * (1 - (q * q)) * (1 - (q * q)));
                } else {
                    rhoLocal[i] += 0.0;
                }
            
            }
            
        }
        
    } else {
        
        // Receiving the corresponding particle allocation
        MPI_Recv(&startN, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&endN, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
        xLocal = new double[2 * localN]();  // Initialising the local copy of the particles' location
        rhoLocal = new double[localN]();    // Initialising the local copy of the particles' density
            
        startN = currRank * localN;         // First particle this processor should iterate
        endN = (currRank + 1) * localN;     // First particle this processor should iterate
        
        // Copying the relevant particle coordinates for this processor    
        for (unsigned int i = 0; i < localN; ++i) {
                
            xLocal[2*i] = x[2*(i + startN)];
            xLocal[2*i+1] = x[2*(i + startN) + 1];
                
        }
            
        double q;
    
        for (unsigned int i = 0; i < localN; ++i) {
            
            for (unsigned int j = 0; j < N; ++j) {
                
                calcRij(i, j);
                q = sqrt(r_ij[2*j]*r_ij[2*j] + r_ij[2*j+1]*r_ij[2*j+1]) / h;
                
                if (q < 1) {
                    rhoLocal[i] += ((4 * m) / (M_PI * h * h)) * ((1 - (q * q)) * (1 - (q * q)) * (1 - (q * q)));
                } else {
                    rhoLocal[i] += 0.0;
                }
            
            }
            
        }
            
    }
    
    // Concatenating the calculated particle densities into the global "rho" array
    MPI_Gather(rhoLocal, localN, MPI_DOUBLE, rho, localN, MPI_DOUBLE, 0, MPI_COMM_WORLD);

}

/**
 * @brief Calculating the particle's pressure based on its relative density with parallelisation
 * @param localN    - Allocates the number of particles iterated per processor
 * @param nProc     - The number of processors available to iterate
 * @param currRank  - The rank of the current processor
 */
void SPH::calcPressureWithMPI(unsigned int localN, unsigned int nProc, unsigned int currRank) {
    
    MPI_Bcast(rho, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
    if (currRank == 0){
            
        for (unsigned int i = 1; i < nProc; ++i) {
                
            startN = i * localN;
            endN = (i+1) * localN;
                
            MPI_Send(&startN, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(&endN, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
             
        }
            
        pLocal = new double[localN]();
            
        startN = currRank * localN;
        endN = (currRank + 1) * localN;
            
        for (unsigned int i = 0; i < localN; ++i) {
      
            pLocal[i] = k * (rho[i+startN] - rho_0);
            
        } 
        
    } else {
        
        MPI_Recv(&startN, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&endN, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
        pLocal = new double[localN]();
            
        startN = currRank * localN;
        endN = (currRank + 1) * localN;
            
        for (unsigned int i = 0; i < localN; ++i) {
            pLocal[i] = k * (rho[i+startN] - rho_0);
        } 
            
    }
    
    // Concatenating the calculated particle pressures into the global "p" array
    MPI_Gather(pLocal, localN, MPI_DOUBLE, p, localN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
}

/**
 * @brief Calculating the pressure forces acting on the particles with parallelisation
 * @param localN    - Allocates the number of particles iterated per processor
 * @param nProc     - The number of processors available to iterate
 * @param currRank  - The rank of the current processor
 */
void SPH::calcPressureForceWithMPI(unsigned int localN, unsigned int nProc, unsigned int currRank) {
    
    MPI_Bcast(x, 2*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(rho, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(p, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
    if (currRank == 0){
            
        for (unsigned int i = 1; i < nProc; ++i) {
                
            startN = i * localN;
            endN = (i+1) * localN;
           
            MPI_Send(&startN, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(&endN, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
            
        }
            
        xLocal = new double[2 * localN]();
        F_pLocal = new double[2 * localN]();
            
        startN = currRank * localN;
        endN = (currRank + 1) * localN;
        
        for (unsigned int i = 0; i < localN; ++i) {
                
            xLocal[2*i] = x[2*(i + startN)];
            xLocal[2*i+1] = x[2*(i + startN) + 1];
                
        }
            
        double q;
    
        for (unsigned int i = 0; i < localN; ++i) {
        
            for (unsigned int j = 0; j < N; ++j) {
            
                calcRij(i, j);
                q = sqrt(r_ij[2*j]*r_ij[2*j] + r_ij[2*j+1]*r_ij[2*j+1]) / h;
                
                if (q < 1 && q != 0) {
                    
                    F_pLocal[2*i] += ((m / rho[j]) * ((p[i+startN] + p[j]) / 2.0)) * ((-30 / (M_PI * h * h * h)) * (r_ij[2*j]) * (((1 - q) * (1 - q)) / q));
                    F_pLocal[2*i+1] += ((m / rho[j]) * ((p[i+startN] + p[j]) / 2.0)) * ((-30 / (M_PI * h * h * h)) * (r_ij[2*j+1]) * (((1 - q) * (1 - q)) / q));
                    
                } else {
                    
                    F_pLocal[2*i] += 0.0;
                    F_pLocal[2*i+1] += 0.0;
                    
                }
            
            }
            
            F_pLocal[2*i] *= -1.0;
            F_pLocal[2*i+1] *= -1.0;
            
        }
        
    } else {
        
        MPI_Recv(&startN, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&endN, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
        xLocal = new double[2 * localN]();
        F_pLocal = new double[2 * localN]();
            
        startN = currRank * localN;
        endN = (currRank + 1) * localN;
            
        for (unsigned int i = 0; i < localN; ++i) {
                
            xLocal[2*i] = x[2*(i + startN)];
            xLocal[2*i+1] = x[2*(i + startN) + 1];
                
        }
            
        double q;
    
        for (unsigned int i = 0; i < localN; ++i) {
        
            for (unsigned int j = 0; j < N; ++j) {
            
                calcRij(i, j);
                q = sqrt(r_ij[2*j]*r_ij[2*j] + r_ij[2*j+1]*r_ij[2*j+1]) / h;
                
                if (q < 1 && q != 0) {
                    
                    F_pLocal[2*i] += ((m / rho[j]) * ((p[i+startN] + p[j]) / 2.0)) * ((-30 / (M_PI * h * h * h)) * (r_ij[2*j]) * (((1 - q) * (1 - q)) / q));
                    F_pLocal[2*i+1] += ((m / rho[j]) * ((p[i+startN] + p[j]) / 2.0)) * ((-30 / (M_PI * h * h * h)) * (r_ij[2*j+1]) * (((1 - q) * (1 - q)) / q));
                    
                } else {
                    
                    F_pLocal[2*i] += 0.0;
                    F_pLocal[2*i+1] += 0.0;
                    
                }
            
            }
            
            F_pLocal[2*i] *= -1.0;
            F_pLocal[2*i+1] *= -1.0;
            
        }
            
    }
    
    // Concatenating the calculated pressure forces into the global "F_p" array
    MPI_Gather(F_pLocal, 2*localN, MPI_DOUBLE, F_p, 2*localN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
}

/**
 * @brief Calculating the viscous forces acting on the particles with parallelisation
 * @param localN    - Allocates the number of particles iterated per processor
 * @param nProc     - The number of processors available to iterate
 * @param currRank  - The rank of the current processor
 */
void SPH::calcViscousForceWithMPI(unsigned int localN, unsigned int nProc, unsigned int currRank) {
    
    MPI_Bcast(x, 2*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(rho, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(v, 2*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
    if (currRank == 0){
            
        for (unsigned int i = 1; i < nProc; ++i) {
                
            startN = i * localN;
            endN = (i+1) * localN;
                
            MPI_Send(&startN, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(&endN, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
             
        }
            
        xLocal = new double[2 * localN]();
        vLocal = new double[2 * localN]();
        F_vLocal = new double[2 * localN]();
            
        startN = currRank * localN;
        endN = (currRank + 1) * localN;
            
        for (unsigned int i = 0; i < localN; ++i) {
                
            xLocal[2*i] = x[2*(i + startN)];
            xLocal[2*i+1] = x[2*(i + startN) + 1];
            
            vLocal[2*i] = v[2*(i + startN)];
            vLocal[2*i+1] = v[2*(i + startN) + 1];
                
        }
            
        double q;
    
        for (unsigned int i = 0; i < localN; ++i) {
        
            for (unsigned int j = 0; j < N; ++j) {
            
                calcVij(i, j);
                calcRij(i, j);
                q = sqrt(r_ij[2*j]*r_ij[2*j] + r_ij[2*j+1]*r_ij[2*j+1]) / h;
                
                if (q < 1 && q != 0) {
                    
                    F_vLocal[2*i] += (m / rho[j]) * (v_ij[2*j]) * ((40 / (M_PI * h * h * h * h)) * (1 - q));
                    F_vLocal[2*i+1] += (m / rho[j]) * (v_ij[2*j+1]) * ((40 / (M_PI * h * h * h * h)) * (1 - q));
                    
                } else {
                    
                    F_vLocal[2*i] += 0.0;
                    F_vLocal[2*i+1] += 0.0;
                    
                }
            
            }
            
            F_vLocal[2*i] *= -mu;
            F_vLocal[2*i+1] *= -mu;
            
        }
        
    } else {
        
        MPI_Recv(&startN, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&endN, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
        xLocal = new double[2 * localN]();
        vLocal = new double[2 * localN]();
        F_vLocal = new double[2 * localN]();
            
        startN = currRank * localN;
        endN = (currRank + 1) * localN;
            
        for (unsigned int i = 0; i < localN; ++i) {
                
            xLocal[2*i] = x[2*(i + startN)];
            xLocal[2*i+1] = x[2*(i + startN) + 1];
            
            vLocal[2*i] = v[2*(i + startN)];
            vLocal[2*i+1] = v[2*(i + startN) + 1];
                
        }
            
        double q;
    
        for (unsigned int i = 0; i < localN; ++i) {
        
            for (unsigned int j = 0; j < N; ++j) {
            
                calcVij(i, j);
                calcRij(i, j);
                q = sqrt(r_ij[2*j]*r_ij[2*j] + r_ij[2*j+1]*r_ij[2*j+1]) / h;
                
                if (q < 1 && q != 0) {
                    
                    F_vLocal[2*i] += (m / rho[j]) * (v_ij[2*j]) * ((40 / (M_PI * h * h * h * h)) * (1 - q));
                    F_vLocal[2*i+1] += (m / rho[j]) * (v_ij[2*j+1]) * ((40 / (M_PI * h * h * h * h)) * (1 - q));
                    
                } else {
                    
                    F_vLocal[2*i] += 0.0;
                    F_vLocal[2*i+1] += 0.0;
                    
                }
            
            }
            
            F_vLocal[2*i] *= -mu;
            F_vLocal[2*i+1] *= -mu;
            
        }
            
    }
    
    MPI_Gather(F_vLocal, 2*localN, MPI_DOUBLE, F_v, 2*localN, MPI_DOUBLE, 0, MPI_COMM_WORLD);

}

/**
 * @brief Calculating the gravitational forces acting on the particles with parallelisation
 * @param localN    - Allocates the number of particles iterated per processor
 * @param nProc     - The number of processors available to iterate
 * @param currRank  - The rank of the current processor
 */
void SPH::calcGravityForceWithMPI(unsigned int localN, unsigned int nProc, unsigned int currRank) {
    
    MPI_Bcast(rho, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
    if (currRank == 0){
            
        for (unsigned int i = 1; i < nProc; ++i) {
                
            startN = i * localN;
            endN = (i+1) * localN;
            
            MPI_Send(&startN, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
            MPI_Send(&endN, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
            
        }
            
        F_gLocal = new double[2 * localN]();
            
        startN = currRank * localN;
        endN = (currRank + 1) * localN;
            
        for (unsigned int i = 0; i < localN; ++i) {
        
            F_gLocal[2*i] = 0.0;
            F_gLocal[2*i+1] = -rho[i+startN] * g;
            
        }
        
    } else {
        
        MPI_Recv(&startN, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&endN, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
        F_gLocal = new double[2 * localN]();
            
        startN = currRank * localN;
        endN = (currRank + 1) * localN;
            
        for (unsigned int i = 0; i < localN; ++i) {
        
            F_gLocal[2*i] = 0.0;
            F_gLocal[2*i+1] = -rho[i+startN] * g;
            
        }
            
    }
    
    MPI_Gather(F_gLocal, 2*localN, MPI_DOUBLE, F_g, 2*localN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
}

/**
 * @brief Calculating the acceleration of each particle due to the forces (relative to the particle density)
 */
void SPH::calcAcceleration() {
    
    for (unsigned int i = 0; i < N; ++i) {
        
        a[2*i] = (F_p[2*i] + F_v[2*i] + F_g[2*i]) / rho[i];
        a[2*i+1] = (F_p[2*i+1] + F_v[2*i+1] + F_g[2*i+1]) / rho[i];
        
    }
    
    
}

/**
 * @brief Calculating the velocity of each particle using BLAS
 */
void SPH::getNextParticleVel() {
    
    if (getCurrT() == 0.0) {                        // Initiate velocity (at half time-step) if yet to do so
        cblas_daxpy(2*N, (0.5*dt), a, 1, v, 1);
    } else {                                        // Calculate velocity based on time-step and previous particle velocity
        cblas_daxpy(2*N, dt, a, 1, v, 1);
    }
    
}

/**
 * @brief Calculating the particles location at the next time-step
 */
void SPH::getNextParticlePos() {
    
    for (unsigned int i = 0; i < N; ++i) {
        
        x[2*i] += v[2*i] * dt;
        x[2*i+1] += v[2*i+1] * dt;
        
        applyBC(i);                                 // Altering particle position and velocity should it hit the boundary of the domain
        
    }
    
}

/**
 * @brief Updates the particle's location and velocity should it hit the boundary of the domain
 * @param i
 */
void SPH::applyBC(unsigned int i) {
    
    if (x[2*i] < h) {                               // boundary x = 0;
        
        x[2*i] = h;                                 // Repositions the particle on the boundary
        v[2*i] *= -e;                               // Reverses the direction of travel for the particle and removing energy due to collision with boundary
        
    } else if (x[2*i] > (1 - h)) {                  // boundary x = 1;
        
        x[2*i] = 1 - h;
        v[2*i] *= -e;
            
    }
        
    if (x[2*i+1] < h) {                             // boundary y = 0;
            
        x[2*i+1] = h;
        v[2*i+1] *= -e; 
            
    } else if (x[2*i+1] > (1 - h)) {                // boundary y = 1;
            
        x[2*i+1] = 1 - h;
        v[2*i+1] *= -e;
            
    }
    
}

/**
 * @brief Wrting the particles' coordinates at final time, T, to "output.txt"
 * @param xCoor - Reference variable to "output.txt"
 */
void SPH::writeParticlePosition(ofstream& xCoor) {
    
    if (xCoor.is_open()) {
        
        for (unsigned int i = 0; i < N; ++i) {
            xCoor << x[2*i] << " " << x[2*i+1] << endl;
        }
        
    } else {
        cout << "Failed to open output.txt" << endl;
    }
    
}

/**
 * @brief Wrting the system's energy breakdown at at each time-step, t, to "energy.txt"
 * @param xCoor - Reference variable to "energy.txt"
 */
void SPH::writeEnergy(ofstream& energyTxt) {
    
    if (energyTxt.is_open()) {
        energyTxt << getCurrT() << " " << calcKineticEnergy() << " " << calcPotentialEnergy() << " " << calcTotalEnergy() << endl;
    } else {
        cout << "Failed to open energy.txt" << endl;
    }

}

/**
 * @brief Closing the text files
 * @param fileName - Reference variable to the file that is to be closed
 */
void SPH::closePPOutputFile(ofstream& fileName) {fileName.close();}

/**
 * @brief Calculating the total kinetic energy of the system at time t
 * @return Total Kinetic Energy of the system at time t
 */
double SPH::calcKineticEnergy() {
    
    double E_k = 0.0;
    
    for (unsigned int i = 0; i < N; ++i) {
        E_k += (v[2*i] * v[2*i]) + v[2*i+1] * v[2*i+1];
    }
    
    return 0.5 * m * E_k;
    
}

/**
 * @brief Calculating the total potential energy of the system at time t
 * @return Total Potential Energy of the system at time t
 */
double SPH::calcPotentialEnergy() {
    
    double E_p = 0.0;
    
    for (unsigned int i = 0; i < N; ++i) {
        E_p += x[2*i+1];
    }
    
    return m * g * E_p;
    
}

/**
 * @brief Summing the total kinetic and potential energy of the system at time t
 * @return Total Energy of the system at time t
 */
double SPH::calcTotalEnergy() {return calcKineticEnergy() + calcPotentialEnergy();}

/**
 * @brief Getter for number of particles, N
 * @return N
 */
int SPH::getN() {return N;}

/**
 * @brief Getter for current time within iteration, t
 * @return t
 */
double SPH::getCurrT() {return t;}

/**
 * @brief Setter for current incremental time-step progression within iteration, t
 */
void SPH::setCurrT() {t += dt;}

/**
 * @brief Getter for time-step, dt
 * @return dt
 */
double SPH::getDeltaT() {return dt;}

/**
 * @brief Getter for final time of iteration, T
 * @return T
 */
double SPH::getTotalIntTime() {return T;}